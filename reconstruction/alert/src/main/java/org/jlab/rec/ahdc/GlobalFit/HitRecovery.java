package org.jlab.rec.ahdc.GlobalFit;

import java.util.Arrays;

/**
 * Hit recovery around the global fit: port of recover.py (the GE9 procedure, frozen windows of ge9_recover.json held by
 * {@link GFitConfig}) for one track. Per round: crossed layers of the helix of the fitted state without a track hit are missing
 * layers (recover.crossings, recover.py:140-166 with fit._layer_point, fit.py:503); the event's pool rows (no finder track) within
 * +-poolW wires of the crossing wire (build_pool, recover.py:182-209) are evaluated as hits of the fitted state with the ATOF term
 * off (probe, recover.py:224-256); the best row per missing layer inside the windows (pass_windows, select, recover.py:274-292)
 * is added and the track refitted (the full search, plan D4 as amended); the refit is accepted when it converged and
 * F1 - F0 <= q n_add (_accept, recover.py:333-338); a rejected refit with two or more added rows is retried once without the
 * row of lowest likelihood in the refit's hit list (one_round, recover.py:359-387). A second round runs on accepted tracks only, on top
 * of the accepted rows and against the new result (recover, recover.py:392-420). Immutable and thread safe: all scratch is in
 * {@link GFitWorkspace#rec}; only the refit track objects are allocated, and only for a track with a selected row.
 */
public final class HitRecovery {

	private static final int NPAR = GFitWorkspace.NPAR;
	private static final int NL = 8;

	private final GlobalTrackFitter fitter;
	private final TrackCost cost;
	private final CellGeometry geo;
	private final GFitConfig cfg;
	private final double bzHelix, beamX, beamY, errScale, tLo, tHi, totMin, adcMin, wfMax, kRes, q;
	private final int rounds, poolW;
	private final boolean dropWorst;

	/** Per-thread scratch of one {@link #recover} call, with the accepted rows of the last call (arrays valid up to nAcc). */
	public static final class Scratch {
		/** accepted rows in Python's order (per round, layers ascending): row in event, layer index, wire, doca, smear class, sigma, round (0, 1). */
		public int nAcc;
		public final int[] accEv = new int[2 * NL], accLayer = new int[2 * NL], accWire = new int[2 * NL], accIsm = new int[2 * NL], accRound = new int[2 * NL];
		public final double[] accDoca = new double[2 * NL], accSigma = new double[2 * NL];
		/** drop-worst retries of the last call: number, and per round the event row dropped (-1 no retry). */
		public int retries;
		public final int[] dropEv = new int[2];

		// crossings of the current state: reaches layer with z inside the wire (crossed), crossing wire (gstar, -1), track has a hit, missing = crossed and no hit
		boolean ok;
		final boolean[] crossed = new boolean[NL], has = new boolean[NL], missing = new boolean[NL];
		final int[] gstar = new int[NL];
		// pool candidates (probe table, in Python's order: layer ascending, event order): index in the event pool, layer, wire offset, wire, ism, doca, sigma, probe l, d, s
		int nc;
		int[] cIdx = new int[64], cLi = new int[64], cDw = new int[64], cG = new int[64], cIsm = new int[64];
		double[] cDoca = new double[64], cSig = new double[64], cL = new double[64], cD = new double[64], cS = new double[64];
		boolean[] cPass = new boolean[64];
		/** selected candidates (best per missing layer, layer ascending) and the kept subset of a retry. */
		int nsel;
		final int[] sel = new int[NL], keep = new int[NL];
		final double[] u = new double[NPAR], xyz = new double[3];
		final GFitResult r1, r2;

		public Scratch(int maxHits) {
			r1 = new GFitResult(maxHits);
			r2 = new GFitResult(maxHits);
		}

		private void grow(int n) {
			if (n <= cIdx.length) return;
			int m = Math.max(n, 2 * cIdx.length);
			cIdx = Arrays.copyOf(cIdx, m); cLi = Arrays.copyOf(cLi, m); cDw = Arrays.copyOf(cDw, m);
			cG = Arrays.copyOf(cG, m); cIsm = Arrays.copyOf(cIsm, m);
			cDoca = Arrays.copyOf(cDoca, m); cSig = Arrays.copyOf(cSig, m); cL = Arrays.copyOf(cL, m);
			cD = Arrays.copyOf(cD, m); cS = Arrays.copyOf(cS, m); cPass = Arrays.copyOf(cPass, m);
		}
	}

	public HitRecovery(GlobalTrackFitter fitter, TrackCost cost, CellGeometry geo, GFitConfig cfg) {
		this.fitter = fitter;
		this.cost = cost;
		this.geo = geo;
		this.cfg = cfg;
		bzHelix = cfg.bzHelix(); beamX = cfg.beamX(); beamY = cfg.beamY(); errScale = cfg.errScale();
		tLo = cfg.recTLo(); tHi = cfg.recTHi(); totMin = cfg.recTotMin(); adcMin = cfg.recAdcMin(); wfMax = cfg.recWfMax(); kRes = cfg.recKRes();
		q = cfg.recQ(); rounds = cfg.recRounds(); poolW = cfg.recPoolW(); dropWorst = cfg.recDropWorst();
	}

	/**
	 * Recovery of one track: base is its fit (t's hits), out receives the result after recovery (a copy of base when nothing is
	 * accepted; base and out may be the same object). Returns the number of rows added in total; the accepted rows are in
	 * {@code ws.rec}. A base with status other than 0 or a non-finite state adds nothing (recover.internal_u).
	 */
	public int recover(GFitTrack t, GFitEvent ev, AtofPositionTable tab, boolean bankTimes, GFitResult base, GFitWorkspace ws, GFitResult out) {
		final Scratch s = ws.rec;
		s.nAcc = 0; s.retries = 0; s.dropEv[0] = -1; s.dropEv[1] = -1;
		if (base != out) out.copyFrom(base);
		GFitTrack cur = t;
		int total = 0;
		for (int rnd = 0; rnd < rounds; rnd++) {
			int ns = prepare(cur, ev, out, ws, s);
			if (ns == 0) break;
			final double f0 = out.cost;
			final int st0 = out.status;
			GFitTrack t1 = extend(cur, ev, s, s.sel, ns);
			fitter.fit(t1, ev, tab, bankTimes, ws, s.r1);
			GFitTrack acc = null;
			GFitResult accRes = null;
			int[] kept = s.sel;
			int nk = ns;
			if (accept(s.r1, f0, st0, ns)) {
				acc = t1; accRes = s.r1;
			} else if (dropWorst && ns >= 2) {
				// likelihood of each added row in the refit's hit list, matched by event row (recover.added_l); a row missing
				// from it (a failed refit keeps no hits) or a non-finite value counts as -1e300
				int worst = 0;
				double lw = 0.0;
				for (int i = 0; i < ns; i++) {
					int row = ev.pEv[s.cIdx[s.sel[i]]];
					double l = Double.NaN;
					for (int k = 0; k < s.r1.nh; k++)
						if (s.r1.hitKind[k] != 2 && s.r1.hitRow[k] == row) { l = s.r1.hitL[k]; break; }
					if (!(l - l == 0.0)) l = -1e300;
					if (i == 0 || l < lw) { lw = l; worst = i; }
				}
				nk = 0;
				for (int i = 0; i < ns; i++) if (i != worst) s.keep[nk++] = s.sel[i];
				kept = s.keep;
				s.retries++;
				s.dropEv[rnd] = ev.pEv[s.cIdx[s.sel[worst]]];
				GFitTrack t2 = extend(cur, ev, s, kept, nk);
				fitter.fit(t2, ev, tab, bankTimes, ws, s.r2);
				if (accept(s.r2, f0, st0, nk)) { acc = t2; accRes = s.r2; }
			}
			if (acc == null) break;
			out.copyFrom(accRes);
			cur = acc;
			for (int i = 0; i < nk; i++) {
				int c = kept[i], a = s.nAcc++;
				s.accEv[a] = ev.pEv[s.cIdx[c]]; s.accLayer[a] = s.cLi[c]; s.accWire[a] = s.cG[c]; s.accIsm[a] = s.cIsm[c];
				s.accDoca[a] = s.cDoca[c]; s.accSigma[a] = s.cSig[c]; s.accRound[a] = rnd;
			}
			total += nk;
		}
		return total;
	}

	/** recover._accept (recover.py:333-338). */
	private boolean accept(GFitResult r, double f0, int st0, int nadd) {
		if (r.status != 0) return false;
		if (st0 != 0 || !(f0 - f0 == 0.0)) return true;
		return (r.cost - f0) <= q * nadd;
	}

	/** the track with the selected candidates cand[0..n) appended to cur's hits (Python's augment: original hits first, then the added rows in layer order). */
	GFitTrack extend(GFitTrack cur, GFitEvent ev, Scratch s, int[] cand, int n) {
		int nb = cur.nHits(), nt = nb + n;
		int[] hg = Arrays.copyOf(cur.hg, nt), hism = Arrays.copyOf(cur.hism, nt), hcls = Arrays.copyOf(cur.hcls, nt),
			hrow = Arrays.copyOf(cur.hrow, nt);
		double[] hm = Arrays.copyOf(cur.hm, nt), hsig = Arrays.copyOf(cur.hsig, nt), hh = Arrays.copyOf(cur.hh, nt);
		for (int i = 0; i < n; i++) {
			int c = cand[i], k = nb + i;
			hg[k] = s.cG[c]; hm[k] = s.cDoca[c]; hsig[k] = s.cSig[c]; hh[k] = geo.half[s.cG[c]]; hism[k] = s.cIsm[c]; hcls[k] = s.cLi[c];
			hrow[k] = ev.pEv[s.cIdx[c]];
		}
		return new GFitTrack(cur.trackid, hg, hm, hsig, hh, hism, hcls, hrow, cur.seedPhi, cur.seedP, cur.seedTheta, cur.seedZ);
	}

	// ================================================================== crossings, pool, probe, selection of one round

	/**
	 * Crossings, pool, probe and selection of one round at the state of res for the hits of cur: fills the scratch (crossed, gstar,
	 * has, missing, the candidate table with probe values and window flags, sel) and returns the number of selected rows.
	 */
	int prepare(GFitTrack cur, GFitEvent ev, GFitResult res, GFitWorkspace ws, Scratch s) {
		crossings(cur, res, ws, s);
		s.nc = 0; s.nsel = 0;
		if (!s.ok) return 0;
		final int np = ev.nPool();
		final int[] lo = geo.layoff;
		for (int li = 0; li < NL; li++) {
			if (!s.missing[li]) continue;
			int nw = lo[li + 1] - lo[li], ic = s.gstar[li] - lo[li];
			for (int k = 0; k < np; k++) {
				if (ev.pLi[k] != li) continue;
				int dw = Math.floorMod(ev.pG[k] - lo[li] - ic + nw / 2, nw) - nw / 2;
				if (Math.abs(dw) > poolW) continue;
				s.grow(s.nc + 1);
				int c = s.nc++;
				double tim = ev.pTime[k];
				s.cIdx[c] = k; s.cLi[c] = li; s.cDw[c] = dw; s.cG[c] = ev.pG[k];
				s.cDoca[c] = ev.pDoca[k];
				s.cSig[c] = errScale * kfError(ev.pAdc[k], tim);
				s.cIsm[c] = tim < -20 ? 2 : (tim < 0 ? 1 : 0);
			}
		}
		if (s.nc == 0) return 0;
		probe(s.u, ws, s);
		// windows (pass_windows) and the best passing row per layer: largest l, the first of equal values (stable lexsort)
		int[] best = s.keep;
		Arrays.fill(best, -1);
		for (int c = 0; c < s.nc; c++) {
			int k = s.cIdx[c];
			double res_ = s.cDoca[c] - Math.abs(s.cD[c]);
			boolean pass = ev.pTime[k] >= tLo && ev.pTime[k] <= tHi && ev.pTot[k] >= totMin && ev.pAdc[k] >= adcMin && ev.pWf[k] <= wfMax
				&& s.cS[c] > 0 && finite(res_) && Math.abs(res_) <= kRes * s.cSig[c] && finite(s.cL[c]);
			s.cPass[c] = pass;
			if (pass && (best[s.cLi[c]] < 0 || s.cL[c] > s.cL[best[s.cLi[c]]])) best[s.cLi[c]] = c;
		}
		for (int li = 0; li < NL; li++) if (best[li] >= 0) s.sel[s.nsel++] = best[li];
		return s.nsel;
	}

	/** recover.crossings for one track: ok (status 0 and finite state), crossed layers, crossing wires, track hits per layer; s.u = the internal state. */
	private void crossings(GFitTrack cur, GFitResult res, GFitWorkspace ws, Scratch s) {
		Arrays.fill(s.crossed, false); Arrays.fill(s.has, false); Arrays.fill(s.missing, false);
		Arrays.fill(s.gstar, -1);
		for (int k = 0; k < cur.nHits(); k++) s.has[geo.layerOf(cur.hg[k])] = true;
		s.ok = loadState(res, s) && res.status == 0;
		if (!s.ok) return;
		final double[] wo = geo.wo, we = geo.we;
		final int[] lo = geo.layoff;
		for (int li = 0; li < NL; li++) {
			if (!layerPoint(s.u, geo.lrad[li], ws.o9, s.xyz)) continue;
			double x = s.xyz[0], y = s.xyz[1], z = s.xyz[2];
			if (!finite(x)) continue;
			int best = -1;
			double bd = 1e300;
			for (int g = lo[li]; g < lo[li + 1]; g++) {
				double tt = (z - wo[3 * g + 2]) / (we[3 * g + 2] - wo[3 * g + 2]);
				double xx = wo[3 * g] + (we[3 * g] - wo[3 * g]) * tt - x, yy = wo[3 * g + 1] + (we[3 * g + 1] - wo[3 * g + 1]) * tt - y;
				double dd = xx * xx + yy * yy;
				if (dd < bd) { bd = dd; best = g; }
			}
			s.gstar[li] = best;
			if (best < 0) continue;
			double zlo = Math.min(wo[3 * best + 2], we[3 * best + 2]), zhi = Math.max(wo[3 * best + 2], we[3 * best + 2]);
			s.crossed[li] = zlo <= z && z <= zhi;
			s.missing[li] = s.crossed[li] && !s.has[li];
		}
	}

	/** s.u = state of res (kinks: nan_to_num); true when (d0, phi0, ln p, z0, theta) are finite (internal_u, recover.py:169-178). */
	private static boolean loadState(GFitResult res, Scratch s) {
		boolean fin = true;
		for (int j = 0; j < NPAR; j++) {
			double v = res.u[j];
			if (j < 5) fin &= finite(v);
			else if (!(v - v == 0.0)) v = 0.0;
			s.u[j] = v;
		}
		return fin;
	}

	/** fit._layer_point (fit.py:503-515): first helix point at transverse radius R (helix with Bz = bzHelix) into xyz; false when none. */
	private boolean layerPoint(double[] u, double R, double[] o, double[] xyz) {
		double prev = -1.0;
		for (int j = 0; j < GFitConfig.NG; j++) {
			HelixModel.pt(u, j * GFitConfig.DS, bzHelix, beamX, beamY, o);
			double r = Math.sqrt(o[0] * o[0] + o[1] * o[1]);
			if (r >= R && j > 0) {
				double f = (R - prev) / (r - prev);
				HelixModel.pt(u, (j - 1 + f) * GFitConfig.DS, bzHelix, beamX, beamY, o);
				xyz[0] = o[0]; xyz[1] = o[1]; xyz[2] = o[2];
				return true;
			}
			prev = r;
		}
		return false;
	}

	/**
	 * recover.probe for one track: (l, signed d, cell margin s) of the candidates 0..nc-1 as hits of the state u, propagation model,
	 * hit term only (no ATOF term). Processed in chunks of the workspace hit capacity.
	 */
	private void probe(double[] u, GFitWorkspace ws, Scratch s) {
		ws.atofOn = false;
		final double[] b3 = ws.box3;
		final int cap = ws.cap;
		for (int a = 0; a < s.nc; a += cap) {
			int m = Math.min(cap, s.nc - a);
			for (int k = 0; k < m; k++) {
				int c = a + k, g = s.cG[c], li = s.cLi[c];
				ws.lg[k] = g; ws.lli[k] = li; ws.lw[k] = g - geo.layoff[li]; ws.lm[k] = s.cDoca[c]; ws.lism[k] = s.cIsm[c];
				ws.lsig[k] = s.cSig[c]; ws.lh[k] = geo.half[g]; ws.lcls[k] = li; ws.lkind[k] = 0; ws.lrow[k] = -1;
			}
			ws.nh = m;
			cost.predict(u, 1, NPAR, false, true, ws.th0, ws);
			for (int k = 0; k < m; k++) {
				double d = ws.D[k];
				HitTerm.dbox(d, ws.lm[k], ws.lsig[k], ws.lh[k], false, cfg, b3);
				s.cL[a + k] = -0.5 * b3[0]; s.cD[a + k] = d; s.cS[a + k] = ws.S[k];
			}
		}
	}

	/** Hit.getMeasurementNoiseMatrix sigma [mm] before errscale (fit.kf_error, fit.py:930-934). */
	public static double kfError(double adc, double t) {
		double ea = (1.15146 * adc + 437.63) / (3.21187 * adc + 878.855);
		double et = (0.4423 * t + 13.7215) / (0.846038 * t + 31.9867);
		return ea * et / 0.471;
	}

	private static boolean finite(double x) { return x - x == 0.0; }
}
