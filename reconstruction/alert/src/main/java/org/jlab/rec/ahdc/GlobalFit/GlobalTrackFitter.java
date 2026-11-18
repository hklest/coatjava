package org.jlab.rec.ahdc.GlobalFit;

/**
 * Global track fit of one track: port of fit._fit_one (fit.py:600-856) for the configuration in use (propagation model,
 * dbox hit term, no missing wires, no added hits, ATOF position term on), with the start grid fit._grid (fit.py:552-596),
 * the outer theta0 passes fit._outer (fit.py:859-878), the near-wire flips, the covariance 2 H^-1 of the exact Hessian by
 * eigen decomposition and the status, flag and counter rules of fit.py:819-856. Search mode FULL is the Python sequence
 * (3 grid starts plus the seed start, helix prefit, outer passes, flips). Search mode FAST (plan D5) runs the seeded fit of
 * the rule in {@link GFitConfig#fastRule()} and falls back to a from-scratch FULL fit when the rule's trigger fires or the seed
 * momentum is not finite or at most 0.03 GeV; rule 6 is the full search for every track. Immutable and thread safe; every
 * scratch array is in the {@link GFitWorkspace} and nothing is allocated inside {@link #fit}. Any exception inside a track
 * gives status 4.
 */
public final class GlobalTrackFitter {

	private static final int NPAR = GFitWorkspace.NPAR, NPL = GFitWorkspace.NPL;

	private final TrackCost cost;
	private final Propagator prop;
	private final CellGeometry geo;
	private final GFitConfig cfg;
	private final LevMar lm;
	private final double bzHelix, beamX, beamY, sigZ, tolCost;
	private final boolean zOn, prefit, atofOn;
	private final int maxit, npass;
	private final double b0Pscale;

	/** The ATOF term (static class AtofTerm) needs the propagator and the grid the wire table, which TrackCost keeps private. */
	public GlobalTrackFitter(TrackCost cost, Propagator prop, CellGeometry geo, GFitConfig cfg) {
		this.cost = cost;
		this.prop = prop;
		this.geo = geo;
		this.cfg = cfg;
		this.lm = new LevMar(cost, cfg);
		bzHelix = cfg.bzHelix(); beamX = cfg.beamX(); beamY = cfg.beamY(); sigZ = cfg.sigZ(); tolCost = cfg.tolCost();
		zOn = cfg.zOn(); prefit = cfg.prefit(); atofOn = cfg.atofOn();
		maxit = cfg.maxit(); npass = cfg.npass(); b0Pscale = cfg.b0Pscale();
	}

	/**
	 * Fits one track into out (cleared first). tab is the ATOF position table of the event's half, bankTimes selects the bank
	 * ATOF times. Never throws: an exception inside the fit gives status 4.
	 */
	public void fit(GFitTrack t, GFitEvent ev, AtofPositionTable tab, boolean bankTimes, GFitWorkspace ws, GFitResult out) {
		try {
			ws.atofOn = atofOn;
			if (atofOn) AtofTerm.attach(t, ev, tab.par(), tab.recalOffset(), tab.recalScale(), bankTimes, prop, ws, ws.atof);
			if (cfg.search() == GFitConfig.Search.FULL) core(t, ev, ws, out, cfg.starts().code, cfg.nFlip());
			else fast(t, ev, ws, out);
		} catch (RuntimeException e) {
			out.clear();
			out.status = 4;
		}
	}

	// ================================================================== fast search (plan D5, D3 ladder)

	private void full(GFitTrack t, GFitEvent ev, GFitWorkspace ws, GFitResult out) {
		core(t, ev, ws, out, GFitConfig.Starts.GRID_B0.code, cfg.nFlip());
	}

	private boolean trigG(GFitResult r) { return r.status != 0; }

	/** G' = G, or flag 2 (grid edge), or cost per hit above c95 (search_measure.py trig_Gp). */
	private boolean trigGp(GFitResult r, int nHits) { return r.status != 0 || (r.flags & 2) != 0 || r.cost / nHits > cfg.c95(); }

	private void fast(GFitTrack t, GFitEvent ev, GFitWorkspace ws, GFitResult out) {
		int rule = cfg.fastRule();
		double p = t.seedP;
		int b0 = GFitConfig.Starts.B0.code, nb = t.nHits();
		if (Double.isNaN(p) || Double.isInfinite(p) || p <= 0.03 || rule < 1 || rule > 5) {
			full(t, ev, ws, out);
			return;
		}
		switch (rule) {
			case 1: case 2:
				core(t, ev, ws, out, b0, 0);
				if (rule == 1 ? trigG(out) : trigGp(out, nb)) full(t, ev, ws, out);
				break;
			case 3:     // B; C where G' fires on B's result (C's flips continue from B's state), the full search where G fires on C's
				core(t, ev, ws, out, b0, 0);
				if (trigGp(out, nb)) {
					core(t, ev, ws, out, b0, 2);
					if (trigG(out)) full(t, ev, ws, out);
				}
				break;
			default:
				core(t, ev, ws, out, b0, 2);
				if (rule == 4 ? trigG(out) : trigGp(out, nb)) full(t, ev, ws, out);
		}
	}

	// ================================================================== one fit (fit._fit_one)

	/** the final state of the fit of one track from the given start modes (bit 1 grid, bit 2 seed) and flip count. */
	private void core(GFitTrack t, GFitEvent ev, GFitWorkspace ws, GFitResult out, int smode, int nFlip) {
		out.clear();
		cost.loadHits(t, ws);
		final int nh = ws.nh;
		final int npar = NPAR;
		final double ze = ev.ze;
		final double[] th0f = ws.th0f;
		java.util.Arrays.fill(th0f, 0.0);
		if (nh < 3) {
			out.status = 2;
			return;
		}
		// ---- starts
		double[] starts = ws.starts;
		java.util.Arrays.fill(starts, 0.0);
		int nst = 0;
		if ((smode & 1) != 0) {
			int ns = grid(ze, t.seedPhi, ws);
			for (int j = 0; j < ns; j++) {
				for (int c = 0; c < 5; c++) starts[nst * NPAR + c] = ws.bestU[j * 5 + c];
				ws.startKind[nst] = 0;
				nst++;
			}
		}
		if ((smode & 2) != 0 && t.seedP > 0.03 && !Double.isInfinite(t.seedP) && !Double.isNaN(t.seedP)) {
			int o = nst * NPAR;
			starts[o] = 0.0;
			starts[o + 1] = t.seedPhi;
			starts[o + 2] = Math.log(t.seedP * b0Pscale);
			starts[o + 3] = t.seedZ;
			starts[o + 4] = finite(t.seedTheta) ? t.seedTheta : Math.PI / 2;
			if (!finite(starts[o + 3])) starts[o + 3] = ze;
			ws.startKind[nst] = 1;
			nst++;
		}
		out.nstart = nst;
		// ---- fit each start (helix prefit first), keep minima
		final double[] minU = ws.minU, minF = ws.minF, u = ws.u;
		java.util.Arrays.fill(minU, 0.0);
		java.util.Arrays.fill(minF, 1e300);
		java.util.Arrays.fill(ws.minC, false);
		int nmin = 0, ntot = 0;
		for (int j = 0; j < nst; j++) {
			System.arraycopy(starts, j * NPAR, u, 0, NPAR);
			if (!LevMar.inb(u)) continue;
			for (int k = 0; k < nh; k++) ws.sarc[k] = -1.0;
			if (prefit) {
				lm.run(u, 0, 5, ze, th0f, maxit, ws);
				ntot += ws.lmIt;
				boolean dup = false;
				for (int m = 0; m < nmin; m++) {     // identical helix minimum already refined
					double mx = Math.abs(minU[m * NPAR] - u[0]) / 0.05;
					mx = Math.max(mx, Math.abs(minU[m * NPAR + 1] - u[1]) / 1e-3);
					mx = Math.max(mx, Math.abs(minU[m * NPAR + 2] - u[2]) / 1e-3);
					mx = Math.max(mx, Math.abs(minU[m * NPAR + 3] - u[3]) / 0.5);
					mx = Math.max(mx, Math.abs(minU[m * NPAR + 4] - u[4]) / 1e-3);
					if (mx < 1.0 && minF[m] < 1e299) dup = true;
				}
				if (dup) continue;
			}
			double F = outer(u, nh, ze, ws);
			ntot += ws.lmIt;
			System.arraycopy(u, 0, minU, nmin * NPAR, NPAR);
			minF[nmin] = F; ws.minC[nmin] = ws.lmConv; ws.minK[nmin] = ws.startKind[j];
			nmin++;
		}
		if (nmin == 0) {
			out.status = 4;
			return;
		}
		int ib = 0;      // numpy argmin: the first NaN wins, otherwise the first smallest
		for (int m = 0; m < nmin; m++) {
			if (Double.isNaN(minF[m])) { ib = m; break; }
			if (minF[m] < minF[ib]) ib = m;
		}
		System.arraycopy(minU, ib * NPAR, u, 0, NPAR);
		double F = minF[ib];
		boolean conv = ws.minC[ib];
		out.startWon = ws.minK[ib];
		int flags = 0;
		// ---- near-wire flips (E5): theta0 frozen at the chosen state
		cost.predict(u, 1, npar, false, false, th0f, ws);
		if (nFlip > 0) {
			double[] hs = ws.flipHs, d0 = ws.flipD0, jd0 = ws.flipJd0, cand = ws.flipCand, jd = ws.flipJd, w = ws.flipW;
			cost.cost(u, 1, npar, true, ze, th0f, ws);
			loadFlipState(ws, nh, npar);
			for (int k = 0; k < nh; k++) {
				cand[k] = 1e300;
				if (ws.lkind[k] != 2 && finite(d0[k]) && Math.abs(d0[k]) < Math.max(0.2, 0.5 * ws.lm[k]) && rowFinite(jd0, k, npar))
					cand[k] = Math.abs(d0[k]);
			}
			int nf = 0;
			double[] u2 = ws.u2;
			while (true) {
				int kk = 0;      // stable order: the first smallest candidate
				for (int k = 1; k < nh; k++) if (cand[k] < cand[kk]) kk = k;
				if (nf >= nFlip || cand[kk] >= 1e299) break;
				cand[kk] = 1e300;
				nf++;
				for (int j = 0; j < npar; j++) jd[j] = jd0[kk * NPAR + j];
				boolean fin = true;
				for (int j = 0; j < npar && fin; j++) for (int m = 0; m < npar; m++) if (!finite(hs[j * NPAR + m])) { fin = false; break; }
				if (!fin) break;
				for (int j = 0; j < npar; j++) for (int m = 0; m < npar; m++) ws.lmA[j * NPAR + m] = hs[j * NPAR + m];
				System.arraycopy(jd, 0, w, 0, npar);
				if (!SmallLinAlg.solve(ws.lmA, NPAR, npar, w, ws.piv)) throw new ArithmeticException("singular matrix in flip step");
				double den = 0.0;
				for (int j = 0; j < npar; j++) den += jd[j] * w[j];
				if (!(den > 0)) continue;
				System.arraycopy(u, 0, u2, 0, NPAR);
				for (int j = 0; j < npar; j++) u2[j] += w[j] * (-2.0 * d0[kk]) / den;
				if (!LevMar.inb(u2)) continue;
				double F2 = outer(u2, nh, ze, ws);
				boolean c2 = ws.lmConv;
				ntot += ws.lmIt;
				boolean won = F2 < F - 1e-6;
				cost.predict(won ? u2 : u, 1, npar, false, false, th0f, ws);      // restore theta0 of the current best state
				if (won) {
					System.arraycopy(u2, 0, u, 0, NPAR);
					F = F2; conv = c2; flags |= 8; out.startWon = 3;
					cost.cost(u, 1, npar, true, ze, th0f, ws);
					loadFlipState(ws, nh, npar);
				}
				if (nmin < 8) {
					System.arraycopy(u2, 0, minU, nmin * NPAR, NPAR);
					minF[nmin] = F2; ws.minK[nmin] = 3; nmin++;
				}
			}
		}
		// ---- final evaluation, covariance, status
		cost.predict(u, 1, npar, false, false, th0f, ws);
		F = cost.cost(u, 1, npar, true, ze, th0f, ws);
		if (!(blockFinite(ws.he, npar) && blockFinite(ws.hc, npar) && finite(F))) {
			System.arraycopy(u, 0, out.u, 0, NPAR);
			out.cost = F;
			out.status = 4;
			return;
		}
		double[] ea = ws.eigA, ev_ = ws.eigV, ew = ws.eigW, iw = ws.eigIw;
		for (int j = 0; j < npar; j++) for (int m = 0; m < npar; m++) ea[j * NPAR + m] = ws.he[j * NPAR + m];
		SmallLinAlg.eigSym(ea, NPAR, npar, ew, ev_);
		double wmin = ew[0], wmax = ew[0];
		for (int j = 1; j < npar; j++) { wmin = Math.min(wmin, ew[j]); wmax = Math.max(wmax, ew[j]); }
		if (!(wmin > 1e-12 * (1e-12 > wmax ? 1e-12 : wmax))) {
			for (int j = 0; j < npar; j++) for (int m = 0; m < npar; m++) ea[j * NPAR + m] = ws.hc[j * NPAR + m];
			SmallLinAlg.eigSym(ea, NPAR, npar, ew, ev_);
			flags |= 4;
			wmax = ew[0];
			for (int j = 1; j < npar; j++) wmax = Math.max(wmax, ew[j]);
		}
		double cut = 1e-12 * (1e-12 > wmax ? 1e-12 : wmax);
		for (int j = 0; j < npar; j++) iw[j] = ew[j] > cut ? 1.0 / ew[j] : 0.0;
		double[] C = out.cov;
		for (int a = 0; a < npar; a++)
			for (int b = 0; b < npar; b++) {
				double s = 0.0;
				for (int j = 0; j < npar; j++) s += ev_[a * npar + j] * iw[j] * ev_[b * npar + j];
				C[a * NPAR + b] = 2.0 * s;
			}
		// distinct minima: |du| > 3 sigma in any state parameter
		double[] sg = ws.eigW;      // eigenvalues are no longer needed
		for (int j = 0; j < 5; j++) {
			double c = C[j * NPAR + j];
			sg[j] = Math.sqrt(1e-30 > c ? 1e-30 : c);
		}
		int nd = 1;
		double d2 = 1e300;
		for (int m = 0; m < nmin; m++) {
			if (minF[m] >= 1e299) continue;
			double mx = Math.abs(minU[m * NPAR] - u[0]) / sg[0];
			for (int j = 1; j < 5; j++) mx = Math.max(mx, Math.abs(minU[m * NPAR + j] - u[j]) / sg[j]);
			if (mx > 3.0) {
				nd++;
				d2 = Math.min(d2, minF[m] - F);
			}
		}
		out.nmin = nd;
		if (d2 < 1.0) flags |= 1;
		double p = Math.exp(u[2]), thd = u[4] * 180 / Math.PI;
		if (p < 0.12 || p > 1.0 || thd < 30 || thd > 150) flags |= 2;
		if (atofOn && ws.atof.fnone == 1) flags |= 32;      // ATOF_NONE
		int st = conv ? 0 : 1;
		boolean usedBad = false;
		for (int k = 0; k < nh; k++) if (ws.lkind[k] != 2 && (ws.PF[k] & 1) == 1) usedBad = true;
		int tfl = ws.tflag;
		if ((tfl & 1) == 1 || ((tfl & 4) == 4 && usedBad)) st = 3;
		boolean uf = true;
		for (int j = 0; j < npar; j++) if (!finite(u[j])) uf = false;
		if (!(finite(F) && uf)) st = 4;
		out.status = st; out.flags = flags; out.niter = ntot;
		System.arraycopy(u, 0, out.u, 0, NPAR);
		out.cost = F;
		out.dcost2 = d2 < 1e299 ? d2 : Double.NaN;
		for (int j = 0; j < NPL; j++) out.theta0[j] = th0f[j];
		int nhh = Math.min(nh, out.hitD.length);
		for (int k = 0; k < nhh; k++) {
			out.hitD[k] = ws.outH[3 * k + 1]; out.hitS[k] = ws.outH[3 * k + 2]; out.hitL[k] = ws.outH[3 * k];
			out.hitKind[k] = ws.lkind[k]; out.hitRow[k] = ws.lkind[k] != 2 ? ws.lrow[k] : -1; out.hitWire[k] = ws.lg[k];
		}
		out.nh = nhh;
	}

	private static boolean finite(double x) { return !(Double.isNaN(x) || Double.isInfinite(x)); }

	private static boolean blockFinite(double[] a, int n) {
		for (int j = 0; j < n; j++) for (int m = 0; m < n; m++) if (!finite(a[j * NPAR + m])) return false;
		return true;
	}

	private static boolean rowFinite(double[] a, int k, int n) {
		for (int j = 0; j < n; j++) if (!finite(a[k * NPAR + j])) return false;
		return true;
	}

	/** Hs = clipped Gauss-Newton matrix + 1e-9 on the diagonal, D0 and JD0 of the last cost(want) (fit.py:683-690). */
	private static void loadFlipState(GFitWorkspace ws, int nh, int npar) {
		for (int j = 0; j < npar; j++) for (int m = 0; m < npar; m++) ws.flipHs[j * NPAR + m] = ws.hc[j * NPAR + m];
		for (int j = 0; j < npar; j++) ws.flipHs[j * NPAR + j] += 1e-9;
		System.arraycopy(ws.D, 0, ws.flipD0, 0, nh);
		System.arraycopy(ws.JD, 0, ws.flipJd0, 0, nh * NPAR);
	}

	// ================================================================== outer passes (fit._outer)

	/** theta0 frozen at the current state, LM, refresh; stops when theta0 changes by less than 1 % and the cost by less than tol_cost. ws.lmConv, ws.lmIt = total iterations. */
	private double outer(double[] u, int nh, double ze, GFitWorkspace ws) {
		double Fp = 1e300, F = 1e300;
		boolean conv = false;
		int ntot = 0;
		double[] th0f = ws.th0f;
		for (int ps = 0; ps < npass; ps++) {
			cost.predict(u, 1, NPAR, false, false, ws.th0, ws);
			double ch = 0.0;
			for (int j = 0; j < NPL; j++) {
				if (ws.th0[j] > 0) ch = Math.max(ch, Math.abs(ws.th0[j] - th0f[j]) / ws.th0[j]);
				th0f[j] = ws.th0[j];
			}
			F = lm.run(u, 1, NPAR, ze, th0f, maxit, ws);
			conv = ws.lmConv;
			ntot += ws.lmIt;
			if (ps > 0 && ch < 0.01 && Math.abs(Fp - F) < tolCost) break;
			Fp = F;
		}
		ws.lmConv = conv;
		ws.lmIt = ntot;
		return F;
	}

	// ================================================================== start grid (fit._grid)

	/**
	 * Coarse grid of 8 momenta x 13 polar angles x 5 azimuth offsets around the seed phi with the helix model and the dbox
	 * term, at z = ze; keeps the lowest costs whose theta differs by at least 10 degrees. Equal costs are not rare (every hit
	 * saturated in the dbox gives a plateau of exactly equal sums), so the order is numpy's default argsort, reproduced by
	 * {@link #argsortNumpy}. Fills ws.bestU (5 per start), ws.bestF and returns the number of starts.
	 */
	private int grid(double ze, double b0phi, GFitWorkspace ws) {
		final int npv = 8, nth = 13;
		final double th0 = 30.0, dth = 10.0;
		final double[] y = ws.gridY, gf = ws.gridF, gu = ws.gridU, o = ws.o9, b3 = ws.box3;
		final int nh = ws.nh;
		int q = 0;
		for (int ip = 0; ip < npv; ip++) {
			double lp = Math.log(0.12) + ip * (Math.log(1.0) - Math.log(0.12)) / (npv - 1);
			for (int it = 0; it < nth; it++) {
				for (int jf = 0; jf < 5; jf++) {
					y[0] = 0.0; y[1] = b0phi + (jf - 2) * Math.PI / 180.0; y[2] = lp; y[3] = ze; y[4] = (th0 + dth * it) * Math.PI / 180.0;
					double F = 0.0;
					for (int k = 0; k < nh; k++) {
						if (ws.lkind[k] != 0) continue;
						double sd = HelixModel.hhit(y, geo.wo, geo.wu, ws.lg[k], bzHelix, beamX, beamY, o, -1.0, false, ws.jx5, ws.j53, ws.sArc1);
						HitTerm.dbox(sd, ws.lm[k], ws.lsig[k], ws.lh[k], false, cfg, b3);
						F += b3[0];
					}
					if (zOn) {
						double r = (y[3] - ze) / sigZ;
						F += r * r;
					}
					gf[q] = F;
					for (int c = 0; c < 5; c++) gu[q * 5 + c] = y[c];
					q++;
				}
			}
		}
		final int ng = q;
		int[] ord = ws.gridOrd;
		for (int i = 0; i < ng; i++) ord[i] = i;
		argsortNumpy(gf, ng, ord, ws.gridTmp);
		final double sep = 10.0 * Math.PI / 180.0 - 1e-9;
		int nsel = 0;
		for (int a = 0; a < ng; a++) {
			if (nsel >= GFitWorkspace.NBEST) break;
			int i = ord[a];
			boolean ok = true;
			for (int j = 0; j < nsel; j++) {
				if (Math.abs(gu[i * 5 + 4] - ws.bestU[j * 5 + 4]) < sep) { ok = false; break; }
			}
			if (ok) {
				for (int c = 0; c < 5; c++) ws.bestU[nsel * 5 + c] = gu[i * 5 + c];
				ws.bestF[nsel] = gf[i];
				nsel++;
			}
		}
		return nsel;
	}

	/** numpy sort order: NaN after every number. */
	private static boolean less(double a, double b) { return a < b || (Double.isNaN(b) && !Double.isNaN(a)); }

	private static void swap(int[] t, int i, int j) { int x = t[i]; t[i] = t[j]; t[j] = x; }

	/**
	 * numpy's default {@code argsort} of v[0..n) (npy_aquicksort of numpy 1.23, the pinned reference: introsort, median-of-3
	 * quicksort while more than 16 elements remain, insertion sort below, heapsort when the depth limit is exhausted) into t, which holds the
	 * identity permutation on entry; equal keys end up in numpy's order. stk is scratch of at least 192 entries.
	 */
	static void argsortNumpy(double[] v, int n, int[] t, int[] stk) {
		int pl = 0, pr = n - 1, sp = 0, dp = 0;
		int cdepth = 2 * (31 - Integer.numberOfLeadingZeros(n));
		for (;;) {
			if (cdepth < 0) {
				aheapsort(v, t, pl, pr - pl + 1);
			} else {
				while (pr - pl > 15) {
					int pm = pl + ((pr - pl) >> 1);
					if (less(v[t[pm]], v[t[pl]])) swap(t, pm, pl);
					if (less(v[t[pr]], v[t[pm]])) swap(t, pr, pm);
					if (less(v[t[pm]], v[t[pl]])) swap(t, pm, pl);
					double vp = v[t[pm]];
					int pi = pl, pj = pr - 1;
					swap(t, pm, pj);
					for (;;) {
						do { ++pi; } while (less(v[t[pi]], vp));
						do { --pj; } while (less(vp, v[t[pj]]));
						if (pi >= pj) break;
						swap(t, pi, pj);
					}
					swap(t, pi, pr - 1);
					if (pi - pl < pr - pi) {     // push the larger partition
						stk[sp++] = pi + 1; stk[sp++] = pr;
						pr = pi - 1;
					} else {
						stk[sp++] = pl; stk[sp++] = pi - 1;
						pl = pi + 1;
					}
					stk[128 + dp++] = --cdepth;
				}
				for (int pi = pl + 1; pi <= pr; ++pi) {
					int vi = t[pi];
					double vp = v[vi];
					int pj = pi;
					while (pj > pl && less(vp, v[t[pj - 1]])) { t[pj] = t[pj - 1]; pj--; }
					t[pj] = vi;
				}
			}
			if (sp == 0) break;
			pr = stk[--sp];
			pl = stk[--sp];
			cdepth = stk[128 + --dp];
		}
	}

	/** numpy's aheapsort of t[off..off+n) (1-based heap indexing as numpy). */
	private static void aheapsort(double[] v, int[] t, int off, int n) {
		final int b = off - 1;
		for (int l = n >> 1; l > 0; --l) {
			int tmp = t[b + l], i = l, j = l << 1;
			while (j <= n) {
				if (j < n && less(v[t[b + j]], v[t[b + j + 1]])) j += 1;
				if (less(v[tmp], v[t[b + j]])) { t[b + i] = t[b + j]; i = j; j += j; } else break;
			}
			t[b + i] = tmp;
		}
		while (n > 1) {
			int tmp = t[b + n];
			t[b + n] = t[b + 1];
			n -= 1;
			int i = 1, j = 2;
			while (j <= n) {
				if (j < n && less(v[t[b + j]], v[t[b + j + 1]])) j++;
				if (less(v[tmp], v[t[b + j]])) { t[b + i] = t[b + j]; i = j; j += j; } else break;
			}
			t[b + i] = tmp;
		}
	}
}
