package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * HitRecovery on the 60-track subset t60 (full search, final ATOF times, per-half 23014 tables) against the pinned Python
 * recovery (recover.py): crossings and probe table at the stored baseline state, accepted rows and the result after recovery
 * (fidelity tiers 1-3: status and hit set, cost and parameters, cost bound), and the drop-worst retry.
 */
public class GFitRecoveryTest {

	/** relative tolerance of the probe values and the floor of the reference magnitude (l, d, s, res are in mm or log units of order 1). */
	static final double REL = 1e-9, FLOOR = 1e-3;

	static GFitReplaySet set;
	static GFitConfig cfg;
	static GlobalTrackFitter ft;
	static HitRecovery hr;
	static int maxHits;
	static Run t60;

	/** results of one pass over a set: base fit, recovery, accepted rows and retry counters per track. */
	static final class Run {
		GFitResult[] base, after;
		int[][] accEv, accLayer, accRound;
		int[] nAdded, retries, dropEv0;
		double recoverMs;
	}

	@BeforeClass
	public static void load() {
		GFitFullSearchTest.load();
		set = GFitFullSearchTest.set;
		maxHits = GFitFullSearchTest.maxHits;
		cfg = GFitConfig.pythonDefault();
		ft = GFitFullSearchTest.fitter(cfg);
		hr = new HitRecovery(ft, new TrackCost(GFitFullSearchTest.prop, GFitFullSearchTest.geo, cfg), GFitFullSearchTest.geo, cfg);
		if (t60 == null) t60 = run(set, new GFitWorkspace(maxHits), maxHits);
	}

	/** base fit then recovery of every track of s; recovery time measured alone (refits included). */
	static Run run(GFitReplaySet s, GFitWorkspace ws, int maxH) {
		Run r = new Run();
		int n = s.n;
		r.base = new GFitResult[n]; r.after = new GFitResult[n];
		r.accEv = new int[n][]; r.accLayer = new int[n][]; r.accRound = new int[n][];
		r.nAdded = new int[n]; r.retries = new int[n]; r.dropEv0 = new int[n];
		long ns = 0;
		for (int k = 0; k < n; k++) {
			r.base[k] = new GFitResult(maxH); r.after[k] = new GFitResult(maxH);
			AtofPositionTable tab = GFitFullSearchTest.consts.atofPosition(GFitFixtureConstants.RUN_PER_HALF, s.events[k].half);
			ft.fit(s.tracks[k], s.events[k], tab, false, ws, r.base[k]);
			long t0 = System.nanoTime();
			r.nAdded[k] = hr.recover(s.tracks[k], s.events[k], tab, false, r.base[k], ws, r.after[k]);
			ns += System.nanoTime() - t0;
			HitRecovery.Scratch sc = ws.rec;
			r.accEv[k] = java.util.Arrays.copyOf(sc.accEv, sc.nAcc); r.accLayer[k] = java.util.Arrays.copyOf(sc.accLayer, sc.nAcc);
			r.accRound[k] = java.util.Arrays.copyOf(sc.accRound, sc.nAcc);
			r.retries[k] = sc.retries; r.dropEv0[k] = sc.dropEv[0];
		}
		r.recoverMs = ns * 1e-6 / n;
		return r;
	}

	private static boolean near(double a, double b) {
		if (Double.isNaN(a) || Double.isNaN(b)) return Double.isNaN(a) && Double.isNaN(b);
		if (a == b) return true;
		return Math.abs(a - b) <= REL * Math.max(Math.abs(b), FLOOR);
	}

	@Test
	public void probeMatchesPython() {
		GFitFixtures.Fx x = set.recover;
		boolean[] xOk = bool(x, "x_ok");
		int[] gstar = x.i("x_gstar"), ptT = x.i("pt_t"), ptLi = x.i("pt_li"), ptEv = x.i("pt_ev"), ptG = x.i("pt_g"), ptDw = x.i("pt_dw"), ptIsm = x.i("pt_ism");
		boolean[] crossed = bool(x, "x_crossed"), has = bool(x, "x_has"), missing = bool(x, "x_missing"), pass = bool(x, "pt_pass"), sel = bool(x, "pt_sel");
		double[] time = x.d("pt_time"), doca = x.d("pt_doca"), sig = x.d("pt_sig"), l = x.d("pt_l"), dpred = x.d("pt_dpred"), sp = x.d("pt_s"), res = x.d("pt_res");
		GFitWorkspace ws = new GFitWorkspace(maxHits);
		HitRecovery.Scratch s = ws.rec;
		int row = 0, nrows = 0, bad = 0;
		double worst = 0;
		List<String> msg = new ArrayList<>();
		for (int k = 0; k < set.n; k++) {
			GFitResult b = set.baselineResult(k, maxHits);
			hr.prepare(set.tracks[k], set.events[k], b, ws, s);
			assertEquals("x_ok " + k, xOk[k], s.ok);
			for (int li = 0; li < 8; li++) {
				assertEquals("gstar " + k + "/" + li, gstar[k * 8 + li], s.gstar[li]);
				assertEquals("crossed " + k + "/" + li, crossed[k * 8 + li], s.crossed[li]);
				assertEquals("has " + k + "/" + li, has[k * 8 + li], s.has[li]);
				assertEquals("missing " + k + "/" + li, missing[k * 8 + li], s.missing[li]);
			}
			int n0 = row;
			while (row < ptT.length && ptT[row] == k) row++;
			assertEquals("pool rows of track " + k, row - n0, s.nc);
			GFitEvent ev = set.events[k];
			for (int c = 0; c < s.nc; c++) {
				int i = n0 + c, p = s.cIdx[c];
				assertEquals(ptLi[i], s.cLi[c]);
				assertEquals(ptEv[i], ev.pEv[p]);
				assertEquals(ptG[i], s.cG[c]);
				assertEquals(ptDw[i], s.cDw[c]);
				assertEquals(ptIsm[i], s.cIsm[c]);
				assertEquals(time[i], ev.pTime[p], 0.0);
				assertEquals(doca[i], s.cDoca[c], 0.0);
				assertTrue("sigma " + sig[i] + " " + s.cSig[c], near(sig[i], s.cSig[c]));
				double[] jv = {s.cL[c], s.cD[c], s.cS[c], s.cDoca[c] - Math.abs(s.cD[c])}, pv = {l[i], dpred[i], sp[i], res[i]};
				for (int q = 0; q < 4; q++) {
					worst = Math.max(worst, Math.abs(jv[q] - pv[q]) / Math.max(Math.abs(pv[q]), FLOOR));
					if (!near(jv[q], pv[q])) { bad++; msg.add("track " + k + " row " + c + " value " + q + " java " + jv[q] + " python " + pv[q]); }
				}
				assertEquals("pass " + k + "/" + c, pass[i], s.cPass[c]);
				boolean isSel = false;
				for (int q = 0; q < s.nsel; q++) if (s.sel[q] == c) isSel = true;
				assertEquals("sel " + k + "/" + c, sel[i], isSel);
			}
			nrows += s.nc;
		}
		assertEquals(ptT.length, nrows);
		System.out.printf("probe: %d tracks, %d pool rows, worst relative deviation %.2e (floor %.0e), %d values outside %.0e%n", set.n, nrows, worst, FLOOR, bad, REL);
		assertTrue(msg.toString(), bad == 0);
	}

	private static boolean[] bool(GFitFixtures.Fx x, String name) {
		long[] a = x.l(name);
		boolean[] b = new boolean[a.length];
		for (int i = 0; i < a.length; i++) b[i] = a[i] != 0;
		return b;
	}

	/** tracks whose accepted rows (event row, layer, round, in order) differ from acc_*; prints them. */
	static int accMismatches(String label, GFitReplaySet s, Run r) {
		GFitFixtures.Fx a = s.after;
		int[] at = a.i("acc_t"), aev = a.i("acc_ev"), arnd = a.i("acc_rnd"), ali = a.i("acc_li");
		int bad = 0;
		for (int k = 0; k < s.n; k++) {
			int[] py = accRows(at, k);       // acc_* lists the rows round by round, so the rows of one track are not contiguous
			boolean eq = py.length == r.accEv[k].length;
			for (int j = 0; eq && j < py.length; j++) eq = aev[py[j]] == r.accEv[k][j] && arnd[py[j]] == r.accRound[k][j] && ali[py[j]] == r.accLayer[k][j];
			if (!eq) {
				bad++;
				int[] pe = new int[py.length];
				for (int j = 0; j < py.length; j++) pe[j] = aev[py[j]];
				System.out.println("[" + label + "] accepted rows differ: trackid " + s.tracks[k].trackid + " row " + k + " java " + r.accEv[k].length
					+ " rows " + java.util.Arrays.toString(r.accEv[k]) + " python " + py.length + " rows " + java.util.Arrays.toString(pe));
			}
		}
		return bad;
	}

	/** indices (file order) of the acc_* rows of track k. */
	private static int[] accRows(int[] at, int k) {
		int n = 0;
		for (int v : at) if (v == k) n++;
		int[] o = new int[n];
		n = 0;
		for (int i = 0; i < at.length; i++) if (at[i] == k) o[n++] = i;
		return o;
	}

	/**
	 * Tiers 1-3 of D5 for the result after recovery against expected_after: tier 1 = status and final hit set (original hits plus the
	 * accepted rows, compared through the accepted rows) identical; tier 2 and 3 as GFitFullSearchTest.checkTiers. Prints the worst
	 * deviations and returns {ok1, n1, ok2, n2, ok3, n3}.
	 */
	static int[] tiersAfter(String label, GFitReplaySet s, Run r) {
		GFitFixtures.Fx ex = s.after;
		int[] st = ex.i("ea_status"), nmin = ex.i("ea_nmin"), fl = ex.i("ea_flags"), nrec = ex.i("n_rec_added");
		double[] cost = ex.d("ea_cost"), state = ex.d("ea_state"), cov = ex.d("ea_cov");
		int n1 = 0, ok1 = 0, n2 = 0, ok2 = 0, n3 = 0, ok3 = 0;
		double wCost = 0, wPar = 0, wCov = 0, wCost3 = Double.NEGATIVE_INFINITY;
		double[] cx = new double[25];
		int[] at = ex.i("acc_t"), aev = ex.i("acc_ev");
		for (int k = 0; k < s.n; k++) {
			GFitResult j = r.after[k];
			String id = "trackid " + s.tracks[k].trackid + " row " + k + " tier " + s.tier[k];
			int[] je = r.accEv[k].clone();
			java.util.Arrays.sort(je);
			int[] py = accRows(at, k), pe = new int[py.length];
			for (int q = 0; q < py.length; q++) pe[q] = aev[py[q]];
			java.util.Arrays.sort(pe);
			boolean hits = java.util.Arrays.equals(pe, je) && nrec[k] == r.nAdded[k];
			n1++;
			if (j.status == st[k] && hits) ok1++;
			else System.out.println("[" + label + "] outside tier 1: " + id + " status java " + j.status + " python " + st[k] + " hit set " + hits);
			boolean t3 = nmin[k] > 1 || (fl[k] & 1) != 0 || (fl[k] & 8) != 0;
			boolean t2 = st[k] == 0 && nmin[k] == 1 && (fl[k] & 8) == 0 && !t3;
			if (t2) {
				n2++;
				double dc = Math.abs(j.cost - cost[k]);
				wCost = Math.max(wCost, dc);
				boolean ok = j.status == 0 && dc <= 1e-3;
				double[] xx = {j.u[0], j.u[1], j.u[3], j.u[4], 1.0 / j.p()};
				j.covX(cx);
				double dp = 0, dcv = 0;
				for (int q = 0; q < 5; q++) {
					double sg = Math.sqrt(cov[k * 25 + q * 6]);
					dp = Math.max(dp, Math.abs(xx[q] - state[k * 5 + q]) / sg);
					dcv = Math.max(dcv, Math.abs(cx[q * 6] - cov[k * 25 + q * 6]) / cov[k * 25 + q * 6]);
				}
				wPar = Math.max(wPar, dp);
				wCov = Math.max(wCov, dcv);
				ok &= dp <= 0.01 && dcv <= 0.01;
				if (ok) ok2++;
				else System.out.println("[" + label + "] outside tier 2: " + id + " cost dev " + dc + " param " + dp + " sigma, cov diag " + dcv);
			}
			if (t3) {
				n3++;
				double d = j.cost - cost[k];
				wCost3 = Math.max(wCost3, d);
				if (d <= 1.0) ok3++;
				else System.out.println("[" + label + "] outside tier 3: " + id + " java cost " + j.cost + " python " + cost[k]);
			}
		}
		System.out.printf("[%s] tier1 %d/%d (status, accepted rows); tier2 %d/%d: worst cost %.2e, params %.2e sigma, cov diag %.2e; tier3 %d/%d: worst Java-Python cost %.3g%n",
			label, ok1, n1, ok2, n2, wCost, wPar, wCov, ok3, n3, wCost3);
		return new int[] {ok1, n1, ok2, n2, ok3, n3};
	}

	@Test
	public void recoverMatchesPython() {
		int bad = accMismatches("t60", set, t60);
		int total = 0, tracks = 0;
		for (int k = 0; k < set.n; k++) { total += t60.nAdded[k]; if (t60.nAdded[k] > 0) tracks++; }
		System.out.printf("recover: %d rows added on %d of %d tracks; %d tracks with accepted rows unlike Python; %.1f ms per track (recovery alone)%n",
			total, tracks, set.n, bad, t60.recoverMs);
		int[] t = tiersAfter("t60", set, t60);
		assertTrue("accepted rows differ on " + bad + " tracks", bad <= 0.005 * set.n);
		assertTrue("tier 1", t[0] >= 0.995 * t[1]);
		assertTrue("tier 2", t[2] >= 0.98 * t[3]);
		assertTrue("tier 3", t[4] >= 0.90 * t[5]);
	}

	/** a track whose first refit was rejected with two or more added rows: one retry without the added row of lowest l (the refit's own hit l). */
	@Test
	public void dropWorst() {
		GFitReplaySet s = set;
		Run r = t60;
		int k = firstRetry(r, s.n);
		GFitWorkspace ws = new GFitWorkspace(maxHits);
		assertTrue("t60 has a track whose first refit was rejected with two or more added rows", k >= 0);
		System.out.println("dropWorst: track " + k + " of t60 (trackid " + s.tracks[k].trackid + "), retries " + r.retries[k] + ", dropped event row " + r.dropEv0[k]);
		// the refit with all selected rows, its own hit l of the added rows: the lowest is the row dropped
		HitRecovery.Scratch sc = ws.rec;
		GFitTrack cur = s.tracks[k];
		GFitEvent ev = s.events[k];
		int ns = hr.prepare(cur, ev, r.base[k], ws, sc);
		assertTrue(ns >= 2);
		GFitTrack t1 = hr.extend(cur, ev, sc, sc.sel, ns);
		int[] selEv = new int[ns];
		for (int i = 0; i < ns; i++) selEv[i] = ev.pEv[sc.cIdx[sc.sel[i]]];
		GFitResult f = new GFitResult(maxHits + 8);
		AtofPositionTable tab = GFitFullSearchTest.consts.atofPosition(GFitFixtureConstants.RUN_PER_HALF, ev.half);
		ft.fit(t1, ev, tab, false, ws, f);
		double q = cfg.recQ();
		assertTrue("first refit must be rejected", f.status != 0 || f.cost - r.base[k].cost > q * ns);
		int worst = -1;
		double lw = 0;
		for (int i = 0; i < ns; i++) {
			double l = Double.NaN;
			for (int h = 0; h < f.nh; h++) if (f.hitRow[h] == selEv[i]) { l = f.hitL[h]; break; }
			if (!Double.isFinite(l)) l = -1e300;
			if (worst < 0 || l < lw) { lw = l; worst = i; }
		}
		assertEquals(selEv[worst], r.dropEv0[k]);
		assertEquals(1, r.retries[k]);
	}

	private static int firstRetry(Run r, int n) {
		for (int k = 0; k < n; k++) if (r.retries[k] > 0 && r.dropEv0[k] >= 0) return k;
		return -1;
	}
}
