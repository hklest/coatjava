package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * GlobalTrackFitter in FULL mode on the 60-track subset t60 (final ATOF times, per-half 23014 tables) against the pinned
 * Python fit: tiers 1-3 of parent plan D5, a 2-hit track and a NaN doca. Also hosts the shared setup of the fit tests.
 */
public class GFitFullSearchTest {

	static GFitReplaySet set;
	static Propagator prop;
	static CellGeometry geo;
	static GFitFixtureConstants consts;
	static int maxHits;

	@BeforeClass
	public static void load() {
		if (set != null) return;
		GFitFixtures.Fx f = GFitFixtures.resource("unit/field.txt.gz");
		double[] fp = f.d("field_par");
		FieldGrid field = FieldGrid.fromArrays(fp[0], fp[1], fp[2], fp[3], f.shape("br")[0], f.shape("br")[1], f.d("br"), f.d("bz"));
		double[] beam = GFitFixtures.resource("unit/stopping.txt.gz").d("beam");
		prop = new Propagator(field, MaterialStack.dataBest(), StoppingTables.proton(), beam[0], beam[1]);
		GFitFixtures.Fx w = GFitFixtures.resource("unit/wires.txt.gz");
		geo = CellGeometry.fromWireEnds(w.d("WO"), w.d("WE"), w.i("layoff"));
		consts = new GFitFixtureConstants(false);
		set = GFitReplaySet.resource("t60");
		maxHits = 0;
		for (GFitTrack t : set.tracks) maxHits = Math.max(maxHits, t.nHits());
	}

	static GlobalTrackFitter fitter(GFitConfig cfg) {
		return new GlobalTrackFitter(new TrackCost(prop, geo, cfg), prop, geo, cfg);
	}

	static GFitWorkspace workspace() { return new GFitWorkspace(maxHits); }

	/** fits every track of the set once (final times, 23014 tables) into fresh results. */
	static GFitResult[] runAll(GlobalTrackFitter ft, GFitWorkspace ws) {
		GFitResult[] r = new GFitResult[set.n];
		for (int k = 0; k < set.n; k++) {
			r[k] = new GFitResult(maxHits);
			ft.fit(set.tracks[k], set.events[k], consts.atofPosition(GFitFixtureConstants.RUN_PER_HALF, set.events[k].half), false, ws, r[k]);
		}
		return r;
	}

	static void fitOne(GlobalTrackFitter ft, GFitWorkspace ws, GFitTrack t, GFitEvent e, GFitResult out) {
		ft.fit(t, e, consts.atofPosition(GFitFixtureConstants.RUN_PER_HALF, e.half), false, ws, out);
	}

	/** every field of two results equal bit for bit (NaN equal to NaN). */
	static String diff(GFitResult a, GFitResult b) {
		if (a.status != b.status || a.flags != b.flags || a.niter != b.niter || a.nstart != b.nstart || a.nmin != b.nmin
			|| a.startWon != b.startWon || a.nh != b.nh) return "counters " + a.status + "/" + b.status + " " + a.flags + "/" + b.flags;
		if (!same(a.u, b.u) || !same(a.cov, b.cov) || !same(a.theta0, b.theta0)) return "u, cov or theta0";
		if (Double.doubleToLongBits(a.cost) != Double.doubleToLongBits(b.cost) || Double.doubleToLongBits(a.dcost2) != Double.doubleToLongBits(b.dcost2))
			return "cost " + a.cost + "/" + b.cost;
		for (int k = 0; k < a.nh; k++) {
			if (Double.doubleToLongBits(a.hitD[k]) != Double.doubleToLongBits(b.hitD[k]) || Double.doubleToLongBits(a.hitS[k]) != Double.doubleToLongBits(b.hitS[k])
				|| Double.doubleToLongBits(a.hitL[k]) != Double.doubleToLongBits(b.hitL[k]) || a.hitKind[k] != b.hitKind[k]
				|| a.hitRow[k] != b.hitRow[k] || a.hitWire[k] != b.hitWire[k]) return "hit " + k;
		}
		return null;
	}

	private static boolean same(double[] a, double[] b) {
		for (int i = 0; i < a.length; i++) if (Double.doubleToLongBits(a[i]) != Double.doubleToLongBits(b[i])) return false;
		return true;
	}

	private static double relDev(double a, double e) { return Math.abs(a - e) / Math.max(1.0, Math.abs(e)); }

	/**
	 * Parent D5 tiers 1-3 of results r against the stored Python fit of the set; prints worst deviations and the ids (trackid, row)
	 * of tracks outside a tolerance; asserts the thresholds (tier 1 99.5 %, tier 2 98 %, tier 3 90 %).
	 */
	static void checkTiers(String label, GFitResult[] r) {
		GFitFixtures.Fx ex = set.before;
		int[] st = ex.i("est_status"), nmin = ex.i("est_nmin"), fl = ex.i("est_flags"), eoff = ex.i("eh_off"), ewire = ex.i("eh_wire"), ekind = ex.i("eh_kind");
		double[] cost = ex.d("est_cost"), state = ex.d("est_state"), cov = ex.d("est_cov");
		List<String> out1 = new ArrayList<>(), out2 = new ArrayList<>(), out3 = new ArrayList<>();
		int n1 = 0, ok1 = 0, n2 = 0, ok2 = 0, n3 = 0, ok3 = 0;
		double wCost = 0, wPar = 0, wCov = 0, wCost3 = Double.NEGATIVE_INFINITY;
		double[] cx = new double[25];
		for (int k = 0; k < set.n; k++) {
			GFitResult j = r[k];
			int tier = set.tier[k];
			String id = "trackid " + set.tracks[k].trackid + " row " + k + " tier " + tier;
			// tier 1 (all tracks): status and hit set
			n1++;
			boolean hits = j.nh == eoff[k + 1] - eoff[k];
			for (int h = 0; hits && h < j.nh; h++) hits = j.hitWire[h] == ewire[eoff[k] + h] && j.hitKind[h] == ekind[eoff[k] + h];
			if (j.status == st[k] && hits) ok1++;
			else out1.add(id + " status java " + j.status + " python " + st[k] + " hits " + hits);
			boolean t3 = nmin[k] > 1 || (fl[k] & 1) != 0 || (fl[k] & 8) != 0;
			boolean t2 = st[k] == 0 && nmin[k] == 1 && (fl[k] & 8) == 0 && !t3;
			if (t2) {
				n2++;
				boolean ok = j.status == 0;
				double dc = relDev(j.cost, cost[k]);
				wCost = Math.max(wCost, dc);
				ok &= dc <= 1e-3;
				double[] x = {j.u[0], j.u[1], j.u[3], j.u[4], 1.0 / j.p()};
				j.covX(cx);
				double dp = 0, dcv = 0;
				for (int q = 0; q < 5; q++) {
					double sg = Math.sqrt(cov[k * 25 + q * 6]);
					dp = Math.max(dp, Math.abs(x[q] - state[k * 5 + q]) / sg);
					dcv = Math.max(dcv, Math.abs(cx[q * 6] - cov[k * 25 + q * 6]) / cov[k * 25 + q * 6]);
				}
				wPar = Math.max(wPar, dp);
				wCov = Math.max(wCov, dcv);
				ok &= dp <= 0.01 && dcv <= 0.01;
				if (ok) ok2++;
				else out2.add(id + " cost dev " + dc + " (java " + j.cost + " python " + cost[k] + ") param " + dp + " sigma, cov diag " + dcv);
			}
			if (t3) {
				n3++;
				double d = j.cost - cost[k];
				wCost3 = Math.max(wCost3, d);
				if (d <= 1.0) ok3++;
				else out3.add(id + " java cost " + j.cost + " python " + cost[k]);
			}
		}
		System.out.printf("[%s] tier1 %d/%d (status and hit set); tier2 %d/%d: worst cost %.2e, params %.2e sigma, cov diag %.2e; tier3 %d/%d: worst Java-Python cost %.3g%n",
			label, ok1, n1, ok2, n2, wCost, wPar, wCov, ok3, n3, wCost3);
		for (String s : out1) System.out.println("[" + label + "] outside tier 1: " + s);
		for (String s : out2) System.out.println("[" + label + "] outside tier 2: " + s);
		for (String s : out3) System.out.println("[" + label + "] outside tier 3: " + s);
		assertTrue("tier 1 " + out1, ok1 >= 0.995 * n1);
		assertTrue("tier 2 " + out2, ok2 >= 0.98 * n2);
		assertTrue("tier 3 " + out3, ok3 >= 0.90 * n3);
	}

	@Test
	public void tiersOnSubset() {
		GlobalTrackFitter ft = fitter(GFitConfig.pythonDefault());
		GFitWorkspace ws = workspace();
		GFitResult[] r = runAll(ft, ws);
		checkTiers("full", r);
		if (Boolean.getBoolean("gfit.timing")) {
			runAll(ft, ws);       // warm-up (the first pass above counts as the first warm-up)
			GFitResult tmp = new GFitResult(maxHits);
			long t0 = System.nanoTime();
			for (int k = 0; k < set.n; k++) fitOne(ft, ws, set.tracks[k], set.events[k], tmp);
			double ms = (System.nanoTime() - t0) * 1e-6 / set.n;
			System.out.printf("timing: FULL single thread %.2f ms per track over %d tracks after warm-up%n", ms, set.n);
		}
	}

	@Test
	public void fewHitsAndBadInput() {
		GlobalTrackFitter ft = fitter(GFitConfig.pythonDefault());
		GFitWorkspace ws = workspace();
		GFitTrack t = set.tracks[0];
		GFitEvent e = set.events[0];
		GFitResult out = new GFitResult(maxHits);
		GFitTrack two = new GFitTrack(t.trackid, Arrays.copyOf(t.hg, 2), Arrays.copyOf(t.hm, 2), Arrays.copyOf(t.hsig, 2), Arrays.copyOf(t.hh, 2),
			Arrays.copyOf(t.hism, 2), Arrays.copyOf(t.hcls, 2), Arrays.copyOf(t.hrow, 2), t.seedPhi, t.seedP, t.seedTheta, t.seedZ);
		fitOne(ft, ws, two, e, out);
		assertEquals(2, out.status);
		// fit._dbox adds the clip (2.0) for a NaN doca (Lc < clip is false), so the faithful port fits such a hit as an outlier: status 0, no exception
		double[] hm = t.hm.clone();
		hm[1] = Double.NaN;
		GFitTrack nanDoca = new GFitTrack(t.trackid, t.hg, hm, t.hsig, t.hh, t.hism, t.hcls, t.hrow, t.seedPhi, t.seedP, t.seedTheta, t.seedZ);
		fitOne(ft, ws, nanDoca, e, out);
		assertEquals(0, out.status);
		assertTrue(Double.isFinite(out.cost));
		// a NaN half pitch makes the cost NaN everywhere: status 4 without an exception
		double[] hh = t.hh.clone();
		hh[1] = Double.NaN;
		GFitTrack bad = new GFitTrack(t.trackid, t.hg, t.hm, t.hsig, hh, t.hism, t.hcls, t.hrow, t.seedPhi, t.seedP, t.seedTheta, t.seedZ);
		fitOne(ft, ws, bad, e, out);
		assertEquals(4, out.status);
		// the workspace is reusable after the failure: the clean track reproduces its normal result
		GFitResult a = new GFitResult(maxHits), b = new GFitResult(maxHits);
		fitOne(ft, ws, t, e, a);
		fitOne(ft, new GFitWorkspace(maxHits), t, e, b);
		assertEquals(null, diff(a, b));
	}
}
