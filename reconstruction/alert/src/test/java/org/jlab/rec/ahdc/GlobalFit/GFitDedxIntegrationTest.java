package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * The switched sequence {@link SwitchedFit} on the t60 tracks: the dE/dx step on the fit's own default-fit state with
 * Python's constants against switch/t60/expected_ge10_dedx (dE/dx alone, no recovery) and expected_ge10_all (time fit, recovery and
 * dE/dx), and the switches-off sequence against the plain fit and recovery.
 */
public class GFitDedxIntegrationTest {

	static final int RUN = GFitFixtureConstants.RUN_PER_HALF;

	private static GFitReplaySet set;
	private static GFitFixtureConstants consts;
	private static AtofPositionTable[] timeTab;
	private static double[] aEnergy, aTime;

	@BeforeClass
	public static void load() {
		GFitFullSearchTest.load();
		consts = GFitFullSearchTest.consts;
		set = GFitFullSearchTest.set;
		timeTab = new AtofPositionTable[2];
		for (int h = 0; h < 2; h++) timeTab[h] = consts.atofPosition(RUN, h).withTime(consts.timeKeys(RUN, h));
		GFitFixtures.Fx hx = GFitFixtures.resource("switch/t60/hits_extra.txt.gz"), ax = GFitFixtures.resource("switch/t60/atof_extra.txt.gz");
		double[] ped = hx.d("h_ped");
		int[] adc = set.tracksFx.i("h_adc"), hoff = set.tracksFx.i("hoff"), aoff = set.atofFx.i("aoff");
		assertEquals(ped.length, adc.length);
		aEnergy = ax.d("a_energy");
		aTime = ax.d("a_t_ge10");
		for (int k = 0; k < set.n; k++) {
			GFitTrack t = set.tracks[k];
			t.hAdc = new double[hoff[k + 1] - hoff[k]];
			for (int j = 0; j < t.hAdc.length; j++) t.hAdc[j] = adc[hoff[k] + j];
			t.hPed = Arrays.copyOfRange(ped, hoff[k], hoff[k + 1]);
			set.events[k].aEnergy = Arrays.copyOfRange(aEnergy, aoff[k], aoff[k + 1]);
			set.events[k].aTime = Arrays.copyOfRange(aTime, aoff[k], aoff[k + 1]);
		}
	}

	private static double rel(double got, double exp) {
		if (Double.isNaN(got) || Double.isNaN(exp)) return Double.isNaN(got) && Double.isNaN(exp) ? 0.0 : Double.POSITIVE_INFINITY;
		if (got == exp) return 0.0;
		return Math.abs(got - exp) / Math.max(Math.abs(exp), 1e-300);
	}

	/** runs the sequence on every track and compares with the expectation file; asserts the share of tracks inside the tolerances. */
	private static void compare(String label, boolean time, boolean recovery, String expected, double minShare) {
		GFitConfig cfg = GFitConfig.pythonDefault().withSwitches(time, true);
		GlobalTrackFitter ft = GFitFullSearchTest.fitter(cfg);
		HitRecovery hr = recovery ? new HitRecovery(ft, new TrackCost(GFitFullSearchTest.prop, GFitFullSearchTest.geo, cfg), GFitFullSearchTest.geo, cfg) : null;
		SwitchedFit sw = new SwitchedFit(ft, hr, GFitFullSearchTest.prop, GFitFullSearchTest.geo, StoppingTables.proton());
		SwitchConstants.DedxConstants dx = consts.dedx(RUN, time);
		GFitFixtures.Fx ex = GFitFixtures.resource("switch/t60/" + expected + ".txt.gz");
		double[] ep = ex.d("est_p"), esg = ex.d("est_sig_p"), elnE = ex.d("dedx_lnE"), elnp = ex.d("dedx_lnp"), esl = ex.d("dedx_sig_lnp");
		int[] en = ex.i("dedx_n"), est = ex.i("est_status"), efl = ex.i("est_flags");
		GFitWorkspace ws = GFitFullSearchTest.workspace();
		GFitResult r = new GFitResult(GFitFullSearchTest.maxHits);
		int nOk = 0, nLnE = 0, nN = 0, nP = 0, nFew = 0, nFewPy = 0, nFlag = 0;
		StringBuilder out = new StringBuilder();
		double wLnE = 0, wP = 0;
		for (int k = 0; k < set.n; k++) {
			GFitEvent e = set.events[k];
			sw.run(set.tracks[k], e, consts.atofPosition(RUN, e.half), timeTab[e.half], false, dx, ws, r);
			double dE = rel(r.dedxLnE, elnE[k]), dP = rel(r.pDedx, ep[k]);
			boolean few = (r.flags & DedxProfile.FLAG_FEWHITS) != 0, fewPy = (efl[k] & 64) != 0;
			boolean cLnE = dE <= 1e-12, cN = r.dedxN == en[k], cP = dP <= 1e-6, cFlag = few == fewPy;
			if (cLnE) nLnE++;
			if (cN) nN++;
			if (cP) nP++;
			if (cFlag) nFlag++;
			if (few) nFew++;
			if (fewPy) nFewPy++;
			if (cLnE && cN && cP && cFlag) nOk++;
			else out.append(String.format(" row %d trackid %d: lnE %.2e n %d/%d p %.2e (java %.6g python %.6g) flag512 %b python64 %b status java %d;", k,
				set.tracks[k].trackid, dE, r.dedxN, en[k], dP, r.pDedx, ep[k], few, fewPy, r.status));
			if (cLnE) wLnE = Math.max(wLnE, dE);
			if (cP) wP = Math.max(wP, dP);
			assertEquals("est_status 0 iff p_dedx finite, row " + k, est[k] == 0, Double.isFinite(r.pDedx));
		}
		System.out.printf("[%s] %d tracks: dedx_lnE equal %d, dedx_n equal %d, p_dedx within 1e-6 %d (worst inside %.2e), flag 512 = Python 64 on %d (Java %d, Python %d), all four %d; outside:%s%n",
			label, set.n, nLnE, nN, nP, wP, nFlag, nFew, nFewPy, nOk, out);
		assertTrue(label + " tracks outside:" + out, nOk >= minShare * set.n);
	}

	@Test
	public void dedxAloneOnDefaultFit() { compare("dedx t60", false, false, "expected_ge10_dedx", 0.95); }

	@Test
	public void dedxWithTimeAndRecovery() { compare("all t60", true, true, "expected_ge10_all", 0.95); }

	/** dE/dx off: the sequence is the plain fit and recovery bit for bit, and the dE/dx outputs stay NaN. */
	@Test
	public void switchOffUnchanged() {
		GFitConfig cfg = GFitConfig.pythonDefault();
		GlobalTrackFitter ft = GFitFullSearchTest.fitter(cfg);
		HitRecovery hr = new HitRecovery(ft, new TrackCost(GFitFullSearchTest.prop, GFitFullSearchTest.geo, cfg), GFitFullSearchTest.geo, cfg);
		SwitchedFit sw = new SwitchedFit(ft, hr, GFitFullSearchTest.prop, GFitFullSearchTest.geo, StoppingTables.proton());
		GFitWorkspace ws = GFitFullSearchTest.workspace();
		GFitResult a = new GFitResult(GFitFullSearchTest.maxHits), b = new GFitResult(GFitFullSearchTest.maxHits);
		for (int k = 0; k < 12; k++) {
			GFitEvent e = set.events[k];
			AtofPositionTable tab = consts.atofPosition(RUN, e.half);
			int na = sw.run(set.tracks[k], e, tab, null, false, null, ws, a);
			ft.fit(set.tracks[k], e, tab, false, ws, b);
			int nb = hr.recover(set.tracks[k], e, tab, false, b, ws, b);
			assertEquals(nb, na);
			assertNull("track " + k, GFitFullSearchTest.diff(a, b));
			assertTrue(Double.isNaN(a.pDedx) && Double.isNaN(a.sigPDedx) && a.dedxN == 0 && !a.prompt);
			assertEquals(0, a.flags & DedxProfile.FLAG_FEWHITS);
		}
	}
}
