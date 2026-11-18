package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * DedxPrediction and DedxProfile against the t60 fixtures of switch/t60 (dedx_pred, dedx_profile, dedx_hits; spec in
 * tools/gr_java/export_switch_dedx.py): pred0, sgc, invert, at_lnp, and the post-fit profile with Python's per-hit ln u as input.
 */
public class GFitDedxTest {

	private static GFitFixtures.Fx hits, pred, prof;
	private static GFitFixtureConstants provider;
	private static int[] hoff;

	@BeforeClass
	public static void load() {
		hits = GFitFixtures.resource("switch/t60/dedx_hits.txt.gz");
		pred = GFitFixtures.resource("switch/t60/dedx_pred.txt.gz");
		prof = GFitFixtures.resource("switch/t60/dedx_profile.txt.gz");
		provider = new GFitFixtureConstants(false);
		hoff = hits.i("d_hoff");
	}

	/** relative difference, 0 for two NaN, infinite for one NaN. */
	private static double rel(double got, double exp) {
		if (Double.isNaN(got) || Double.isNaN(exp)) return Double.isNaN(got) && Double.isNaN(exp) ? 0.0 : Double.POSITIVE_INFINITY;
		if (got == exp) return 0.0;
		return Math.abs(got - exp) / Math.max(Math.abs(exp), 1e-300);
	}

	/** kept hits of track t: x, lnu, use compacted to the front of the arrays, returns the count. */
	private static int kept(int t, double[] x, double[] lnu, boolean[] use) {
		double[] dx = hits.d("dh_X"), dl = hits.d("dh_lnu"), du = hits.d("dh_use"), dk = hits.d("dh_kept");
		int n = 0;
		for (int i = hoff[t]; i < hoff[t + 1]; i++)
			if (dk[i] != 0) {
				x[n] = dx[i];
				lnu[n] = dl[i];
				use[n] = du[i] != 0;
				n++;
			}
		return n;
	}

	@Test
	public void dedxPredictionMatchesPython() {
		SwitchConstants.DedxConstants c = provider.dedx(23014, false);
		DedxPrediction dp = new DedxPrediction(StoppingTables.proton());
		int[] track = pred.i("pr_track"), nUsed = pred.i("pr_n_used"), gi = pred.i("pr_grid_idx"), si = pred.i("pr_sgc_idx");
		int[] half = hits.i("d_half");
		double aRun = pred.d0("pr_a_run"), g = pred.d0("pr_g");
		assertEquals(0.0, rel(c.aRun, aRun), 1e-12);
		assertEquals(0.0, rel(c.g, g), 1e-12);
		double[] p0 = pred.d("pr_pred0"), sg = pred.d("pr_sgc"), qE = pred.d("pr_q_lnE"), il = pred.d("pr_inv_lnp"), isl = pred.d("pr_inv_slope"),
			qp = pred.d("pr_q_lnp"), a0 = pred.d("pr_at_pred0"), as = pred.d("pr_at_sgc"), at = pred.d("pr_at_predt");
		int[] nanFirst = pred.i("pr_nan_first");
		double[] x = new double[hoff[hoff.length - 1] + 1], lnu = new double[x.length];
		boolean[] use = new boolean[x.length];
		double mPred = 0, mSg = 0, mInv = 0, mSlope = 0, mAt = 0;
		int nNan = 0, nCmp = 0, nBad = 0;
		for (int r = 0; r < track.length; r++) {
			int t = track[r];
			int n = kept(t, x, lnu, use);
			dp.compute(x, n, nUsed[r], half[t], c);
			double[] predt = dp.response(c);
			int nan = 0;
			for (int j = 0; j < DedxPrediction.NG; j++) if (Double.isNaN(dp.pred0[j])) nan++;
			assertEquals("NaN grid points of track " + r, nanFirst[r], nan);
			for (int k = 0; k < gi.length; k++) {
				double d = rel(dp.pred0[gi[k]], p0[r * gi.length + k]);
				mPred = Math.max(mPred, d);
				if (d > 1e-9) nBad++;
				if (Double.isNaN(dp.pred0[gi[k]])) nNan++;
				nCmp++;
			}
			for (int k = 0; k < si.length; k++) mSg = Math.max(mSg, rel(dp.sgc[si[k]], sg[r * si.length + k]));
			for (int q = 0; q < 7; q++) {
				int o = r * 7 + q;
				double ln = dp.invert(qE[o], predt);
				mInv = Math.max(mInv, rel(ln, il[o]));
				mSlope = Math.max(mSlope, rel(dp.slope, isl[o]));
				mAt = Math.max(mAt, rel(DedxPrediction.atLnp(dp.pred0, qp[q]), a0[o]));
				mAt = Math.max(mAt, rel(DedxPrediction.atLnp(dp.sgc, qp[q]), as[o]));
				mAt = Math.max(mAt, rel(DedxPrediction.atLnp(predt, qp[q]), at[o]));
			}
		}
		System.out.println("GFitDedxTest pred: max rel pred0 " + mPred + " sgc " + mSg + " invert " + mInv + " slope " + mSlope + " at_lnp " + mAt
			+ "; NaN points " + nNan + " of " + nCmp);
		assertEquals("pred0 points outside 1e-9", 0, nBad);
		assertTrue("sgc " + mSg, mSg < 1e-9);
		assertTrue("invert " + mInv, mInv < 1e-9);
		assertTrue("slope " + mSlope, mSlope < 1e-9);
		assertTrue("at_lnp " + mAt, mAt < 1e-9);
	}

	private static void profileVariant(String v, boolean all) {
		SwitchConstants.DedxConstants c = provider.dedx(23014, all);
		DedxProfile dp = new DedxProfile(StoppingTables.proton());
		int nt = prof.shape("p_half")[0];
		int[] half = prof.i("p_half"), pn = prof.i("p_n");
		double[] lnE = prof.d("p_lnE"), dl = prof.d("p_dedx_lnp"), ds = prof.d("p_dedx_sig_lnp"), bp = prof.d(v + "_base_p"), bs = prof.d(v + "_base_sig_p"),
			cv = prof.d(v + "_class_var"), pl = prof.d(v + "_prof_lnp"), pe = prof.d(v + "_prof_err"), ep = prof.d(v + "_est_p"),
			es = prof.d(v + "_est_sig_p");
		int[] bst = prof.i(v + "_base_status"), est = prof.i(v + "_est_status"), efl = prof.i(v + "_est_flags");
		boolean[] prompt = new boolean[nt];
		for (int k = 0; k < nt; k++) prompt[k] = prof.d("p_prompt")[k] != 0;
		int[] sel = prof.i("p_sel"), dsel = hits.i("d_sel");
		double[] x = new double[hoff[hoff.length - 1] + 1], lnu = new double[x.length];
		boolean[] use = new boolean[x.length];
		double mE = 0, mLnp = 0, mSigLnp = 0, mP = 0, mSig = 0, mCv = 0, mPl = 0, mPe = 0;
		StringBuilder out = new StringBuilder();
		int nFew = 0, nNanSig = 0;
		for (int k = 0; k < nt; k++) {
			assertEquals(sel[k], dsel[k]);
			int n = kept(k, x, lnu, use);
			dp.evaluate(x, lnu, use, n, bp[k], bs[k], bst[k], prompt[k], half[k], c);
			assertEquals("n_used " + k, pn[k], dp.n);
			assertEquals("status " + k, est[k], dp.status);
			assertEquals("flags " + k, efl[k] == 64 ? DedxProfile.FLAG_FEWHITS : efl[k], dp.flags);
			if (dp.flags != 0) nFew++;
			if (Double.isNaN(es[k])) nNanSig++;
			double[] d = {rel(dp.lnE, lnE[k]), rel(dp.dedxLnp, dl[k]), rel(dp.dedxSigLnp, ds[k]), rel(dp.pDedx, ep[k]), rel(dp.sigPDedx, es[k]),
				rel(dp.classVar, cv[k]), rel(dp.profLnp, pl[k]), rel(dp.profErr, pe[k])};
			mE = Math.max(mE, d[0]);
			mLnp = Math.max(mLnp, d[1]);
			mSigLnp = Math.max(mSigLnp, d[2]);
			mP = Math.max(mP, d[3]);
			mSig = Math.max(mSig, d[4]);
			mCv = Math.max(mCv, d[5]);
			mPl = Math.max(mPl, d[6]);
			mPe = Math.max(mPe, d[7]);
			if (d[0] > 1e-12 || d[1] > 1e-6 || d[2] > 1e-6 || d[3] > 1e-6 || d[4] > 1e-6 || d[6] > 1e-6 || d[7] > 1e-6)
				out.append(" track ").append(k).append(String.format(" lnE %.2e lnp %.2e sig_lnp %.2e p %.2e sig %.2e prof %.2e err %.2e;", d[0], d[1], d[2], d[3], d[4], d[6], d[7]));
		}
		System.out.println("GFitDedxTest profile " + v + ": max rel lnE " + mE + " dedx_lnp " + mLnp + " dedx_sig_lnp " + mSigLnp + " p " + mP + " sig_p " + mSig
			+ " class_var " + mCv + " prof_lnp " + mPl + " prof_err " + mPe + "; few " + nFew + ", NaN sig " + nNanSig + (out.length() > 0 ? "; outside:" + out : ""));
		assertTrue("tracks outside tolerance:" + out, out.length() == 0);
	}

	@Test
	public void dedxProfileMatchesPython() {
		profileVariant("dd", false);
		profileVariant("al", true);
	}

	@Test
	public void dedxFewHits() {
		SwitchConstants.DedxConstants c = provider.dedx(23014, true);
		DedxProfile dp = new DedxProfile(StoppingTables.proton());
		double[] x = {0.02, 0.05, 0.08}, lnu = {1.0, 1.1, 5.0};
		boolean[] use = {true, true, false};
		dp.evaluate(x, lnu, use, 3, 0.45, 0.03, 0, false, 0, c);
		assertEquals(2, dp.n);
		assertEquals(0, dp.status);
		assertEquals(DedxProfile.FLAG_FEWHITS, dp.flags);
		assertEquals(0.45, dp.pDedx, 0.0);
		assertEquals(0.03, dp.sigPDedx, 0.0);
		assertTrue(Double.isNaN(dp.dedxLnp));
		// failed base: nothing, status 1
		dp.evaluate(x, lnu, use, 3, 0.45, 0.03, 1, false, 0, c);
		assertEquals(1, dp.status);
		assertEquals(0, dp.flags);
		assertTrue(Double.isNaN(dp.pDedx) && Double.isNaN(dp.sigPDedx));
	}
}
