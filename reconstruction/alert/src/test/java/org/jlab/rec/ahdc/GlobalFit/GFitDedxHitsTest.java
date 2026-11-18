package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/**
 * DedxHits against switch/t60/dedx_hits.txt.gz (export_switch_dedx.py: hit table of the default fit, dedx.py hit_table and
 * dedx_cal.py calibrate) at the stored default-fit state of each track, with the constants of switch/dedx_const_23014.txt.gz.
 */
public class GFitDedxHitsTest {

	private static final double TOL = 1e-9;

	private static double rel(double a, double e, double floor) {
		if (Double.isNaN(a) && Double.isNaN(e)) return 0.0;
		if (a == e) return 0.0;
		return Math.abs(a - e) / Math.max(floor, Math.abs(e));
	}

	@Test
	public void dedxHitsMatchPython() {
		GFitFixtures.Fx f = GFitFixtures.resource("unit/field.txt.gz");
		double[] fp = f.d("field_par");
		FieldGrid field = FieldGrid.fromArrays(fp[0], fp[1], fp[2], fp[3], f.shape("br")[0], f.shape("br")[1], f.d("br"), f.d("bz"));
		double[] beam = GFitFixtures.resource("unit/stopping.txt.gz").d("beam");
		StoppingTables tab = StoppingTables.proton();
		Propagator prop = new Propagator(field, MaterialStack.dataBest(), tab, beam[0], beam[1]);
		GFitFixtures.Fx w = GFitFixtures.resource("unit/wires.txt.gz");
		CellGeometry geo = CellGeometry.fromWireEnds(w.d("WO"), w.d("WE"), w.i("layoff"));
		GFitFixtures.Fx c = GFitFixtures.resource("switch/dedx_const_23014.txt.gz");
		double[] offset = c.d("dx_offset"), f3All = c.d("dx_f3"), f3Norm = c.d("dx_f3_norm");
		GFitFixtures.Fx x = GFitFixtures.resource("switch/t60/dedx_hits.txt.gz");

		int n = x.shape("d_sel")[0];
		int[] half = x.i("d_half"), status = x.i("d_status"), tstat = x.i("d_tstat"), nOn = x.i("d_n_on"), nUsed = x.i("d_n_used"),
			hoff = x.i("d_hoff"), hwire = x.i("dh_wire"), hkind = x.i("dh_kind"), hkept = x.i("dh_kept"), huse = x.i("dh_use");
		double[] st = x.d("d_state"), kk = x.d("d_kinks"), p0 = x.d("d_p_traj0"), gm = x.d("d_gm"), lnE = x.d("d_lnE");
		double[] hadc = x.d("dh_adc"), hped = x.d("dh_ped"), hdoca = x.d("dh_doca");
		double[] ePath = x.d("dh_path"), eS = x.d("dh_s"), eX = x.d("dh_X"), eLnu = x.d("dh_lnu");

		DedxHits dh = new DedxHits(prop, geo, tab, 64);
		GFitResult r = new GFitResult(64);
		double[] ad = new double[72], pd = new double[72], dd = new double[72], f3 = new double[10];
		double wPath = 0, wS = 0, wX = 0, wLnu = 0, wP0 = 0, wGm = 0, wLnE = 0;
		int bad = 0, nhits = 0, nkept = 0;
		StringBuilder msg = new StringBuilder();
		for (int t = 0; t < n; t++) {
			r.clear();
			r.status = status[t];
			r.u[0] = st[t * 5]; r.u[1] = st[t * 5 + 1]; r.u[2] = -Math.log(Math.abs(st[t * 5 + 4])); r.u[3] = st[t * 5 + 2]; r.u[4] = st[t * 5 + 3];
			for (int j = 0; j < 10; j++) r.u[5 + j] = kk[t * 10 + j];
			int a = hoff[t], nh = hoff[t + 1] - a;
			r.nh = nh;
			for (int k = 0; k < nh; k++) {
				r.hitWire[k] = hwire[a + k]; r.hitKind[k] = hkind[a + k];
				r.hitRow[k] = Double.isNaN(hadc[a + k]) ? -1 : k;
				ad[k] = hadc[a + k]; pd[k] = hped[a + k]; dd[k] = hdoca[a + k];
			}
			System.arraycopy(f3All, half[t] * 10, f3, 0, 10);
			int s = dh.compute(r, ad, pd, dd, offset, f3, f3Norm[half[t]]);
			assertEquals("status track " + t, tstat[t], s);
			double tw = 0;
			if (s == 0) {
				assertEquals("n_on track " + t, nOn[t], dh.nOn);
				tw = Math.max(rel(dh.pTraj0, p0[t], 1e-300), Math.max(rel(dh.gm, gm[t], 1.0), rel(dh.lnE, lnE[t], 1.0)));
				wP0 = Math.max(wP0, rel(dh.pTraj0, p0[t], 1e-300));
				wGm = Math.max(wGm, rel(dh.gm, gm[t], 1.0));
				wLnE = Math.max(wLnE, rel(dh.lnE, lnE[t], 1.0));
				assertEquals("n_used track " + t, nUsed[t], dh.nUsed);
				for (int k = 0; k < nh; k++) {
					nhits++;
					assertEquals("kept track " + t + " hit " + k, hkept[a + k] != 0, dh.kept[k]);
					assertEquals("use track " + t + " hit " + k, huse[a + k] != 0, dh.use[k]);
					if (!dh.kept[k]) continue;
					nkept++;
					double dp = rel(dh.path[k], ePath[a + k], 1e-300), ds = rel(dh.s[k], eS[a + k], 1e-300), dx = rel(dh.X[k], eX[a + k], 1e-300),
						dl = rel(dh.lnu[k], eLnu[a + k], 1.0);
					wPath = Math.max(wPath, dp); wS = Math.max(wS, ds); wX = Math.max(wX, dx); wLnu = Math.max(wLnu, dl);
					tw = Math.max(tw, Math.max(Math.max(dp, ds), Math.max(dx, dl)));
				}
			}
			if (!(tw <= TOL)) {
				bad++;
				msg.append(" track ").append(t).append(" (").append(tw).append(")");
			}
		}
		System.out.printf("DedxHits: %d tracks, %d hits, %d kept; worst relative: path %.3e s %.3e X %.3e lnu(floor 1) %.3e p_traj0 %.3e gm %.3e lnE %.3e (limit %.0e)%n",
			n, nhits, nkept, wPath, wS, wX, wLnu, wP0, wGm, wLnE, TOL);
		assertTrue("tracks outside " + TOL + ":" + msg, bad == 0);
	}
}
