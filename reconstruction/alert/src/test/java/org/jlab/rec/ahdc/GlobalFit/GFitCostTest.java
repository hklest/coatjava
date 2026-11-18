package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.BeforeClass;
import org.junit.Test;

/** HitTerm and TrackCost against unit/prop_states.txt (wire crossings) and unit/cost_states.txt (cost, gradient, Hessians). */
public class GFitCostTest {

	private static Propagator prop;
	private static CellGeometry geo;
	private static TrackCost tc;
	private static GFitFixtures.Fx ps, cs, ct;
	private static GFitTrack[] tracks;
	private static final int NS = 90, NPAR = 15;
	// results of cost(want = true, ATOF off) at the 90 states
	private static double[] fF, fHit, fZ, fD0, fK, fG, fHc, fHe, fOut;
	private static int[] fhoff;

	@BeforeClass
	public static void load() {
		GFitFixtures.Fx f = GFitFixtures.resource("unit/field.txt.gz");
		double[] fp = f.d("field_par");
		FieldGrid field = FieldGrid.fromArrays(fp[0], fp[1], fp[2], fp[3], f.shape("br")[0], f.shape("br")[1], f.d("br"), f.d("bz"));
		double[] beam = GFitFixtures.resource("unit/stopping.txt.gz").d("beam");
		prop = new Propagator(field, MaterialStack.dataBest(), StoppingTables.proton(), beam[0], beam[1]);
		GFitFixtures.Fx w = GFitFixtures.resource("unit/wires.txt.gz");
		geo = CellGeometry.fromWireEnds(w.d("WO"), w.d("WE"), w.i("layoff"));
		tc = new TrackCost(prop, geo, GFitConfig.pythonDefault());
		ps = GFitFixtures.resource("unit/prop_states.txt.gz");
		cs = GFitFixtures.resource("unit/cost_states.txt.gz");
		ct = GFitFixtures.resource("unit/cost_tracks.txt.gz");
		int[] hoff = ct.i("ct_hoff"), hg = ct.i("ct_hg"), hism = ct.i("ct_hism"), hcls = ct.i("ct_hcls");
		double[] hm = ct.d("ct_hm"), hsig = ct.d("ct_hsig"), hh = ct.d("ct_hh"), tcr = ct.d("ct_TC");
		int n = ct.shape("ct_track")[0];
		tracks = new GFitTrack[n];
		for (int j = 0; j < n; j++) {
			int a = hoff[j], b = hoff[j + 1];
			tracks[j] = new GFitTrack(j, java.util.Arrays.copyOfRange(hg, a, b), java.util.Arrays.copyOfRange(hm, a, b),
				java.util.Arrays.copyOfRange(hsig, a, b), java.util.Arrays.copyOfRange(hh, a, b), java.util.Arrays.copyOfRange(hism, a, b),
				java.util.Arrays.copyOfRange(hcls, a, b), new int[b - a], tcr[j * 7 + 1], tcr[j * 7 + 2], tcr[j * 7 + 3], tcr[j * 7 + 4]);
		}
		fhoff = cs.i("cs_hoff");
		fF = new double[NS]; fHit = new double[NS]; fZ = new double[NS]; fD0 = new double[NS]; fK = new double[NS];
		fG = new double[NS * NPAR]; fHc = new double[NS * 225]; fHe = new double[NS * 225]; fOut = new double[fhoff[NS] * 3];
		GFitWorkspace ws = new GFitWorkspace(64);
		for (int s = 0; s < NS; s++) {
			run(s, ws, true);
			fF[s] = lastF;
			fHit[s] = ws.fHit; fZ[s] = ws.fZ; fD0[s] = ws.fD0; fK[s] = ws.fK;
			System.arraycopy(ws.g, 0, fG, s * NPAR, NPAR);
			System.arraycopy(ws.hc, 0, fHc, s * 225, 225);
			System.arraycopy(ws.he, 0, fHe, s * 225, 225);
			System.arraycopy(ws.outH, 0, fOut, fhoff[s] * 3, ws.nh * 3);
		}
	}

	private static double lastF;

	/** loads the hits of state s's track and evaluates cost (ATOF off, th0f = cs_th0). */
	private static void run(int s, GFitWorkspace ws, boolean want) {
		int j = java.util.Arrays.binarySearch(ct.i("ct_track"), cs.i("cs_track")[s]);
		tc.loadHits(tracks[j], ws);
		assertEquals(fhoff[s + 1] - fhoff[s], ws.nh);
		lastF = tc.cost(cs.row("cs_u", s), 1, NPAR, want, ct.d("ct_TC")[j * 7], cs.row("cs_th0", s), ws);
	}

	private static double dev(double a, double e) {
		if (Double.isNaN(a) && Double.isNaN(e)) return 0.0;
		if (a == e) return 0.0;
		return Math.abs(a - e) / Math.max(1.0, Math.abs(e));
	}

	private static double rel(double a, double e) {
		if (a == e) return 0.0;
		return Math.abs(a - e) / Math.abs(e);
	}

	@Test
	public void wireDistanceMatchesPython() {
		double[] ru = ps.d("ps_u"), pd = ps.d("pw_d"), px = ps.d("pw_x");
		int[] off = ps.i("pw_off"), pg = ps.i("pw_g"), pf = ps.i("pw_f");
		int n = ps.shape("ps_u")[0];
		GFitWorkspace ws = new GFitWorkspace(0);
		double[] u = new double[NPAR], gs = new double[GFitWorkspace.NMAX], ds = new double[GFitWorkspace.NMAX], r = new double[12];
		double wd = 0, wx = 0;
		int count = 0;
		for (int s = 0; s < n; s++) {
			System.arraycopy(ru, s * NPAR, u, 0, NPAR);
			int nn = prop.run(u, ws);
			for (int i = off[s]; i < off[s + 1]; i++) {
				int f = HitTerm.pca1(ws.out, nn, geo.wo, geo.we, pg[i], gs, ds, r);
				assertEquals("flag crossing " + i, pf[i], f);
				wd = Math.max(wd, dev(r[0], pd[i]));
				for (int c = 0; c < 3; c++) wx = Math.max(wx, dev(r[1 + c], px[3 * i + c]));
				count++;
			}
		}
		System.out.printf("pca1: %d crossings, worst d %.3e, worst point %.3e (limit 1e-9, floor 1 mm)%n", count, wd, wx);
		assertEquals(1152, count);
		assertTrue(wd <= 1e-9);
		assertTrue(wx <= 1e-9);
	}

	@Test
	public void costTermsMatchPython() {
		double[] eF = cs.d("cs_F"), eFa = cs.d("cs_Fatof"), eHit = cs.d("cs_Fhit"), eZ = cs.d("cs_Fz"), eD0 = cs.d("cs_Fd0"), eK = cs.d("cs_Fk");
		double[] el = cs.d("cs_l"), ed = cs.d("cs_d"), es = cs.d("cs_s");
		double wF = 0, wHit = 0, wZ = 0, wD0 = 0, wK = 0, wL = 0, wD = 0, wS = 0;
		for (int s = 0; s < NS; s++) {
			wF = Math.max(wF, rel(fF[s], eF[s] - eFa[s]));
			wHit = Math.max(wHit, rel(fHit[s], eHit[s]));
			wZ = Math.max(wZ, rel(fZ[s], eZ[s]));
			wD0 = Math.max(wD0, rel(fD0[s], eD0[s]));
			wK = Math.max(wK, rel(fK[s], eK[s]));
			for (int k = fhoff[s]; k < fhoff[s + 1]; k++) {
				wL = Math.max(wL, dev(fOut[3 * k], el[k]));
				wD = Math.max(wD, dev(fOut[3 * k + 1], ed[k]));
				wS = Math.max(wS, dev(fOut[3 * k + 2], es[k]));
			}
		}
		System.out.printf("cost terms worst relative: F %.3e Fhit %.3e Fz %.3e Fd0 %.3e Fk %.3e; per hit l %.3e d %.3e s %.3e (limit 1e-9)%n",
			wF, wHit, wZ, wD0, wK, wL, wD, wS);
		for (double w : new double[] {wF, wHit, wZ, wD0, wK, wL, wD, wS}) assertTrue(w <= 1e-9);
	}

	@Test
	public void gradientMatchesPython() {
		double[] eg = cs.d("cs_g"), ega = cs.d("cs_gatof");
		double worst = 0;
		for (int s = 0; s < NS; s++) {
			double nd = 0, ne = 0;
			for (int j = 0; j < NPAR; j++) {
				double e = eg[s * NPAR + j] - ega[s * NPAR + j], d = fG[s * NPAR + j] - e;
				nd += d * d;
				ne += e * e;
			}
			// converged states have |g| ~ 1e-7: the scale is floored at 1 (cost units)
			worst = Math.max(worst, Math.sqrt(nd) / Math.max(Math.sqrt(ne), 1.0));
		}
		int[] sub = cs.i("cs_sub");
		double[] hc = cs.d("cs_Hc_tri"), he = cs.d("cs_He_tri"), hac = cs.d("cs_Hatof_c_tri"), hae = cs.d("cs_Hatof_e_tri");
		double wHc = 0, wHe = 0;
		for (int q = 0; q < sub.length; q++) {
			int s = sub[q];
			double dc = 0, nc = 0, de = 0, ne = 0;
			int t = 0;
			for (int j = 0; j < NPAR; j++)
				for (int m = j; m < NPAR; m++, t++) {
					double w = j == m ? 1.0 : 2.0;
					double ec = hc[q * 120 + t] - hac[q * 120 + t], ee = he[q * 120 + t] - hae[q * 120 + t];
					double d1 = fHc[s * 225 + j * NPAR + m] - ec, d2 = fHe[s * 225 + j * NPAR + m] - ee;
					assertEquals("symmetric hc", fHc[s * 225 + j * NPAR + m], fHc[s * 225 + m * NPAR + j], 0.0);
					assertEquals("symmetric he", fHe[s * 225 + j * NPAR + m], fHe[s * 225 + m * NPAR + j], 0.0);
					dc += w * d1 * d1; nc += w * ec * ec; de += w * d2 * d2; ne += w * ee * ee;
				}
			wHc = Math.max(wHc, Math.sqrt(dc / nc));
			wHe = Math.max(wHe, Math.sqrt(de / ne));
		}
		System.out.printf("gradient worst |dg|/max(|g|,1) %.3e; Hessians worst Frobenius ratio clipped %.3e exact %.3e (limit 1e-6)%n", worst, wHc, wHe);
		assertEquals(18, sub.length);
		assertTrue(worst <= 1e-6);
		assertTrue(wHc <= 1e-6);
		assertTrue(wHe <= 1e-6);
	}

	@Test
	public void restartGradientEqualsFull() {
		GFitWorkspace a = new GFitWorkspace(64), b = new GFitWorkspace(64);
		TrackCost full = new TrackCost(prop, geo, GFitConfig.pythonDefault());
		full.useRestart = false;
		for (int s = 0; s < NS; s++) {
			run(s, a, true);
			int j = java.util.Arrays.binarySearch(ct.i("ct_track"), cs.i("cs_track")[s]);
			full.loadHits(tracks[j], b);
			full.cost(cs.row("cs_u", s), 1, NPAR, true, ct.d("ct_TC")[j * 7], cs.row("cs_th0", s), b);
			for (int i = 0; i < a.nh * NPAR; i++) {
				assertEquals("JD state " + s, Double.doubleToLongBits(b.JD[i]), Double.doubleToLongBits(a.JD[i]));
				assertEquals("JS state " + s, Double.doubleToLongBits(b.JS[i]), Double.doubleToLongBits(a.JS[i]));
			}
			for (int i = 0; i < NPAR; i++) assertEquals(Double.doubleToLongBits(b.g[i]), Double.doubleToLongBits(a.g[i]));
			for (int i = 0; i < 225; i++) assertEquals(Double.doubleToLongBits(b.he[i]), Double.doubleToLongBits(a.he[i]));
		}
	}

	@Test
	public void helixCost() {
		TrackCost hx = new TrackCost(prop, geo, GFitConfig.pythonDefault());
		GFitWorkspace ws = new GFitWorkspace(64);
		double[] steps = {1e-4, 1e-6, 1e-6, 1e-4, 1e-6};
		double worst = 0;
		int used = 0;
		for (int s = 0; s < NS; s += 7) {
			if (cs.i("cs_kind")[s] == 2) continue;
			int j = java.util.Arrays.binarySearch(ct.i("ct_track"), cs.i("cs_track")[s]);
			double ze = ct.d("ct_TC")[j * 7];
			double[] u = cs.row("cs_u", s);
			hx.loadHits(tracks[j], ws);
			double F = hx.cost(u, 0, 5, true, ze, new double[6], ws);
			assertTrue("finite F state " + s, !Double.isNaN(F) && !Double.isInfinite(F));
			double[] g = java.util.Arrays.copyOf(ws.g, 5);
			double ng = 0, nd = 0;
			for (int k = 0; k < 5; k++) {
				double[] up = u.clone(), um = u.clone();
				up[k] += steps[k]; um[k] -= steps[k];
				hx.loadHits(tracks[j], ws);
				double fp = hx.cost(up, 0, 5, false, ze, new double[6], ws);
				hx.loadHits(tracks[j], ws);
				double fm = hx.cost(um, 0, 5, false, ze, new double[6], ws);
				double fd = (fp - fm) / (2 * steps[k]);
				nd += (fd - g[k]) * (fd - g[k]);
				ng += g[k] * g[k];
			}
			worst = Math.max(worst, Math.sqrt(nd) / Math.max(Math.sqrt(ng), 1.0));
			used++;
		}
		System.out.printf("helix cost: %d states, gradient vs central difference worst %.3e%n", used, worst);
		assertTrue(used > 5);
		assertTrue(worst <= 1e-4);
	}

	@Test
	public void timing() {
		if (!Boolean.getBoolean("gfit.timing")) return;
		GFitWorkspace ws = new GFitWorkspace(64);
		int s = 2;
		for (int i = 0; i < 2000; i++) run(s, ws, true);
		int reps = 5000;
		long t0 = System.nanoTime();
		for (int i = 0; i < reps; i++) run(s, ws, true);
		double us = (System.nanoTime() - t0) / 1e3 / reps;
		System.out.printf("cost(want) at a 15-parameter state with %d hits: %.1f us (Python one gradient = 16 propagations, 500 us)%n", ws.nh, us);
	}
}
