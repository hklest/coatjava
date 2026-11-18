package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.BeforeClass;
import org.junit.Test;

/** AtofTerm (atof.py attach, _layer_eval, term_cost) against t60 candidates and the unit/cost_states.txt ATOF values. */
public class GFitAtofTest {

	private static Propagator prop;
	private static GFitFixtures.Fx tabs;

	@BeforeClass
	public static void load() {
		GFitFixtures.Fx f = GFitFixtures.resource("unit/field.txt.gz");
		double[] fp = f.d("field_par");
		int[] sh = f.shape("br");
		FieldGrid field = FieldGrid.fromArrays(fp[0], fp[1], fp[2], fp[3], sh[0], sh[1], f.d("br"), f.d("bz"));
		double[] beam = GFitFixtures.resource("unit/stopping.txt.gz").d("beam");
		prop = new Propagator(field, MaterialStack.dataBest(), StoppingTables.proton(), beam[0], beam[1]);
		tabs = GFitFixtures.resource("unit/atof_tables.txt.gz");
	}

	private static double[] par(int half) { return tabs.row("par_23014", half); }

	/** recal_23014 layout [half][offset, scale][60]. */
	private static double[] recal(int half, int k) {
		double[] a = tabs.d("recal_23014"), o = new double[60];
		System.arraycopy(a, (half * 2 + k) * 60, o, 0, 60);
		return o;
	}

	private static double rel(double got, double exp) { return Math.abs(got - exp) / Math.max(Math.abs(exp), 1e-300); }

	/** wedge time pdf of atof._tpdf written out (value only). */
	private static double tpdfValue(double r, double[] p) {
		double[] s = {p[AtofTerm.P_TS1], p[AtofTerm.P_TS2], p[AtofTerm.P_TSL]};
		double[] c = {r - p[AtofTerm.P_TOFF], r - p[AtofTerm.P_TOFF], r - p[AtofTerm.P_TOFF] - p[AtofTerm.P_TML]};
		double[] n = new double[3];
		for (int i = 0; i < 3; i++) n[i] = Math.exp(-0.5 * (c[i] / s[i]) * (c[i] / s[i])) / Math.sqrt(2 * Math.PI) / s[i];
		double f2 = p[AtofTerm.P_TF2], fl = p[AtofTerm.P_TFL], fb = p[AtofTerm.P_TFB];
		return (1.0 - fl - fb) * ((1 - f2) * n[0] + f2 * n[1]) + fl * n[2] + fb / (2.0 * p[AtofTerm.P_TWIN]);
	}

	@Test
	public void attachMatchesPython() {
		GFitReplaySet set = GFitReplaySet.resource("t60");
		GFitWorkspace ws = new GFitWorkspace(0);
		AtofTerm.TrackState st = new AtofTerm.TrackState();
		AtofTerm.TrackState sb = new AtofTerm.TrackState();
		double wSeed = 0, wCz = 0, wPhi = 0, wCt = 0, wBank = 0;
		int nBarCand = 0, nWedCand = 0;
		for (int k = 0; k < set.n; k++) {
			GFitTrack t = set.tracks[k];
			GFitEvent e = set.events[k];
			int h = e.half;
			AtofTerm.attach(t, e, par(h), recal(h, 0), recal(h, 1), false, prop, ws, st);
			for (int q = 0; q < 20; q++) {
				double d = Math.abs(st.record[q] - e.atofSeed[q]);
				double r = d / Math.max(Math.abs(e.atofSeed[q]), 1e-9);
				wSeed = Math.max(wSeed, Math.min(r, d / 1e-9));
				assertTrue("track " + k + " seed record " + q + " java " + st.record[q] + " python " + e.atofSeed[q],
					d <= 1e-9 * Math.max(Math.abs(e.atofSeed[q]), 1.0) || r <= 1e-9);
			}
			assertEquals("track " + k + " fnone", e.fnone, st.fnone);
			assertEquals("track " + k + " ncand", e.cRow.length, st.ncand);
			for (int c = 0; c < st.ncand; c++) {
				String m = "track " + k + " cand " + c;
				assertEquals(m + " row", e.cRow[c], st.crow[c]);
				assertEquals(m + " bar", e.cBar[c], st.cbar[c]);
				assertEquals(m + " sec", e.cSec[c], st.csec[c]);
				assertEquals(m + " comp", e.cComp[c], st.ccomp[c]);
				wCz = Math.max(wCz, rel(st.cz[c], e.cz[c]));
				wPhi = Math.max(wPhi, Math.abs(st.cphi[c] - e.cphi[c]) / Math.max(Math.abs(e.cphi[c]), 1e-9));
				wCt = Math.max(wCt, rel(st.ct[c], e.ct[c]));
				assertTrue(m + " cz java " + st.cz[c] + " python " + e.cz[c], rel(st.cz[c], e.cz[c]) <= 1e-9 || Math.abs(st.cz[c] - e.cz[c]) <= 1e-9);
				assertTrue(m + " cphi", Math.abs(st.cphi[c] - e.cphi[c]) <= 1e-9 * Math.max(Math.abs(e.cphi[c]), 1e-3));
				assertTrue(m + " weight java " + st.ct[c] + " python " + e.ct[c], rel(st.ct[c], e.ct[c]) <= 1e-9);
				// Assumption 4: bar index b = 4 * sector + layer, sector 0-14, layer 0-3, component 10; wedge rows carry their own sector
				int row = e.cRow[c];
				assertEquals(m + " c_sec vs a_sector", e.sector[row], e.cSec[c]);
				assertEquals(m + " c_comp vs a_comp", e.comp[row], e.cComp[c]);
				if (e.cBar[c] == 1) {
					nBarCand++;
					assertTrue(m, e.sector[row] >= 0 && e.sector[row] <= 14 && e.layer[row] >= 0 && e.layer[row] <= 3);
					assertEquals(m + " bar comp", 10, e.comp[row]);
					int b = 4 * e.sector[row] + e.layer[row];
					assertTrue(m + " b " + b, b >= 0 && b < 60);
				} else {
					nWedCand++;
				}
			}
			// bank times: the weight recomputed from c_time_bank with the _tpdf formula
			AtofTerm.attach(t, e, par(h), recal(h, 0), recal(h, 1), true, prop, ws, sb);
			assertEquals(e.cRow.length, sb.ncand);
			boolean reached = (((int) e.atofSeed[15]) & 1) != 0;
			double tp = par(h)[AtofTerm.P_TMODE] > 0.5 ? e.atofSeed[11] : e.atofSeed[10];
			for (int c = 0; c < sb.ncand; c++) {
				double w = 1.0;
				if (e.cBar[c] == 0 && reached && Double.isFinite(e.cTimeBank[c]))
					w = tpdfValue(e.cTimeBank[c] - tp, par(h)) * 2.0 * par(h)[AtofTerm.P_TWIN];
				wBank = Math.max(wBank, rel(sb.ct[c], w));
				assertTrue("track " + k + " cand " + c + " bank weight java " + sb.ct[c] + " test " + w, rel(sb.ct[c], w) <= 1e-9);
			}
		}
		System.out.printf("attach worst: seed %.2e cz %.2e cphi %.2e weight %.2e bank weight %.2e; bar cands %d wedge %d%n",
			wSeed, wCz, wPhi, wCt, wBank, nBarCand, nWedCand);
	}

	private static void fill(AtofTerm.TrackState st, GFitFixtures.Fx ca, int tr, double[] par) {
		int[] coff = ca.i("ca_coff"), cbar = ca.i("ca_cbar");
		double[] cz = ca.d("ca_cz"), cphi = ca.d("ca_cphi"), ct = ca.d("ca_ct");
		int c0 = coff[tr], c1 = coff[tr + 1];
		st.grow(c1 - c0);
		for (int c = c0; c < c1; c++) {
			st.cz[c - c0] = cz[c]; st.cphi[c - c0] = cphi[c]; st.ct[c - c0] = ct[c]; st.cbar[c - c0] = cbar[c];
		}
		st.ncand = c1 - c0;
		st.par = par;
	}

	private static double norm(double[] a, int n) {
		double s = 0;
		for (int i = 0; i < n; i++) s += a[i] * a[i];
		return Math.sqrt(s);
	}

	@Test
	public void atofTermMatchesPython() {
		GFitFixtures.Fx cs = GFitFixtures.resource("unit/cost_states.txt.gz");
		GFitFixtures.Fx ca = GFitFixtures.resource("unit/cost_atof.txt.gz");
		int ns = cs.shape("cs_u")[0];
		int[] trk = cs.i("cs_track"), sub = cs.i("cs_sub"), half = ca.i("ca_half");
		double[] fa = cs.d("cs_Fatof"), ga = cs.d("cs_gatof"), ia0 = cs.d("cs_ia0"), iafd = cs.d("cs_ia_fd");
		double[] hcT = cs.d("cs_Hatof_c_tri"), heT = cs.d("cs_Hatof_e_tri"), lay = cs.d("cs_layer");
		AtofTerm.TrackState st = new AtofTerm.TrackState();
		double[] g = new double[15], hc = new double[225], he = new double[225], ia = new double[16 * 20];
		double wF = 0, wG = 0, wHc = 0, wHe = 0, wLay = 0;
		int[] subIdx = new int[ns];
		java.util.Arrays.fill(subIdx, -1);
		for (int k = 0; k < sub.length; k++) subIdx[sub[k]] = k;
		for (int s = 0; s < ns; s++) {
			fill(st, ca, trk[s], par(half[trk[s]]));
			double[] u = cs.row("cs_u", s);
			int k = subIdx[s];
			java.util.Arrays.fill(g, 0); java.util.Arrays.fill(hc, 0); java.util.Arrays.fill(he, 0);
			boolean want = k >= 0;
			if (want) System.arraycopy(iafd, k * 320, ia, 0, 320);
			else {
				java.util.Arrays.fill(ia, 0);
				System.arraycopy(ia0, s * 20, ia, 0, 20);
			}
			double f = AtofTerm.termCost(u, ia, want, 15, 0.0, st, g, hc, he);
			double r = Math.abs(f - fa[s]) / Math.max(Math.abs(fa[s]), 1e-9);
			wF = Math.max(wF, r);
			assertTrue("state " + s + " Fatof java " + f + " python " + fa[s], r <= 1e-9);
			if (!want) continue;
			double[] gE = new double[15];
			System.arraycopy(ga, s * 15, gE, 0, 15);
			double d = 0;
			for (int j = 0; j < 15; j++) d += (g[j] - gE[j]) * (g[j] - gE[j]);
			d = Math.sqrt(d);
			double nrm = Math.max(norm(gE, 15), 1e-300);
			wG = Math.max(wG, d / nrm);
			assertTrue("state " + s + " gradient err " + d + " norm " + nrm, d <= 1e-6 * nrm);
			double[] w = new double[2];
			for (int q = 0; q < 2; q++) {
				double[] tri = q == 0 ? hcT : heT, h = q == 0 ? hc : he;
				double dd = 0, nn = 0;
				int idx = k * 120;
				for (int j = 0; j < 15; j++) {
					for (int m = j; m < 15; m++) {
						double e = tri[idx++];
						dd += (h[j * 15 + m] - e) * (h[j * 15 + m] - e);
						nn += e * e;
					}
				}
				dd = Math.sqrt(dd); nn = Math.sqrt(nn);
				w[q] = dd / Math.max(nn, 1e-300);
				assertTrue("state " + s + (q == 0 ? " Hc" : " He") + " err " + dd + " norm " + nn, dd <= 1e-6 * nn);
			}
			wHc = Math.max(wHc, w[0]); wHe = Math.max(wHe, w[1]);
			// per layer evaluation at the stored intercept variables
			for (int L = 0; L < 2; L++) {
				int o = (k * 2 + L) * 32;
				double[] x = {lay[o], lay[o + 1], lay[o + 2], lay[o + 3]};
				double c = AtofTerm.layerEval(L, x, st, true);
				String m = "state " + s + " layer " + L;
				assertTrue(m + " C java " + c + " python " + lay[o + 4], tol(c, lay[o + 4]));
				for (int q = 0; q < 4; q++) assertTrue(m + " Cx" + q, tol(st.Cx[q], lay[o + 5 + q]));
				for (int q = 0; q < 16; q++) assertTrue(m + " Cxx" + q + " java " + st.Cxx[q] + " python " + lay[o + 9 + q], tol(st.Cxx[q], lay[o + 9 + q]));
				for (int q = 0; q < 7; q++) assertTrue(m + " info" + q + " java " + st.info[q] + " python " + lay[o + 25 + q], tol(st.info[q], lay[o + 25 + q]));
				wLay = Math.max(wLay, relTol(c, lay[o + 4]));
			}
		}
		System.out.printf("atof term worst: Fatof rel %.2e, gradient %.2e of norm, Hc %.2e, He %.2e of Frobenius norm, layer %.2e rel%n",
			wF, wG, wHc, wHe, wLay);
	}

	private static double relTol(double a, double e) { return Math.abs(a - e) / Math.max(Math.abs(e), 1e-12); }

	private static boolean tol(double a, double e) { return relTol(a, e) <= 1e-9; }
}
