package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * The ATOF time term (plan P4 S2): the ge10_t configuration of {@link AtofTerm} (cost, gradient, Hessians) against the pin's
 * atof.term_cost on the states of switch/cost_time, the time fit and the time fit plus recovery on t60 against ge10_t and ge10_trec
 * (parent plan D5 tiers 1-3), and a missing start time. Flag 256 is carried in GFitResult.flags; no comparison here or in the
 * harness reads that bit (tiers use bits 1 and 8), so Python's flags, which have no 256, need no masking.
 */
public class GFitTimeTest {

	static final int RUN = GFitFixtureConstants.RUN_PER_HALF;

	private static GFitFixtureConstants consts;
	private static AtofPositionTable[] timeTab;

	@BeforeClass
	public static void load() {
		GFitFullSearchTest.load();
		consts = GFitFullSearchTest.consts;
		timeTab = new AtofPositionTable[2];
		for (int h = 0; h < 2; h++) timeTab[h] = consts.atofPosition(RUN, h).withTime(consts.timeKeys(RUN, h));
	}

	private static double rel(double a, double e) { return Math.abs(a - e) / Math.max(Math.abs(e), 1e-12); }

	private static double norm(double[] a, int n) {
		double s = 0;
		for (int i = 0; i < n; i++) s += a[i] * a[i];
		return Math.sqrt(s);
	}

	// ================================================================== unit states

	@Test
	public void timeTermMatchesPython() {
		GFitFixtures.Fx cs = GFitFixtures.resource("unit/cost_states.txt.gz");
		GFitFixtures.Fx ct = GFitFixtures.resource("switch/cost_time.txt.gz");
		int ns = ct.shape("cs_orig")[0];
		assertTrue("states in cost_time", ns > 0);
		int[] orig = ct.i("cs_orig"), subk = ct.i("cs_subk"), trk = ct.i("cs_track"), coff = ct.i("ca_coff"), cbar = ct.i("ca_cbar"), half = ct.i("ca_half");
		double[] fa = ct.d("cs_Fatof"), fTot = ct.d("cs_F"), ga = ct.d("cs_gatof"), hcT = ct.d("cs_Hatof_c_tri"), heT = ct.d("cs_Hatof_e_tri"), lay = ct.d("cs_layer");
		double[] cz = ct.d("ca_cz"), cphi = ct.d("ca_cphi"), cct = ct.d("ca_ct");
		double[] fhit = cs.d("cs_Fhit"), fz = cs.d("cs_Fz"), fd0 = cs.d("cs_Fd0"), fk = cs.d("cs_Fk"), iafd = cs.d("cs_ia_fd"), ia0 = cs.d("cs_ia0");
		// withTime reproduces the pin's ge10_t parameter vector of each half slot for slot
		for (int h = 0; h < 2; h++) {
			double[] py = ct.row("ct_par", h), jv = timeTab[h].par();
			for (int q = 0; q < AtofTerm.NP; q++) assertEquals("ct_par half " + h + " slot " + q, py[q], jv[q], 0.0);
		}
		AtofTerm.TrackState st = new AtofTerm.TrackState();
		double[] g = new double[15], hc = new double[225], he = new double[225], ia = new double[16 * 20];
		double wF = 0, wG = 0, wHc = 0, wHe = 0, wLay = 0, wTot = 0;
		for (int s = 0; s < ns; s++) {
			int tr = trk[s], c0 = coff[tr], c1 = coff[tr + 1];
			st.grow(c1 - c0);
			for (int c = c0; c < c1; c++) {
				st.cz[c - c0] = cz[c]; st.cphi[c - c0] = cphi[c]; st.ct[c - c0] = cct[c]; st.cbar[c - c0] = cbar[c];
			}
			st.ncand = c1 - c0;
			st.par = timeTab[half[tr]].par();
			double[] u = cs.row("cs_u", orig[s]);
			Arrays.fill(g, 0); Arrays.fill(hc, 0); Arrays.fill(he, 0);
			System.arraycopy(iafd, subk[s] * 320, ia, 0, 320);
			double f = AtofTerm.termCost(u, ia, true, 15, 0.0, st, g, hc, he);
			double r = rel(f, fa[s]);
			wF = Math.max(wF, r);
			assertTrue("state " + s + " Fatof java " + f + " python " + fa[s], r <= 1e-9);
			double tot = fhit[orig[s]] + fz[orig[s]] + fd0[orig[s]] + fk[orig[s]] + f;
			wTot = Math.max(wTot, rel(tot, fTot[s]));
			assertTrue("state " + s + " total cost java " + tot + " python " + fTot[s], rel(tot, fTot[s]) <= 1e-9);
			double d = 0;
			for (int j = 0; j < 15; j++) d += (g[j] - ga[s * 15 + j]) * (g[j] - ga[s * 15 + j]);
			d = Math.sqrt(d);
			double nrm = Math.max(norm(Arrays.copyOfRange(ga, s * 15, s * 15 + 15), 15), 1e-300);
			wG = Math.max(wG, d / nrm);
			assertTrue("state " + s + " gradient err " + d + " norm " + nrm, d <= 1e-6 * nrm);
			for (int q = 0; q < 2; q++) {
				double[] tri = q == 0 ? hcT : heT, h = q == 0 ? hc : he;
				double dd = 0, nn = 0;
				int idx = s * 120;
				for (int j = 0; j < 15; j++)
					for (int m = j; m < 15; m++) {
						double e = tri[idx++];
						dd += (h[j * 15 + m] - e) * (h[j * 15 + m] - e);
						nn += e * e;
					}
				dd = Math.sqrt(dd); nn = Math.sqrt(nn);
				double w = dd / Math.max(nn, 1e-300);
				if (q == 0) wHc = Math.max(wHc, w); else wHe = Math.max(wHe, w);
				assertTrue("state " + s + (q == 0 ? " Hc" : " He") + " err " + dd + " norm " + nn, dd <= 1e-6 * nn);
			}
			for (int L = 0; L < 2; L++) {
				int o = (s * 2 + L) * 32;
				double[] x = {lay[o], lay[o + 1], lay[o + 2], lay[o + 3]};
				double c = AtofTerm.layerEval(L, x, st, true);
				String m = "state " + s + " layer " + L;
				assertTrue(m + " C java " + c + " python " + lay[o + 4], rel(c, lay[o + 4]) <= 1e-9);
				for (int q = 0; q < 4; q++) assertTrue(m + " Cx" + q, rel(st.Cx[q], lay[o + 5 + q]) <= 1e-9);
				for (int q = 0; q < 16; q++) assertTrue(m + " Cxx" + q, rel(st.Cxx[q], lay[o + 9 + q]) <= 1e-9);
				for (int q = 0; q < 7; q++) assertTrue(m + " info" + q, rel(st.info[q], lay[o + 25 + q]) <= 1e-9);
				wLay = Math.max(wLay, rel(c, lay[o + 4]));
			}
		}
		System.out.printf("time term worst over %d states: Fatof rel %.2e, total cost rel %.2e, gradient %.2e of norm, Hc %.2e, He %.2e of Frobenius norm, layer %.2e rel%n",
			ns, wF, wTot, wG, wHc, wHe, wLay);
	}

	// ================================================================== time fit and recovery on t60

	/** a fresh t60 set with the corrected times (a_t_ge10) and energies of switch/t60/atof_extra, and the ge10 expectations in before/after. */
	private static GFitReplaySet timeSet() {
		GFitReplaySet s = GFitReplaySet.resource("t60");
		GFitFixtures.Fx x = GFitFixtures.resource("switch/t60/atof_extra.txt.gz");
		double[] en = x.d("a_energy"), tm = x.d("a_t_ge10");
		int[] aoff = s.atofFx.i("aoff");
		assertEquals(aoff[s.n], tm.length);
		for (int k = 0; k < s.n; k++) {
			s.events[k].aEnergy = Arrays.copyOfRange(en, aoff[k], aoff[k + 1]);
			s.events[k].aTime = Arrays.copyOfRange(tm, aoff[k], aoff[k + 1]);
		}
		s.before = GFitFixtures.resource("switch/t60/expected_ge10_t.txt.gz");
		s.after = GFitFixtures.resource("switch/t60/expected_ge10_trec.txt.gz");
		return s;
	}

	/** parent D5 tiers 1-3 of the time fit r against s.before (GFitFullSearchTest.checkTiers on a set of its own); asserts 99.5 / 98 / 90 %. */
	static void tiersBefore(String label, GFitReplaySet s, GFitResult[] r) {
		GFitFixtures.Fx ex = s.before;
		int[] st = ex.i("est_status"), nmin = ex.i("est_nmin"), fl = ex.i("est_flags"), eoff = ex.i("eh_off"), ewire = ex.i("eh_wire"), ekind = ex.i("eh_kind");
		double[] cost = ex.d("est_cost"), state = ex.d("est_state"), cov = ex.d("est_cov");
		List<String> out1 = new ArrayList<>(), out2 = new ArrayList<>(), out3 = new ArrayList<>();
		int n1 = 0, ok1 = 0, n2 = 0, ok2 = 0, n3 = 0, ok3 = 0;
		double wCost = 0, wPar = 0, wCov = 0, wCost3 = Double.NEGATIVE_INFINITY;
		double[] cx = new double[25];
		for (int k = 0; k < s.n; k++) {
			GFitResult j = r[k];
			String id = "trackid " + s.tracks[k].trackid + " row " + k;
			n1++;
			boolean hits = j.nh == eoff[k + 1] - eoff[k];
			for (int h = 0; hits && h < j.nh; h++) hits = j.hitWire[h] == ewire[eoff[k] + h] && j.hitKind[h] == ekind[eoff[k] + h];
			if (j.status == st[k] && hits) ok1++;
			else out1.add(id + " status java " + j.status + " python " + st[k] + " hits " + hits);
			boolean t3 = nmin[k] > 1 || (fl[k] & 1) != 0 || (fl[k] & 8) != 0;
			boolean t2 = st[k] == 0 && nmin[k] == 1 && (fl[k] & 8) == 0 && !t3;
			if (t2) {
				n2++;
				double dc = Math.abs(j.cost - cost[k]);
				wCost = Math.max(wCost, dc);
				boolean ok = j.status == 0 && dc <= 1e-3;
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
		for (String m : out1) System.out.println("[" + label + "] outside tier 1: " + m);
		for (String m : out2) System.out.println("[" + label + "] outside tier 2: " + m);
		for (String m : out3) System.out.println("[" + label + "] outside tier 3: " + m);
		assertTrue("tier 1 " + out1, ok1 >= 0.995 * n1);
		assertTrue("tier 2 " + out2, ok2 >= 0.98 * n2);
		assertTrue("tier 3 " + out3, ok3 >= 0.90 * n3);
	}

	@Test
	public void timeTiersOnSubset() {
		GFitReplaySet s = timeSet();
		int maxHits = GFitFullSearchTest.maxHits;
		GFitConfig cfg = GFitConfig.pythonDefault().withSwitches(true, false);
		GlobalTrackFitter ft = GFitFullSearchTest.fitter(cfg);
		HitRecovery hr = new HitRecovery(ft, new TrackCost(GFitFullSearchTest.prop, GFitFullSearchTest.geo, cfg), GFitFullSearchTest.geo, cfg);
		GFitWorkspace ws = new GFitWorkspace(maxHits);
		GFitRecoveryTest.Run run = new GFitRecoveryTest.Run();
		GFitResult[] base = new GFitResult[s.n];
		run.base = base;
		run.after = new GFitResult[s.n];
		run.accEv = new int[s.n][]; run.accLayer = new int[s.n][]; run.accRound = new int[s.n][];
		run.nAdded = new int[s.n]; run.retries = new int[s.n]; run.dropEv0 = new int[s.n];
		int live = 0;
		for (int k = 0; k < s.n; k++) {
			base[k] = new GFitResult(maxHits); run.after[k] = new GFitResult(maxHits);
			AtofPositionTable tab = timeTab[s.events[k].half];
			ft.fit(s.tracks[k], s.events[k], tab, false, ws, base[k]);
			if ((base[k].flags & 256) != 0) live++;
			run.nAdded[k] = hr.recover(s.tracks[k], s.events[k], tab, false, base[k], ws, run.after[k]);
			HitRecovery.Scratch sc = ws.rec;
			run.accEv[k] = Arrays.copyOf(sc.accEv, sc.nAcc); run.accLayer[k] = Arrays.copyOf(sc.accLayer, sc.nAcc);
			run.accRound[k] = Arrays.copyOf(sc.accRound, sc.nAcc);
			run.retries[k] = sc.retries; run.dropEv0[k] = sc.dropEv[0];
		}
		System.out.println("time fit: flag 256 on " + live + " of " + s.n + " tracks");
		assertTrue("flag 256 set on tracks with a finite wedge time", live > 0);
		tiersBefore("time t60 before recovery", s, base);
		int bad = GFitRecoveryTest.accMismatches("time t60", s, run);
		int[] t = GFitRecoveryTest.tiersAfter("time t60 after recovery", s, run);
		int rows = 0, tracks = 0;
		for (int k = 0; k < s.n; k++) { rows += run.nAdded[k]; if (run.nAdded[k] > 0) tracks++; }
		System.out.printf("time recover: %d rows added on %d of %d tracks; %d tracks with accepted rows unlike Python%n", rows, tracks, s.n, bad);
		assertTrue("accepted rows differ on " + bad + " tracks", bad <= 0.005 * s.n);
		assertTrue("tier 1 after", t[0] >= 0.995 * t[1]);
		assertTrue("tier 2 after", t[2] >= 0.98 * t[3]);
		assertTrue("tier 3 after", t[4] >= 0.90 * t[5]);
	}

	/** a NaN start time makes every ct NaN (wedge factor 1 / (2 t_win)): the fit stays finite and flag 256 is clear. */
	@Test
	public void timeNoStartTime() {
		GFitReplaySet s = timeSet();
		GFitConfig cfg = GFitConfig.pythonDefault().withSwitches(true, false);
		GlobalTrackFitter ft = GFitFullSearchTest.fitter(cfg);
		GFitWorkspace ws = new GFitWorkspace(GFitFullSearchTest.maxHits);
		GFitResult r = new GFitResult(GFitFullSearchTest.maxHits);
		int fin = 0, n = Math.min(s.n, 12);
		for (int k = 0; k < n; k++) {
			GFitEvent e = s.events[k];
			e.startTime = Double.NaN;
			ft.fit(s.tracks[k], e, timeTab[e.half], false, ws, r);
			assertTrue("track " + k + " status " + r.status, r.status != 4);
			assertEquals("track " + k + " flag 256", 0, r.flags & 256);
			if (Double.isFinite(r.cost) && r.status == 0) fin++;
		}
		System.out.println("no start time: " + fin + " of " + n + " fits converged with finite cost");
		assertTrue(fin >= n - 1);
	}
}
