package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Assume;
import org.junit.BeforeClass;
import org.junit.Test;

/** FieldGrid and Propagator against the fixtures unit/field.txt and unit/prop_states.txt of the pinned prop.py / fit.py. */
public class GFitPropagatorTest {

	private static FieldGrid field;
	private static Propagator prop;
	private static GFitFixtures.Fx ps;
	private static double[] beam;

	@BeforeClass
	public static void load() {
		GFitFixtures.Fx f = GFitFixtures.resource("unit/field.txt.gz");
		double[] fp = f.d("field_par");
		int[] sh = f.shape("br");
		field = FieldGrid.fromArrays(fp[0], fp[1], fp[2], fp[3], sh[0], sh[1], f.d("br"), f.d("bz"));
		beam = GFitFixtures.resource("unit/stopping.txt.gz").d("beam");
		prop = new Propagator(field, MaterialStack.dataBest(), StoppingTables.proton(), beam[0], beam[1]);
		ps = GFitFixtures.resource("unit/prop_states.txt.gz");
	}

	/**
	 * Expected values: computed once by the pinned Python (prop._bfield.py_func on prop.load_field() of the pin, the
	 * grid of unit/field.txt) at these points and hard-coded; the points: x, y, z in mm.
	 */
	private static final double[][] EDGE = {
		// rho > 120
		{150.0, 0.0, 0.0}, {0.0, 200.0, 10.0}, {100.0, 100.0, 50.0}, {-130.0, 5.0, -100.0}, {121.0, 0.0, 399.5}, {90.0, -90.0, -370.0},
		// z beyond +-400 after the shift (z_map = z + 30)
		{10.0, 5.0, 450.0}, {20.0, -30.0, 371.0}, {50.0, 50.0, 1000.0}, {-60.0, 20.0, -500.0},
		// below zmin
		{5.0, 5.0, -431.0}, {0.0, 70.0, -440.0}, {-100.0, 0.0, -2000.0},
		// rho = 0 (and tiny)
		{0.0, 0.0, 0.0}, {0.0, 0.0, 500.0}, {0.0, 0.0, -500.0}, {0.0, 0.0, 100.0}, {1e-10, 0.0, 30.0}, {0.0, -1e-12, -450.0},
		// both
		{300.0, 300.0, 900.0}};
	private static final double[][] EDGE_B = {
		{-0.0037329068804085254, -0.0, -5.023385939666749},
		{-0.0, -0.004104512302428484, -5.023607725044251},
		{-0.007022751787187457, -0.007022751787187457, -5.025035731758118},
		{-0.0003116975669683011, 1.1988367960319274e-05, -5.0339986955719},
		{-0.3673007890534401, -0.0, -4.18789207698822},
		{0.20520718838929228, -0.20520718838929228, -4.558354097709656},
		{-0.029667680533108767, -0.014833840266554384, -4.125344188578955},
		{-0.059483280203295665, 0.0892249203049435, -4.130337021354346},
		{-0.1498928774076938, -0.1498928774076938, -4.146238533990424},
		{-0.17960383047339, 0.059867943491130005, -4.171553087943867},
		{0.014839752234101304, 0.014839752234101304, -4.154313938949709},
		{0.0, 0.20994267516279222, -4.17549929232788},
		{-0.30332944952726365, 0.0, -4.198324033760072},
		{0.0, 0.0, -5.0178367086410525},
		{0.0, 0.0, -4.1247950947113035},
		{0.0, 0.0, -4.154071147590638},
		{0.0, 0.0, -4.9812325451278685},
		{0.0, 0.0, -5.0125035172576915},
		{0.0, 0.0, -4.154071147590638},
		{-0.2597208786748571, -0.2597208786748571, -4.18789207698822}};

	@Test
	public void fieldClampsAtEdges() {
		assertEquals(20, EDGE.length);
		assertEquals(20, EDGE_B.length);
		double[] b = new double[3];
		double worst = 0;
		for (int i = 0; i < EDGE.length; i++) {
			field.field(EDGE[i][0], EDGE[i][1], EDGE[i][2], b);
			for (int c = 0; c < 3; c++) {
				double d = Math.abs(b[c] - EDGE_B[i][c]) / Math.max(1.0, Math.abs(EDGE_B[i][c]));
				worst = Math.max(worst, d);
				assertTrue("point " + i + " comp " + c + " java " + b[c] + " python " + EDGE_B[i][c], d <= 1e-15);
			}
		}
		System.out.printf("field edge clamp worst deviation %.3e (limit 1e-15)%n", worst);
	}

	private static double dev(double a, double e) {
		if (Double.isNaN(a) && Double.isNaN(e)) return 0.0;
		if (a == e) return 0.0;
		return Math.abs(a - e) / Math.max(1.0, Math.abs(e));
	}

	@Test
	public void propMatchesPython() {
		double[] ru = ps.d("ps_u");
		int n = ps.shape("ps_u")[0];
		int[] nnE = ps.i("ps_nn"), flE = ps.i("ps_flags"), regE = ps.i("ps_regnode");
		double[] auxE = ps.d("ps_aux"), th0E = ps.d("ps_theta0"), iaE = ps.d("ps_ia");
		int nregp1 = ps.shape("ps_regnode")[1];
		assertEquals(MaterialStack.dataBest().nreg + 1, nregp1);
		int[] pnState = ps.i("pn_state"), pnN = ps.i("pn_n"), knCnt = ps.i("pn_kn_cnt"), knIdx = ps.i("pn_kn_idx");
		double[] pnNodes = ps.d("pn_nodes");
		GFitWorkspace ws = new GFitWorkspace(0);
		double[] u = new double[15], rec = new double[Propagator.NI];
		double wTh = 0, wAux = 0, wIa = 0, wNode = 0;
		int off = 0, kOff = 0;
		MaterialStack stk = MaterialStack.dataBest();
		for (int r = 0; r < n; r++) {
			System.arraycopy(ru, r * 15, u, 0, 15);
			int nn = prop.run(u, ws);
			assertEquals("nn state " + r, nnE[r], nn);
			assertEquals("flags state " + r, flE[r], ws.flags);
			for (int k = 0; k < nregp1; k++) assertEquals("regnode state " + r + " k " + k, regE[r * nregp1 + k], ws.regnode[k]);
			for (int j = 0; j < 6; j++) wTh = Math.max(wTh, dev(ws.theta0[j], th0E[r * 6 + j]));
			for (int j = 0; j < 4; j++) wAux = Math.max(wAux, dev(ws.aux[j], auxE[r * 4 + j]));
			Propagator.atofIntercepts(ws.out, nn, ws.regnode, ws.aux, stk.rOut, rec, 0);
			for (int j = 0; j < Propagator.NI; j++) wIa = Math.max(wIa, dev(rec[j], iaE[r * Propagator.NI + j]));
			int m = -1;
			for (int q = 0; q < pnState.length; q++) if (pnState[q] == r) m = q;
			if (m >= 0) {
				assertEquals("pn_n state " + r, pnN[m], nn);
				for (int i = 0; i < nn * GFitWorkspace.NC; i++) wNode = Math.max(wNode, dev(ws.out[i], pnNodes[off + i]));
				int c = 0;
				for (int i = 1; i < nn; i++) {
					if (ws.out[i * GFitWorkspace.NC] == ws.out[(i - 1) * GFitWorkspace.NC]) {
						assertTrue("kink node count state " + r, c < knCnt[m]);
						assertEquals("kink node index state " + r, knIdx[kOff + c], i);
						c++;
					}
				}
				assertEquals("kink node count state " + r, knCnt[m], c);
				off += nn * GFitWorkspace.NC;
				kOff += knCnt[m];
			} else {
				// nodes are only stored for pn_state; states in pn_state are in increasing order
			}
		}
		System.out.printf("propagation worst deviation (relative, floor 1 for absolute): theta0 %.3e  aux %.3e  atof record %.3e  nodes %.3e (limit 1e-9)%n",
			wTh, wAux, wIa, wNode);
		assertTrue(wTh <= 1e-9);
		assertTrue(wAux <= 1e-9);
		assertTrue(wIa <= 1e-9);
		assertTrue(wNode <= 1e-9);
	}

	private static void assertBits(String what, double[] a, double[] b, int n) {
		for (int i = 0; i < n; i++)
			if (Double.doubleToLongBits(a[i]) != Double.doubleToLongBits(b[i])) throw new AssertionError(what + " differs at " + i + ": " + a[i] + " vs " + b[i]);
	}

	@Test
	public void restartEqualsFull() {
		double[] ru = ps.d("ps_u");
		int n = ps.shape("ps_u")[0];
		GFitWorkspace wsN = new GFitWorkspace(0), wsF = new GFitWorkspace(0);
		wsN.saveSnapshots = true;
		MaterialStack stk = MaterialStack.dataBest();
		double[] u = new double[15], up = new double[15], recF = new double[Propagator.NI], recR = new double[Propagator.NI];
		int unreached = 0, checked = 0;
		for (int r = 0; r < n; r++) {
			System.arraycopy(ru, r * 15, u, 0, 15);
			prop.run(u, wsN);
			for (int j = 0; j < 10; j++) {
				System.arraycopy(u, 0, up, 0, 15);
				up[5 + j] += 1e-5;
				int nnF = prop.run(up, wsF);
				int nnR = prop.runRestart(up, j / 2, wsN);
				String id = "state " + r + " kink " + j;
				assertEquals(id + " nn", nnF, nnR);
				assertEquals(id + " ws.nn", nnF, wsN.nn);
				assertEquals(id + " flags", wsF.flags, wsN.flags);
				assertBits(id + " nodes", wsF.out, wsN.out, nnF * GFitWorkspace.NC);
				for (int k = 0; k <= stk.nreg; k++) assertEquals(id + " regnode " + k, wsF.regnode[k], wsN.regnode[k]);
				assertBits(id + " aux", wsF.aux, wsN.aux, 4);
				assertBits(id + " theta0", wsF.theta0, wsN.theta0, 6);
				Propagator.atofIntercepts(wsF.out, nnF, wsF.regnode, wsF.aux, stk.rOut, recF, 0);
				Propagator.atofIntercepts(wsN.out, nnR, wsN.regnode, wsN.aux, stk.rOut, recR, 0);
				assertBits(id + " atof", recF, recR, Propagator.NI);
				if (!wsN.snap[j / 2].valid) unreached++;
				checked++;
			}
		}
		System.out.println("restartEqualsFull: " + checked + " restarts bit-identical (" + unreached + " on planes the nominal run did not reach)");
	}

	@Test
	public void timing() {
		Assume.assumeTrue("timing only with -Dgfit.timing=true", Boolean.getBoolean("gfit.timing"));
		double[] ru = ps.d("ps_u");
		int n = ps.shape("ps_u")[0];
		GFitWorkspace ws = new GFitWorkspace(0);
		double[] u = new double[15];
		long sink = 0;
		for (int rep = 0; rep < 300; rep++)
			for (int r = 0; r < n; r++) {
				System.arraycopy(ru, r * 15, u, 0, 15);
				sink += prop.run(u, ws);
			}
		int reps = 2000;
		long t0 = System.nanoTime();
		for (int rep = 0; rep < reps; rep++)
			for (int r = 0; r < n; r++) {
				System.arraycopy(ru, r * 15, u, 0, 15);
				sink += prop.run(u, ws);
			}
		double perRun = (System.nanoTime() - t0) * 1e-3 / ((double) reps * n);
		ws.saveSnapshots = true;
		long t1 = System.nanoTime();
		for (int rep = 0; rep < reps; rep++)
			for (int r = 0; r < n; r++) {
				System.arraycopy(ru, r * 15, u, 0, 15);
				sink += prop.run(u, ws);
			}
		double perSave = (System.nanoTime() - t1) * 1e-3 / ((double) reps * n);
		double[] up = new double[15];
		GFitWorkspace ws2 = new GFitWorkspace(0);
		long full = 0;
		double[] tp = new double[5];
		long[] cp = new long[5];
		for (int rep = 0; rep < 700; rep++)
			for (int r = 0; r < n; r++) {
				System.arraycopy(ru, r * 15, u, 0, 15);
				prop.run(u, ws);
				System.arraycopy(u, 0, up, 0, 15);
				for (int j = 0; j < 10; j++) {
					up[5 + j] += 1e-5;
					long a = System.nanoTime();
					sink += prop.runRestart(up, j / 2, ws);
					long d = System.nanoTime() - a;
					a = System.nanoTime();
					sink += prop.run(up, ws2);
					if (rep >= 300) full += System.nanoTime() - a;
					up[5 + j] = u[5 + j];
					if (rep >= 300) {
						tp[j / 2] += d;
						cp[j / 2]++;
					}
				}
			}
		System.out.printf("  full run of perturbed states: %.1f us%n", full * 1e-3 / (400.0 * n * 10));
		double t2 = 0;
		long cnt = 0;
		for (int q = 0; q < 5; q++) {
			System.out.printf("  runRestart plane %d: %.1f us%n", q, tp[q] * 1e-3 / cp[q]);
			t2 += tp[q];
			cnt += cp[q];
		}
		System.out.printf("timing (single thread, after warm-up, %d states): run %.1f us, run with snapshots %.1f us, runRestart %.1f us (sink %d)%n",
			n, perRun, perSave, t2 * 1e-3 / cnt, sink);
	}
}
