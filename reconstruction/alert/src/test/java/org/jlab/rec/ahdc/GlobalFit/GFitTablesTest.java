package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/** StoppingTables and MaterialStack against the fixture unit/stopping.txt of the pinned prop.py (prop.Tables, make_stack). */
public class GFitTablesTest {

	@Test
	public void stoppingMatchesPython() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/stopping.txt.gz");
		StoppingTables t = StoppingTables.proton();
		double[] par = fx.d("tab_par");
		System.out.printf("tab_par java lnT0 %s dl %s invd %s ; python %s %s %s%n", t.lnT0, t.dl, t.invd, par[0], par[1], par[3]);
		assertEquals(par[0], t.lnT0, 0.0);
		assertEquals(par[1], t.dl, 0.0);
		assertEquals(par[2], StoppingTables.NG, 0.0);
		assertEquals(par[3], t.invd, 0.0);
		int[] idx = fx.i("tab_idx");
		assertEquals(401, idx.length);
		double[] s0 = fx.d("tab_S0"), r0 = fx.d("tab_R0");
		double worstS = 0, worstR = 0;
		for (int m = 0; m < StoppingTables.NMAT; m++) {
			for (int j = 0; j < idx.length; j++) {
				double es = s0[m * idx.length + j], er = r0[m * idx.length + j];
				double ds = Math.abs(t.S[m * StoppingTables.NG + idx[j]] - es) / Math.abs(es);
				double dr = Math.abs(t.R[m * StoppingTables.NG + idx[j]] - er) / Math.abs(er);
				worstS = Math.max(worstS, ds);
				worstR = Math.max(worstR, dr);
			}
		}
		System.out.printf("stopping worst relative deviation: S %.3e  R %.3e (limit 1e-9)%n", worstS, worstR);
		assertTrue("S " + worstS, worstS <= 1e-9);
		assertTrue("R " + worstR, worstR <= 1e-9);
		double wx = 0, wz = 0;
		double[] x0 = fx.d("tab_x0g"), zoa = fx.d("tab_zoa");
		for (int m = 0; m < StoppingTables.NMAT; m++) {
			wx = Math.max(wx, Math.abs(t.x0g[m] - x0[m]) / x0[m]);
			wz = Math.max(wz, Math.abs(t.zoa[m] - zoa[m]) / zoa[m]);
		}
		System.out.printf("x0g worst relative %.3e, zoa %.3e%n", wx, wz);
		assertArrayEquals(x0, t.x0g, 0.0);
		assertArrayEquals(zoa, t.zoa, 0.0);
		double[] misc = fx.d("const_misc");
		assertEquals(misc[0], StoppingTables.RHO_WIRE_EQ, 0.0);
		assertEquals(misc[1], StoppingTables.W_HE_GAS, 0.0);
		assertEquals(misc[2], Propagator.KAPPA, 0.0);
		assertEquals(misc[3], Propagator.CMMNS, 0.0);
		assertEquals(misc[4], 0.51099895e-3, 0.0);
		double[] hyp = fx.d("hyp_p");
		assertEquals(hyp[0], Propagator.MASS, 0.0);
		assertEquals(hyp[1], Propagator.ZP, 0.0);
		assertEquals(hyp[3], Propagator.TSC, 0.0);
		assertEquals(hyp[4], Propagator.ZR2, 0.0);
		assertEquals(fx.d0("tfloor"), Propagator.T_FLOOR, 0.0);
	}

	@Test
	public void stackMatchesPython() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/stopping.txt.gz");
		MaterialStack k = MaterialStack.dataBest();
		assertEquals(15, k.nreg);
		assertEquals(6, k.npl);
		assertArrayEquals(fx.d("stk_r_out"), k.rOut, 0.0);
		assertArrayEquals(fx.d("stk_rho"), k.rho, 0.0);
		assertArrayEquals(fx.i("stk_mat"), k.mat);
		assertArrayEquals(fx.i("stk_nst"), k.nst);
		assertArrayEquals(fx.i("stk_seg"), k.seg);
		assertArrayEquals(fx.i("stk_pla"), k.pla);
		assertArrayEquals(fx.d("stk_xk"), k.xk, 0.0);
		assertArrayEquals(fx.d("stk_zoarho"), k.zoarho, 0.0);
		assertArrayEquals(fx.d("stk_plane_r"), k.planeR, 0.0);
	}
}
