package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/** Cell polygons, wire table and helix model against the pinned Python (fit.wire_tables, cell._eval, fit._hhit). */
public class GFitGeometryTest {

	private static CellGeometry geometry() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/wires.txt.gz");
		return CellGeometry.fromWireEnds(fx.d("WO"), fx.d("WE"), fx.i("layoff"));
	}

	private static double rel(double a, double b) {
		return Math.abs(a - b) / Math.max(Math.abs(b), 1e-12);
	}

	@Test
	public void cellMatchesPython() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/wires.txt.gz");
		CellGeometry c = geometry();
		double[] wu = fx.d("WU"), lrad = fx.d("lrad"), half = fx.d("half"), cell = fx.d("cellarr_flat");
		assertEquals(wu.length, c.wu.length);
		assertEquals(cell.length, c.cellarr.length);
		double wWu = 0, wLrad = 0, wHalf = 0, wCell = 0;
		for (int k = 0; k < wu.length; k++) wWu = Math.max(wWu, Math.abs(wu[k] - c.wu[k]));
		for (int k = 0; k < 8; k++) wLrad = Math.max(wLrad, rel(c.lrad[k], lrad[k]));
		for (int k = 0; k < half.length; k++) wHalf = Math.max(wHalf, rel(c.half[k], half[k]));
		for (int k = 0; k < cell.length; k++) wCell = Math.max(wCell, Math.abs(cell[k] - c.cellarr[k]));
		System.out.printf("cellMatchesPython worst: wu %.3g lrad %.3g half %.3g cellarr %.3g mm%n", wWu, wLrad, wHalf, wCell);
		assertTrue(wWu < 1e-15);
		assertTrue(wLrad < 1e-12);
		assertTrue(wHalf < 1e-12);
		assertTrue(wCell < 1e-12);
		int[] lo = fx.i("layoff");
		for (int g = 0; g < c.nWires(); g++) {
			int li = c.layerOf(g);
			assertTrue(g >= lo[li] && g < lo[li + 1]);
		}
	}

	@Test
	public void cellEvalMatchesPython() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/cell_eval.txt.gz");
		CellGeometry c = geometry();
		// eval alone: use Python's own polygons (the 3e-14 mm rebuild difference flips the gradient at a vertex, where d = 0)
		double[] pyCell = GFitFixtures.resource("unit/wires.txt.gz").d("cellarr_flat");
		System.arraycopy(pyCell, 0, c.cellarr, 0, pyCell.length);
		int[] li = fx.i("ce_li"), w = fx.i("ce_w");
		double[] p = fx.d("ce_p"), r = fx.d("ce_r");
		double[] out = new double[4];
		double worst = 0;
		for (int k = 0; k < li.length; k++) {
			c.eval(li[k], w[k], p[3 * k], p[3 * k + 1], p[3 * k + 2], out);
			for (int j = 0; j < 4; j++) worst = Math.max(worst, Math.abs(out[j] - r[4 * k + j]));
		}
		System.out.printf("cellEvalMatchesPython worst abs %.3g over %d points%n", worst, li.length);
		assertTrue(worst < 1e-12);
	}

	@Test
	public void helixJacobian() {
		double bz = 5.0207, bx = -1.05, by = 0.0;
		java.util.Random rnd = new java.util.Random(7);
		double[] o = new double[9], jac = new double[15], xp = new double[5], xm = new double[5];
		double[] a = new double[9], b = new double[9];
		double worst = 0;
		for (int n = 0; n < 20; n++) {
			double p = 0.15 + 0.45 * rnd.nextDouble(), th = Math.toRadians(40 + 100 * rnd.nextDouble());
			double[] x = {4 * rnd.nextDouble() - 2, 2 * Math.PI * rnd.nextDouble() - Math.PI, Math.log(p), 200 * rnd.nextDouble() - 100, th};
			double s = 30 + 90 * rnd.nextDouble();
			HelixModel.dpdx(x, s, bz, jac);
			double[] h = {1e-4, 1e-6, 1e-6, 1e-4, 1e-6};
			for (int j = 0; j < 5; j++) {
				System.arraycopy(x, 0, xp, 0, 5);
				System.arraycopy(x, 0, xm, 0, 5);
				xp[j] += h[j];
				xm[j] -= h[j];
				HelixModel.pt(xp, s, bz, bx, by, a);
				HelixModel.pt(xm, s, bz, bx, by, b);
				for (int k = 0; k < 3; k++) {
					double fd = (a[k] - b[k]) / (2 * h[j]);
					double an = jac[3 * j + k];
					// relative to the row scale (an entry can be exactly 0)
					double scale = Math.max(Math.abs(an), 1.0);
					worst = Math.max(worst, Math.abs(fd - an) / scale);
				}
			}
			assertEquals(HelixModel.pt(x, s, bz, bx, by, o), 1000.0 * p * Math.sin(th) / (0.299792458 * bz), 1e-9);
		}
		System.out.printf("helixJacobian worst scaled deviation %.3g%n", worst);
		assertTrue(worst < 1e-6);
	}

	@Test
	public void helixMatchesPython() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/helix_hhit.txt.gz");
		CellGeometry c = geometry();
		double[] cf = fx.d("hh_cf"), xs = fx.d("hh_x"), s0 = fx.d("hh_s0"), dref = fx.d("hh_d"), sref = fx.d("hh_s"), jref = fx.d("hh_jx");
		int[] g = fx.i("hh_g");
		double[] o = new double[9], jx = new double[5], j53 = new double[15], sa = new double[1], x = new double[5];
		double wd = 0, ws = 0, wj = 0;
		for (int k = 0; k < g.length; k++) {
			System.arraycopy(xs, 5 * k, x, 0, 5);
			double d = HelixModel.hhit(x, c.wo, c.wu, g[k], cf[0], cf[1], cf[2], o, s0[k], true, jx, j53, sa);
			if (Double.isNaN(dref[k])) {
				assertTrue(Double.isNaN(d));
				continue;
			}
			wd = Math.max(wd, Math.abs(d - dref[k]) / Math.max(Math.abs(dref[k]), 1e-12));
			ws = Math.max(ws, Math.abs(sa[0] - sref[k]) / Math.max(Math.abs(sref[k]), 1e-12));
			for (int j = 0; j < 5; j++)
				wj = Math.max(wj, Math.abs(jx[j] - jref[5 * k + j]) / Math.max(Math.abs(jref[5 * k + j]), 1e-12));
		}
		System.out.printf("helixMatchesPython worst rel: d %.3g s %.3g Jx %.3g (%d cases)%n", wd, ws, wj, g.length);
		assertTrue(wd < 1e-12);
		assertTrue(ws < 1e-12);
		assertTrue(wj < 1e-12);
	}
}
