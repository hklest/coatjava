package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

/** GFitMath.erfc against unit/erfc_ref.txt.gz, Python math.erfc (glibc libm) at dense points, region boundaries and special values. */
public class GFitMathTest {

	private static long ulps(double a, double b) {
		if (Double.isNaN(a) || Double.isNaN(b)) {
			return Double.isNaN(a) && Double.isNaN(b) ? 0 : Long.MAX_VALUE;
		}
		if (a == b) {
			return 0;
		}
		if (Double.isInfinite(a) || Double.isInfinite(b)) {
			return Long.MAX_VALUE;
		}
		return Math.abs(Double.doubleToLongBits(a) - Double.doubleToLongBits(b));
	}

	@Test
	public void erfcMatchesPython() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/erfc_ref.txt.gz");
		double[] x = fx.d("x"), y = fx.d("erfc");
		long worst = 0;
		int exact = 0;
		double wx = 0;
		for (int i = 0; i < x.length; i++) {
			long u = ulps(GFitMath.erfc(x[i]), y[i]);
			if (u == 0) {
				exact++;
			}
			if (u > worst) {
				worst = u;
				wx = x[i];
			}
		}
		System.out.printf("erfc: %d points, %d bit-identical, worst %d ulp (x = %s, limit 2)%n", x.length, exact, worst, wx);
		assertTrue("worst " + worst + " ulp at x = " + wx, worst <= 2);
	}
}
