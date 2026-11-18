package org.jlab.rec.ahdc.GlobalFit;

/**
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * <p>Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this software is freely granted, provided that this notice is preserved.
 *
 * <p>Allocation-free complementary error function. Java port of erfc from the Sun fdlibm file s_erf.c
 * (the algorithm and coefficients of glibc libm erfc, which Python math.erfc calls): piecewise rational
 * approximations on |x| &lt; 0.84375, &lt; 1.25, &lt; 1/0.35 and &lt; 28, with
 * exp(-z*z-0.5625)*exp((z-x)(z+x)+R/S) where z is x with the low 32 bits cleared. The polynomials are grouped as in glibc
 * sysdeps/ieee754/dbl-64/s_erf.c (not fdlibm's Horner form), which reproduces glibc to the bit in the rational branches. Replaces commons-math3
 * Erf.erfc, whose regularizedGammaQ path allocates.
 */
final class GFitMath {

	private static final double TINY = 1e-300;
	private static final double ERX = 8.45062911510467529297e-01;
	// |x| < 0.84375
	private static final double PP0 = 1.28379167095512558561e-01;
	private static final double PP1 = -3.25042107247001499370e-01;
	private static final double PP2 = -2.84817495755985104766e-02;
	private static final double PP3 = -5.77027029648944159157e-03;
	private static final double PP4 = -2.37630166566501626084e-05;
	private static final double QQ1 = 3.97917223959155352819e-01;
	private static final double QQ2 = 6.50222499887672944485e-02;
	private static final double QQ3 = 5.08130628187576562776e-03;
	private static final double QQ4 = 1.32494738004321644526e-04;
	private static final double QQ5 = -3.96022827877536812320e-06;
	// 0.84375 <= |x| < 1.25
	private static final double PA0 = -2.36211856075265944077e-03;
	private static final double PA1 = 4.14856118683748331666e-01;
	private static final double PA2 = -3.72207876035701323847e-01;
	private static final double PA3 = 3.18346619901161753674e-01;
	private static final double PA4 = -1.10894694282396677476e-01;
	private static final double PA5 = 3.54783043256182359371e-02;
	private static final double PA6 = -2.16637559486879084300e-03;
	private static final double QA1 = 1.06420880400844228286e-01;
	private static final double QA2 = 5.40397917702171048937e-01;
	private static final double QA3 = 7.18286544141962662868e-02;
	private static final double QA4 = 1.26171219808761642112e-01;
	private static final double QA5 = 1.36370839120290507362e-02;
	private static final double QA6 = 1.19844998467991074170e-02;
	// 1.25 <= |x| < 1/0.35
	private static final double RA0 = -9.86494403484714822705e-03;
	private static final double RA1 = -6.93858572707181764372e-01;
	private static final double RA2 = -1.05586262253232909814e+01;
	private static final double RA3 = -6.23753324503260060396e+01;
	private static final double RA4 = -1.62396669462573470355e+02;
	private static final double RA5 = -1.84605092906711035994e+02;
	private static final double RA6 = -8.12874355063065934246e+01;
	private static final double RA7 = -9.81432934416914548592e+00;
	private static final double SA1 = 1.96512716674392571292e+01;
	private static final double SA2 = 1.37657754143519042600e+02;
	private static final double SA3 = 4.34565877475229228821e+02;
	private static final double SA4 = 6.45387271733267880336e+02;
	private static final double SA5 = 4.29008140027567833386e+02;
	private static final double SA6 = 1.08635005541779435134e+02;
	private static final double SA7 = 6.57024977031928170135e+00;
	private static final double SA8 = -6.04244152148580987438e-02;
	// 1/0.35 <= |x| < 28
	private static final double RB0 = -9.86494292470009928597e-03;
	private static final double RB1 = -7.99283237680523006574e-01;
	private static final double RB2 = -1.77579549177547519889e+01;
	private static final double RB3 = -1.60636384855821916062e+02;
	private static final double RB4 = -6.37566443368389627722e+02;
	private static final double RB5 = -1.02509513161107724954e+03;
	private static final double RB6 = -4.83519191608651397019e+02;
	private static final double SB1 = 3.03380607434824582924e+01;
	private static final double SB2 = 3.25792512996573918826e+02;
	private static final double SB3 = 1.53672958608443695994e+03;
	private static final double SB4 = 3.19985821950859553908e+03;
	private static final double SB5 = 2.55305040643316442583e+03;
	private static final double SB6 = 4.74528541206955367215e+02;
	private static final double SB7 = -2.24409524465858183362e+01;

	private GFitMath() {
	}

	/** Complementary error function erfc(x) = 1 - erf(x), accurate to about 1 ulp (fdlibm s_erf.c). */
	static double erfc(double x) {
		final int hx = (int) (Double.doubleToRawLongBits(x) >> 32);
		final int ix = hx & 0x7fffffff;
		if (ix >= 0x7ff00000) { // nan, +-inf
			return (hx >>> 31) * 2 + 1.0 / x;
		}
		if (ix < 0x3feb0000) { // |x| < 0.84375
			if (ix < 0x3c700000) {
				return 1.0 - x; // |x| < 2^-56
			}
			final double z = x * x;
			final double z2 = z * z;
			final double z4 = z2 * z2;
			final double r = (PP0 + z * PP1) + z2 * (PP2 + z * PP3) + z4 * PP4;
			final double s = (1.0 + z * QQ1) + z2 * (QQ2 + z * QQ3) + z4 * (QQ4 + z * QQ5);
			final double y = r / s;
			if (hx < 0x3fd00000) { // x < 1/4
				return 1.0 - (x + x * y);
			}
			return 0.5 - (x * y + (x - 0.5));
		}
		if (ix < 0x3ff40000) { // 0.84375 <= |x| < 1.25
			final double s = Math.abs(x) - 1.0;
			final double s2 = s * s;
			final double s4 = s2 * s2;
			final double s6 = s4 * s2;
			final double p = (PA0 + s * PA1) + s2 * (PA2 + s * PA3) + s4 * (PA4 + s * PA5) + s6 * PA6;
			final double q = (1.0 + s * QA1) + s2 * (QA2 + s * QA3) + s4 * (QA4 + s * QA5) + s6 * QA6;
			if (hx >= 0) {
				return 1.0 - ERX - p / q;
			}
			return 1.0 + (ERX + p / q);
		}
		if (ix < 0x403c0000) { // |x| < 28
			final double ax = Math.abs(x);
			final double s = 1.0 / (ax * ax);
			final double s2 = s * s;
			final double s4 = s2 * s2;
			final double s6 = s4 * s2;
			final double r;
			final double ss;
			if (ix < 0x4006db6d) { // |x| < 1/0.35
				r = (RA0 + s * RA1) + s2 * (RA2 + s * RA3) + s4 * (RA4 + s * RA5) + s6 * (RA6 + s * RA7);
				ss = (1.0 + s * SA1) + s2 * (SA2 + s * SA3) + s4 * (SA4 + s * SA5) + s6 * (SA6 + s * SA7)
						+ (s4 * s4) * SA8;
			} else {
				if (hx < 0 && ix >= 0x40180000) {
					return 2.0 - TINY; // x < -6
				}
				r = (RB0 + s * RB1) + s2 * (RB2 + s * RB3) + s4 * (RB4 + s * RB5) + s6 * RB6;
				ss = (1.0 + s * SB1) + s2 * (SB2 + s * SB3) + s4 * (SB4 + s * SB5) + s6 * (SB6 + s * SB7);
			}
			final double z = Double.longBitsToDouble(Double.doubleToRawLongBits(ax) & 0xffffffff00000000L);
			final double e = Math.exp(-z * z - 0.5625) * Math.exp((z - ax) * (z + ax) + r / ss);
			if (hx > 0) {
				return e / ax;
			}
			return 2.0 - e / ax;
		}
		return hx > 0 ? TINY * TINY : 2.0 - TINY;
	}
}
