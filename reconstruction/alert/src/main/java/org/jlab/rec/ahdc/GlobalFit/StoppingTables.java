package org.jlab.rec.ahdc.GlobalFit;

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;

/**
 * Stopping power and CSDA range of the nine propagation materials on a uniform ln T grid, proton family only (PSTAR).
 * Port of prop.Tables (pinned prop.py:126-170, 134-171) with the lookups prop._stop, _range_g and _t_from_range
 * (prop.py:267-302). Each base table is a natural cubic spline of ln S against ln T (commons-math3 SplineInterpolator,
 * the same boundary condition as scipy CubicSpline(bc_type="natural")) on the clipped grid; beyond the last table point
 * S falls as exp(-0.8 dlnT); mixtures add by mass fraction (Bragg); the range is the cumulative trapezoid of T/S dlnT plus
 * 2 T0/S0. Immutable and thread safe. Material ids: H2 0, He 1, CO2 2, KAPTON 3, AL 4, VT 5, HECO2 6, GAS 7, GASW 8.
 */
public final class StoppingTables {

	public static final int H2 = 0, HE = 1, CO2 = 2, KAPTON = 3, AL = 4, VT = 5, HECO2 = 6, GAS = 7, GASW = 8;
	public static final int NMAT = 9;
	public static final int NG = 4001;

	// prop.py:122-125: X0 g/cm2 and Z/A of the bases (order H2 He CO2 Kapton Al VT)
	private static final double[] BASE_X0 = {63.04, 94.32, 36.2, 40.56, 24.01, 43.79};
	private static final double[] BASE_ZOA = {0.99212, 0.49968, 0.49989, 0.51264, 0.48181, 0.54141};
	static final double W_HE_GAS = 0.8 * 1.664e-4 / (0.8 * 1.664e-4 + 0.2 * 1.84e-3);
	static final double RHO_WIRE_EQ = 3026 * Math.PI * (0.015 * 0.015) / (Math.PI * (69.0 * 69.0 - 32.0 * 32.0)) * 2.64;

	/** prop.py:126-131: per material, base indices and mass fractions in the order of the Python dicts. */
	private static final int[][] COMP_BASE = {
		{0}, {1}, {2}, {3}, {4}, {5},
		{1, 2},
		{1, 2},
		{1, 2, 4}};
	private static final double[][] COMP_W = {
		{1.0}, {1.0}, {1.0}, {1.0}, {1.0}, {1.0},
		{4.0026 / 48.0126, 44.0095 / 48.0126},
		{W_HE_GAS, 1 - W_HE_GAS},
		{W_HE_GAS * 5.015e-4 / (5.015e-4 + RHO_WIRE_EQ), (1 - W_HE_GAS) * 5.015e-4 / (5.015e-4 + RHO_WIRE_EQ),
			RHO_WIRE_EQ / (5.015e-4 + RHO_WIRE_EQ)}};

	private static StoppingTables cached;

	public final double lnT0, dl, invd;
	/** stopping power MeV cm2/g, flat [mat * NG + g]. */
	public final double[] S;
	/** CSDA range g/cm2, flat [mat * NG + g]. */
	public final double[] R;
	/** radiation length g/cm2 and Z/A of the nine materials. */
	public final double[] x0g, zoa;

	/** the proton-family tables (cached, immutable). */
	public static synchronized StoppingTables proton() {
		if (cached == null) cached = new StoppingTables();
		return cached;
	}

	private StoppingTables() {
		lnT0 = Math.log(1e-3);
		double lnT1 = Math.log(1e4);
		dl = (lnT1 - lnT0) / (NG - 1);
		invd = 1.0 / dl;
		double[] lnT = new double[NG];
		for (int g = 0; g < NG; g++) lnT[g] = lnT0 + dl * g;
		double[][] base = new double[BASE_X0.length][NG];
		for (int b = 0; b < BASE_X0.length; b++) {
			double[] t = NistPstar.T[b], sv = NistPstar.S[b];
			double[] lt = new double[t.length], ls = new double[t.length];
			for (int k = 0; k < t.length; k++) {
				lt[k] = Math.log(t[k]);
				ls[k] = Math.log(sv[k]);
			}
			double ltLast = lt[lt.length - 1];
			PolynomialSplineFunction sp = new SplineInterpolator().interpolate(lt, ls);
			double[] s = base[b];
			int last = -1;
			for (int g = 0; g < NG; g++) {
				double x = Math.min(Math.max(lnT[g], lt[0]), ltLast);
				s[g] = Math.exp(sp.value(x));
				if (lnT[g] <= ltLast) last = g;
			}
			double sl = s[last];
			for (int g = 0; g < NG; g++) if (lnT[g] > ltLast) s[g] = sl * Math.exp(-0.8 * (lnT[g] - ltLast));
		}
		S = new double[NMAT * NG];
		R = new double[NMAT * NG];
		x0g = new double[NMAT];
		zoa = new double[NMAT];
		double t0 = Math.exp(lnT0);
		double[] dr = new double[NG];
		for (int m = 0; m < NMAT; m++) {
			int o = m * NG;
			for (int g = 0; g < NG; g++) {
				double acc = 0.0;
				for (int c = 0; c < COMP_BASE[m].length; c++) acc += COMP_W[m][c] * base[COMP_BASE[m][c]][g];
				S[o + g] = acc;
				dr[g] = Math.exp(lnT[g]) / acc * dl;
			}
			double c0 = 2 * t0 / S[o];
			double cum = 0.0;
			R[o] = c0 + cum;
			for (int g = 1; g < NG; g++) {
				cum += 0.5 * (dr[g] + dr[g - 1]);
				R[o + g] = c0 + cum;
			}
			double ax = 0.0, az = 0.0;
			for (int c = 0; c < COMP_BASE[m].length; c++) {
				ax += COMP_W[m][c] / BASE_X0[COMP_BASE[m][c]];
				az += COMP_W[m][c] * BASE_ZOA[COMP_BASE[m][c]];
			}
			x0g[m] = 1.0 / ax;
			zoa[m] = az;
		}
	}

	/** prop._stop (fam 0): stopping power MeV cm2/g at kinetic energy tmev (MeV), linear in the ln T grid. */
	public double stop(int mat, double tmev) {
		double u = (Math.log(tmev) - lnT0) * invd;
		if (u < 0.0) u = 0.0;
		int i = (int) u;
		if (i > NG - 2) {
			i = NG - 2;
			u = (double) (i + 1);
		}
		double f = u - i;
		int o = mat * NG + i;
		return S[o] * (1.0 - f) + S[o + 1] * f;
	}

	/** prop._range_g (fam 0): CSDA range g/cm2 at tmev (MeV). */
	public double rangeG(int mat, double tmev) {
		double u = (Math.log(Math.max(tmev, 1.0001e-3)) - lnT0) * invd;
		int i = Math.min((int) u, NG - 2);
		double f = u - i;
		int o = mat * NG + i;
		return R[o] * (1.0 - f) + R[o + 1] * f;
	}

	/** prop._t_from_range (fam 0): kinetic energy (MeV) whose CSDA range is rg (g/cm2), binary search on the range table. */
	public double tFromRange(int mat, double rg) {
		int o = mat * NG;
		int lo = 0, hi = NG - 1;
		while (hi - lo > 1) {
			int mid = (lo + hi) >> 1;
			if (R[o + mid] < rg) lo = mid;
			else hi = mid;
		}
		double f = (rg - R[o + lo]) / (R[o + hi] - R[o + lo]);
		return Math.exp(lnT0 + (lo + Math.min(Math.max(f, 0.0), 1.0)) / invd);
	}
}
