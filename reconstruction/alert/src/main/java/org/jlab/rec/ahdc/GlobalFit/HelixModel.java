package org.jlab.rec.ahdc.GlobalFit;

/**
 * Closed-form helix in a uniform solenoid field and its closest approach to a wire. Ports fit._pt, _dpdx, _newton,
 * _closest and _hhit (fit.py:100-191, copied there from alfit.py). State x = (d0 mm, phi0, ln p [GeV], z0 mm, theta);
 * wire g has origin wo[3g..3g+2] and unit direction wu[3g..3g+2]; bz in T, (bx, by) the beam position in mm.
 * Static and allocation free.
 */
public final class HelixModel {

	/** arc-length grid of the closest-approach search (fit.py:49) */
	public static final int NG = 131;
	public static final double SMAX = 129.75;
	public static final double DS = SMAX / (NG - 1);

	private HelixModel() {}

	/** helix point, tangent and curvature vector at arc length s into o9 (fit._pt, fit.py:101-109); returns the radius R. */
	public static double pt(double[] x, double s, double bz, double bx, double by, double[] o9) {
		double p = Math.exp(x[2]), st = Math.sin(x[4]), ct = Math.cos(x[4]), k = 1.0;
		double R = p * st / (0.299792458 * bz) * 1000.0;
		double a = x[1] + k * s / R;
		double sa = Math.sin(a), ca = Math.cos(a), s0 = Math.sin(x[1]), c0 = Math.cos(x[1]);
		o9[0] = bx - x[0] * s0 + k * R * (sa - s0);
		o9[1] = by + x[0] * c0 - k * R * (ca - c0);
		o9[2] = x[3] + s * ct / st;
		o9[3] = ca;
		o9[4] = sa;
		o9[5] = ct / st;
		o9[6] = -sa * k / R;
		o9[7] = ca * k / R;
		o9[8] = 0.0;
		return R;
	}

	/** d(point)/d(x0..x4) at arc length s into the flat 5x3 j53 (fit._dpdx, fit.py:112-124). */
	public static void dpdx(double[] x, double s, double bz, double[] j53) {
		double p = Math.exp(x[2]), st = Math.sin(x[4]), ct = Math.cos(x[4]), k = 1.0;
		double R = p * st / (0.299792458 * bz) * 1000.0;
		double a = x[1] + k * s / R;
		double sa = Math.sin(a), ca = Math.cos(a), s0 = Math.sin(x[1]), c0 = Math.cos(x[1]);
		j53[0] = -s0; j53[1] = c0; j53[2] = 0.0;
		j53[3] = -x[0] * c0 + k * R * (ca - c0); j53[4] = -x[0] * s0 + k * R * (sa - s0); j53[5] = 0.0;
		double lx = k * R * (sa - s0) - s * ca, ly = -k * R * (ca - c0) - s * sa;
		j53[6] = lx; j53[7] = ly; j53[8] = 0.0;
		j53[9] = 0.0; j53[10] = 0.0; j53[11] = 1.0;
		double cot = ct / st;
		j53[12] = cot * lx; j53[13] = cot * ly; j53[14] = -s / (st * st);
	}

	private static double newton(double[] x, double[] wo, double[] wu, int g, double bz, double bx, double by, double[] o, double s,
			double lo, double hi) {
		int b = 3 * g;
		double ox = wo[b], oy = wo[b + 1], oz = wo[b + 2], ux = wu[b], uy = wu[b + 1], uz = wu[b + 2];
		for (int it = 0; it < 40; it++) {
			pt(x, s, bz, bx, by, o);
			double dx = o[0] - ox, dy = o[1] - oy, dz = o[2] - oz;
			double t = dx * ux + dy * uy + dz * uz;
			double px = dx - t * ux, py = dy - t * uy, pz = dz - t * uz;
			double tu = o[3] * ux + o[4] * uy + o[5] * uz;
			double f1 = px * o[3] + py * o[4] + pz * o[5];
			double g0 = o[3] * o[3] + o[4] * o[4] + o[5] * o[5] - tu * tu;
			double f2 = g0 + px * o[6] + py * o[7] + pz * o[8];
			if (f2 <= 1e-12) f2 = g0;
			double sn = Math.min(Math.max(s - f1 / f2, lo), hi);
			if (Math.abs(sn - s) < 1e-11) return sn;
			s = sn;
		}
		return s;
	}

	/**
	 * arc length of the closest approach to wire g: warm Newton from s0 > 0 (bracket +-6 mm), else grid argmin (DS) plus
	 * Newton; NaN at a grid end (fit._closest, fit.py:152-177).
	 */
	public static double closest(double[] x, double[] wo, double[] wu, int g, double bz, double bx, double by, double[] o9, double s0) {
		if (s0 > 0.0) {
			double lo = Math.max(s0 - 6.0, 1e-3), hi = Math.min(s0 + 6.0, SMAX);
			double s = newton(x, wo, wu, g, bz, bx, by, o9, s0, lo, hi);
			if (s > lo + 1e-9 && s < hi - 1e-9) return s;
		}
		int b = 3 * g;
		double ox = wo[b], oy = wo[b + 1], oz = wo[b + 2], ux = wu[b], uy = wu[b + 1], uz = wu[b + 2];
		int best = -1;
		double bd = 1e300;
		double rw = Math.sqrt(ox * ox + oy * oy);
		int j0 = Math.max((int) ((rw - 6.0) / DS), 0), j1 = Math.min((int) ((1.4 * rw + 6.0) / DS) + 2, NG);
		for (int j = j0; j < j1; j++) {
			pt(x, j * DS, bz, bx, by, o9);
			double dx = o9[0] - ox, dy = o9[1] - oy, dz = o9[2] - oz;
			double t = dx * ux + dy * uy + dz * uz;
			double p2 = dx * dx + dy * dy + dz * dz - t * t;
			if (p2 < bd) {
				bd = p2;
				best = j;
			}
		}
		if (best <= 0 || best >= NG - 1 || (best == j0 && j0 > 0) || best == j1 - 1) return Double.NaN;
		double s = newton(x, wo, wu, g, bz, bx, by, o9, best * DS, best * DS - DS, best * DS + DS);
		if (s <= 0.0 || s >= SMAX) return Double.NaN;
		return s;
	}

	/**
	 * signed helix-to-wire distance (legacy sign D . (T x U)/|T x U|) or NaN; sArc[0] = arc length or NaN; with want also
	 * d(d)/d(x0..x4) (envelope) into jx5 and the helix Jacobian into j53 (fit._hhit, fit.py:180-191).
	 */
	public static double hhit(double[] x, double[] wo, double[] wu, int g, double bz, double bx, double by, double[] o9, double s0,
			boolean want, double[] jx5, double[] j53, double[] sArc) {
		double s = closest(x, wo, wu, g, bz, bx, by, o9, s0);
		sArc[0] = s;
		if (Double.isNaN(s) || Double.isInfinite(s)) {
			sArc[0] = Double.NaN;
			return Double.NaN;
		}
		pt(x, s, bz, bx, by, o9);
		int b = 3 * g;
		double ux = wu[b], uy = wu[b + 1], uz = wu[b + 2];
		double dx = o9[0] - wo[b], dy = o9[1] - wo[b + 1], dz = o9[2] - wo[b + 2];
		double nx = o9[4] * uz - o9[5] * uy, ny = o9[5] * ux - o9[3] * uz, nz = o9[3] * uy - o9[4] * ux;
		double nn = Math.sqrt(nx * nx + ny * ny + nz * nz);
		if (nn < 1e-12) {
			sArc[0] = Double.NaN;
			return Double.NaN;
		}
		nx /= nn;
		ny /= nn;
		nz /= nn;
		if (want) {
			dpdx(x, s, bz, j53);
			for (int j = 0; j < 5; j++) jx5[j] = nx * j53[3 * j] + ny * j53[3 * j + 1] + nz * j53[3 * j + 2];
		}
		return dx * nx + dy * ny + dz * nz;
	}
}
