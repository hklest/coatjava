package org.jlab.rec.ahdc.GlobalFit;

/**
 * Hit likelihood pieces of the pinned fit: the box-times-Cauchy cost of one drift-distance hit (fit._dbox, fit.py:195-226,
 * l = -c/2 as in fit._hit_l kinds 0 and 1, fit.py:229-246) and the signed closest approach of a propagated node table to a
 * wire (fit._pca1, fit.py:256-306, with prop._hs and prop._gfun, prop.py:768-788). Static and allocation free.
 */
public final class HitTerm {

	private static final double SQ2 = Math.sqrt(2.0), ISQ2PI = 1.0 / Math.sqrt(2.0 * Math.PI);
	private static final int NC = GFitWorkspace.NC;

	private HitTerm() {}

	/**
	 * dbox cost c = -2 l and its derivatives dc/dd, d2c/dd2 at signed distance d into out3 (fit._dbox). doca = measured
	 * drift distance, sig its sigma, h the half pitch. Legacy: raw doca and hard |d|; default: doca clamped at 0 and smooth
	 * |d| = sqrt(d^2 + delta^2) - delta. A non-finite d is replaced by 3 h and has zero derivatives.
	 */
	public static void dbox(double d, double doca, double sig, double h, boolean legacy, GFitConfig c, double[] out3) {
		double e = c.boxEdge(), eps = c.boxEps(), clip = c.dboxClip(), dl = c.delta();
		if (!legacy && doca < 0.0) doca = 0.0;
		boolean fin = !(Double.isNaN(d) || Double.isInfinite(d));
		if (!fin) d = 3.0 * h;
		double z0 = 0.5 * GFitMath.erfc(-(h / e) / SQ2) - 0.5 * GFitMath.erfc((h / e) / SQ2);
		double u1 = (h - d) / e, u2 = (-h - d) / e;
		double pb = (0.5 * GFitMath.erfc(-u1 / SQ2) - 0.5 * GFitMath.erfc(-u2 / SQ2)) / z0;
		double q = pb * (1.0 - eps) + eps;
		double cc = -2.0 * Math.log(q);
		double a, a1, a2;
		if (legacy) {
			a = Math.abs(d);
			a1 = d >= 0 ? 1.0 : -1.0;
			a2 = 0.0;
		} else {
			double rr = Math.sqrt(d * d + dl * dl);
			a = rr - dl;
			a1 = d / rr;
			a2 = dl * dl / (rr * rr * rr);
		}
		double r = doca - a;
		double rs = r / sig;
		double lc = Math.log1p(rs * rs);
		cc += lc < clip ? lc : clip;
		if (!fin) {
			out3[0] = cc;
			out3[1] = 0.0;
			out3[2] = 0.0;
			return;
		}
		double f1 = Math.exp(-0.5 * u1 * u1) * ISQ2PI, f2 = Math.exp(-0.5 * u2 * u2) * ISQ2PI;
		double pb1 = (-f1 + f2) / (e * z0), pb2 = (-u1 * f1 + u2 * f2) / (e * e * z0);
		double c1 = -2.0 * (1.0 - eps) * pb1 / q;
		double c2 = -2.0 * (1.0 - eps) * (pb2 / q - (1.0 - eps) * pb1 * pb1 / (q * q));
		if (lc < clip) {
			double den = sig * sig + r * r;
			c1 += 2.0 * r / den * (-a1);
			c2 += 2.0 * (sig * sig - r * r) / (den * den) * a1 * a1 + 2.0 * r / den * (-a2);
		}
		out3[0] = cc;
		out3[1] = c1;
		out3[2] = c2;
	}

	/** cubic Hermite position (r[o..o+2]) and dr/ds (r[o+3..o+5]) of segment i at tau (prop._hs). */
	private static void hs(double[] a, int i, double tau, double[] r, int o) {
		int p = i * NC, q = (i + 1) * NC;
		double d = a[q] - a[p];
		double t2 = tau * tau, t3 = tau * tau * tau;
		double h00 = 2 * t3 - 3 * t2 + 1, h10 = (t3 - 2 * t2 + tau) * d, h01 = -2 * t3 + 3 * t2, h11 = (t3 - t2) * d;
		double g00 = (6 * t2 - 6 * tau) / d, g10 = 3 * t2 - 4 * tau + 1, g01 = (-6 * t2 + 6 * tau) / d, g11 = 3 * t2 - 2 * tau;
		r[o] = h00 * a[p + 1] + h10 * a[p + 4] + h01 * a[q + 1] + h11 * a[q + 4];
		r[o + 1] = h00 * a[p + 2] + h10 * a[p + 5] + h01 * a[q + 2] + h11 * a[q + 5];
		r[o + 2] = h00 * a[p + 3] + h10 * a[p + 6] + h01 * a[q + 3] + h11 * a[q + 6];
		r[o + 3] = g00 * a[p + 1] + g10 * a[p + 4] + g01 * a[q + 1] + g11 * a[q + 4];
		r[o + 4] = g00 * a[p + 2] + g10 * a[p + 5] + g01 * a[q + 2] + g11 * a[q + 5];
		r[o + 5] = g00 * a[p + 3] + g10 * a[p + 6] + g01 * a[q + 3] + g11 * a[q + 6];
	}

	/** first component of prop._gfun at the point and direction r[o..o+5] of the wire (a0, a1, a2, w). */
	private static double gfun(double[] r, int o, double a0, double a1, double a2, double wx, double wy, double wz) {
		double dx = r[o] - a0, dy = r[o + 1] - a1, dz = r[o + 2] - a2;
		double dw = dx * wx + dy * wy + dz * wz;
		double px = dx - dw * wx, py = dy - dw * wy, pz = dz - dw * wz;
		double vx = r[o + 3], vy = r[o + 4], vz = r[o + 5];
		double vw = vx * wx + vy * wy + vz * wz;
		return px * (vx - vw * wx) + py * (vy - vw * wy) + pz * (vz - vw * wz);
	}

	/**
	 * signed closest approach of the stepped trajectory (node table a, flat [i * NC + c], nn nodes) to wire g (ends wo/we, flat
	 * n x 3): fit._pca1. res4 = (signed d, PCA x, y, z), all NaN with flag 1 when no node pair has positive length; returns
	 * the flag (1 no root, 2 PCA outside the wire). gs and ds are scratch of length >= nn; res4 must have length >= 10
	 * (slots 4..9 are scratch of the Hermite evaluation). With res4.length > 10 slot 10 receives the arc length s of the PCA
	 * (prop._pca1 return value s = a[i0].s + tau (a[i0+1].s - a[i0].s), prop.py:838; NaN with flag 1 and no valid segment).
	 * A trajectory of fewer than two nodes (the track stopped at its first step) has no segment either: the reference reads
	 * outside its node table there (a[-1], numpy wrap-around), this port returns the NaN result instead of throwing.
	 */
	public static int pca1(double[] a, int nn, double[] wo, double[] we, int g, double[] gs, double[] ds, double[] res4) {
		if (nn < 2) {
			res4[0] = res4[1] = res4[2] = res4[3] = Double.NaN;
			if (res4.length > 10) res4[10] = Double.NaN;
			return 1;
		}
		int b = 3 * g;
		double ax = wo[b], ay = wo[b + 1], az = wo[b + 2];
		double l0 = Math.sqrt((we[b] - ax) * (we[b] - ax) + (we[b + 1] - ay) * (we[b + 1] - ay) + (we[b + 2] - az) * (we[b + 2] - az));
		double wx = (we[b] - ax) / l0, wy = (we[b + 1] - ay) / l0, wz = (we[b + 2] - az) / l0;
		double best = 1e30;
		int ibest = 0;
		for (int i = 0; i < nn; i++) {
			int o = i * NC;
			double dx = a[o + 1] - ax, dy = a[o + 2] - ay, dz = a[o + 3] - az;
			double dw = dx * wx + dy * wy + dz * wz;
			double px = dx - dw * wx, py = dy - dw * wy, pz = dz - dw * wz;
			double vx = a[o + 4], vy = a[o + 5], vz = a[o + 6];
			double vw = vx * wx + vy * wy + vz * wz;
			gs[i] = px * (vx - vw * wx) + py * (vy - vw * wy) + pz * (vz - vw * wz);
			ds[i] = px * px + py * py + pz * pz;
			if (ds[i] < best) {
				best = ds[i];
				ibest = i;
			}
		}
		int root = -1;
		for (int i = 0; i < nn - 1; i++) {
			if (a[(i + 1) * NC] - a[i * NC] > 1e-12 && gs[i] <= 0.0 && gs[i + 1] > 0.0
					&& (root < 0 || Math.min(ds[i], ds[i + 1]) < Math.min(ds[root], ds[root + 1])))
				root = i;
		}
		int flag = 0, i0;
		double tau;
		if (root < 0) {
			flag = 1;
			i0 = Math.min(ibest, nn - 2);
			tau = ibest < nn - 1 ? 0.0 : 1.0;
			while (i0 > 0 && !(a[(i0 + 1) * NC] - a[i0 * NC] > 1e-12)) {
				i0--;
				tau = 1.0;
			}
			if (!(a[(i0 + 1) * NC] - a[i0 * NC] > 1e-12)) {
				res4[0] = res4[1] = res4[2] = res4[3] = Double.NaN;
				if (res4.length > 10) res4[10] = Double.NaN;
				return 1;
			}
		} else {
			i0 = root;
			double lo = 0.0, hi = 1.0, glo = gs[root], ghi = gs[root + 1];
			tau = 0.5;
			for (int it = 0; it < 40; it++) {
				tau = ghi != glo ? hi - ghi * (hi - lo) / (ghi - glo) : 0.5 * (lo + hi);
				hs(a, i0, tau, res4, 4);
				double gtv = gfun(res4, 4, ax, ay, az, wx, wy, wz);
				if (Math.abs(gtv) < 1e-10 || hi - lo < 1e-10) break;
				if (gtv < 0.0) {
					lo = tau;
					glo = gtv;
					ghi *= 0.5;
				} else {
					hi = tau;
					ghi = gtv;
					glo *= 0.5;
				}
			}
		}
		hs(a, i0, tau, res4, 4);
		double rx = res4[4], ry = res4[5], rz = res4[6], vx = res4[7], vy = res4[8], vz = res4[9];
		double nv = Math.sqrt(vx * vx + vy * vy + vz * vz);
		double dx = rx - ax, dy = ry - ay, dz = rz - az;
		double dw = dx * wx + dy * wy + dz * wz;
		double dist = Math.sqrt((dx - dw * wx) * (dx - dw * wx) + (dy - dw * wy) * (dy - dw * wy) + (dz - dw * wz) * (dz - dw * wz));
		double nx = vy * wz - vz * wy, ny = vz * wx - vx * wz, nz = vx * wy - vy * wx;
		double sign = (Math.sqrt(nx * nx + ny * ny + nz * nz) > 1e-9 * nv && (dx * nx + dy * ny + dz * nz) < 0.0) ? -1.0 : 1.0;
		if (dw < -1e-6 || dw > l0 + 1e-6) flag |= 2;
		if (res4.length > 10) res4[10] = a[i0 * NC] + tau * (a[(i0 + 1) * NC] - a[i0 * NC]);
		res4[0] = sign * dist;
		res4[1] = rx;
		res4[2] = ry;
		res4[3] = rz;
		return flag;
	}
}
