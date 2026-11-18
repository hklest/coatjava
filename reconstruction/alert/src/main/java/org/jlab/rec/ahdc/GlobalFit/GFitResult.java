package org.jlab.rec.ahdc.GlobalFit;

import java.util.Arrays;

/**
 * Result of one track fit: the fields of fit.pack_est (fit.py:1114-1144) in the internal u-space (state u = d0, phi0,
 * ln p, z0, theta, 10 kinks; covariance C = 2 H^-1, 15 x 15 row-major), cost, status, flags, counters and the per-hit
 * (d, s, l). Mutable and reusable: {@link #clear()} resets it; arrays have room for maxHits + 8 hits.
 */
public final class GFitResult {
	public final double[] u = new double[GFitConfig.NPAR];
	public final double[] cov = new double[GFitConfig.NPAR * GFitConfig.NPAR];
	public final double[] theta0 = new double[GFitConfig.NPL];
	public double cost, dcost2;
	public int status, flags, niter, nstart, nmin, startWon;
	/** number of hits of the final hit set and per hit: d, s, l, kind, row in event (-1 missing wire), wire. */
	public int nh;
	public final double[] hitD, hitS, hitL;
	public final int[] hitKind, hitRow, hitWire;

	public GFitResult(int maxHits) {
		int n = maxHits + GFitConfig.NADD;
		hitD = new double[n]; hitS = new double[n]; hitL = new double[n];
		hitKind = new int[n]; hitRow = new int[n]; hitWire = new int[n];
	}

	public void clear() {
		Arrays.fill(u, 0.0); Arrays.fill(cov, 0.0); Arrays.fill(theta0, 0.0);
		cost = 0.0; dcost2 = 0.0;
		status = 0; flags = 0; niter = 0; nstart = 0; nmin = 0; startWon = 0; nh = 0;
		Arrays.fill(hitD, 0.0); Arrays.fill(hitS, 0.0); Arrays.fill(hitL, 0.0);
		Arrays.fill(hitKind, 0); Arrays.fill(hitRow, 0); Arrays.fill(hitWire, 0);
	}

	/** momentum [GeV]. */
	public double p() { return Math.exp(u[2]); }

	public double thetaDeg() { return Math.toDegrees(u[4]); }

	/** phi [deg] wrapped to [-180, 180) as pack_est: (degrees + 180) % 360 - 180 with the Python (floored) modulo. */
	public double phiDeg() {
		double m = (Math.toDegrees(u[1]) + 180) % 360;
		if (m != 0 && m < 0) m += 360;
		return m - 180;
	}

	public double zv() { return u[3]; }

	/** 5 x 5 covariance of pack_est (est_cov), state (d0, phi0, z0, theta, q/p), row-major into out[0..24]. */
	public void covX(double[] out) {
		final int[] idx = {0, 1, 3, 4, 2};
		double[] s = {1, 1, 1, 1, -1.0 / p()};
		int n = GFitConfig.NPAR;
		for (int i = 0; i < 5; i++)
			for (int l = 0; l < 5; l++)
				out[i * 5 + l] = s[i] * cov[idx[i] * n + idx[l]] * s[l];
	}
}
