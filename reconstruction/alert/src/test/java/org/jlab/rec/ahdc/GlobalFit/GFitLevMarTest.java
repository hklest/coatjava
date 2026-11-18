package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;

/** SmallLinAlg (LU solve with partial pivoting, symmetric eigen decomposition) against commons-math3. */
public class GFitLevMarTest {

	private static final int N = 15;

	/** random symmetric matrix Q diag(lambda) Q^T: SPD with condition about 1e3 (even k), indefinite (odd k). */
	private static double[] matrix(Random rnd, int k) {
		RealMatrix g = new Array2DRowRealMatrix(N, N);
		for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) g.setEntry(i, j, rnd.nextGaussian());
		RealMatrix q = new org.apache.commons.math3.linear.QRDecomposition(g).getQ();
		double[] lam = new double[N];
		for (int i = 0; i < N; i++) {
			double mag = Math.pow(10.0, 3.0 * rnd.nextDouble() - 1.0);      // 0.1 .. 100
			lam[i] = (k % 2 == 0 || rnd.nextBoolean()) ? mag : -mag;
		}
		RealMatrix a = q.multiply(new Array2DRowRealMatrix(diag(lam))).multiply(q.transpose());
		double[] o = new double[N * N];
		for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) o[i * N + j] = 0.5 * (a.getEntry(i, j) + a.getEntry(j, i));
		return o;
	}

	private static double[][] diag(double[] d) {
		double[][] m = new double[N][N];
		for (int i = 0; i < N; i++) m[i][i] = d[i];
		return m;
	}

	@Test
	public void linAlgAgainstCommonsMath() {
		Random rnd = new Random(20261118L);
		double worstSolve = 0, worstEig = 0;
		int[] piv = new int[N];
		double[] w = new double[N], v = new double[N * N];
		for (int k = 0; k < 200; k++) {
			double[] a = matrix(rnd, k);
			double[] b = new double[N];
			for (int i = 0; i < N; i++) b[i] = rnd.nextGaussian();
			RealMatrix am = new Array2DRowRealMatrix(N, N);
			for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) am.setEntry(i, j, a[i * N + j]);
			// solve
			double[] x = b.clone(), ac = a.clone();
			assertTrue(SmallLinAlg.solve(ac, N, N, x, piv));
			double[] xc = new LUDecomposition(am).getSolver().solve(new ArrayRealVector(b)).toArray();
			double dn = 0, xn = 0;
			for (int i = 0; i < N; i++) { dn += (x[i] - xc[i]) * (x[i] - xc[i]); xn += xc[i] * xc[i]; }
			double rs = Math.sqrt(dn / xn);
			worstSolve = Math.max(worstSolve, rs);
			assertTrue("system " + k + " solve rel " + rs, rs <= 1e-12);
			// eigenvalues (relative to the largest |eigenvalue|, as both are accurate to that), eigenvector residual and orthogonality
			double[] ae = a.clone();
			SmallLinAlg.eigSym(ae, N, N, w, v);
			double[] ec = new EigenDecomposition(am).getRealEigenvalues();
			Arrays.sort(ec);
			double scale = Math.max(Math.abs(ec[0]), Math.abs(ec[N - 1]));
			for (int i = 0; i < N; i++) {
				double r = Math.abs(w[i] - ec[i]) / scale;
				worstEig = Math.max(worstEig, r);
				assertTrue("system " + k + " eigenvalue " + i + " java " + w[i] + " commons " + ec[i], r <= 1e-12);
			}
			for (int j = 0; j < N; j++) {
				double res = 0, nrm = 0;
				for (int i = 0; i < N; i++) {
					double s = -w[j] * v[i * N + j];
					for (int m = 0; m < N; m++) s += a[i * N + m] * v[m * N + j];
					res = Math.max(res, Math.abs(s));
					nrm += v[i * N + j] * v[i * N + j];
				}
				assertTrue("system " + k + " eigenvector " + j + " residual " + res, res <= 1e-12 * scale);
				assertTrue("system " + k + " eigenvector " + j + " norm " + nrm, Math.abs(nrm - 1.0) <= 1e-12);
			}
		}
		System.out.printf("linalg worst: solve %.2e relative, eigenvalue %.2e of the largest |eigenvalue|%n", worstSolve, worstEig);
	}

	@Test
	public void singularSolveReportsFailure() {
		double[] a = new double[N * N];
		double[] b = new double[N];
		assertTrue(!SmallLinAlg.solve(a, N, N, b, new int[N]));
	}
}
