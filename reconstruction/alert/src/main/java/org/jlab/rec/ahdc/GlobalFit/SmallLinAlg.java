package org.jlab.rec.ahdc.GlobalFit;

/**
 * Fixed-size dense linear algebra on flat row-major arrays without allocation, for the Levenberg-Marquardt loop and the
 * covariance of the track fit: LU solve with partial pivoting (what numpy.linalg.solve does in fit._lm, fit.py:463, and in
 * the flip step, fit.py:703) and the symmetric eigen decomposition (numpy.linalg.eigh in fit.py:806-813) by cyclic Jacobi
 * rotations. Eigenvalues and vectors agree with LAPACK to rounding.
 */
public final class SmallLinAlg {

	private SmallLinAlg() {}

	/**
	 * Solves A x = b for the leading n x n block of the row-major matrix a (leading dimension lda), overwriting a with its LU
	 * factors and b with the solution; piv is scratch of length n. Returns false when a pivot is exactly zero (LAPACK info &gt; 0,
	 * where numpy raises LinAlgError).
	 */
	public static boolean solve(double[] a, int lda, int n, double[] b, int[] piv) {
		for (int k = 0; k < n; k++) {
			int p = k;
			double big = Math.abs(a[k * lda + k]);
			for (int i = k + 1; i < n; i++) {
				double v = Math.abs(a[i * lda + k]);
				if (v > big) { big = v; p = i; }
			}
			piv[k] = p;
			if (big == 0.0) return false;
			if (p != k) {
				for (int j = 0; j < n; j++) {
					double t = a[k * lda + j]; a[k * lda + j] = a[p * lda + j]; a[p * lda + j] = t;
				}
				double t = b[k]; b[k] = b[p]; b[p] = t;
			}
			double d = a[k * lda + k];
			for (int i = k + 1; i < n; i++) {
				double f = a[i * lda + k] / d;
				a[i * lda + k] = f;
				if (f != 0.0) {
					for (int j = k + 1; j < n; j++) a[i * lda + j] -= f * a[k * lda + j];
					b[i] -= f * b[k];
				}
			}
		}
		for (int i = n - 1; i >= 0; i--) {
			double s = b[i];
			for (int j = i + 1; j < n; j++) s -= a[i * lda + j] * b[j];
			b[i] = s / a[i * lda + i];
		}
		return true;
	}

	/**
	 * Eigen decomposition of the symmetric n x n block of a (leading dimension lda, both triangles used, destroyed): eigenvalues
	 * ascending into w[0..n), eigenvector j in column j of v (row-major n x n, v[i * n + j]).
	 */
	public static void eigSym(double[] a, int lda, int n, double[] w, double[] v) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) v[i * n + j] = i == j ? 1.0 : 0.0;
		}
		for (int sweep = 0; sweep < 60; sweep++) {
			double off = 0.0;
			for (int p = 0; p < n - 1; p++)
				for (int q = p + 1; q < n; q++) off += a[p * lda + q] * a[p * lda + q];
			if (off == 0.0) break;
			double diag = 0.0;
			for (int p = 0; p < n; p++) diag += a[p * lda + p] * a[p * lda + p];
			if (off <= 1e-34 * diag) break;
			for (int p = 0; p < n - 1; p++) {
				for (int q = p + 1; q < n; q++) {
					double apq = a[p * lda + q];
					if (apq == 0.0) continue;
					double app = a[p * lda + p], aqq = a[q * lda + q];
					double theta = (aqq - app) / (2.0 * apq);
					double t = Math.copySign(1.0, theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1.0));
					double c = 1.0 / Math.sqrt(t * t + 1.0), s = t * c;
					for (int k = 0; k < n; k++) {
						double akp = a[k * lda + p], akq = a[k * lda + q];
						a[k * lda + p] = c * akp - s * akq;
						a[k * lda + q] = s * akp + c * akq;
					}
					for (int k = 0; k < n; k++) {
						double apk = a[p * lda + k], aqk = a[q * lda + k];
						a[p * lda + k] = c * apk - s * aqk;
						a[q * lda + k] = s * apk + c * aqk;
					}
					a[p * lda + q] = 0.0;
					a[q * lda + p] = 0.0;
					for (int k = 0; k < n; k++) {
						double vkp = v[k * n + p], vkq = v[k * n + q];
						v[k * n + p] = c * vkp - s * vkq;
						v[k * n + q] = s * vkp + c * vkq;
					}
				}
			}
		}
		for (int i = 0; i < n; i++) w[i] = a[i * lda + i];
		// ascending order (selection sort of eigenvalues, columns swapped with them)
		for (int i = 0; i < n - 1; i++) {
			int m = i;
			for (int j = i + 1; j < n; j++) if (w[j] < w[m]) m = j;
			if (m != i) {
				double t = w[i]; w[i] = w[m]; w[m] = t;
				for (int k = 0; k < n; k++) {
					double x = v[k * n + i]; v[k * n + i] = v[k * n + m]; v[k * n + m] = x;
				}
			}
		}
	}
}
