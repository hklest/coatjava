package org.jlab.rec.ahdc.GlobalFit;

import org.apache.commons.math3.special.Erf;

/**
 * ATOF position term of the global fit: per-track candidate attach from the event's ATOF rows and the B0 seed (port of
 * atof.seed_intercepts, atof.attach and the frozen wedge association weight, atof.py:102-123, 387-450), the likelihood
 * pieces (_boxg, _gau, _wrap, _tpdf, _layer_eval, atof.py:158-230) and the cost with its gradient and both Hessians chained
 * through the finite-difference intercept Jacobian (term_cost, atof.py:233-293). The parameter vector is the 34 slots P_*
 * of atof.py:41-45 as loaded by load_par (atof.py:330-356). All per-track state lives in {@link TrackState}; nothing is
 * allocated by termCost.
 */
public final class AtofTerm {

	/** intercept record length, window constants (atof.py:29-33). */
	static final int NI = Propagator.NI;
	static final double BAR_HALF = 139.85, WEDGE_HALF = 5 * 27.7, ZACC = 5.0, SRCH = 2.0;
	static final double ROW = Math.toRadians(6.0), WIN_ROWS = 3.0, NONE_ROWS = 1.5;
	private static final double SQ2 = Math.sqrt(2.0), ISQ2PI = 1.0 / Math.sqrt(2.0 * Math.PI), TWOPI = 2.0 * Math.PI;
	private static final double CEIL_STH = Math.toRadians(0.05), CEIL_SZ = 0.5;
	private static final int NPAR = GFitWorkspace.NPAR;
	private static final int NC = GFitWorkspace.NC;

	/** indices of the parameter vector (atof.py:41-43). */
	public static final int P_WZLEN = 0, P_WZSIG = 1, P_WZOFF = 2, P_ROWW = 3, P_WPSIG = 4, P_WPOFF = 5, P_BZOFF = 6, P_BZS1 = 7,
		P_BZS2 = 8, P_BZF2 = 9, P_BPSIG = 10, P_BPOFF = 11, P_ETAW = 12, P_ETAB = 13, P_RHOW = 14, P_RHOB = 15, P_TOFF = 16,
		P_TS1 = 17, P_TS2 = 18, P_TF2 = 19, P_TML = 20, P_TSL = 21, P_TFL = 22, P_TFB = 23, P_TWIN = 24, P_TON = 25, P_TMODE = 26,
		P_BARON = 27, P_FORCEW = 28, P_CEIL = 29, P_ZMID = 30, P_TASSOC = 31, P_TRA = 32, P_TRB = 33;
	public static final int NP = 34;

	private AtofTerm() {}

	/** Per-track ATOF state: candidates, parameter vector, seed record and the fixed scratch of termCost. One per thread. */
	public static final class TrackState {
		/** candidates: z [mm], row-centre phi [rad], frozen weight (or time - start), bar flag, sector, component, row in event. */
		public double[] cz = new double[64], cphi = new double[64], ct = new double[64];
		public int[] cbar = new int[64], csec = new int[64], ccomp = new int[64], crow = new int[64];
		public int ncand;
		/** ATOF_NONE flag (no candidate near the seed intercept). */
		public int fnone;
		/** parameter vector of the track's half (reference, not copied). */
		public double[] par;
		/** intercept record of the seed (atof_seed). */
		public final double[] record = new double[NI];
		/** electron-vz prior fallback for an ATOF_NONE track (ci[14]) and its width (cf[4]); off in the default variant. */
		public int fb;
		public double sigz = 1.0;

		// scratch of the layer evaluation and term_cost
		final double[] x = new double[4], Cx = new double[4], Cxx = new double[16], info = new double[7];
		final double[] J = new double[4 * NPAR], Cc = new double[16];
		final double[] Sx = new double[3], Sxx = new double[9];
		final double[] gz3 = new double[3], gp3 = new double[3], g1 = new double[3], g2 = new double[3], g3 = new double[3];
		final double[] tp3 = new double[3];
		final double[] ea = new double[16], ev = new double[16], ew = new double[4];

		void grow(int n) {
			if (n <= cz.length) return;
			int m = Math.max(n, 2 * cz.length);
			cz = java.util.Arrays.copyOf(cz, m); cphi = java.util.Arrays.copyOf(cphi, m); ct = java.util.Arrays.copyOf(ct, m);
			cbar = java.util.Arrays.copyOf(cbar, m); csec = java.util.Arrays.copyOf(csec, m);
			ccomp = java.util.Arrays.copyOf(ccomp, m); crow = java.util.Arrays.copyOf(crow, m);
		}
	}

	// ================================================================== likelihood pieces

	static double sstep(double x) {
		if (x <= 0.0) return 0.0;
		if (x >= 1.0) return 1.0;
		return x * x * (3.0 - 2.0 * x);
	}

	/** box(ell) convolved with N(0, sig) at d: value and first, second derivative in d (atof._boxg). */
	static void boxg(double d, double ell, double sig, double[] o) {
		double a = (d + 0.5 * ell) / sig, b = (d - 0.5 * ell) / sig;
		o[0] = (0.5 * Erf.erfc(-a / SQ2) - 0.5 * Erf.erfc(-b / SQ2)) / ell;
		double pa = Math.exp(-0.5 * a * a) * ISQ2PI, pb = Math.exp(-0.5 * b * b) * ISQ2PI;
		o[1] = (pa - pb) / (sig * ell);
		o[2] = (-a * pa + b * pb) / (sig * sig * ell);
	}

	/** Gaussian and its derivatives in r (atof._gau). */
	static void gau(double r, double s, double[] o) {
		double q = r / s;
		double n = Math.exp(-0.5 * (q * q)) * ISQ2PI / s;
		o[0] = n;
		o[1] = -r / (s * s) * n;
		o[2] = (r * r / Math.pow(s, 4) - 1.0 / (s * s)) * n;
	}

	static double wrap(double x) { return x - TWOPI * Math.floor((x + Math.PI) / TWOPI); }

	/** wedge time pdf at the residual r and its first and second derivative (atof._tpdf); scratch g1, g2, g3 of st. */
	static void tpdf(double r, double[] par, TrackState st, double[] o) {
		gau(r - par[P_TOFF], par[P_TS1], st.g1);
		gau(r - par[P_TOFF], par[P_TS2], st.g2);
		gau(r - par[P_TOFF] - par[P_TML], par[P_TSL], st.g3);
		double f2 = par[P_TF2], fl = par[P_TFL], fb = par[P_TFB];
		double c = 1.0 - fl - fb;
		o[0] = c * ((1 - f2) * st.g1[0] + f2 * st.g2[0]) + fl * st.g3[0] + fb / (2.0 * par[P_TWIN]);
		o[1] = c * ((1 - f2) * st.g1[1] + f2 * st.g2[1]) + fl * st.g3[1];
		o[2] = c * ((1 - f2) * st.g1[2] + f2 * st.g2[2]) + fl * st.g3[2];
	}

	/**
	 * Cost C = -2 log LR of one layer (L 0 bar, 1 wedge) at the intercept x = (phi, z, t, w); with want the gradient st.Cx (4)
	 * and exact Hessian st.Cxx (4 x 4 flat) in x. st.info: posterior hit probability, index of the most probable candidate
	 * (-1 none), its residuals dz, dphi, dt. x[3] &lt;= 0: layer not reached, C = 0 (atof._layer_eval).
	 */
	static double layerEval(int L, double[] x, TrackState st, boolean want) {
		final double[] par = st.par, info = st.info, Cx = st.Cx, Cxx = st.Cxx;
		for (int i = 0; i < 7; i++) info[i] = 0.0;
		info[1] = -1.0;
		if (want) {
			for (int i = 0; i < 4; i++) Cx[i] = 0.0;
			for (int i = 0; i < 16; i++) Cxx[i] = 0.0;
		}
		double w = x[3];
		if (!(w > 0.0)) return 0.0;
		boolean bar = L == 0;
		double eta = bar ? par[P_ETAB] : par[P_ETAW];
		boolean ton = !bar && par[P_TON] > 0.5;
		double rho = bar ? par[P_RHOB] : (ton ? par[P_RHOW] / (2.0 * par[P_TWIN]) : par[P_RHOW]);
		double S = 0.0;
		final double[] Sx = st.Sx, Sxx = st.Sxx;
		for (int i = 0; i < 3; i++) Sx[i] = 0.0;
		for (int i = 0; i < 9; i++) Sxx[i] = 0.0;
		double best = -1.0;
		final double[] fz = st.gz3, fp = st.gp3;
		for (int c = 0; c < st.ncand; c++) {
			if (st.cbar[c] >= 2 || (st.cbar[c] == 1) != bar) continue;
			double dz, dp;
			if (bar) {
				dz = st.cz[c] - x[1] - par[P_BZOFF];
				gau(dz, par[P_BZS1], st.g1);
				gau(dz, par[P_BZS2], st.g2);
				double f2 = par[P_BZF2];
				fz[0] = (1 - f2) * st.g1[0] + f2 * st.g2[0];
				fz[1] = (1 - f2) * st.g1[1] + f2 * st.g2[1];
				fz[2] = (1 - f2) * st.g1[2] + f2 * st.g2[2];
				dp = wrap(st.cphi[c] - x[0] - par[P_BPOFF]);
				boxg(dp, par[P_ROWW], par[P_BPSIG], fp);
			} else {
				dz = st.cz[c] - x[1] - par[P_WZOFF];
				boxg(dz, par[P_WZLEN], par[P_WZSIG], fz);
				dp = wrap(st.cphi[c] - x[0] - par[P_WPOFF]);
				boxg(dp, par[P_ROWW], par[P_WPSIG], fp);
			}
			double ft = 1.0, ft1 = 0.0, ft2 = 0.0, dt = 0.0;
			if (!bar && par[P_TASSOC] > 0.5 && Double.isFinite(st.ct[c])) ft = st.ct[c];
			if (ton) {
				if (Double.isFinite(st.ct[c])) {
					dt = st.ct[c] - x[2];
					tpdf(dt, par, st, st.tp3);
					ft = st.tp3[0]; ft1 = st.tp3[1]; ft2 = st.tp3[2];
				} else {
					ft = 1.0 / (2.0 * par[P_TWIN]);
				}
			}
			double gz = -fz[1], gp = -fp[1], gt = -ft1;
			double Lc = fz[0] * fp[0] * ft;
			S += Lc;
			if (Lc > best) {
				best = Lc; info[1] = c; info[2] = dz; info[3] = dp; info[4] = dt;
			}
			if (want) {
				Sx[0] += gp * fz[0] * ft; Sx[1] += gz * fp[0] * ft; Sx[2] += fz[0] * fp[0] * gt;
				Sxx[0] += fp[2] * fz[0] * ft; Sxx[4] += fz[2] * fp[0] * ft; Sxx[8] += fz[0] * fp[0] * ft2;
				Sxx[1] += gp * gz * ft; Sxx[2] += gp * fz[0] * gt; Sxx[5] += fp[0] * gz * gt;
			}
		}
		S /= rho;
		double LR = (1.0 - eta * w) + eta * w * S;
		if (1e-300 > LR) LR = 1e-300;
		info[0] = eta * w * S / LR;
		if (want) {
			for (int i = 0; i < 3; i++) Sx[i] /= rho;
			for (int i = 0; i < 9; i++) Sxx[i] /= rho;
			Sxx[3] = Sxx[1]; Sxx[6] = Sxx[2]; Sxx[7] = Sxx[5];
			double[] LRx = st.ea, LRxx = st.ev;  // 4 and 16 of the eigen scratch, free here
			for (int i = 0; i < 16; i++) LRxx[i] = 0.0;
			for (int a = 0; a < 3; a++) {
				LRx[a] = eta * w * Sx[a];
				LRxx[a * 4 + 3] = eta * Sx[a];
				LRxx[3 * 4 + a] = eta * Sx[a];
				for (int b = 0; b < 3; b++) LRxx[a * 4 + b] = eta * w * Sxx[a * 3 + b];
			}
			LRx[3] = eta * (S - 1.0);
			for (int a = 0; a < 4; a++) {
				Cx[a] = -2.0 * LRx[a] / LR;
				for (int b = 0; b < 4; b++) Cxx[a * 4 + b] = -2.0 * LRxx[a * 4 + b] / LR + 2.0 * LRx[a] * LRx[b] / (LR * LR);
			}
		}
		return -2.0 * Math.log(LR);
	}

	/** eigen decomposition of the symmetric 4 x 4 matrix st.Cxx by cyclic Jacobi rotations: eigenvalues st.ew, vectors in the columns of st.ev. */
	static void eigh4(TrackState st) {
		final double[] a = st.ea, v = st.ev, w = st.ew;
		System.arraycopy(st.Cxx, 0, a, 0, 16);
		for (int i = 0; i < 16; i++) v[i] = 0.0;
		for (int i = 0; i < 4; i++) v[i * 5] = 1.0;
		double fro = 0.0;
		for (int i = 0; i < 16; i++) fro += a[i] * a[i];
		for (int sweep = 0; sweep < 60; sweep++) {
			double off = 0.0;
			for (int p = 0; p < 3; p++) for (int q = p + 1; q < 4; q++) off += 2.0 * a[p * 4 + q] * a[p * 4 + q];
			if (!(off > 1e-34 * fro)) break;
			for (int p = 0; p < 3; p++) {
				for (int q = p + 1; q < 4; q++) {
					double apq = a[p * 4 + q];
					if (apq == 0.0) continue;
					double theta = (a[q * 4 + q] - a[p * 4 + p]) / (2.0 * apq);
					double t = (theta >= 0 ? 1.0 : -1.0) / (Math.abs(theta) + Math.sqrt(theta * theta + 1.0));
					double c = 1.0 / Math.sqrt(t * t + 1.0), s = t * c;
					for (int k = 0; k < 4; k++) {  // A <- A J
						double akp = a[k * 4 + p], akq = a[k * 4 + q];
						a[k * 4 + p] = c * akp - s * akq;
						a[k * 4 + q] = s * akp + c * akq;
					}
					for (int k = 0; k < 4; k++) {  // A <- J^T A
						double apk = a[p * 4 + k], aqk = a[q * 4 + k];
						a[p * 4 + k] = c * apk - s * aqk;
						a[q * 4 + k] = s * apk + c * aqk;
					}
					for (int k = 0; k < 4; k++) {  // V <- V J
						double vkp = v[k * 4 + p], vkq = v[k * 4 + q];
						v[k * 4 + p] = c * vkp - s * vkq;
						v[k * 4 + q] = s * vkp + c * vkq;
					}
				}
			}
		}
		for (int i = 0; i < 4; i++) w[i] = a[i * 5];
	}

	/**
	 * ATOF cost at u from the intercept array ia (16 rows x 20: row 0 nominal, row 1 + j the propagation at u + step_j e_j).
	 * With want adds the gradient to g and the Gauss-Newton (PSD-clipped I-space Hessian) and exact I-space Hessians to hc,
	 * he (upper triangle j &lt;= m, flat [j * 15 + m]); the caller symmetrises (atof.term_cost, atof.py:233-293).
	 */
	public static double termCost(double[] u, double[] ia, boolean want, int npar, double ze, TrackState st, double[] g,
			double[] hc, double[] he) {
		final double[] par = st.par, x = st.x, Cx = st.Cx, Cxx = st.Cxx, J = st.J, Cc = st.Cc;
		double fv = 0.0;
		for (int L = 0; L < 2; L++) {
			if (L == 0 && par[P_BARON] < 0.5) continue;
			int b = 8 * L;
			boolean zmid = par[P_ZMID] > 0.5;
			int ip = zmid ? 16 + 2 * L : b;
			int iz = zmid ? ip + 1 : b + 1;
			boolean tmode = par[P_TMODE] > 0.5;
			x[0] = ia[ip]; x[1] = ia[iz]; x[3] = ia[b + 5];
			x[2] = par[P_TRA] + par[P_TRB] * (tmode ? ia[b + 3] : ia[b + 2]);
			if (par[P_FORCEW] > 0.5) x[3] = (((int) ia[b + 7]) & 1) != 0 ? 1.0 : 0.0;
			fv += layerEval(L, x, st, want);
			if (want && x[3] > 0.0) {
				int fl0 = ((int) ia[b + 7]) & 1;
				int ti = tmode ? b + 3 : b + 2;
				for (int j = 0; j < npar; j++) {
					int r = (1 + j) * NI;
					if ((((int) ia[r + b + 7]) & 1) != fl0) {
						for (int a = 0; a < 4; a++) J[a * NPAR + j] = 0.0;
						continue;
					}
					double hj = j < 5 ? GFitConfig.stepU(j) : 1e-5;
					J[j] = wrap(ia[r + ip] - ia[ip]) / hj;
					J[NPAR + j] = (ia[r + iz] - ia[iz]) / hj;
					J[2 * NPAR + j] = par[P_TRB] * (ia[r + ti] - ia[ti]) / hj;
					J[3 * NPAR + j] = par[P_FORCEW] > 0.5 ? 0.0 : (ia[r + b + 5] - ia[b + 5]) / hj;
				}
				eigh4(st);
				for (int i = 0; i < 16; i++) Cc[i] = 0.0;
				for (int q = 0; q < 4; q++) {
					double e = st.ew[q];
					if (e > 0.0) {
						for (int a = 0; a < 4; a++) for (int c = 0; c < 4; c++) Cc[a * 4 + c] += e * st.ev[a * 4 + q] * st.ev[c * 4 + q];
					}
				}
				for (int j = 0; j < npar; j++) {
					for (int a = 0; a < 4; a++) g[j] += Cx[a] * J[a * NPAR + j];
					for (int m = j; m < npar; m++) {
						for (int a = 0; a < 4; a++) {
							for (int c = 0; c < 4; c++) {
								hc[j * NPAR + m] += J[a * NPAR + j] * Cc[a * 4 + c] * J[c * NPAR + m];
								he[j * NPAR + m] += J[a * NPAR + j] * Cxx[a * 4 + c] * J[c * NPAR + m];
							}
						}
					}
				}
			}
		}
		if (par[P_CEIL] > 0.5) {  // CEILING diagnostic: theta and z0 pinned to the target carried by the cbar == 2 pseudo-candidate
			for (int c = 0; c < st.ncand; c++) {
				if (st.cbar[c] == 2) {
					double r1 = (u[4] - st.cz[c]) / CEIL_STH, r2 = (u[3] - st.cphi[c]) / CEIL_SZ;
					fv += r1 * r1 + r2 * r2;
					if (want) {
						g[4] += 2.0 * r1 / CEIL_STH; hc[4 * NPAR + 4] += 2.0 / (CEIL_STH * CEIL_STH); he[4 * NPAR + 4] += 2.0 / (CEIL_STH * CEIL_STH);
						g[3] += 2.0 * r2 / CEIL_SZ; hc[3 * NPAR + 3] += 2.0 / (CEIL_SZ * CEIL_SZ); he[3 * NPAR + 3] += 2.0 / (CEIL_SZ * CEIL_SZ);
					}
				}
			}
		}
		if (st.fb == 1 && st.fnone == 1) {  // electron-vz prior fallback for an ATOF_NONE track
			double sz = st.sigz;
			double r = (u[3] - ze) / sz;
			fv += r * r;
			if (want) {
				g[3] += 2.0 * r / sz; hc[3 * NPAR + 3] += 2.0 / (sz * sz); he[3 * NPAR + 3] += 2.0 / (sz * sz);
			}
		}
		return fv;
	}

	// ================================================================== attach

	/** numpy.mod for a positive divisor. */
	private static double npMod(double a, double b) {
		double m = a % b;
		if (m != 0.0 && m < 0.0) m += b;
		return m;
	}

	/**
	 * Candidate attach of one track (atof.attach, atof.py:387-450; seed_intercepts atof.py:102-123): propagates the B0 seed
	 * (st.record = intercept record), takes every ATOF row of the event within the phi window of the seed intercept (the last
	 * node's phi for a layer the seed does not reach), the ATOF_NONE flag, the bar z recalibration z_cal = (z - offset) /
	 * scale (recalOffset/recalScale: 60 values of the track's half, null = none) and, when par[P_TASSOC] is on, the frozen
	 * wedge weight 2 t_win pdf_t(t_c - t_pred(seed) - off) (1 for bar rows, an unreached seed layer or a non-finite time);
	 * otherwise ct is the time - start. bankTimes selects the bank times instead of the final ones.
	 */
	public static void attach(GFitTrack t, GFitEvent ev, double[] par, double[] recalOffset, double[] recalScale,
			boolean bankTimes, Propagator prop, GFitWorkspace ws, TrackState st) {
		st.par = par;
		final double[] rec = st.record;
		for (int q = 0; q < NI; q++) rec[q] = 0.0;
		double pl = 0.0;
		double p = t.seedP;
		if (p > 0.03 && Double.isFinite(t.seedPhi) && Double.isFinite(t.seedTheta) && Double.isFinite(t.seedZ)) {
			int nn = prop.runSeed(t.seedPhi, p, t.seedTheta, t.seedZ, ws);
			Propagator.atofIntercepts(ws.out, nn, ws.regnode, ws.aux, MaterialStack.dataBest().rOut, rec, 0);
			pl = Math.atan2(ws.out[(nn - 1) * NC + 2], ws.out[(nn - 1) * NC + 1]);
		}
		final int na = ev.nAtof();
		final double[] times = bankTimes ? ev.tBank : ev.tFinal;
		final double bwin = WIN_ROWS * ROW, wwin = WIN_ROWS * ROW;
		final double c0bar = (((int) rec[7]) & 1) != 0 ? rec[0] : pl;
		final double c0wed = (((int) rec[15]) & 1) != 0 ? rec[8] : pl;
		final boolean tassoc = par[P_TASSOC] > 0.5;
		final boolean reached = (((int) rec[15]) & 1) != 0;
		final double tp = par[P_TMODE] > 0.5 ? rec[11] : rec[10];
		int none = 1, n = 0;
		for (int r = 0; r < na; r++) {
			double ph = Math.atan2(ev.ay[r], ev.ax[r]) + 0.0;
			ph = npMod(ph + Math.PI, TWOPI) - Math.PI;
			boolean bar = ev.isBar[r];
			double cen = bar ? c0bar : c0wed;
			double dph = Math.abs(npMod(ph - cen + Math.PI, TWOPI) - Math.PI);
			if (dph <= NONE_ROWS * ROW + 1e-9) none = 0;
			if (!(dph <= (bar ? bwin : wwin) + 1e-9)) continue;
			st.grow(n + 1);
			double zz = ev.az[r] + 0.0;
			if (recalOffset != null && bar) {
				int bi = Math.min(Math.max(4 * ev.sector[r] + ev.layer[r], 0), 59);
				zz = (zz - recalOffset[bi]) / recalScale[bi];
			}
			double time = times[r] - ev.startTime;
			double w = time;
			if (tassoc) {
				if (bar) w = 1.0;
				else if (reached) {
					w = Double.isFinite(time) ? pdfTime(time - tp, par, st) * 2.0 * par[P_TWIN] : 1.0;
				} else w = 1.0;
			}
			st.cz[n] = zz; st.cphi[n] = ph; st.ct[n] = w; st.cbar[n] = bar ? 1 : 0;
			st.csec[n] = ev.sector[r]; st.ccomp[n] = ev.comp[r]; st.crow[n] = r;
			n++;
		}
		st.ncand = n;
		st.fnone = none;
	}

	private static double pdfTime(double r, double[] par, TrackState st) {
		tpdf(r, par, st, st.tp3);
		return st.tp3[0];
	}
}
