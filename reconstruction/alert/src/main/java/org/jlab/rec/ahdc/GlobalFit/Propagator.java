package org.jlab.rec.ahdc.GlobalFit;

/**
 * Proton propagation through the region stack in the solenoid map: port of prop._core with the proton hypothesis of
 * prop._hy (mass 0.938272088 GeV, z 1, tsc 1, zr2 1) and the propagation part of fit._prop_run (pinned prop.py:305-513,
 * 539-541; fit.py:310-321), plus the ATOF intercepts of atof.intercepts_nb (atof.py:45-99). RK4 steps with mean energy
 * loss, CSDA range for solids, a Newton correction onto the last cylinder of each region, Highland variances per scattering
 * plane and direction kinks at the planes. The arithmetic follows the Python operation by operation.
 *
 * <p>Snapshot and restart (plan D8): with {@code ws.saveSnapshots} set, {@link #run} saves, at each kink plane 0-4, the
 * complete internal state just before the kink rotation, and copies the nominal node table to {@code ws.outNom} with the
 * nominal results. {@link #runRestart} restores the snapshot of a plane, copies the nominal node prefix and executes the
 * same loop code from the kink on, so for a state that differs from the nominal one only in the kinks of that plane (and
 * later ones) it gives bit-identical results to {@link #run}. Buffers: {@code ws.out}, {@code ws.nn}, {@code ws.flags},
 * {@code ws.regnode}, {@code ws.aux}, {@code ws.theta0} hold the result of the last run or restart; {@code ws.outNom}
 * and the private nominal copies hold the last saving run.
 *
 * <p>Immutable; thread safe with one {@link GFitWorkspace} per thread.
 */
public final class Propagator {

	static final double KAPPA = 0.299792458e-3;  // prop.py:32 (GeV/c per T mm)
	static final double CMMNS = 299.792458;  // prop.py:33
	static final double MASS = 0.938272088, ZP = 1.0, TSC = 0.938272088 / 0.938272088, ZR2 = 1.0;  // prop._hy, prop.py:539-541
	static final double T_FLOOR = 0.3;  // MeV, make_cfg (prop.py:257)
	private static final double M2 = MASS * MASS;
	private static final int NC = GFitWorkspace.NC;

	/** Internal state of the core loop at a kink plane (the done state is always false there). */
	public static final class State {
		boolean valid;
		double x, y, z, ux, uy, uz, p, t, s, xx0, th2, varp;
		final double[] thseg = new double[GFitWorkspace.NPL], xxseg = new double[GFitWorkspace.NPL], betaseg = new double[GFitWorkspace.NPL];
		/** region index k of the kink, node count before the kink rotation, flags. */
		int k, nn, flags;
		final int[] regnode = new int[GFitWorkspace.NREG_MAX];
	}

	private final FieldGrid fld;
	private final MaterialStack stk;
	private final StoppingTables tab;
	private final double beamX, beamY;

	public Propagator(FieldGrid f, MaterialStack s, StoppingTables t, double beamX, double beamY) {
		this.fld = f;
		this.stk = s;
		this.tab = t;
		this.beamX = beamX;
		this.beamY = beamY;
	}

	/**
	 * Propagation of state u = (d0, phi0, ln p, z0, theta, k0..k9) with kinks u[5..14] (kink j on plane j/2, component j%2;
	 * plane 5 has none), as fit._prop_run (fit.py:310-321). Fills the nodes in ws.out, ws.nn, ws.flags, ws.regnode, ws.aux,
	 * ws.theta0; returns the node count. With ws.saveSnapshots it also saves the D8 snapshots and the nominal results.
	 */
	public int run(double[] u, GFitWorkspace ws) {
		final double[] kin = ws.kin;
		for (int j = 0; j < GFitWorkspace.NPL; j++) {
			kin[2 * j] = 0.0;
			kin[2 * j + 1] = 0.0;
		}
		for (int j = 0; j < GFitWorkspace.NK; j++) kin[(j / 2) * 2 + j % 2] = u[5 + j];
		double d0 = u[0], ph = u[1], z0 = u[3], th = u[4];
		double x = beamX - d0 * Math.sin(ph), y = beamY + d0 * Math.cos(ph), z = z0;
		double ux = Math.sin(th) * Math.cos(ph), uy = Math.sin(th) * Math.sin(ph), uz = Math.cos(th);
		double p = Math.exp(u[2]);
		// prop._core prologue (prop.py:366-391)
		final int nreg = stk.nreg;
		final int[] regnode = ws.regnode;
		for (int k = 0; k <= nreg; k++) regnode[k] = -1;
		double e = Math.sqrt(p * p + M2);
		final double[] out = ws.out;
		out[0] = 0.0; out[1] = x; out[2] = y; out[3] = z; out[4] = ux; out[5] = uy; out[6] = uz;
		out[7] = p; out[8] = p / e; out[9] = 0.0; out[10] = -1.0; out[11] = 0.0; out[12] = 0.0; out[13] = 0.0;
		double rs = Math.sqrt(x * x + y * y);
		int k0 = 0;
		while (k0 < nreg - 1 && stk.rOut[k0] <= rs) k0++;
		ws.aux[0] = 0.0; ws.aux[1] = -1.0; ws.aux[2] = -1.0; ws.aux[3] = 0.0;
		State c = ws.cur;
		c.x = x; c.y = y; c.z = z; c.ux = ux; c.uy = uy; c.uz = uz; c.p = p;
		c.t = 0.0; c.s = 0.0; c.xx0 = 0.0; c.th2 = 0.0; c.varp = 0.0;
		for (int j = 0; j < GFitWorkspace.NPL; j++) {
			c.thseg[j] = 0.0;
			c.xxseg[j] = 0.0;
			c.betaseg[j] = 0.0;
		}
		c.k = k0;
		c.nn = 1;
		c.flags = 0;
		boolean save = ws.saveSnapshots;
		if (save) for (State sn : ws.snap) sn.valid = false;
		core(ws, c, false, save);
		if (save) {
			System.arraycopy(out, 0, ws.outNom, 0, ws.nn * NC);
			ws.nomNn = ws.nn;
			ws.nomFlags = ws.flags;
			System.arraycopy(regnode, 0, ws.nomRegnode, 0, regnode.length);
			System.arraycopy(ws.aux, 0, ws.nomAux, 0, 4);
			System.arraycopy(ws.theta0, 0, ws.nomTheta0, 0, ws.theta0.length);
		}
		return ws.nn;
	}

	/**
	 * Restart from the snapshot of kink plane {@code plane} (0-4) of the last saving run: same result as run(u, ws)
	 * (bit-identical) when u differs from the snapshotted state only in the kinks of that plane and later planes. A
	 * plane the nominal run never reached has no snapshot; the nominal result is returned then. Copies the nominal node
	 * prefix into ws.out. Needs a run(..) with ws.saveSnapshots first.
	 */
	public int runRestart(double[] u, int plane, GFitWorkspace ws) {
		State sn = ws.snap[plane];
		final double[] out = ws.out;
		if (!sn.valid) {
			System.arraycopy(ws.outNom, 0, out, 0, ws.nomNn * NC);
			ws.nn = ws.nomNn;
			ws.flags = ws.nomFlags;
			System.arraycopy(ws.nomRegnode, 0, ws.regnode, 0, ws.regnode.length);
			System.arraycopy(ws.nomAux, 0, ws.aux, 0, 4);
			System.arraycopy(ws.nomTheta0, 0, ws.theta0, 0, ws.theta0.length);
			return ws.nn;
		}
		final double[] kin = ws.kin;
		for (int j = 0; j < GFitWorkspace.NPL; j++) {
			kin[2 * j] = 0.0;
			kin[2 * j + 1] = 0.0;
		}
		for (int j = 0; j < GFitWorkspace.NK; j++) kin[(j / 2) * 2 + j % 2] = u[5 + j];
		System.arraycopy(ws.outNom, 0, out, 0, sn.nn * NC);
		System.arraycopy(sn.regnode, 0, ws.regnode, 0, ws.regnode.length);
		ws.aux[0] = 0.0; ws.aux[1] = -1.0; ws.aux[2] = -1.0; ws.aux[3] = 0.0;
		core(ws, sn, true, false);
		return ws.nn;
	}

	/**
	 * Propagation of the B0 seed state as atof.seed_intercepts (atof.py:110-122): position (beam, z), direction from
	 * (theta, phi), momentum p exactly, d0 = 0 and all kinks zero (duplicate kink nodes are appended). Never saves snapshots.
	 * Fills the same workspace fields as run(); returns the node count.
	 */
	/** the material stack this propagator steps through. */
	public MaterialStack stack() { return stk; }

	public int runSeed(double phi, double p, double theta, double z, GFitWorkspace ws) {
		final double[] kin = ws.kin;
		for (int j = 0; j < 2 * GFitWorkspace.NPL; j++) kin[j] = 0.0;
		double x = beamX, y = beamY;
		double ux = Math.sin(theta) * Math.cos(phi), uy = Math.sin(theta) * Math.sin(phi), uz = Math.cos(theta);
		final int nreg = stk.nreg;
		final int[] regnode = ws.regnode;
		for (int k = 0; k <= nreg; k++) regnode[k] = -1;
		double e = Math.sqrt(p * p + M2);
		final double[] out = ws.out;
		out[0] = 0.0; out[1] = x; out[2] = y; out[3] = z; out[4] = ux; out[5] = uy; out[6] = uz;
		out[7] = p; out[8] = p / e; out[9] = 0.0; out[10] = -1.0; out[11] = 0.0; out[12] = 0.0; out[13] = 0.0;
		double rs = Math.sqrt(x * x + y * y);
		int k0 = 0;
		while (k0 < nreg - 1 && stk.rOut[k0] <= rs) k0++;
		ws.aux[0] = 0.0; ws.aux[1] = -1.0; ws.aux[2] = -1.0; ws.aux[3] = 0.0;
		State c = ws.cur;
		c.x = x; c.y = y; c.z = z; c.ux = ux; c.uy = uy; c.uz = uz; c.p = p;
		c.t = 0.0; c.s = 0.0; c.xx0 = 0.0; c.th2 = 0.0; c.varp = 0.0;
		for (int j = 0; j < GFitWorkspace.NPL; j++) {
			c.thseg[j] = 0.0;
			c.xxseg[j] = 0.0;
			c.betaseg[j] = 0.0;
		}
		c.k = k0;
		c.nn = 1;
		c.flags = 0;
		core(ws, c, false, false);
		return ws.nn;
	}

	/**
	 * The loop of prop._core (prop.py:392-513) from state {@code from} (region from.k; with resume, directly at the kink of
	 * region from.k). Reads ws.kin, writes ws.out/nn/flags/regnode/aux/theta0.
	 */
	private void core(GFitWorkspace ws, State from, boolean resume, boolean save) {
		final int nreg = stk.nreg;
		final double[] rOut = stk.rOut, rhoR = stk.rho, xk = stk.xk, zoarho = stk.zoarho;
		final int[] matR = stk.mat, nstR = stk.nst, segR = stk.seg, pla = stk.pla;
		final double[] out = ws.out, kin = ws.kin, b = ws.b3, st = ws.st, ns = ws.ns, aux = ws.aux;
		final double[] thseg = ws.thseg, xxseg = ws.xxseg, betaseg = ws.betaseg;
		final int[] regnode = ws.regnode;
		double x = from.x, y = from.y, z = from.z, ux = from.ux, uy = from.uy, uz = from.uz, p = from.p, t = from.t, s = from.s;
		double xx0 = from.xx0, th2 = from.th2, varp = from.varp;
		int nn = from.nn, flags = from.flags;
		System.arraycopy(from.thseg, 0, thseg, 0, thseg.length);
		System.arraycopy(from.xxseg, 0, xxseg, 0, xxseg.length);
		System.arraycopy(from.betaseg, 0, betaseg, 0, betaseg.length);
		final int kStart = from.k;
		final int nmax = GFitWorkspace.NMAX - 2;
		final double pfl = Math.sqrt((MASS + T_FLOOR * 1e-3) * (MASS + T_FLOOR * 1e-3) - M2);
		boolean done = false;
		double e;
		boolean skip = resume;
		outer:
		for (int k = kStart; k < nreg; k++) {
			if (!skip) {
				regnode[k] = nn - 1;
				final double rho = rhoR[k];
				final int mat = matR[k], sg = segR[k];
				fld.field(x, y, z, b);
				double sest = sToCyl(x, y, ux, uy, uz, KAPPA / p * b[2], rOut[k]);
				if (sest < 0.0) {
					flags |= 2;
					break outer;
				}
				final int n = nstR[k];
				final double h = sest / n;
				for (int i = 0; i < n; i++) {
					if (nn >= nmax) {
						flags |= 2;
						done = true;
						break;
					}
					double tmv = (Math.sqrt(p * p + M2) - MASS) * 1e3;
					double lr = rho > 0.0 ? 10.0 * tab.rangeG(mat, tmv * TSC) / (TSC * rho) : 1e30;
					if (lr <= h) {  // stops inside this step: place the stop point at the range, no RK4 step across the Bragg peak
						double xn = x + ux * lr, yn = y + uy * lr, zn = z + uz * lr;
						double bet0 = p / Math.sqrt(p * p + M2);
						double tn = t + lr / (0.5 * bet0 * CMMNS);
						double dx = xk[k] * lr;
						double bp = bet0 * p;
						xx0 += dx;
						th2 += dx / (bp * bp);
						thseg[sg] += dx / (bp * bp);
						xxseg[sg] += dx;
						betaseg[sg] = bet0;
						s += lr;
						x = xn; y = yn; z = zn; p = pfl; t = tn;
						int o = nn * NC;
						out[o] = s; out[o + 1] = x; out[o + 2] = y; out[o + 3] = z; out[o + 4] = ux; out[o + 5] = uy; out[o + 6] = uz;
						out[o + 7] = p; out[o + 8] = pfl / Math.sqrt(pfl * pfl + M2); out[o + 9] = t; out[o + 10] = (double) k;
						out[o + 11] = xx0; out[o + 12] = th2; out[o + 13] = varp;
						nn++;
						aux[0] = 1.0; aux[1] = s; aux[2] = Math.sqrt(x * x + y * y);
						flags |= 4;
						done = true;
						break;
					}
					st[0] = x; st[1] = y; st[2] = z; st[3] = ux; st[4] = uy; st[5] = uz; st[6] = p; st[7] = t;
					rk4(ws, st, h, mat, rho, ns);
					double xn = ns[0], yn = ns[1], zn = ns[2], uxn = ns[3], uyn = ns[4], uzn = ns[5], pn = ns[6], tn = ns[7];
					double hEff;
					if (rho > 0.1) {  // solids: energy from the CSDA table over the step
						double lremS = 10.0 * tab.rangeG(mat, (Math.sqrt(p * p + M2) - MASS) * 1e3 * TSC) / (TSC * rho);
						double tref = tab.tFromRange(mat, Math.max(lremS - h, 1e-9) * TSC * rho / 10.0);
						double enN = MASS + tref / TSC * 1e-3;
						pn = Math.sqrt(Math.max(enN * enN - M2, 1e-30));
					}
					if (i == n - 1) {  // one Newton correction onto the exact cylinder
						double rr = Math.sqrt(xn * xn + yn * yn);
						double drds = (xn * uxn + yn * uyn) / rr;
						if (drds < 1e-6) {
							flags |= 2;
							done = true;
							break;
						}
						double ds = (rOut[k] - rr) / drds;
						fld.field(xn, yn, zn, b);
						double bx = b[0], by = b[1], bz = b[2];
						double c = KAPPA / pn;
						double ax = c * (uyn * bz - uzn * by), ay = c * (uzn * bx - uxn * bz), az = c * (uxn * by - uyn * bx);
						xn += uxn * ds + 0.5 * ax * ds * ds;
						yn += uyn * ds + 0.5 * ay * ds * ds;
						zn += uzn * ds + 0.5 * az * ds * ds;
						uxn += ax * ds;
						uyn += ay * ds;
						uzn += az * ds;
						double en = Math.sqrt(pn * pn + M2);
						pn -= rho * tab.stop(mat, Math.max((en - MASS) * 1e3 * TSC, 1.0001e-3)) * ZR2 * 1e-4 * en / pn * ds;
						tn += en / (pn * CMMNS) * ds;
						hEff = h + ds;
					} else {
						hEff = h;
					}
					double un = Math.sqrt(uxn * uxn + uyn * uyn + uzn * uzn);
					uxn = uxn / un; uyn = uyn / un; uzn = uzn / un;
					double ts = (p + pn) * 0.5;
					double en2 = Math.sqrt(ts * ts + M2);
					double bet = ts / en2;
					if (pn <= pfl || !(pn > 0.0)) {  // stopped within this step: remaining CSDA range
						double tmev = (Math.sqrt(p * p + M2) - MASS) * 1e3 * TSC;
						double lrem = 10.0 * tab.rangeG(mat, tmev) / (TSC * Math.max(rho, 1e-12));
						lrem = Math.min(lrem, hEff);
						xn = x + ux * lrem; yn = y + uy * lrem; zn = z + uz * lrem;
						uxn = ux; uyn = uy; uzn = uz; pn = pfl;
						tn = t + lrem / (0.5 * p / Math.sqrt(p * p + M2) * CMMNS);
						hEff = lrem;
						aux[0] = 1.0; aux[1] = s + lrem; aux[2] = Math.sqrt(xn * xn + yn * yn);
						flags |= 4;
						done = true;
					}
					double dx = xk[k] * hEff;
					double bs = bet * ts;
					xx0 += dx;
					th2 += dx / (bs * bs);
					thseg[sg] += dx / (bs * bs);
					xxseg[sg] += dx;
					betaseg[sg] = bet;
					varp += 0.1569e-6 * ZP * ZP * zoarho[k] * 0.1 * hEff * (1.0 - 0.5 * bet * bet) / (1.0 - bet * bet) / (bet * bet);
					x = xn; y = yn; z = zn; ux = uxn; uy = uyn; uz = uzn; p = pn; t = tn;
					s += hEff;
					e = Math.sqrt(p * p + M2);
					int o = nn * NC;
					out[o] = s; out[o + 1] = x; out[o + 2] = y; out[o + 3] = z; out[o + 4] = ux; out[o + 5] = uy; out[o + 6] = uz;
					out[o + 7] = p; out[o + 8] = p / e; out[o + 9] = t; out[o + 10] = (double) k;
					out[o + 11] = xx0; out[o + 12] = th2; out[o + 13] = varp;
					nn++;
					if (done) break;
				}
				if (done) break;
			}
			skip = false;
			if (pla[k] >= 0) {  // kink: rotate the direction by two projected angles (e1 azimuthal, e2 = u x e1)
				int pl = pla[k];
				if (save && pl < 5) {
					State sn = ws.snap[pl];
					sn.valid = true;
					sn.x = x; sn.y = y; sn.z = z; sn.ux = ux; sn.uy = uy; sn.uz = uz; sn.p = p; sn.t = t; sn.s = s;
					sn.xx0 = xx0; sn.th2 = th2; sn.varp = varp;
					System.arraycopy(thseg, 0, sn.thseg, 0, thseg.length);
					System.arraycopy(xxseg, 0, sn.xxseg, 0, xxseg.length);
					System.arraycopy(betaseg, 0, sn.betaseg, 0, betaseg.length);
					sn.k = k;
					sn.nn = nn;
					sn.flags = flags;
					System.arraycopy(regnode, 0, sn.regnode, 0, regnode.length);
				}
				double k1 = kin[2 * pl], k2 = kin[2 * pl + 1];
				double ut = Math.sqrt(ux * ux + uy * uy);
				if (ut > 1e-9) {
					double e1x = -uy / ut, e1y = ux / ut, e1z = 0.0;
					double e2x = uy * e1z - uz * e1y, e2y = uz * e1x - ux * e1z, e2z = ux * e1y - uy * e1x;
					double nx = ux + k1 * e1x + k2 * e2x, ny = uy + k1 * e1y + k2 * e2y, nz = uz + k1 * e1z + k2 * e2z;
					double un = Math.sqrt(nx * nx + ny * ny + nz * nz);
					ux = nx / un; uy = ny / un; uz = nz / un;
					int o = nn * NC;
					System.arraycopy(out, o - NC, out, o, NC);
					out[o + 4] = ux; out[o + 5] = uy; out[o + 6] = uz;
					nn++;
				}
			}
		}
		if (!done) regnode[nreg] = nn - 1;
		final double[] theta0 = ws.theta0;
		for (int j = 0; j < theta0.length; j++) {
			if (thseg[j] > 0.0) {
				double bj = betaseg[j];
				theta0[j] = 0.0136 * ZP * Math.sqrt(thseg[j]) * (1.0 + 0.038 * Math.log(Math.max(xxseg[j] * ZP * ZP / (bj * bj), 1e-12)));
			} else {
				theta0[j] = 0.0;
			}
		}
		if (x * x + y * y > 1e6 || Math.abs(z + fld.zsh) > 400.0) flags |= 1;
		ws.nn = nn;
		ws.flags = flags;
	}

	/** prop._deriv (prop.py:305-312): right-hand side of the equations of motion for st = (x y z ux uy uz p t). */
	private void deriv(double[] st, double[] k, int mat, double rho, double[] b) {
		fld.field(st[0], st[1], st[2], b);
		double bx = b[0], by = b[1], bz = b[2];
		double ux = st[3], uy = st[4], uz = st[5], p = st[6];
		double c = KAPPA / p;
		double e = Math.sqrt(p * p + M2);
		double s = tab.stop(mat, Math.max((e - MASS) * 1e3 * TSC, 1.0001e-3)) * ZR2;
		k[0] = ux; k[1] = uy; k[2] = uz;
		k[3] = c * (uy * bz - uz * by);
		k[4] = c * (uz * bx - ux * bz);
		k[5] = c * (ux * by - uy * bx);
		k[6] = -rho * s * 1e-4 * e / p;
		k[7] = e / (p * CMMNS);
	}

	/** prop._rk4 (prop.py:315-331): one classical RK4 step of length h into r. */
	private void rk4(GFitWorkspace ws, double[] st, double h, int mat, double rho, double[] r) {
		final double[] k1 = ws.k1, k2 = ws.k2, k3 = ws.k3, k4 = ws.k4, tmp = ws.tmp, b = ws.b3;
		deriv(st, k1, mat, rho, b);
		for (int j = 0; j < 8; j++) tmp[j] = st[j] + 0.5 * h * k1[j];
		deriv(tmp, k2, mat, rho, b);
		for (int j = 0; j < 8; j++) tmp[j] = st[j] + 0.5 * h * k2[j];
		deriv(tmp, k3, mat, rho, b);
		for (int j = 0; j < 8; j++) tmp[j] = st[j] + h * k3[j];
		deriv(tmp, k4, mat, rho, b);
		double n = h / 6.0;
		for (int j = 0; j < 8; j++) r[j] = st[j] + n * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
	}

	/** prop._s_to_cyl (prop.py:334-360): path length (mm) along a transverse circle to the cylinder r = rr; -1 if never reached. */
	private static double sToCyl(double x, double y, double ux, double uy, double uz, double cBz, double rr) {
		double ut = Math.sqrt(ux * ux + uy * uy);
		if (ut < 1e-9) return -1.0;
		double phi0 = Math.atan2(uy, ux);
		double om = -cBz / ut;
		double tx = ux / ut, ty = uy / ut;
		double b = x * tx + y * ty;
		double rad = b * b - (x * x + y * y - rr * rr);
		if (rad < 0) return -1.0;
		double sg = -b + Math.sqrt(rad);
		for (int it = 0; it < 4; it++) {
			double px, py, dx, dy;
			if (Math.abs(om * sg) < 1e-8) {
				px = x + tx * sg; py = y + ty * sg; dx = tx; dy = ty;
			} else {
				double a = phi0 + om * sg;
				px = x + (Math.sin(a) - Math.sin(phi0)) / om;
				py = y - (Math.cos(a) - Math.cos(phi0)) / om;
				dx = Math.cos(a);
				dy = Math.sin(a);
			}
			double f = px * px + py * py - rr * rr;
			double fp = 2.0 * (px * dx + py * dy);
			if (fp <= 1e-12) return -1.0;
			sg -= f / fp;
		}
		return sg / ut;
	}

	// ---- ATOF intercepts ----

	static final int NI = 20;
	private static final double BAR_HALF = 139.85, WEDGE_HALF = 5 * 27.7, ZACC = 5.0, SRCH = 2.0;

	private static double sstep(double x) {
		if (x <= 0.0) return 0.0;
		if (x >= 1.0) return 1.0;
		return x * x * (3.0 - 2.0 * x);
	}

	/**
	 * atof.intercepts_nb (atof.py:54-99): rec[recOff + 0..8) bar, +8..16 wedge = (phi, z, t_entry, t_mid, beta, w,
	 * s_entry, flag), +16..19 = mid-depth phi, z of bar and wedge. out is the flat node table [i * NC + c]. All NI = 20 slots
	 * are cleared first (Python receives a zeroed record). regnode has the entries [0, nreg].
	 */
	public static void atofIntercepts(double[] out, int nn, int[] regnode, double[] aux, double[] rOut, double[] rec, int recOff) {
		final int nreg = rOut.length;
		for (int q = 0; q < NI; q++) rec[recOff + q] = 0.0;
		int iin = 1;
		for (int i = 0; i < nreg; i++) {
			if (Math.abs(rOut[i] - 77.0) < 1e-9) {
				iin = i + 1;
				break;
			}
		}
		boolean stopped = aux[0] > 0.5;
		double sStop = aux[1], rStop = aux[2];
		for (int L = 0; L < 2; L++) {
			int k0 = L == 0 ? iin : iin + 2;
			int k1 = L == 0 ? iin + 1 : iin + 2;
			int b = recOff + 8 * L;
			for (int q = 0; q < 8; q++) rec[b + q] = 0.0;
			if (k1 >= nreg) continue;
			int n0 = regnode[k0];
			int n1 = k1 + 1 < nreg + 1 ? regnode[k1 + 1] : -1;  // regnode.shape[0] = nreg + 1
			if (n0 < 0 || n0 >= nn) continue;
			if (n1 < 0 || n1 <= n0 || n1 >= nn) n1 = nn - 1;
			double x = out[n0 * NC + 1], y = out[n0 * NC + 2], z = out[n0 * NC + 3];
			if (!(Double.isFinite(x) && Double.isFinite(y) && Double.isFinite(z) && Double.isFinite(out[n0 * NC + 9])
				&& Double.isFinite(out[n1 * NC + 9]))) continue;
			rec[b] = Math.atan2(y, x);
			rec[b + 1] = z;
			rec[b + 2] = out[n0 * NC + 9];
			rec[b + 3] = 0.5 * (out[n0 * NC + 9] + out[n1 * NC + 9]);
			rec[b + 4] = out[n0 * NC + 8];
			rec[b + 6] = out[n0 * NC];
			double fl = 1.0;
			if (stopped && rStop >= rOut[k0 - 1] - 1e-6 && rStop <= rOut[k1] + 1e-6) fl += 2.0;
			rec[b + 7] = fl;
			double sm = 0.5 * (out[n0 * NC] + out[n1 * NC]);  // mid-depth point (arc-length midpoint of entry and exit/stop node)
			int im = n0;
			while (im < n1 - 1 && out[(im + 1) * NC] < sm) im++;
			double ds = n1 > n0 ? out[(im + 1) * NC] - out[im * NC] : 0.0;
			double fr = ds > 1e-12 ? Math.min(Math.max((sm - out[im * NC]) / ds, 0.0), 1.0) : 0.0;
			int i2 = Math.min(im + 1, nn - 1);
			rec[recOff + 16 + 2 * L] = Math.atan2(out[im * NC + 2] + fr * (out[i2 * NC + 2] - out[im * NC + 2]),
				out[im * NC + 1] + fr * (out[i2 * NC + 1] - out[im * NC + 1]));
			rec[recOff + 17 + 2 * L] = out[im * NC + 3] + fr * (out[i2 * NC + 3] - out[im * NC + 3]);
			double wr = stopped ? sstep((sStop - out[n0 * NC]) / SRCH) : 1.0;
			double half = L == 0 ? BAR_HALF : WEDGE_HALF;
			rec[b + 5] = wr * sstep((half - Math.abs(z)) / ZACC);
		}
	}
}
