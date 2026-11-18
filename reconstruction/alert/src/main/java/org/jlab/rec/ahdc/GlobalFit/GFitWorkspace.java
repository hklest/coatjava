package org.jlab.rec.ahdc.GlobalFit;

/**
 * Scratch arrays of one fit thread; nothing is allocated inside a fit. One instance per thread, never shared. Port of the
 * workspace tuple of the pinned fit (fit.py:609-633), grown unit by unit: each block names the class that uses it.
 */
public final class GFitWorkspace {

	/** node columns (s x y z ux uy uz p beta t reg xx0 th2 varp), maximum nodes of one propagation, planes, kinks, parameters. */
	public static final int NC = 14, NMAX = 600, NPL = 6, NK = 10, NPAR = 15;
	/** capacity of per-region arrays (the data_best stack has 15 regions). */
	static final int NREG_MAX = 32;

	// ---- Propagator ----

	/** when true, run() saves the D8 snapshots at the kink planes 0-4 and the nominal results (see Propagator). */
	public boolean saveSnapshots;
	/** node table of the last run() or runRestart(), flat [i * NC + c], valid rows [0, nn). */
	public final double[] out = new double[NMAX * NC];
	/** node table of the last run() made with saveSnapshots (the nominal run); runRestart() copies its prefix into out. */
	public final double[] outNom = new double[NMAX * NC];
	/** node count and trajectory flags (1 left the map, 2 no region exit, 4 stopped) of the last run or restart. */
	public int nn, flags;
	/** node index at entry of each region (-1 never entered), entry nreg = exit node; length NREG_MAX, valid [0, nreg]. */
	public final int[] regnode = new int[NREG_MAX];
	/** (stopped, s_stop, r_stop, 0). */
	public final double[] aux = new double[4];
	/** Highland theta0 per scattering plane. */
	public final double[] theta0 = new double[NPL];

	/** snapshots of the nominal run at kink planes 0-4, and its final results. */
	final Propagator.State[] snap = new Propagator.State[5];
	int nomNn, nomFlags;
	final int[] nomRegnode = new int[NREG_MAX];
	final double[] nomAux = new double[4];
	final double[] nomTheta0 = new double[NPL];
	/** working state and scratch of the core loop. */
	final Propagator.State cur = new Propagator.State();
	final double[] kin = new double[2 * NPL];
	final double[] b3 = new double[3];
	final double[] st = new double[8], ns = new double[8], k1 = new double[8], k2 = new double[8], k3 = new double[8], k4 = new double[8],
		tmp = new double[8];
	final double[] thseg = new double[NPL], xxseg = new double[NPL], betaseg = new double[NPL];

	public GFitWorkspace(int maxHits) {
		for (int j = 0; j < snap.length; j++) snap[j] = new Propagator.State();
	}
}
