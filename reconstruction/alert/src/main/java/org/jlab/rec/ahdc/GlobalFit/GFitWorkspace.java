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


	// ---- TrackCost (hit list, prediction and cost scratch; users: TrackCost, later GlobalTrackFitter) ----

	/** capacity of the hit arrays: maxHits + NADD + NMISS (fit._fit_one maxn, fit.py:603). */
	public final int cap;
	/** active hit list: count, global wire, layer index, wire in layer, doca, smear class, sigma, half pitch, layer class, kind (0 hit), warm arc length, row in event. */
	public int nh;
	public final int[] lg, lli, lw, lism, lcls, lkind, lrow;
	public final double[] lm, lsig, lh, sarc;
	/** per hit prediction: signed d, cell margin s, flag (1 = no valid prediction), Jacobians [k * NPAR + j]. */
	public final double[] D, S, JD, JS;
	public final int[] PF;
	/** per hit (l, d, s) of the last cost(): outH[k * 3 + {0, 1, 2}]. */
	public final double[] outH;
	/** Highland theta0 of the nominal propagation of the last predict (fit TH0); trajectory flags of the last predict. */
	public final double[] th0 = new double[NPL];
	public int tflag;
	/** gradient, clipped Gauss-Newton and exact Hessians (full symmetric, NPAR x NPAR) of the last cost(want). */
	public final double[] g = new double[NPAR], hc = new double[NPAR * NPAR], he = new double[NPAR * NPAR];
	/** terms of the last cost(): hit, z prior, d0 prior, kink prior. */
	public double fHit, fZ, fD0, fK;
	/** when true, predict() records the ATOF intercepts of the nominal and the 15 perturbed runs in ia[16 * NI] (set by the ATOF user). */
	public boolean atofOn;
	public final double[] ia = new double[16 * 20];
	/** ATOF candidates and scratch of the current track (AtofTerm.attach) and the ATOF cost of the last cost(). */
	public final AtofTerm.TrackState atof = new AtofTerm.TrackState();
	public double fAtof;
	/** scratch of predict and cost. */
	final int[] PF1;
	final double[] pcaG, pcaD, pcaR = new double[12], nomX, pertD, pertX, gradG, upert = new double[NPAR];
	final double[] box3 = new double[3], cell4 = new double[4], o9 = new double[9], o9b = new double[9], jx5 = new double[5],
		j53 = new double[15], sArc1 = new double[1];

	// ---- LevMar and GlobalTrackFitter ----

	/** number of grid points of the start search (8 x 13 x 5) and of the best starts kept. */
	static final int NGRID = 8 * 13 * 5, NBEST = 3;
	/** LM outputs of the last LevMar.run: converged flag and iterations. */
	public boolean lmConv;
	public int lmIt;
	/** LM scratch: gradient, matrix, damped matrix, step, trial and probe states, pivots. */
	final double[] lmGs = new double[NPAR], lmHs = new double[NPAR * NPAR], lmA = new double[NPAR * NPAR], lmDx = new double[NPAR],
		lmUn = new double[NPAR], lmUu = new double[NPAR];
	final int[] piv = new int[NPAR];
	/** grid scratch: costs, states (5 each), sort order and merge buffer, best starts. */
	final double[] gridF = new double[NGRID], gridU = new double[NGRID * 5], gridY = new double[NPAR], bestU = new double[NBEST * 5],
		bestF = new double[NBEST];
	final int[] gridOrd = new int[NGRID], gridTmp = new int[NGRID];
	/** starts (6 x NPAR, kind), minima found (8 x NPAR, cost, converged, kind), frozen theta0 and working states. */
	final double[] starts = new double[6 * NPAR], minU = new double[8 * NPAR], minF = new double[8], u = new double[NPAR],
		u2 = new double[NPAR], th0f = new double[NPL];
	final int[] startKind = new int[6], minK = new int[8];
	final boolean[] minC = new boolean[8];
	/** flip scratch: candidate |D|, D and JD at the minimum, clipped matrix, JD row, solve vector. */
	final double[] flipCand, flipD0, flipJd0, flipHs = new double[NPAR * NPAR], flipJd = new double[NPAR], flipW = new double[NPAR];
	/** covariance scratch: matrix copy, eigenvectors, eigenvalues, inverse eigenvalues. */
	final double[] eigA = new double[NPAR * NPAR], eigV = new double[NPAR * NPAR], eigW = new double[NPAR], eigIw = new double[NPAR];

	// ---- HitRecovery ----

	/** recovery scratch and the accepted rows of the last HitRecovery.recover call. */
	public final HitRecovery.Scratch rec;

	public GFitWorkspace(int maxHits) {
		cap = maxHits + 8 + 40;
		rec = new HitRecovery.Scratch(maxHits);
		flipCand = new double[cap]; flipD0 = new double[cap]; flipJd0 = new double[cap * NPAR];
		lg = new int[cap]; lli = new int[cap]; lw = new int[cap]; lism = new int[cap]; lcls = new int[cap]; lkind = new int[cap]; lrow = new int[cap];
		lm = new double[cap]; lsig = new double[cap]; lh = new double[cap]; sarc = new double[cap];
		D = new double[cap]; S = new double[cap]; JD = new double[cap * NPAR]; JS = new double[cap * NPAR];
		PF = new int[cap]; PF1 = new int[cap]; outH = new double[cap * 3];
		pcaG = new double[NMAX]; pcaD = new double[NMAX];
		nomX = new double[cap * 3]; pertD = new double[cap]; pertX = new double[cap * 3]; gradG = new double[cap * 3];
		for (int j = 0; j < snap.length; j++) snap[j] = new Propagator.State();
	}
}
