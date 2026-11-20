package org.jlab.rec.ahdc.GlobalFit;

/**
 * Immutable holders of the constants of the global fit's switchable terms (S1: the ATOF time correction of the pinned Python
 * tcal.py / atof_cal.py). Built by a {@link GFitConstants} implementation ({@code CcdbConstants} from CCDB tables, the test
 * provider from fixtures); the arrays are shared, not copied: callers must not modify them.
 */
public final class SwitchConstants {

	/** channels (15 sectors x 4 layers x 10 wedge components, tcal.py:292), gain groups and walk nodes per group (tcal.py NG_NODE = 12). */
	public static final int NCH = 600, NGROUP = 3, NNODE = 12;
	/** gain group limits of the p50 ToT (tcal.py:30, G_LO and G_HI). */
	public static final double G_LO = 1200.0, G_HI = 2800.0;

	private SwitchConstants() {}

	/** channel index (4 sector + layer) * 10 + component (tcal.py:37 {@code _atof_index}). */
	public static int channel(int sector, int layer, int component) { return (4 * sector + layer) * 10 + component; }

	/** gain group of a channel from its p50 ToT: 0 below G_LO, 2 above G_HI, else 1 (tcal.py:102-106 {@code group_of}). */
	public static int groupOf(double p50) { return p50 < G_LO ? 0 : p50 > G_HI ? 2 : 1; }

	/** the 8 ATOF time keys of one track half in the column order off s1 s2 f2 mlate slate flate fbg (P_TOFF..P_TFB). */
	public static final class TimeKeys {
		private final double[] v;

		public TimeKeys(double[] v) {
			if (v.length != 8) throw new IllegalArgumentException("time keys need 8 values");
			this.v = v;
		}

		/** the 8 values (shared, do not modify). */
		public double[] values() { return v; }
	}

	/**
	 * Per-channel time constants: T0 and p50 ToT of the walk fit (tcal.py), gain group, the walk nodes (x = ln(max(ToT, 1) / p50),
	 * W [ns]) per group, and the deployed upstream T0, slope and intercept of /calibration/alert/atof/time_offsets that
	 * {@link AtofTimeCorrection} adds back to the bank time.
	 */
	public static final class TimeChannels {
		public final double[] t0, p50, depT0, depSlope, depIcpt;
		public final int[] group;
		public final double[][] walkX, walkW;

		public TimeChannels(double[] t0, double[] p50, int[] group, double[][] walkX, double[][] walkW, double[] depT0,
			double[] depSlope, double[] depIcpt) {
			if (t0.length != NCH || p50.length != NCH || group.length != NCH || depT0.length != NCH || depSlope.length != NCH
				|| depIcpt.length != NCH || walkX.length != NGROUP || walkW.length != NGROUP)
				throw new IllegalArgumentException("time channels need " + NCH + " channels and " + NGROUP + " groups");
			for (int g = 0; g < NGROUP; g++)
				if (walkX[g].length != NNODE || walkW[g].length != NNODE)
					throw new IllegalArgumentException("walk needs " + NNODE + " nodes per group");
			this.t0 = t0;
			this.p50 = p50;
			this.group = group;
			this.walkX = walkX;
			this.walkW = walkW;
			this.depT0 = depT0;
			this.depSlope = depSlope;
			this.depIcpt = depIcpt;
		}
	}

	/** dE/dx constants: hit counts n = 0..16 (n = 0 unused), toy nodes in beta*gamma, z grid points (dedx.py:55 NMAX, toy_table: geomspace(0.12, 1.4, 36), linspace(-6, 10, 321)). */
	public static final int DX_NN = 17, DX_NBG = 36, DX_NZ = 321, DX_NLAYER = 8, DX_NF3 = 10;

	/**
	 * Constants of the dE/dx term of one run: the response line a, g, ln c, the ADC offsets,
	 * the F3 corrections, the toy tables (tau, sigma per half, hit count and beta*gamma node, and the standardised pdf per half and hit
	 * count) and the class variance [half][prompt][hit-count class]. Per-half rows are indexed by the half that uses them. The arrays
	 * are shared, not copied: callers must not modify them.
	 */
	public static final class DedxConstants {
		/** response a, g, ln c of the reference run and of this run; a_run = a + (lnc - lnc_ref) (comb_run.py:278-279). */
		public final double a, g, lnCRef, lnC, aRun;
		/** additive ADC offsets by layer index (8), F3 corrections [half][10] (gain_21 .. path_exp), F3 normalisation [half]. */
		public final double[] offset, f3Norm;
		public final double[][] f3;
		/** instrumental smear [half] and efficiency parameters [half][3] (documentation: the toy tables replace them). */
		public final double[] sigInst;
		public final double[][] effPar;
		/** toy nodes beta*gamma (36), their logarithms, tau and sigma [half][n][node]. */
		public final double[] toyBg, lnBg;
		public final double[][][] toyTau, toySig;
		/** z grid (321) and standardised pdf [half][n][bin 0..319]. */
		public final double[] zgrid;
		public final double[][][] pdf;
		/** class variance [half][prompt][hit-count class 0..2]. */
		public final double[][][] classVar;

		public DedxConstants(double a, double g, double lnCRef, double lnC, double[] offset, double[][] f3, double[] f3Norm, double[] sigInst,
			double[][] effPar, double[] toyBg, double[][][] toyTau, double[][][] toySig, double[] zgrid, double[][][] pdf, double[][][] classVar) {
			if (offset.length != DX_NLAYER || f3.length != 2 || f3[0].length != DX_NF3 || f3[1].length != DX_NF3 || f3Norm.length != 2
				|| sigInst.length != 2 || effPar.length != 2 || toyBg.length != DX_NBG || zgrid.length != DX_NZ)
				throw new IllegalArgumentException("dedx constants: bad array length");
			if (toyTau.length != 2 || toySig.length != 2 || pdf.length != 2 || classVar.length != 2 || classVar[0].length != 2
				|| classVar[1].length != 2 || classVar[0][0].length != 3)
				throw new IllegalArgumentException("dedx constants: need 2 halves");
			for (int h = 0; h < 2; h++) {
				if (toyTau[h].length != DX_NN || toySig[h].length != DX_NN || pdf[h].length != DX_NN)
					throw new IllegalArgumentException("dedx constants: toy tables need " + DX_NN + " hit counts");
				for (int n = 0; n < DX_NN; n++)
					if (toyTau[h][n].length != DX_NBG || toySig[h][n].length != DX_NBG || pdf[h][n].length != DX_NZ - 1)
						throw new IllegalArgumentException("dedx constants: bad toy table row");
			}
			this.a = a;
			this.g = g;
			this.lnCRef = lnCRef;
			this.lnC = lnC;
			this.aRun = a + (lnC - lnCRef);
			this.offset = offset;
			this.f3 = f3;
			this.f3Norm = f3Norm;
			this.sigInst = sigInst;
			this.effPar = effPar;
			this.toyBg = toyBg;
			this.lnBg = new double[DX_NBG];
			for (int k = 0; k < DX_NBG; k++) lnBg[k] = Math.log(toyBg[k]);
			this.toyTau = toyTau;
			this.toySig = toySig;
			this.zgrid = zgrid;
			this.pdf = pdf;
			this.classVar = classVar;
		}
	}
}
