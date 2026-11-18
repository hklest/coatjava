package org.jlab.rec.ahdc.GlobalFit;

import static org.jlab.rec.ahdc.GlobalFit.AtofTerm.*;

import java.util.function.ToDoubleFunction;

/**
 * Immutable ATOF position parameters of one track half: the vector {@code PAR[h]} of 34 values in the P_* order of
 * atof.py:36-37 and the v2 bar recalibration (offset and scale per bar, 60 each). {@link #fromKeys} ports
 * {@code atof.load_par} (atof.py:330-356) for one table with the acfg of variant ge5_vz_zt (use_bar, tassoc, zr_mode mid,
 * time off, tmode entry, tresp (0, 1), force_w and ceiling off) and {@code load_recal} (atof.py:359-364). The arrays are
 * shared, not copied: callers must not modify them.
 */
public final class AtofPositionTable {

	/** parameter count and P_* indices are those of {@link AtofTerm}. */
	public static final int NP = AtofTerm.NP, NBAR = 60;

	private static final String[] SPATIAL = {"wz_len", "wz_sig", "wz_off", "row_w", "wphi_sig", "wphi_off", "bz_off", "bz_sig1",
		"bz_sig2", "bz_f2", "bphi_sig", "bphi_off", "eta_w", "eta_b", "rho_w", "rho_b"};
	private static final String[] TIME = {"off", "s1", "s2", "f2", "mlate", "slate", "flate", "fbg"};

	private final double[] par, offset, scale;

	public AtofPositionTable(double[] par, double[] offset, double[] scale) {
		if (par.length != NP || offset.length != NBAR || scale.length != NBAR)
			throw new IllegalArgumentException("ATOF table needs " + NP + " parameters and " + NBAR + " bars");
		this.par = par;
		this.offset = offset;
		this.scale = scale;
	}

	/** the 34 parameters in P_* order. */
	public double[] par() { return par; }

	/** bar recalibration z_cal = (z - offset) / scale, per bar b = 4 sector + layer. */
	public double[] recalOffset() { return offset; }

	public double[] recalScale() { return scale; }

	/**
	 * Builds the table from the named values of one calibration table. {@code key} returns the value of a key of the
	 * Python npz ({@code wz_len ... rho_b}, {@code bz_off_wW}, {@code t_win}, {@code t_<calset>_entry_<off|s1|s2|f2|mlate|slate|flate|fbg>}).
	 * {@code calset} is "final" (faithful tests) or "deployed" (bank times, engine). The recalibration comes from
	 * {@code bz_recal} (2 x 60 x 2, row-major) as {@link #recalFromBzRecal}.
	 */
	public static AtofPositionTable fromKeys(ToDoubleFunction<String> key, double[] bzRecal, String calset) {
		double[] p = new double[NP];
		for (int i = 0; i < SPATIAL.length; i++) p[i] = key.applyAsDouble(SPATIAL[i]);
		p[P_BZOFF] = key.applyAsDouble("bz_off_wW");      // cal v2
		p[P_TWIN] = key.applyAsDouble("t_win");
		p[P_BARON] = 1.0;
		p[P_TRA] = 0.0;
		p[P_TRB] = 1.0;
		p[P_ZMID] = 1.0;
		for (int i = 0; i < TIME.length; i++) p[P_TOFF + i] = key.applyAsDouble("t_" + calset + "_entry_" + TIME[i]);
		p[P_TON] = 0.0;
		p[P_TMODE] = 0.0;
		p[P_TASSOC] = 1.0;
		double[][] r = recalFromBzRecal(bzRecal);
		return new AtofPositionTable(p, r[0], r[1]);
	}

	/**
	 * The ge10_t configuration of this table (comb_time.py ge10_t: time on, tassoc off): a copy with P_TON 1, P_TASSOC 0 and the
	 * eight time keys {@code off s1 s2 f2 mlate slate flate fbg} of {@code keys} in P_TOFF..P_TFB. Built once per run and half
	 * by the caller, never per track; the fit then reads the wedge time from {@link GFitEvent#aTime} (see {@link AtofTerm#attach}).
	 */
	public AtofPositionTable withTime(SwitchConstants.TimeKeys keys) {
		double[] p = par.clone(), k = keys.values();
		for (int i = 0; i < TIME.length; i++) p[P_TOFF + i] = k[i];
		p[P_TON] = 1.0;
		p[P_TASSOC] = 0.0;
		return new AtofPositionTable(p, offset, scale);
	}

	/** {offset, scale} = bz_recal[0, :, 0], bz_recal[1, :, 0] of a flat 2 x 60 x 2 array. */
	public static double[][] recalFromBzRecal(double[] bzRecal) {
		double[] o = new double[NBAR], s = new double[NBAR];
		for (int b = 0; b < NBAR; b++) {
			o[b] = bzRecal[2 * b];
			s[b] = bzRecal[2 * NBAR + 2 * b];
		}
		return new double[][] {o, s};
	}
}
