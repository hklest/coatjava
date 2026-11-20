package org.jlab.rec.ahdc.GlobalFit;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Per-run constants of the global track fit: the beam position and the ATOF position table (the constants that
 * {@code atof.load_par}/{@code load_recal} read, atof.py:330-366 of the pin). The engine implementation reads CCDB
 * ({@code CcdbConstants}); the test implementation serves the fixture tables.
 */
public interface GFitConstants {

	/** CCDB index depth of every table of the fit (as AHDCEngine and ATOFEngine request their tables). */
	int DEPTH = 3;

	/** beam position (x, y) in mm, AHDC frame. */
	double[] beamXY(int run);

	/** ATOF position parameters and bar recalibration for tracks of {@code half} (0 or 1) of run {@code run}. */
	AtofPositionTable atofPosition(int run, int half);

	/** tables of the ATOF time correction (switch on): keys and channels depth {@link #DEPTH}, time_offsets depth 4 as ATOFEngine requests it. */
	String TIME_KEYS = "/calibration/alert/atof/gfit_time_keys", TIME_CHANNEL = "/calibration/alert/atof/gfit_time_channel",
		TIME_WALK = "/calibration/alert/atof/gfit_time_walk", TIME_OFFSETS = "/calibration/alert/atof/time_offsets";
	/** depth of {@link #TIME_OFFSETS} (sector, layer, component, order). */
	int TIME_OFFSETS_DEPTH = 4;
	/** tables of the dE/dx term (switch on), depth {@link #DEPTH}. */
	String[] DEDX_TABLES = {"/calibration/alert/ahdc/gfit_dedx", "/calibration/alert/ahdc/gfit_dedx_offset",
		"/calibration/alert/ahdc/gfit_dedx_toy", "/calibration/alert/ahdc/gfit_dedx_pdf", "/calibration/alert/ahdc/gfit_dedx_classvar"};

	/** time keys of {@code run} for tracks of {@code half}, row indexed by the using half (never swapped); needs the time switch. */
	default SwitchConstants.TimeKeys timeKeys(int run, int half) { throw new UnsupportedOperationException("time keys"); }

	/** per-channel time constants of {@code run}; needs the time switch. */
	default SwitchConstants.TimeChannels timeChannels(int run) { throw new UnsupportedOperationException("time channels"); }

	/** dE/dx constants of {@code run} (class variance of the all set, as the CCDB table holds it); needs the dE/dx switch. */
	default SwitchConstants.DedxConstants dedx(int run) { throw new UnsupportedOperationException("dedx constants"); }

	/** tables of the default fit (no recovery), name to index depth, as ReconstructionEngine.requireConstants takes them. */
	static Map<String, Integer> requiredTables(GFitConfig cfg) {
		return requiredTables(cfg, false);
	}

	/**
	 * tables of the default fit; {@code recovery} adds the two AHDC time tables of the hit recovery (D7, D11).
	 * GFitConfig has no recovery switch, so it is a parameter.
	 */
	static Map<String, Integer> requiredTables(GFitConfig cfg, boolean recovery) {
		Map<String, Integer> m = new LinkedHashMap<>();
		m.put("/geometry/beam/position", DEPTH);
		m.put("/geometry/alert/ahdc/beam_offset", DEPTH);
		m.put("/calibration/alert/atof/gfit_position", DEPTH);
		m.put("/calibration/alert/atof/gfit_bar_z", DEPTH);
		if (recovery) {
			m.put("/calibration/alert/ahdc/time_offsets", DEPTH);
			m.put("/calibration/alert/ahdc/time_to_distance_wire", DEPTH);
		}
		if (cfg.timeOn()) {
			m.put(TIME_KEYS, DEPTH);
			m.put(TIME_CHANNEL, DEPTH);
			m.put(TIME_WALK, DEPTH);
			m.put(TIME_OFFSETS, TIME_OFFSETS_DEPTH);
		}
		if (cfg.dedxOn()) for (String t : DEDX_TABLES) m.put(t, DEPTH);
		return m;
	}
}
