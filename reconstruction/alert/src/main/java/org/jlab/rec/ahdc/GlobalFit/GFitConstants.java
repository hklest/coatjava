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

	/** beam position (x, y) in mm, AHDC frame (plan D8). */
	double[] beamXY(int run);

	/** ATOF position parameters and bar recalibration for tracks of {@code half} (0 or 1) of run {@code run}. */
	AtofPositionTable atofPosition(int run, int half);

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
		return m;
	}
}
