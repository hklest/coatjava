package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jlab.utils.groups.IndexedTable;
import org.junit.Test;

public class GFitCcdbConstantsTest {

	private static final String[] SPATIAL = {"wz_len", "wz_sig", "wz_off", "row_w", "wphi_sig", "wphi_off", "bz_off", "bz_sig1",
		"bz_sig2", "bz_f2", "bphi_sig", "bphi_off", "eta_w", "eta_b", "rho_w", "rho_b"};
	private static final String[] TIME = {"off", "s1", "s2", "f2", "mlate", "slate", "flate", "fbg"};

	/** in-memory tables of the pooled fixture values, beam (0.052, -0.023) cm and offset (-0.157, 0.023) cm. */
	private static Map<String, IndexedTable> tables() {
		GFitFixtures.Fx fx = GFitFixtures.resource("unit/atof_tables.txt.gz");
		Map<String, IndexedTable> m = new HashMap<>();
		IndexedTable beam = new IndexedTable(3, "x_offset:y_offset");
		beam.addEntry(0, 0, 0);
		beam.setDoubleValue(0.052, "x_offset", 0, 0, 0);
		beam.setDoubleValue(-0.023, "y_offset", 0, 0, 0);
		m.put(CcdbConstants.BEAM, beam);
		IndexedTable off = new IndexedTable(3, "x_offset:y_offset");
		off.addEntry(0, 0, 0);
		off.setDoubleValue(-0.157, "x_offset", 0, 0, 0);
		off.setDoubleValue(0.023, "y_offset", 0, 0, 0);
		m.put(CcdbConstants.BEAM_OFFSET, off);

		String[] names = new String[SPATIAL.length + 1 + TIME.length];
		for (int i = 0; i < SPATIAL.length; i++) names[i] = SPATIAL[i];
		names[SPATIAL.length] = "t_win";
		for (int i = 0; i < TIME.length; i++) names[SPATIAL.length + 1 + i] = "t_" + TIME[i];
		IndexedTable pos = new IndexedTable(3, String.join(":", names));
		pos.addEntry(0, 0, 0);
		for (String n : names) {
			String key = n.equals("bz_off") ? "bz_off_wW" : n.startsWith("t_") && !n.equals("t_win") ? "t_deployed_entry_" + n.substring(2) : n;
			pos.setDoubleValue(fx.d0("hall." + key), n, 0, 0, 0);
		}
		m.put(CcdbConstants.ATOF_POSITION, pos);

		double[] rc = fx.d("hall.bz_recal");
		IndexedTable bar = new IndexedTable(3, "offset:scale");
		for (int b = 0; b < 60; b++) {
			bar.addEntry(b / 4, b % 4, 10);
			bar.setDoubleValue(rc[2 * b], "offset", b / 4, b % 4, 10);
			bar.setDoubleValue(rc[120 + 2 * b], "scale", b / 4, b % 4, 10);
		}
		m.put(CcdbConstants.ATOF_BAR_Z, bar);
		return m;
	}

	@Test
	public void ccdbConstantsFromTables() {
		Map<String, IndexedTable> m = tables();
		CcdbConstants c = new CcdbConstants((run, name) -> m.get(name));
		assertTrue(c.missing(23014).isEmpty());
		GFitFixtureConstants fix = new GFitFixtureConstants(true);
		for (int h = 0; h < 2; h++) {
			AtofPositionTable want = fix.atofPosition(23015, h), got = c.atofPosition(23014, h);
			assertArrayEquals("par " + h, want.par(), got.par(), 0.0);
			assertArrayEquals("offset " + h, want.recalOffset(), got.recalOffset(), 0.0);
			assertArrayEquals("scale " + h, want.recalScale(), got.recalScale(), 0.0);
		}
		assertArrayEquals(new double[] {-1.05, 0.0}, c.beamXY(23014), 1e-12);
	}

	@Test
	public void missingTableIsReported() {
		for (String gone : new String[] {CcdbConstants.BEAM, CcdbConstants.BEAM_OFFSET, CcdbConstants.ATOF_POSITION, CcdbConstants.ATOF_BAR_Z}) {
			Map<String, IndexedTable> m = tables();
			m.remove(gone);
			CcdbConstants c = new CcdbConstants((run, name) -> m.get(name));
			assertEquals(Arrays.asList(gone), c.missing(23014));
			try {
				if (gone.startsWith("/geometry")) c.beamXY(23014); else c.atofPosition(23014, 0);
				fail("no exception for missing " + gone);
			} catch (CcdbConstants.MissingTableException e) {
				assertEquals(List.of(gone), e.tables());
				assertTrue(e.getMessage().contains(gone));
			}
		}
	}
}
