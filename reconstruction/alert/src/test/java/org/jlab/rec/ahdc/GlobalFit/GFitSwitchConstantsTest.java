package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.jlab.utils.groups.IndexedTable;
import org.junit.Test;

public class GFitSwitchConstantsTest {

	private static final String[] TIME = {"off", "s1", "s2", "f2", "mlate", "slate", "flate", "fbg"};
	private static final String[] FOUR = {"/geometry/beam/position", "/geometry/alert/ahdc/beam_offset",
		"/calibration/alert/atof/gfit_position", "/calibration/alert/atof/gfit_bar_z"};
	private static final String[] TIME_TABLES = {"/calibration/alert/atof/gfit_time_keys", "/calibration/alert/atof/gfit_time_channel",
		"/calibration/alert/atof/gfit_time_walk", "/calibration/alert/atof/time_offsets"};
	private static final String[] DEDX = {"/calibration/alert/ahdc/gfit_dedx", "/calibration/alert/ahdc/gfit_dedx_offset",
		"/calibration/alert/ahdc/gfit_dedx_toy", "/calibration/alert/ahdc/gfit_dedx_pdf", "/calibration/alert/ahdc/gfit_dedx_classvar"};

	/** in-memory time tables of the constants {@code c} for the CcdbConstants of run {@code run} (keys by half, never swapped). */
	static Map<String, IndexedTable> timeTables(GFitConstants c, int run) {
		Map<String, IndexedTable> m = new HashMap<>();
		IndexedTable keys = new IndexedTable(3, String.join(":", TIME));
		for (int h = 0; h < 2; h++) {
			keys.addEntry(0, 0, h);
			double[] v = c.timeKeys(run, h).values();
			for (int i = 0; i < TIME.length; i++) keys.setDoubleValue(v[i], TIME[i], 0, 0, h);
		}
		m.put(GFitConstants.TIME_KEYS, keys);
		SwitchConstants.TimeChannels tc = c.timeChannels(run);
		IndexedTable chan = new IndexedTable(3, "t0:p50"), off = new IndexedTable(4, "t0:extra1:extra2");
		for (int s = 0; s < 15; s++)
			for (int l = 0; l < 4; l++)
				for (int k = 0; k < 10; k++) {
					int ch = SwitchConstants.channel(s, l, k);
					chan.addEntry(s, l, k);
					chan.setDoubleValue(tc.t0[ch], "t0", s, l, k);
					chan.setDoubleValue(tc.p50[ch], "p50", s, l, k);
					off.addEntry(s, l, k, 0);
					off.setDoubleValue(tc.depT0[ch], "t0", s, l, k, 0);
					off.setDoubleValue(tc.depSlope[ch], "extra1", s, l, k, 0);
					off.setDoubleValue(tc.depIcpt[ch], "extra2", s, l, k, 0);
				}
		m.put(GFitConstants.TIME_CHANNEL, chan);
		m.put(GFitConstants.TIME_OFFSETS, off);
		StringBuilder cols = new StringBuilder();
		for (int k = 0; k < SwitchConstants.NNODE; k++) cols.append(k == 0 ? "" : ":").append("x").append(k);
		for (int k = 0; k < SwitchConstants.NNODE; k++) cols.append(":w").append(k);
		IndexedTable walk = new IndexedTable(3, cols.toString());
		for (int g = 0; g < SwitchConstants.NGROUP; g++) {
			walk.addEntry(0, 0, g);
			for (int k = 0; k < SwitchConstants.NNODE; k++) {
				walk.setDoubleValue(tc.walkX[g][k], "x" + k, 0, 0, g);
				walk.setDoubleValue(tc.walkW[g][k], "w" + k, 0, 0, g);
			}
		}
		m.put(GFitConstants.TIME_WALK, walk);
		return m;
	}

	private static final String[] F3_COLUMNS = {"gain_21", "gain_22", "gain_31", "gain_32", "gain_41", "gain_42", "gain_51", "doca_per_mm", "ped_per_50", "path_exp"};

	/** in-memory dE/dx tables of the constants {@code c} for {@code run}, filled with the layout of the CCDB files (half rows not swapped). */
	static Map<String, IndexedTable> dedxTables(GFitConstants c, int run) {
		SwitchConstants.DedxConstants d = c.dedx(run);
		Map<String, IndexedTable> m = new HashMap<>();
		IndexedTable main = new IndexedTable(3, "a:g:lnc_ref:lnc:sig_inst:eff_lnthr:eff_width:eff_plateau:" + String.join(":", F3_COLUMNS) + ":f3_norm");
		for (int h = 0; h < 2; h++) {
			main.addEntry(0, 0, h);
			main.setDoubleValue(d.a, "a", 0, 0, h);
			main.setDoubleValue(d.g, "g", 0, 0, h);
			main.setDoubleValue(d.lnCRef, "lnc_ref", 0, 0, h);
			main.setDoubleValue(d.lnC, "lnc", 0, 0, h);
			main.setDoubleValue(d.sigInst[h], "sig_inst", 0, 0, h);
			main.setDoubleValue(d.effPar[h][0], "eff_lnthr", 0, 0, h);
			main.setDoubleValue(d.effPar[h][1], "eff_width", 0, 0, h);
			main.setDoubleValue(d.effPar[h][2], "eff_plateau", 0, 0, h);
			for (int i = 0; i < F3_COLUMNS.length; i++) main.setDoubleValue(d.f3[h][i], F3_COLUMNS[i], 0, 0, h);
			main.setDoubleValue(d.f3Norm[h], "f3_norm", 0, 0, h);
		}
		m.put(GFitConstants.DEDX_TABLES[0], main);
		IndexedTable off = new IndexedTable(3, "offset"), toy = new IndexedTable(3, "bg:tau:sig"), pdf = new IndexedTable(3, "pdf"), cv = new IndexedTable(3, "var");
		for (int li = 0; li < d.offset.length; li++) {
			off.addEntry(0, li, 0);
			off.setDoubleValue(d.offset[li], "offset", 0, li, 0);
		}
		for (int h = 0; h < 2; h++) {
			for (int n = 0; n < SwitchConstants.DX_NN; n++) {
				for (int k = 0; k < SwitchConstants.DX_NBG; k++) {
					toy.addEntry(h, n, k);
					toy.setDoubleValue(d.toyBg[k], "bg", h, n, k);
					toy.setDoubleValue(d.toyTau[h][n][k], "tau", h, n, k);
					toy.setDoubleValue(d.toySig[h][n][k], "sig", h, n, k);
				}
				for (int iz = 0; iz < SwitchConstants.DX_NZ - 1; iz++) {
					pdf.addEntry(h, n, iz);
					pdf.setDoubleValue(d.pdf[h][n][iz], "pdf", h, n, iz);
				}
			}
			for (int w = 0; w < 2; w++)
				for (int k = 0; k < 3; k++) {
					cv.addEntry(h, w, k);
					cv.setDoubleValue(d.classVar[h][w][k], "var", h, w, k);
				}
		}
		m.put(GFitConstants.DEDX_TABLES[1], off);
		m.put(GFitConstants.DEDX_TABLES[2], toy);
		m.put(GFitConstants.DEDX_TABLES[3], pdf);
		m.put(GFitConstants.DEDX_TABLES[4], cv);
		return m;
	}

	@Test
	public void requiredTablesSwitches() {
		GFitConfig base = GFitConfig.pythonDefault();
		assertTrue(!base.timeOn() && !base.dedxOn());
		Map<String, Integer> off = GFitConstants.requiredTables(base.withSwitches(false, false), false);
		assertEquals(Arrays.asList(FOUR), List.copyOf(off.keySet()));
		for (String t : FOUR) assertEquals(Integer.valueOf(3), off.get(t));
		assertEquals(GFitConstants.requiredTables(base, true), GFitConstants.requiredTables(base.withSwitches(false, false), true));
		Map<String, Integer> time = GFitConstants.requiredTables(base.withSwitches(true, false), false);
		assertEquals(8, time.size());
		for (String t : FOUR) assertEquals(Integer.valueOf(3), time.get(t));
		for (String t : TIME_TABLES) assertEquals(t.endsWith("time_offsets") ? Integer.valueOf(4) : Integer.valueOf(3), time.get(t));
		Map<String, Integer> dedx = GFitConstants.requiredTables(base.withSwitches(false, true), false);
		assertEquals(9, dedx.size());
		for (String t : DEDX) assertEquals(Integer.valueOf(3), dedx.get(t));
		Map<String, Integer> both = new LinkedHashMap<>(GFitConstants.requiredTables(base.withSwitches(true, true), true));
		assertEquals(4 + 2 + 4 + 5, both.size());
	}

	@Test
	public void switchConstantsFromTables() {
		GFitFixtureConstants fix = new GFitFixtureConstants(true);
		for (int run : new int[] {23014, 23015}) {
			Map<String, IndexedTable> m = timeTables(fix, run);
			CcdbConstants c = new CcdbConstants((r, name) -> m.get(name), GFitConfig.pythonDefault().withSwitches(true, false));
			assertEquals("only the position tables are not in the map", Arrays.asList(FOUR), c.missing(run));
			for (String t : FOUR) m.put(t, new IndexedTable(3, "x"));
			assertTrue(c.missing(run).isEmpty());
			for (int h = 0; h < 2; h++)
				assertArrayEquals("keys " + run + " half " + h, fix.timeKeys(run, h).values(), c.timeKeys(run, h).values(), 0.0);
			assertTrue("half rows differ, so a swap would show", !Arrays.equals(c.timeKeys(run, 0).values(), c.timeKeys(run, 1).values()));
			SwitchConstants.TimeChannels w = fix.timeChannels(run), g = c.timeChannels(run);
			assertArrayEquals(w.t0, g.t0, 0.0);
			assertArrayEquals(w.p50, g.p50, 0.0);
			assertArrayEquals(w.group, g.group);
			assertArrayEquals(w.depT0, g.depT0, 0.0);
			assertArrayEquals(w.depSlope, g.depSlope, 0.0);
			assertArrayEquals(w.depIcpt, g.depIcpt, 0.0);
			for (int k = 0; k < SwitchConstants.NGROUP; k++) {
				assertArrayEquals(w.walkX[k], g.walkX[k], 0.0);
				assertArrayEquals(w.walkW[k], g.walkW[k], 0.0);
			}
		}
	}

	@Test
	public void switchDedxConstantsFromTables() {
		GFitFixtureConstants fix = new GFitFixtureConstants(true);
		for (int run : new int[] {23014, 23015}) {
			Map<String, IndexedTable> m = dedxTables(fix, run);
			CcdbConstants c = new CcdbConstants((r, name) -> m.get(name), GFitConfig.pythonDefault().withSwitches(false, true));
			assertTrue(c.missing(run).containsAll(Arrays.asList(FOUR)) && c.missing(run).size() == FOUR.length);
			SwitchConstants.DedxConstants w = fix.dedx(run), g = c.dedx(run);
			assertEquals(w.a, g.a, 0.0);
			assertEquals(w.g, g.g, 0.0);
			assertEquals(w.lnCRef, g.lnCRef, 0.0);
			assertEquals(w.lnC, g.lnC, 0.0);
			assertEquals(w.aRun, g.aRun, 0.0);
			assertArrayEquals(w.offset, g.offset, 0.0);
			assertArrayEquals(w.f3Norm, g.f3Norm, 0.0);
			assertArrayEquals(w.sigInst, g.sigInst, 0.0);
			assertArrayEquals(w.toyBg, g.toyBg, 0.0);
			assertArrayEquals(w.lnBg, g.lnBg, 0.0);
			assertArrayEquals(w.zgrid, g.zgrid, 0.0);
			for (int h = 0; h < 2; h++) {
				assertArrayEquals(w.f3[h], g.f3[h], 0.0);
				assertArrayEquals(w.effPar[h], g.effPar[h], 0.0);
				for (int n = 0; n < SwitchConstants.DX_NN; n++) {
					assertArrayEquals(w.toyTau[h][n], g.toyTau[h][n], 0.0);
					assertArrayEquals(w.toySig[h][n], g.toySig[h][n], 0.0);
					assertArrayEquals(w.pdf[h][n], g.pdf[h][n], 0.0);
				}
				for (int p = 0; p < 2; p++) assertArrayEquals(w.classVar[h][p], g.classVar[h][p], 0.0);
			}
			assertTrue("half rows differ, so a swap would show", !Arrays.equals(g.f3[0], g.f3[1]) && !Arrays.equals(g.f3Norm, new double[] {g.f3Norm[1], g.f3Norm[0]}));
		}
	}

	@Test
	public void missingDedxTable() {
		GFitFixtureConstants fix = new GFitFixtureConstants(true);
		Map<String, IndexedTable> m = dedxTables(fix, 23014);
		for (String t : FOUR) m.put(t, new IndexedTable(3, "x"));
		m.remove(DEDX[3]);
		CcdbConstants c = new CcdbConstants((r, name) -> m.get(name), GFitConfig.pythonDefault().withSwitches(false, true));
		assertEquals(List.of(DEDX[3]), c.missing(23014));
		try {
			c.dedx(23014);
			fail("no exception for the missing pdf table");
		} catch (CcdbConstants.MissingTableException e) {
			assertEquals(List.of(DEDX[3]), e.tables());
		}
		assertTrue(new CcdbConstants((r, name) -> m.get(name)).missing(23014).isEmpty());
		assertTrue(new CcdbConstants((r, name) -> m.get(name), GFitConfig.pythonDefault().withSwitches(true, false)).missing(23014).containsAll(Arrays.asList(TIME_TABLES)));
	}

	@Test
	public void missingSwitchTable() {
		GFitFixtureConstants fix = new GFitFixtureConstants(true);
		Map<String, IndexedTable> m = timeTables(fix, 23014);
		for (String t : FOUR) m.put(t, new IndexedTable(3, "x"));
		m.remove(GFitConstants.TIME_WALK);
		GFitConfig on = GFitConfig.pythonDefault().withSwitches(true, false);
		CcdbConstants c = new CcdbConstants((r, name) -> m.get(name), on);
		assertEquals(List.of(GFitConstants.TIME_WALK), c.missing(23014));
		try {
			c.timeChannels(23014);
			fail("no exception for the missing walk table");
		} catch (CcdbConstants.MissingTableException e) {
			assertEquals(List.of(GFitConstants.TIME_WALK), e.tables());
		}
		// switch off: the table is not requested
		assertTrue(new CcdbConstants((r, name) -> m.get(name)).missing(23014).isEmpty());
		assertTrue(new CcdbConstants((r, name) -> m.get(name), GFitConfig.pythonDefault()).missing(23014).isEmpty());
		// dE/dx on: its five tables are requested
		assertEquals(Arrays.asList(DEDX), new CcdbConstants((r, name) -> m.get(name), GFitConfig.pythonDefault().withSwitches(false, true)).missing(23014));
	}
}
