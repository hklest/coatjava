package org.jlab.service.alert;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.rec.ahdc.GlobalFit.CcdbConstants;
import org.jlab.rec.ahdc.GlobalFit.CellGeometry;
import org.jlab.rec.ahdc.GlobalFit.FieldGrid;
import org.jlab.rec.ahdc.GlobalFit.GFitConfig;
import org.jlab.rec.ahdc.GlobalFit.GFitEvent;
import org.jlab.rec.ahdc.GlobalFit.GFitFixtureConstants;
import org.jlab.rec.ahdc.GlobalFit.GFitFixtures;
import org.jlab.rec.ahdc.GlobalFit.GFitReplaySet;
import org.jlab.rec.ahdc.GlobalFit.GFitResult;
import org.jlab.rec.ahdc.GlobalFit.GFitTrack;
import org.jlab.rec.ahdc.GlobalFit.GFitWorkspace;
import org.jlab.rec.ahdc.GlobalFit.HitRecovery;
import org.jlab.rec.ahdc.Hit.Hit;
import org.jlab.rec.ahdc.Hit.HitReader;
import org.jlab.utils.groups.IndexedTable;
import org.junit.Test;

/** GFitEngineAdapter input building, recovery pool, rows and per-track failure handling, without engine init(). */
public class GFitEngineAdapterTest {

    private static CellGeometry geometry() {
        GFitFixtures.Fx w = GFitFixtures.resource("unit/wires.txt.gz");
        return CellGeometry.fromWireEnds(w.d("WO"), w.d("WE"), w.i("layoff"));
    }

    /** an event whose AHDC::adc rows are (layer byte, component) pairs, raw ADC 100 + row, leading edge 300 + row. */
    private static GFitEvent event(CellGeometry geo, int[] layer, int[] comp, double start) {
        GFitEvent ev = new GFitEvent();
        int n = layer.length;
        ev.startTime = start;
        ev.rAdc = new int[n]; ev.rWf = new int[n]; ev.rLet = new double[n]; ev.rTot = new double[n]; ev.rPed = new double[n];
        for (int r = 0; r < n; r++) { ev.rAdc[r] = 100 + r; ev.rLet[r] = 300 + r; ev.rTot[r] = 50; ev.rPed[r] = 200; }
        GFitEngineAdapter.wires(layer, comp, geo.layoff, ev);
        return ev;
    }

    @Test
    public void layerAndWireFromLayerByte() {
        CellGeometry geo = geometry();
        GFitEvent ev = event(geo, new int[] {11, 21, 51, 51, 12, 22}, new int[] {1, 3, 99, 1, 1, 0}, 30);
        int[] lo = geo.layoff;
        int nw51 = lo[8] - lo[7];
        assertArrayEquals(new int[] {0, 1, 7, 7, -1, -1}, ev.rLi);
        assertArrayEquals(new int[] {0, 2, nw51 < 99 ? -1 : 98, 0, -1, -1}, ev.rWire);
        assertArrayEquals(new int[] {lo[0], lo[1] + 2, nw51 < 99 ? -1 : lo[7] + 98, lo[7], -1, -1}, ev.eg);
    }

    @Test
    public void seedFromMomentum() {
        double[] s = GFitEngineAdapter.seed(100.0, 100.0, 0.0, 12.5);
        assertEquals(Math.PI / 4, s[0], 1e-15);
        assertEquals(Math.sqrt(2) * 0.1, s[1], 1e-15);
        assertEquals(Math.PI / 2, s[2], 1e-15);
        assertEquals(12.5, s[3], 0);
        assertTrue(Double.isNaN(GFitEngineAdapter.seed(Double.NaN, 1, 1, 0)[1]));
        assertEquals(67.0 + 10 * 0.5, GFitEngineAdapter.zElectron(0.5, 1.0, false), 1e-12);
        assertEquals(36.0 + 10 * 0.5, GFitEngineAdapter.zElectron(0.5, -1.0, false), 1e-12);
        assertEquals(5.0, GFitEngineAdapter.zElectron(0.5, 1.0, true), 0);
    }

    @Test
    public void trackInputs() {
        CellGeometry geo = geometry();
        GFitEvent ev = event(geo, new int[] {11, 21, 22, 31, 32}, new int[] {5, 6, 7, 8, 9}, 30);
        List<Hit> hits = new ArrayList<>();
        hits.add(new Hit(4, 3, 1, 8, 1.2, 0, 150.0));      // row 3
        hits.add(new Hit(1, 1, 1, 5, 0.4, 0, -5.0));       // row 0
        hits.add(new Hit(2, 2, 1, 6, Double.NaN, 0, 80));  // row 1, NaN doca: left out
        hits.add(new Hit(4, 3, 1, 8, 1.2, 0, 150.0));      // row 3 again: taken once
        hits.add(new Hit(3, 2, 2, 7, 2.0, 0, -30.0));      // row 2
        boolean[] in = new boolean[5];
        GFitTrack t = GFitEngineAdapter.track(7, hits, ev, geo, 2.0, new double[] {0.1, 0.4, 1.2, 3.0, 0}, in);
        assertEquals(7, t.trackid);
        assertArrayEquals(new int[] {0, 2, 3}, t.hrow);
        assertArrayEquals(new int[] {ev.eg[0], ev.eg[2], ev.eg[3]}, t.hg);
        assertArrayEquals(new double[] {0.4, 2.0, 1.2}, t.hm, 0);
        assertArrayEquals(new int[] {1, 2, 0}, t.hism);
        assertArrayEquals(new int[] {0, 2, 3}, t.hcls);
        assertEquals(2.0 * HitRecovery.kfError(100, -5.0), t.hsig[0], 0);
        assertEquals(2.0 * HitRecovery.kfError(103, 150.0), t.hsig[2], 0);
        assertEquals(geo.half[ev.eg[2]], t.hh[1], 0);
        assertEquals(0.4, t.seedP, 0);
        assertArrayEquals(new boolean[] {true, true, true, true, false}, in);      // the NaN-doca row is on the track: not a candidate
    }

    private static IndexedTable t0Table(double t0) {
        IndexedTable t = new IndexedTable(3, "t0");
        for (int l : GFitEngineAdapter.LAYERS) for (int w = 1; w <= 99; w++) {
            t.addEntry(1, l, w);
            t.setDoubleValue(t0, "t0", 1, l, w);
        }
        return t;
    }

    private static IndexedTable t2dTable() {
        String[] c = new String[16];
        for (int i = 0; i < 16; i++) c[i] = "c" + i;
        IndexedTable t = new IndexedTable(3, String.join(":", c));
        double[] v = {0.1, 0.02, 0.3, 0.015, 0.5, 0.01, 60, 10, 160, 20, 0, 0, 0, 0, 0, 0};
        for (int l : GFitEngineAdapter.LAYERS) for (int w = 1; w <= 99; w++) {
            t.addEntry(1, l, w);
            for (int i = 0; i < 16; i++) t.setDoubleValue(v[i], c[i], 1, l, w);
        }
        return t;
    }

    @Test
    public void poolFromRowsOnNoTrack() {
        CellGeometry geo = geometry();
        GFitEvent ev = event(geo, new int[] {11, 21, 22, 12}, new int[] {5, 6, 7, 1}, 30);
        IndexedTable t2d = t2dTable();
        boolean[] in = {true, false, false, false};      // row 0 on a track, rows 1, 2 pruned or never on one, row 3 unknown layer
        GFitEngineAdapter.pool(ev, in, GFitConfig.pythonDefault(), t0Table(200), t2d);
        assertArrayEquals(new int[] {1, 2}, ev.pEv);
        assertArrayEquals(new int[] {ev.eg[1], ev.eg[2]}, ev.pG);
        assertEquals(301 - 200 - 30, ev.pTime[0], 0);
        assertEquals(HitReader.timeToDistance(1, 21, 6, 71.0, t2d), ev.pDoca[0], 0);
        assertEquals(50.0, ev.pTot[1], 0);
        GFitEvent nostart = event(geo, new int[] {11, 21}, new int[] {5, 6}, Double.NaN);
        GFitEngineAdapter.pool(nostart, new boolean[2], GFitConfig.pythonDefault(), t0Table(200), t2d);
        assertEquals(0, nostart.nPool());
    }

    /** a fit stub: status 0 at a fixed state, throws for the track with id 2. */
    private static final class Stub implements GFitEngineAdapter.TrackFit {
        @Override
        public int fit(GFitTrack t, GFitEvent ev, boolean zPrior, GFitWorkspace ws, GFitResult out) {
            if (t.trackid == 2) throw new IllegalStateException("stub");
            out.status = 0; out.u[0] = 0.5; out.u[1] = 0.3; out.u[2] = Math.log(0.25); out.u[3] = -12; out.u[4] = 1.1;
            out.cov[0] = 0.04; out.cov[1] = 0.001; out.cov[2 * GFitConfig.NPAR + 2] = 0.002; out.cost = 7.5; out.nh = t.nHits();
            return 1;
        }

        @Override
        public int[] atof(GFitResult r, GFitEvent ev, GFitEngineAdapter.Scratch s) { return new int[] {3, 1, 4}; }

        @Override
        public double[] beam() { return new double[] {-1.05, 0.0}; }
    }

    private static GFitTrack stubTrack(int id) {
        return new GFitTrack(id, new int[] {1, 2, 3}, new double[3], new double[3], new double[3], new int[3], new int[3],
            new int[3], 0, 0.3, 1, 0);
    }

    @Test
    public void fitThrows() {
        GFitTrack[] t = {stubTrack(1), stubTrack(2), stubTrack(3)};
        GFitEngineAdapter.Row[] rows = new GFitEngineAdapter.Row[3];
        GFitEvent ev = new GFitEvent();
        ev.ze = 70;
        GFitEngineAdapter.fitTracks(new Stub(), t, new int[] {0, 128, 0}, ev, true, new GFitEngineAdapter.Scratch(), rows);
        for (int k : new int[] {0, 2}) {
            assertEquals(k + 1, rows[k].trackid);
            assertEquals(0, rows[k].status);
            assertEquals(250.0, rows[k].p, 1e-12);
            assertEquals(1, rows[k].nAdded);
            assertEquals(3, rows[k].atofSector);
        }
        assertEquals(2, rows[1].trackid);
        assertEquals(4, rows[1].status);
        assertEquals(128, rows[1].flags);
        assertTrue(Double.isNaN(rows[1].p) && Double.isNaN(rows[1].cost));
        assertEquals(-1, rows[1].atofSector);
        assertEquals(70, ev.ze, 0);
    }

    @Test
    public void rowFromResult() {
        GFitResult r = new GFitResult(10);
        new Stub().fit(stubTrack(1), null, true, null, r);
        r.flags = 9;
        GFitEngineAdapter.Row w = GFitEngineAdapter.Row.of(5, r, 64, 2, new int[] {-1, -1, -1}, new double[] {-1.05, 0.0});
        assertEquals(9 | 64, w.flags);
        assertEquals(-1.05 - 0.5 * Math.sin(0.3), w.x, 1e-12);
        assertEquals(0.5 * Math.cos(0.3), w.y, 1e-12);
        assertEquals(-12, w.z, 0);
        assertEquals(250 * Math.sin(1.1) * Math.cos(0.3), w.px, 1e-9);
        assertEquals(250 * Math.cos(1.1), w.pz, 1e-9);
        assertEquals(0.04, w.cov[0], 0);       // c00
        assertEquals(0.001, w.cov[1], 0);      // c01
        assertEquals(0.002, w.cov[9], 0);      // c22 (ln p)
        assertTrue(Double.isNaN(w.pDedx) && Double.isNaN(w.sigPDedx));
        r.u[1] = 3 * Math.PI / 2;
        assertEquals(-Math.PI / 2, GFitEngineAdapter.Row.of(5, r, 0, 0, new int[] {-1, -1, -1}, new double[2]).phi0, 1e-12);
        r.status = 1;
        GFitEngineAdapter.Row bad = GFitEngineAdapter.Row.of(5, r, 0, 0, new int[] {-1, -1, -1}, new double[2]);
        assertTrue(Double.isNaN(bad.p) && Double.isNaN(bad.x) && Double.isNaN(bad.cov[0]));
        assertEquals(7.5, bad.cost, 0);
    }

    private static IndexedTable beam(double x, double y) {
        IndexedTable t = new IndexedTable(3, "x_offset:y_offset");
        t.addEntry(0, 0, 0);
        t.setDoubleValue(x, "x_offset", 0, 0, 0);
        t.setDoubleValue(y, "y_offset", 0, 0, 0);
        return t;
    }

    private static GFitEngineAdapter.Fits fits(double[] beamXY) {
        GFitFixtures.Fx f = GFitFixtures.resource("unit/field.txt.gz");
        double[] fp = f.d("field_par");
        FieldGrid field = FieldGrid.fromArrays(fp[0], fp[1], fp[2], fp[3], f.shape("br")[0], f.shape("br")[1], f.d("br"), f.d("bz"));
        return new GFitEngineAdapter.Fits(GFitConfig.pythonDefault().withSearch(GFitConfig.Search.FAST), beamXY, field, geometry(),
            new org.jlab.rec.ahdc.GlobalFit.AtofPositionTable[] {new GFitFixtureConstants(true).atofPosition(23015, 0), new GFitFixtureConstants(true).atofPosition(23015, 1)}, true);
    }

    @Test
    public void beamIsSum() {
        Map<String, IndexedTable> m = new HashMap<>();
        m.put(CcdbConstants.BEAM, beam(0.052, -0.023));
        m.put(CcdbConstants.BEAM_OFFSET, beam(-0.157, 0.023));
        double[] b = new CcdbConstants((run, name) -> m.get(name)).beamXY(23014);
        GFitEngineAdapter.Fits f = fits(b);
        assertEquals(-1.05, f.cfg.beamX(), 1e-12);
        assertEquals(0.0, f.cfg.beamY(), 1e-12);
        assertArrayEquals(b, f.beam(), 0);
    }

    @Test
    public void noElectron() {
        GFitReplaySet set = GFitReplaySet.resource("t60");
        GFitEngineAdapter.Fits f = fits(new double[] {-1.05, 0.0});
        GFitEvent ev = set.events[0];
        double ze = ev.ze;
        GFitTrack[] t = {set.tracks[0]};
        GFitEngineAdapter.Row[] with = new GFitEngineAdapter.Row[1], without = new GFitEngineAdapter.Row[1];
        GFitEngineAdapter.Scratch s = new GFitEngineAdapter.Scratch();
        GFitEngineAdapter.fitTracks(f, t, new int[] {0}, ev, true, s, with);
        GFitEngineAdapter.fitTracks(f, t, new int[] {GFitEngineAdapter.FLAG_NO_ELECTRON}, ev, false, s, without);
        assertEquals(ze, ev.ze, 0);
        assertEquals(0, with[0].status);
        assertEquals(0, without[0].status);
        assertTrue((without[0].flags & GFitEngineAdapter.FLAG_NO_ELECTRON) != 0);
        assertFalse((with[0].flags & GFitEngineAdapter.FLAG_NO_ELECTRON) != 0);
        assertTrue(Double.isFinite(without[0].p) && Double.isFinite(without[0].z0) && Double.isFinite(without[0].cov[0]));
        assertTrue(Arrays.toString(new double[] {with[0].p, without[0].p}), with[0].p > 100 && with[0].p < 1000);
    }

    /** Fits with the given switches and the fixture constants of run 23014 (time tables per half, dE/dx class_var_all). */
    private static GFitEngineAdapter.Fits switchedFits(boolean time, boolean dedx) {
        GFitFixtures.Fx f = GFitFixtures.resource("unit/field.txt.gz");
        double[] fp = f.d("field_par");
        FieldGrid field = FieldGrid.fromArrays(fp[0], fp[1], fp[2], fp[3], f.shape("br")[0], f.shape("br")[1], f.d("br"), f.d("bz"));
        GFitFixtureConstants fc = new GFitFixtureConstants(true);
        org.jlab.rec.ahdc.GlobalFit.AtofPositionTable[] tab = {fc.atofPosition(23014, 0), fc.atofPosition(23014, 1)};
        org.jlab.rec.ahdc.GlobalFit.AtofPositionTable[] tt = time
            ? new org.jlab.rec.ahdc.GlobalFit.AtofPositionTable[] {tab[0].withTime(fc.timeKeys(23014, 0)), tab[1].withTime(fc.timeKeys(23014, 1))}
            : null;
        GFitConfig cfg = GFitEngineAdapter.config(GFitConfig.Search.FAST, time, dedx);
        return new GFitEngineAdapter.Fits(cfg, new double[] {-1.05, 0.0}, field, geometry(), tab, tt, dedx ? fc.dedx(23014) : null, true);
    }

    /** t60 tracks as the engine would hand them over with the switches on: raw AHDC::adc ADC and ped per event row, fixture ATOF times. */
    private static GFitReplaySet switchSet() {
        GFitReplaySet set = GFitReplaySet.resource("t60");
        GFitFixtures.Fx hx = GFitFixtures.resource("switch/t60/hits_extra.txt.gz"), ax = GFitFixtures.resource("switch/t60/atof_extra.txt.gz");
        int[] hoff = set.tracksFx.i("hoff"), hadc = set.tracksFx.i("h_adc"), aoff = set.atofFx.i("aoff");
        double[] hped = hx.d("h_ped"), aE = ax.d("a_energy"), aT = ax.d("a_t_ge10");
        for (int k = 0; k < set.n; k++) {
            GFitEvent ev = set.events[k];
            GFitTrack t = set.tracks[k];
            int nrow = ev.eg.length;
            ev.rAdc = new int[nrow];
            ev.rPed = new double[nrow];
            for (int i = 0; i < t.nHits(); i++) {
                ev.rAdc[t.hrow[i]] = hadc[hoff[k] + i];
                ev.rPed[t.hrow[i]] = hped[hoff[k] + i];
            }
            ev.aEnergy = Arrays.copyOfRange(aE, aoff[k], aoff[k + 1]);
            ev.aTime = Arrays.copyOfRange(aT, aoff[k], aoff[k + 1]);
        }
        return set;
    }

    @Test
    public void switchInputs() {
        assertEquals(1, GFitEngineAdapter.half(23015));
        assertEquals(0, GFitEngineAdapter.half(-4));
        GFitReplaySet set = switchSet();
        int[] hoff = set.tracksFx.i("hoff"), hadc = set.tracksFx.i("h_adc");
        double[] hped = GFitFixtures.resource("switch/t60/hits_extra.txt.gz").d("h_ped");
        GFitEngineAdapter.Fits f = switchedFits(true, true);
        GFitEngineAdapter.Scratch s = new GFitEngineAdapter.Scratch();
        int n256 = 0, nfin = 0, n512 = 0;
        for (int k = 0; k < 20; k++) {
            GFitTrack t = set.tracks[k];
            GFitEngineAdapter.dedxInputs(t, set.events[k]);
            for (int i = 0; i < t.nHits(); i++) {
                assertEquals(hadc[hoff[k] + i], t.hAdc[i], 0);
                assertEquals(hped[hoff[k] + i], t.hPed[i], 0);
            }
            GFitEngineAdapter.Row[] row = new GFitEngineAdapter.Row[1];
            GFitEngineAdapter.fitTracks(f, new GFitTrack[] {t}, new int[] {0}, set.events[k], true, s, row);
            if ((row[0].flags & 256) != 0) n256++;
            if ((row[0].flags & 512) != 0) n512++;
            if (Double.isFinite(row[0].pDedx)) {
                nfin++;
                assertTrue("p_dedx in MeV " + row[0].pDedx, row[0].pDedx > 50 && row[0].pDedx < 2000);
                assertTrue(Double.isFinite(row[0].sigPDedx) && row[0].sigPDedx > 0);
            }
        }
        System.out.println("adapter switches on (20 t60 tracks): flag 256 on " + n256 + ", finite p_dedx " + nfin + ", flag 512 " + n512);
        assertTrue(n256 >= 15);
        assertTrue(nfin >= 15);
    }

    @Test
    public void switchesOffBank() {
        GFitReplaySet set = switchSet();
        GFitEngineAdapter.Fits off = fits(new double[] {-1.05, 0.0});
        GFitEngineAdapter.Fits offTables = switchedFits(false, false);
        assertTrue(off.sw == null && offTables.sw == null);
        GFitEngineAdapter.Scratch s = new GFitEngineAdapter.Scratch();
        for (int k = 0; k < 10; k++) {
            GFitTrack t = set.tracks[k];
            GFitEngineAdapter.Row[] a = new GFitEngineAdapter.Row[1], b = new GFitEngineAdapter.Row[1];
            GFitEngineAdapter.fitTracks(off, new GFitTrack[] {t}, new int[] {0}, set.events[k], true, s, a);
            GFitEngineAdapter.fitTracks(off, new GFitTrack[] {t}, new int[] {0}, set.events[k], true, s, b);
            assertTrue(Double.isNaN(a[0].pDedx) && Double.isNaN(a[0].sigPDedx));
            assertEquals(0, a[0].flags & (256 | 512));
            assertEquals(a[0].p, b[0].p, 0);
            assertEquals(a[0].cost, b[0].cost, 0);
        }
        assertEquals(GFitEngineAdapter.config(GFitConfig.Search.FAST).cf()[1], GFitEngineAdapter.config(GFitConfig.Search.FAST, false, false).cf()[1], 0);
        assertTrue(!GFitEngineAdapter.config(GFitConfig.Search.FAST, false, false).timeOn());
        assertEquals(new GFitEngineAdapter(GFitConfig.Search.FAST, true).requiredTables(),
            new GFitEngineAdapter(GFitConfig.Search.FAST, true, false, false).requiredTables());
        assertEquals(new GFitEngineAdapter(GFitConfig.Search.FAST, true).requiredTables().size() + 9,
            new GFitEngineAdapter(GFitConfig.Search.FAST, true, true, true).requiredTables().size());
    }
}
