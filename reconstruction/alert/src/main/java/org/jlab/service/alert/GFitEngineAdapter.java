package org.jlab.service.alert;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jlab.clas.swimtools.Swim;
import org.jlab.clas.swimtools.Swimmer;
import org.jlab.geom.detector.alert.AHDC.AlertDCDetector;
import org.jlab.geom.detector.alert.AHDC.AlertDCLayer;
import org.jlab.geom.prim.Line3D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.ahdc.GlobalFit.AtofPositionTable;
import org.jlab.rec.ahdc.GlobalFit.AtofTerm;
import org.jlab.rec.ahdc.GlobalFit.AtofTimeCorrection;
import org.jlab.rec.ahdc.GlobalFit.CcdbConstants;
import org.jlab.rec.ahdc.GlobalFit.CellGeometry;
import org.jlab.rec.ahdc.GlobalFit.FieldGrid;
import org.jlab.rec.ahdc.GlobalFit.GFitConfig;
import org.jlab.rec.ahdc.GlobalFit.GFitConstants;
import org.jlab.rec.ahdc.GlobalFit.GFitEvent;
import org.jlab.rec.ahdc.GlobalFit.GFitResult;
import org.jlab.rec.ahdc.GlobalFit.GFitTrack;
import org.jlab.rec.ahdc.GlobalFit.GFitWorkspace;
import org.jlab.rec.ahdc.GlobalFit.GlobalTrackFitter;
import org.jlab.rec.ahdc.GlobalFit.HitRecovery;
import org.jlab.rec.ahdc.GlobalFit.MaterialStack;
import org.jlab.rec.ahdc.GlobalFit.Propagator;
import org.jlab.rec.ahdc.GlobalFit.RecoveryPool;
import org.jlab.rec.ahdc.GlobalFit.StoppingTables;
import org.jlab.rec.ahdc.GlobalFit.SwitchConstants;
import org.jlab.rec.ahdc.GlobalFit.SwitchedFit;
import org.jlab.rec.ahdc.GlobalFit.TrackCost;
import org.jlab.rec.ahdc.Hit.Hit;
import org.jlab.rec.ahdc.Hit.HitReader;
import org.jlab.rec.alert.Track.Track;
import org.jlab.utils.groups.IndexedTable;

/**
 * Engine side of the global track fit (package org.jlab.rec.ahdc.GlobalFit): builds the fit inputs of one event from the
 * engine's objects and banks, runs the fit (and the hit recovery) for every Kalman-filter track and fills AHDC::gtrack.
 * Constructed by ALERTEngine only when the {@code gfit} key is on.
 * <p>
 * Inputs follow the reference dump and loader (gr_dump.py, fit.load_tracks): the hits of a track are its hits after the
 * second Kalman pass, in AHDC::adc row order, with the AHDC::hits doca and time, sigma = 2 x the Kalman-filter error of the
 * row's raw AHDC::adc ADC; the seed is the AHDC::kftrack state, or the AHDC::track helix fit where that is not finite (flag
 * 128); the electron vertex is the first pid-11 row of REC::Particle (none: z prior off, flag 64); the start time is
 * REC::Event startTime (a value &le; 0 counts as none); the ATOF rows are ATOF::hits with their bank times. The recovery pool
 * is every AHDC::adc row on no track, converted with the CCDB time offsets and {@link HitReader#timeToDistance}.
 * <p>
 * Per run: field grid sampled from the solenoid map as MAGFIELDS configured it, wire geometry from the engine's detector,
 * beam position and ATOF tables from CCDB. A missing table skips the fit for the run (no AHDC::gtrack bank), logged once.
 * <p>
 * Switches (plan P4, off by default): with {@code gfitTime} the per-half time tables (gfit_time_keys) and the corrected wedge
 * times of {@link AtofTimeCorrection} (ATOF::hits time and energy, upstream time_offsets, gfit_time_channel, gfit_time_walk);
 * with {@code gfitDedx} the per-hit AHDC::adc ADC and ped and the dE/dx constants, written as p_dedx, sig_p_dedx and flag 512.
 * With both off none of this is read or built and the fit runs exactly as without the switches.
 */
final class GFitEngineAdapter {

    static final Logger LOGGER = Logger.getLogger(GFitEngineAdapter.class.getName());

    /** AHDC layers (10 superlayer + layer) in the layer order of the fit's wire table. */
    static final int[] LAYERS = {11, 21, 22, 31, 32, 41, 42, 51};
    static final String T0_TABLE = "/calibration/alert/ahdc/time_offsets", T2D_TABLE = "/calibration/alert/ahdc/time_to_distance_wire";
    static final int FLAG_NO_ELECTRON = 64, FLAG_HELIX_SEED = 128;
    /** field grid of the reference (prop.load_field): rho 0-120 mm, z_map -400-400 mm, 5 mm steps. */
    static final int FIELD_NR = 25, FIELD_NZ = 161;
    static final double FIELD_STEP = 5.0, FIELD_ZMIN = -400.0;

    private final GFitConfig cfg;
    private final boolean recovery;
    private volatile RunContext ctx;
    /** runs already reported (missing table or failure), to log once per run. */
    private final Set<Integer> reported = java.util.Collections.synchronizedSet(new HashSet<>());
    private final ThreadLocal<Scratch> scratch = ThreadLocal.withInitial(Scratch::new);

    GFitEngineAdapter(GFitConfig.Search search, boolean recovery) {
        this(search, recovery, false, false);
    }

    GFitEngineAdapter(GFitConfig.Search search, boolean recovery, boolean time, boolean dedx) {
        this.cfg = config(search, time, dedx);
        this.recovery = recovery;
    }

    /** the fit configuration of the engine's switches, built here only (the per-run beam is added in {@link Fits}). */
    static GFitConfig config(GFitConfig.Search search) {
        return config(search, false, false);
    }

    static GFitConfig config(GFitConfig.Search search, boolean time, boolean dedx) {
        GFitConfig c = GFitConfig.pythonDefault().withSearch(search);
        return time || dedx ? c.withSwitches(time, dedx) : c;
    }

    /** CCDB tables the fit needs, as ReconstructionEngine.requireConstants takes them; the per-run check uses the same list. */
    Map<String, Integer> requiredTables() {
        return GFitConstants.requiredTables(cfg, recovery);
    }

    // ================================================================== per run

    /** the fit of one track as fitTracks uses it; {@link Fits} in the engine, a stub in tests. */
    interface TrackFit {
        /** fits t (and recovers hits) into out; returns the rows added by recovery. */
        int fit(GFitTrack t, GFitEvent ev, boolean zPrior, GFitWorkspace ws, GFitResult out);

        /** (sector, layer, component) of the ATOF hit associated with the fitted track, -1 when none. */
        int[] atof(GFitResult r, GFitEvent ev, Scratch s);

        /** beam position (x, y) in mm. */
        double[] beam();
    }

    /** fitters, recovery and constants of one run; immutable. */
    static final class Fits implements TrackFit {
        final GFitConfig cfg;
        final Propagator prop;
        final CellGeometry geo;
        /** ATOF position table per event half (RUN::config event % 2). */
        final AtofPositionTable[] table;
        /** index 0 with the electron-vertex z prior, 1 without. */
        final GlobalTrackFitter[] fitter = new GlobalTrackFitter[2];
        final HitRecovery[] rec = new HitRecovery[2];
        /** switched terms: time table per half (null with gfitTime off), dE/dx constants (null with gfitDedx off), the
         * switched sequence per z-prior index (null with both switches off). */
        final AtofPositionTable[] tableTime;
        final SwitchConstants.DedxConstants dx;
        final SwitchedFit[] sw;

        Fits(GFitConfig base, double[] beam, FieldGrid field, CellGeometry geo, AtofPositionTable[] table, boolean recovery) {
            this(base, beam, field, geo, table, null, null, recovery);
        }

        Fits(GFitConfig base, double[] beam, FieldGrid field, CellGeometry geo, AtofPositionTable[] table, AtofPositionTable[] tableTime,
                SwitchConstants.DedxConstants dx, boolean recovery) {
            this.cfg = base.withBeam(beam[0], beam[1]);
            this.geo = geo;
            this.table = table;
            this.tableTime = tableTime;
            this.dx = dx;
            this.prop = new Propagator(field, MaterialStack.dataBest(), StoppingTables.proton(), beam[0], beam[1]);
            boolean switched = cfg.timeOn() || cfg.dedxOn();
            this.sw = switched ? new SwitchedFit[2] : null;
            for (int z = 0; z < 2; z++) {
                GFitConfig c = cfg.withZPrior(z == 0);
                TrackCost tc = new TrackCost(prop, geo, c);
                fitter[z] = new GlobalTrackFitter(tc, prop, geo, c);
                if (recovery) rec[z] = new HitRecovery(fitter[z], tc, geo, c);
                if (switched) sw[z] = new SwitchedFit(fitter[z], rec[z], prop, geo, StoppingTables.proton());
            }
        }

        @Override
        public int fit(GFitTrack t, GFitEvent ev, boolean zPrior, GFitWorkspace ws, GFitResult out) {
            int z = zPrior ? 0 : 1;
            AtofPositionTable tab = table[ev.half];
            if (sw != null) return sw[z].run(t, ev, tab, tableTime == null ? null : tableTime[ev.half], true, dx, ws, out);
            fitter[z].fit(t, ev, tab, true, ws, out);
            return rec[z] != null ? rec[z].recover(t, ev, tab, true, out, ws, out) : 0;
        }

        /**
         * ATOF wedge hit associated at the fitted state: the most probable wedge candidate of atof.diag_nb (atof_run.py
         * atof_comp, INF[:, 1, 1]), -1 when there is none or the fit did not converge (status other than 0).
         */
        @Override
        public int[] atof(GFitResult r, GFitEvent ev, Scratch s) {
            if (r.status != 0) return new int[] {-1, -1, -1};
            AtofTerm.TrackState st = s.ws.atof;
            AtofTerm.associate(r.u, prop, s.ws, st, s.rec, s.assoc);
            int c = (int) s.assoc[3];
            if (c < 0) return new int[] {-1, -1, -1};
            int row = st.crow[c];
            return new int[] {ev.sector[row], ev.layer[row], ev.comp[row]};
        }

        @Override
        public double[] beam() { return new double[] {cfg.beamX(), cfg.beamY()}; }
    }

    static final class RunContext {
        final int run;
        final AlertDCDetector det;
        /** missing CCDB tables; the fit is skipped when not empty. */
        final List<String> missing;
        final Fits fits;
        final IndexedTable t0, t2d;
        /** wedge time constants of the run (null with gfitTime off). */
        final SwitchConstants.TimeChannels tc;

        RunContext(int run, AlertDCDetector det, List<String> missing, Fits fits, IndexedTable t0, IndexedTable t2d,
                SwitchConstants.TimeChannels tc) {
            this.run = run; this.det = det; this.missing = missing; this.fits = fits; this.t0 = t0; this.t2d = t2d; this.tc = tc;
        }
    }

    private RunContext context(int run, AlertDCDetector det, CcdbConstants.TableSource src) {
        RunContext c = ctx;
        if (c != null && c.run == run && c.det == det) return c;
        synchronized (this) {
            c = ctx;
            if (c != null && c.run == run && c.det == det) return c;
            List<String> missing = new ArrayList<>();
            for (String name : requiredTables().keySet()) if (src.table(run, name) == null) missing.add(name);
            Fits fits = null;
            IndexedTable t0 = null, t2d = null;
            SwitchConstants.TimeChannels tc = null;
            if (missing.isEmpty()) {
                CcdbConstants consts = new CcdbConstants(src, cfg);
                AtofPositionTable[] tab = {consts.atofPosition(run, 0), consts.atofPosition(run, 1)};
                AtofPositionTable[] tabTime = null;
                if (cfg.timeOn()) {
                    tabTime = new AtofPositionTable[] {tab[0].withTime(consts.timeKeys(run, 0)), tab[1].withTime(consts.timeKeys(run, 1))};
                    tc = consts.timeChannels(run);
                }
                SwitchConstants.DedxConstants dx = cfg.dedxOn() ? consts.dedx(run) : null;
                fits = new Fits(cfg, consts.beamXY(run), fieldFromSwim(new Swim()), geometry(det), tab, tabTime, dx, recovery);
                if (recovery) {
                    t0 = src.table(run, T0_TABLE);
                    t2d = src.table(run, T2D_TABLE);
                }
            } else if (reported.add(run)) {
                LOGGER.severe("ALERT global fit skipped for run " + run + " (no AHDC::gtrack): CCDB table(s) missing: "
                        + String.join(", ", missing));
            }
            c = new RunContext(run, det, missing, fits, t0, t2d, tc);
            ctx = c;
            return c;
        }
    }

    /** the solenoid map of MAGFIELDS (its scale and z shift) on the reference grid; z_map = z + zsh with zsh = -10 x shift [cm]. */
    static FieldGrid fieldFromSwim(Swim swim) {
        double zsh = -10.0 * Swimmer.get_zShift();
        double[] br = new double[FIELD_NR * FIELD_NZ], bz = new double[FIELD_NR * FIELD_NZ];
        float[] b = new float[3];
        for (int i = 0; i < FIELD_NR; i++) {
            for (int j = 0; j < FIELD_NZ; j++) {
                double zmap = FIELD_ZMIN + j * FIELD_STEP;
                swim.BfieldLab(i * FIELD_STEP / 10.0, 0.0, (zmap - zsh) / 10.0, b);
                br[i * FIELD_NZ + j] = b[0];
                bz[i * FIELD_NZ + j] = b[2];
            }
        }
        return FieldGrid.fromArrays(FIELD_STEP, FIELD_STEP, FIELD_ZMIN, zsh, FIELD_NR, FIELD_NZ, br, bz);
    }

    /** wire table of the fit from the engine's AHDC geometry: layers in {@link #LAYERS} order, wires 1..n of each layer. */
    static CellGeometry geometry(AlertDCDetector det) {
        int[] layoff = new int[LAYERS.length + 1];
        for (int li = 0; li < LAYERS.length; li++) layoff[li + 1] = layoff[li] + layer(det, li).getNumComponents();
        int n = layoff[LAYERS.length];
        double[] wo = new double[3 * n], we = new double[3 * n];
        for (int li = 0; li < LAYERS.length; li++) {
            AlertDCLayer layer = layer(det, li);
            for (int w = 0; w < layoff[li + 1] - layoff[li]; w++) {
                Line3D line = layer.getComponent(w + 1).getLine();
                int g = layoff[li] + w;
                wo[3 * g] = line.origin().x(); wo[3 * g + 1] = line.origin().y(); wo[3 * g + 2] = line.origin().z();
                we[3 * g] = line.end().x(); we[3 * g + 1] = line.end().y(); we[3 * g + 2] = line.end().z();
            }
        }
        return CellGeometry.fromWireEnds(wo, we, layoff);
    }

    private static AlertDCLayer layer(AlertDCDetector det, int li) {
        return det.getSector(1).getSuperlayer(LAYERS[li] / 10).getLayer(LAYERS[li] % 10);
    }

    // ================================================================== per event

    /**
     * Fits every Kalman-filter track of the event and appends AHDC::gtrack (one row per kftrack row, same order); appends
     * nothing when a table of the run is missing.
     */
    void process(DataEvent event, int run, AlertDCDetector det, CcdbConstants.TableSource src, List<Track> tracks,
            DataBank kfBank, DataBank trackBank, DataBank runBank, boolean mc) {
        try {
            fill(event, run, det, src, tracks, kfBank, trackBank, runBank, mc);
        } catch (RuntimeException e) {      // the fit must never change the other banks (plan D14)
            if (reported.add(run)) LOGGER.log(Level.SEVERE, "ALERT global fit failed for an event of run " + run + "; no AHDC::gtrack for it", e);
        }
    }

    private void fill(DataEvent event, int run, AlertDCDetector det, CcdbConstants.TableSource src, List<Track> tracks,
            DataBank kfBank, DataBank trackBank, DataBank runBank, boolean mc) {
        RunContext c = context(run, det, src);
        if (!c.missing.isEmpty()) return;
        Row[] rows = new Row[tracks.size()];
        if (!tracks.isEmpty()) {
            double[] el = electron(event);
            double ze = el == null ? Double.NaN : zElectron(el[0], runBank.getFloat("torus", 0), mc);
            double start = mc ? 0.0 : startTime(event);
            GFitEvent ev = event(event, ze, start, c.fits.geo);
            ev.half = half(runBank.getInt("event", 0));
            if (cfg.timeOn()) {
                ev.aEnergy = atofEnergy(event);
                AtofTimeCorrection.correct(ev, c.tc, run);
            }
            GFitTrack[] in = new GFitTrack[tracks.size()];
            int[] flags = new int[tracks.size()];
            boolean[] inTrack = new boolean[ev.eg.length];
            for (int k = 0; k < tracks.size(); k++) {
                Track t = tracks.get(k);
                double[] seed = seed(kfBank, k, trackBank, t.get_trackId());
                flags[k] = (seed[4] != 0 ? FLAG_HELIX_SEED : 0) | (el == null ? FLAG_NO_ELECTRON : 0);
                in[k] = track(t.get_trackId(), t.getHits(), ev, c.fits.geo, cfg.errScale(), seed, inTrack);
                if (cfg.dedxOn()) dedxInputs(in[k], ev);
            }
            if (recovery) pool(ev, inTrack, cfg, c.t0, c.t2d);
            fitTracks(c.fits, in, flags, ev, el != null, scratch.get(), rows);
        }
        DataBank bank = event.createBank("AHDC::gtrack", rows.length);
        for (int k = 0; k < rows.length; k++) rows[k].write(bank, k);
        event.appendBank(bank);
    }

    /** ATOF::hits energy per row, in the row order of {@link #event} (input of the corrected wedge time, gfitTime only). */
    static double[] atofEnergy(DataEvent event) {
        DataBank a = event.hasBank("ATOF::hits") ? event.getBank("ATOF::hits") : null;
        int m = a == null ? 0 : a.rows();
        double[] e = new double[m];
        for (int r = 0; r < m; r++) e[r] = a.getFloat("energy", r);
        return e;
    }

    /** per-hit raw AHDC::adc ADC and ped of the track's hits through their AHDC::adc row (gfitDedx only). */
    static void dedxInputs(GFitTrack t, GFitEvent ev) {
        int n = t.nHits();
        t.hAdc = new double[n];
        t.hPed = new double[n];
        for (int i = 0; i < n; i++) {
            t.hAdc[i] = ev.rAdc[t.hrow[i]];
            t.hPed[i] = ev.rPed[t.hrow[i]];
        }
    }

    /** vz [cm] of the first pid-11 row of REC::Particle, null when there is none (as the Kalman filter's vertex). */
    static double[] electron(DataEvent event) {
        if (!event.hasBank("REC::Particle")) return null;
        DataBank b = event.getBank("REC::Particle");
        for (int r = 0; r < b.rows(); r++) if (b.getInt("pid", r) == 11) return new double[] {b.getFloat("vz", r)};
        return null;
    }

    /** electron vertex in the AHDC frame [mm]: 10 vz + 10 z_shift, z_shift 6.7 cm (torus &gt; 0) or 3.6 cm, 0 in simulation (gr_dump.py tag_zv). */
    static double zElectron(double vzCm, double torus, boolean mc) {
        return 10 * vzCm + 10 * (mc ? 0.0 : torus > 0 ? 6.7 : 3.6);
    }

    /** table half of an event: event number parity (the half split of the half-indexed constants). */
    static int half(int eventNumber) {
        return Math.floorMod(eventNumber, 2);
    }

    /** REC::Event startTime, NaN when absent or not positive. */
    static double startTime(DataEvent event) {
        if (!event.hasBank("REC::Event")) return Double.NaN;
        DataBank b = event.getBank("REC::Event");
        if (b.rows() == 0) return Double.NaN;
        double t = b.getFloat("startTime", 0);
        return t > 0 ? t : Double.NaN;
    }

    /** event inputs: AHDC::adc rows (global wire, raw columns), ATOF::hits rows, ze, start time. */
    static GFitEvent event(DataEvent event, double ze, double start, CellGeometry geo) {
        GFitEvent ev = new GFitEvent();
        ev.ze = ze;
        ev.startTime = start;
        DataBank adc = event.getBank("AHDC::adc");
        int n = adc == null ? 0 : adc.rows();
        int[] layer = new int[n], comp = new int[n];
        ev.rAdc = new int[n]; ev.rWf = new int[n]; ev.rLet = new double[n]; ev.rTot = new double[n]; ev.rPed = new double[n];
        for (int r = 0; r < n; r++) {
            layer[r] = adc.getByte("layer", r);
            comp[r] = adc.getShort("component", r);
            ev.rAdc[r] = adc.getInt("ADC", r);
            ev.rWf[r] = adc.getShort("wfType", r);
            ev.rLet[r] = adc.getFloat("leadingEdgeTime", r);
            ev.rTot[r] = adc.getFloat("timeOverThreshold", r);
            ev.rPed[r] = adc.getFloat("ped", r);
        }
        wires(layer, comp, geo.layoff, ev);
        DataBank a = event.hasBank("ATOF::hits") ? event.getBank("ATOF::hits") : null;
        int m = a == null ? 0 : a.rows();
        ev.ax = new double[m]; ev.ay = new double[m]; ev.az = new double[m]; ev.tBank = new double[m];
        ev.isBar = new boolean[m]; ev.sector = new int[m]; ev.layer = new int[m]; ev.comp = new int[m];
        for (int r = 0; r < m; r++) {
            ev.ax[r] = a.getFloat("x", r); ev.ay[r] = a.getFloat("y", r); ev.az[r] = a.getFloat("z", r);
            ev.tBank[r] = a.getFloat("time", r);
            ev.sector[r] = a.getInt("sector", r); ev.layer[r] = a.getInt("layer", r); ev.comp[r] = a.getInt("component", r);
            ev.isBar[r] = ev.comp[r] == 10;
        }
        ev.tFinal = ev.tBank;
        return ev;
    }

    /**
     * Layer index, wire index and global wire of AHDC::adc rows from the bank's layer byte (10 superlayer + layer) and
     * component (wire 1..n): rLi, rWire (component - 1) and eg; a row whose layer or wire is not in the table gets -1 in all three.
     */
    static void wires(int[] layer, int[] comp, int[] layoff, GFitEvent ev) {
        int n = layer.length;
        ev.eg = new int[n]; ev.rLi = new int[n]; ev.rWire = new int[n];
        for (int r = 0; r < n; r++) {
            int li = -1;
            for (int k = 0; k < LAYERS.length; k++) if (LAYERS[k] == layer[r]) li = k;
            int w = comp[r] - 1;
            if (li < 0 || w < 0 || w >= layoff[li + 1] - layoff[li]) {
                ev.eg[r] = -1; ev.rLi[r] = -1; ev.rWire[r] = -1;
            } else {
                ev.eg[r] = layoff[li] + w; ev.rLi[r] = li; ev.rWire[r] = w;
            }
        }
    }

    /**
     * Seed (phi rad, p GeV, theta rad, z mm, helix flag) of kftrack row k, or of the AHDC::track row of trackid when the
     * kftrack state is not finite; angles through degrees as fit.load_tracks reads them from the dump.
     */
    static double[] seed(DataBank kf, int k, DataBank helix, int trackid) {
        double[] s = seed(kf.getFloat("px", k), kf.getFloat("py", k), kf.getFloat("pz", k), kf.getFloat("z", k));
        if (finite(s)) return new double[] {s[0], s[1], s[2], s[3], 0};
        if (helix != null) {
            for (int r = 0; r < helix.rows(); r++) {
                if (helix.getInt("trackid", r) != trackid) continue;
                double[] h = seed(helix.getFloat("px", r), helix.getFloat("py", r), helix.getFloat("pz", r), helix.getFloat("z", r));
                return new double[] {h[0], h[1], h[2], h[3], 1};
            }
        }
        return new double[] {s[0], s[1], s[2], s[3], 0};
    }

    /** (phi, p, theta, z) from a momentum in MeV and z in mm (gr_dump.py b0 rows, fit.load_tracks TC). */
    static double[] seed(double px, double py, double pz, double z) {
        double phi = Math.toRadians(Math.toDegrees(Math.atan2(py, px)));
        double theta = Math.toRadians(Math.toDegrees(Math.atan2(Math.hypot(px, py), pz)));
        return new double[] {phi, Math.sqrt(px * px + py * py + pz * pz) / 1000.0, theta, z};
    }

    private static boolean finite(double[] s) {
        for (int i = 0; i < 4; i++) if (!Double.isFinite(s[i])) return false;
        return true;
    }

    /**
     * Fit inputs of one track (fit.load_tracks): its hits in AHDC::adc row order (row = hit id - 1; a row taken once), rows
     * with a non-finite doca or time or an unknown wire left out; doca, smear class (time &lt; -20: 2, &lt; 0: 1), sigma =
     * errScale x the Kalman-filter error of the row's raw ADC, half pitch, layer class. Marks the rows in inTrack.
     */
    static GFitTrack track(int trackid, List<Hit> hits, GFitEvent ev, CellGeometry geo, double errScale, double[] seed,
            boolean[] inTrack) {
        int nh = hits.size();
        long[] o = new long[nh];      // (id, index in the list), sorted by id
        for (int i = 0; i < nh; i++) o[i] = ((long) hits.get(i).getId() << 32) | i;
        Arrays.sort(o);
        int[] hg = new int[nh], hism = new int[nh], hcls = new int[nh], hrow = new int[nh];
        double[] hm = new double[nh], hsig = new double[nh], hh = new double[nh];
        int n = 0, last = -1;
        for (int i = 0; i < nh; i++) {
            Hit h = hits.get((int) o[i]);
            int row = h.getId() - 1;
            if (row == last) continue;
            last = row;
            if (row < 0 || row >= ev.eg.length || ev.eg[row] < 0) continue;
            inTrack[row] = true;      // a track row is never a recovery candidate, even when the fit leaves it out
            double doca = h.getDoca(), time = h.getTime();
            if (!Double.isFinite(doca) || !Double.isFinite(time)) continue;
            int g = ev.eg[row];
            hg[n] = g; hm[n] = doca; hism[n] = time < -20 ? 2 : time < 0 ? 1 : 0;
            hsig[n] = errScale * HitRecovery.kfError(ev.rAdc[row], time);
            hh[n] = geo.half[g]; hcls[n] = ev.rLi[row]; hrow[n] = row;
            n++;
        }
        return new GFitTrack(trackid, Arrays.copyOf(hg, n), Arrays.copyOf(hm, n), Arrays.copyOf(hsig, n), Arrays.copyOf(hh, n),
                Arrays.copyOf(hism, n), Arrays.copyOf(hcls, n), Arrays.copyOf(hrow, n), seed[0], seed[1], seed[2], seed[3]);
    }

    /**
     * Recovery pool of the event (plan D11): AHDC::adc rows on no track and with a known wire, time = leading edge - t0 - start
     * through the CCDB time offsets and {@link HitReader#timeToDistance}; empty without a start time.
     */
    static void pool(GFitEvent ev, boolean[] inTrack, GFitConfig cfg, IndexedTable t0, IndexedTable t2d) {
        boolean[] skip = inTrack.clone();
        for (int r = 0; r < skip.length; r++) if (ev.eg[r] < 0 || !Double.isFinite(ev.startTime)) skip[r] = true;
        RecoveryPool.fromRawRows(ev, cfg, skip,
                (li, w) -> t0.getDoubleValue("t0", 1, LAYERS[li], w + 1),
                (li, w, time) -> HitReader.timeToDistance(1, LAYERS[li], w + 1, time, t2d));
    }

    // ================================================================== fit and rows

    /** per-thread workspace and result, grown with the largest track. */
    static final class Scratch {
        int cap = -1;
        GFitWorkspace ws;
        GFitResult r;
        final double[] rec = new double[20], assoc = new double[4];

        /** room for maxHits hits plus the rows recovery can add (two rounds of at most one per layer). */
        void ensure(int maxHits) {
            int need = maxHits + 2 * LAYERS.length;
            if (need <= cap) return;
            cap = Math.max(need, 64);
            ws = new GFitWorkspace(cap);
            r = new GFitResult(cap);
        }
    }

    /**
     * Fits tracks[k] (z prior when the event has an electron; without one the grid starts at the seed z, or 0) into rows[k];
     * an exception inside one track gives status 4 for it and leaves the others alone (plan D14).
     */
    static void fitTracks(TrackFit fit, GFitTrack[] tracks, int[] flags, GFitEvent ev, boolean zPrior, Scratch s, Row[] rows) {
        int maxHits = 0;
        for (GFitTrack t : tracks) maxHits = Math.max(maxHits, t.nHits());
        s.ensure(maxHits);
        double ze = ev.ze;
        for (int k = 0; k < tracks.length; k++) {
            GFitTrack t = tracks[k];
            ev.ze = zPrior ? ze : Double.isFinite(t.seedZ) ? t.seedZ : 0.0;
            try {
                s.r.clear();
                int added = fit.fit(t, ev, zPrior, s.ws, s.r);
                rows[k] = Row.of(t.trackid, s.r, flags[k], added, fit.atof(s.r, ev, s), fit.beam());
            } catch (RuntimeException e) {
                LOGGER.log(Level.FINE, "ALERT global fit failed for track " + t.trackid, e);
                GFitResult bad = s.r;
                bad.clear();
                bad.status = 4;
                rows[k] = Row.of(t.trackid, bad, flags[k], 0, new int[] {-1, -1, -1}, fit.beam());
            }
        }
        ev.ze = ze;
    }

    /** one AHDC::gtrack row (bank definition in etc/bankdefs/hipo4/alert.json). */
    static final class Row {
        int trackid, status, flags, nHits, nAdded, nIter, nMin, atofSector, atofLayer, atofComp;
        double x, y, z, px, py, pz, d0, phi0, p, z0, theta, cost, pDedx = Double.NaN, sigPDedx = Double.NaN;
        final double[] cov = new double[15];

        /** row of a fit result: state, momentum and covariance only for status 0 (NaN otherwise, as fit.pack_est). */
        static Row of(int trackid, GFitResult r, int flags, int added, int[] atof, double[] beam) {
            Row w = new Row();
            w.trackid = trackid;
            w.status = r.status;
            w.flags = r.flags | flags;
            w.cost = r.status == 2 || r.status == 4 ? Double.NaN : r.cost;
            w.nHits = r.nh; w.nAdded = added; w.nIter = r.niter; w.nMin = r.nmin;
            w.pDedx = 1000.0 * r.pDedx; w.sigPDedx = 1000.0 * r.sigPDedx;      // GeV to MeV as p; NaN unless gfitDedx
            w.atofSector = atof[0]; w.atofLayer = atof[1]; w.atofComp = atof[2];
            if (r.status != 0) {
                w.x = w.y = w.z = w.px = w.py = w.pz = w.d0 = w.phi0 = w.p = w.z0 = w.theta = Double.NaN;
                Arrays.fill(w.cov, Double.NaN);
                return w;
            }
            double[] u = r.u;
            double ph = u[1], th = u[4], pm = 1000.0 * Math.exp(u[2]);
            w.d0 = u[0];
            w.phi0 = ph - 2 * Math.PI * Math.floor((ph + Math.PI) / (2 * Math.PI));
            w.p = pm; w.z0 = u[3]; w.theta = th;
            w.x = beam[0] - u[0] * Math.sin(ph); w.y = beam[1] + u[0] * Math.cos(ph); w.z = u[3];
            w.px = pm * Math.sin(th) * Math.cos(ph); w.py = pm * Math.sin(th) * Math.sin(ph); w.pz = pm * Math.cos(th);
            int q = 0;
            for (int i = 0; i < 5; i++) for (int j = i; j < 5; j++) w.cov[q++] = r.cov[i * GFitConfig.NPAR + j];
            return w;
        }

        static final String[] COV = {"c00", "c01", "c02", "c03", "c04", "c11", "c12", "c13", "c14", "c22", "c23", "c24", "c33", "c34", "c44"};

        void write(DataBank b, int k) {
            b.setInt("trackid", k, trackid);
            b.setByte("status", k, (byte) status);
            b.setShort("flags", k, (short) flags);
            b.setFloat("x", k, (float) x); b.setFloat("y", k, (float) y); b.setFloat("z", k, (float) z);
            b.setFloat("px", k, (float) px); b.setFloat("py", k, (float) py); b.setFloat("pz", k, (float) pz);
            b.setFloat("d0", k, (float) d0); b.setFloat("phi0", k, (float) phi0); b.setFloat("p", k, (float) p);
            b.setFloat("z0", k, (float) z0); b.setFloat("theta", k, (float) theta);
            for (int i = 0; i < 15; i++) b.setFloat(COV[i], k, (float) cov[i]);
            b.setFloat("cost", k, (float) cost);
            b.setShort("n_hits", k, (short) nHits); b.setShort("n_added", k, (short) nAdded);
            b.setShort("n_iter", k, (short) nIter); b.setShort("n_min", k, (short) nMin);
            b.setByte("atof_sector", k, (byte) atofSector); b.setByte("atof_layer", k, (byte) atofLayer);
            b.setByte("atof_component", k, (byte) atofComp);
            b.setFloat("p_dedx", k, (float) pDedx); b.setFloat("sig_p_dedx", k, (float) sigPDedx);
        }
    }
}
