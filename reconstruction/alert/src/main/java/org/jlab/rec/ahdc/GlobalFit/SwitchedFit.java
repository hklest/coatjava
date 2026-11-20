package org.jlab.rec.ahdc.GlobalFit;

import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The whole fit sequence of one track with the two switches, as one call for the replay and the engine adapter:
 * the fit (the default fit, or the time fit with the time table), the hit recovery when a {@link HitRecovery} is given, and with the
 * dE/dx switch on the post-fit dE/dx step.
 *
 * <p>The dE/dx step takes the hit geometry and the prompt wedge class from the default fit before recovery (the default table, as
 * pinned comb_run.py: DX:366-442, DR:35-39), and as base momentum the fit the sequence ends with (comb_run.py:289-291: the time fit
 * plus recovery for ge10_all, the default fit for ge10_dedx without recovery). Time off: the default fit is the fit; its hit
 * quantities and prompt class are taken straight after it, before recovery overwrites nothing of them. Time on: the default fit
 * runs first into a scratch result (about one more fit per track, only with both switches on) and the time fit follows. The
 * prompt class is {@code wedge posterior > 0.5, |time residual| < 1 ns, default-fit status 0}, the posterior and best candidate
 * from {@link AtofTerm#associate} at the default-fit state, the residual {@code candidate time - start - predicted entry time -
 * off} of the default table with the bank or final times the caller selects. The outputs go to {@link GFitResult#pDedx} and
 * related fields and flag {@link DedxProfile#FLAG_FEWHITS} (512) into {@code flags}; state, covariance and p are never changed.
 * With the dE/dx switch off none of this runs and nothing is allocated. One instance may serve many threads: all scratch is in the
 * {@link GFitWorkspace} ({@link Scratch}).
 */
public final class SwitchedFit {

	private static final Logger LOGGER = Logger.getLogger(SwitchedFit.class.getName());

	private final GlobalTrackFitter fitter;
	private final HitRecovery rec;
	private final Propagator prop;
	private final CellGeometry geo;
	private final StoppingTables st;
	private final boolean timeOn, dedxOn;

	/** Per-thread dE/dx scratch, created on the first track with the switch on. */
	public static final class Scratch {
		final DedxHits hits;
		final DedxProfile prof;
		/** propagator and wire table the hit geometry was built with: the engine builds new ones per run. */
		final Propagator prop;
		final CellGeometry geo;
		final GFitResult dflt;
		final double[] x, lnu, rec = new double[AtofTerm.NI], assoc = new double[4];
		final boolean[] use;
		int nKept;

		Scratch(Propagator prop, CellGeometry geo, StoppingTables st, int maxHits) {
			this.prop = prop;
			this.geo = geo;
			hits = new DedxHits(prop, geo, st, maxHits);
			prof = new DedxProfile(st);
			dflt = new GFitResult(maxHits);
			int n = maxHits + GFitConfig.NADD;
			x = new double[n]; lnu = new double[n]; use = new boolean[n];
		}
	}

	/**
	 * @param rec the hit recovery of the sequence, null for none
	 * @param st the stopping tables of the propagator
	 */
	public SwitchedFit(GlobalTrackFitter fitter, HitRecovery rec, Propagator prop, CellGeometry geo, StoppingTables st) {
		this.fitter = fitter;
		this.rec = rec;
		this.prop = prop;
		this.geo = geo;
		this.st = st;
		timeOn = fitter.timeOn();
		dedxOn = fitter.dedxOn();
	}

	/**
	 * Runs the sequence for one track into out (cleared by the fit). The caller prepares: {@code tabDefault}, the ATOF table of the
	 * event's half of the default configuration (bank or final time keys as {@code bankTimes} says); {@code tabTime}, the time table
	 * of that half ({@link AtofPositionTable#withTime}) and {@code ev.aTime}, only with the time switch; {@code dx} and
	 * {@code t.hAdc}, {@code t.hPed} (and {@code ev.half}), only with the dE/dx switch (the doca is {@code t.hm}). Never throws for
	 * a failed fit (status 4 as {@link GlobalTrackFitter#fit}); an exception in the dE/dx step leaves its outputs NaN.
	 *
	 * @return rows added by the recovery (0 without one)
	 */
	public int run(GFitTrack t, GFitEvent ev, AtofPositionTable tabDefault, AtofPositionTable tabTime, boolean bankTimes,
			SwitchConstants.DedxConstants dx, GFitWorkspace ws, GFitResult out) {
		final AtofPositionTable tab = timeOn ? tabTime : tabDefault;
		if (!dedxOn) {
			fitter.fit(t, ev, tab, bankTimes, ws, out);
			return rec != null ? rec.recover(t, ev, tab, bankTimes, out, ws, out) : 0;
		}
		if (t.hAdc == null || t.hPed == null) throw new IllegalArgumentException("dE/dx switch on: track without hAdc/hPed");
		Scratch s = ws.sw;
		if (s == null || s.prop != prop || s.geo != geo) s = ws.sw = new Scratch(prop, geo, st, ws.cap - GFitConfig.NADD - GFitConfig.NMISS);
		final GFitResult d = timeOn ? s.dflt : out;
		fitter.fit(t, ev, tabDefault, bankTimes, ws, d);
		boolean prompt = false, ok = true;
		try {
			prompt = hitsAndPrompt(t, ev, bankTimes, dx, ws, s, d);
		} catch (RuntimeException e) {
			ok = false;
			LOGGER.log(Level.FINE, "dE/dx hit step failed for track " + t.trackid, e);
		}
		if (timeOn) fitter.fit(t, ev, tabTime, bankTimes, ws, out);
		final int added = rec != null ? rec.recover(t, ev, tab, bankTimes, out, ws, out) : 0;
		if (ok) {
			try {
				profile(ev, dx, s, out, prompt);
			} catch (RuntimeException e) {
				LOGGER.log(Level.FINE, "dE/dx profile failed for track " + t.trackid, e);
				out.pDedx = Double.NaN; out.sigPDedx = Double.NaN; out.dedxLnE = Double.NaN; out.dedxLnp = Double.NaN;
				out.dedxSigLnp = Double.NaN; out.dedxN = 0; out.prompt = false;
			}
		}
		return added;
	}

	/** hit quantities of the default fit d (compacted to the kept hits in s) and its prompt class. */
	private boolean hitsAndPrompt(GFitTrack t, GFitEvent ev, boolean bankTimes, SwitchConstants.DedxConstants dx, GFitWorkspace ws, Scratch s,
			GFitResult d) {
		final DedxHits h = s.hits;
		h.compute(d, t.hAdc, t.hPed, t.hm, dx.offset, dx.f3[ev.half], dx.f3Norm[ev.half]);
		int n = 0;
		for (int k = 0; k < d.nh; k++)
			if (h.kept[k]) {
				s.x[n] = h.X[k]; s.lnu[n] = h.lnu[k]; s.use[n] = h.use[k];
				n++;
			}
		s.nKept = n;
		if (d.status != 0 || !ws.atofOn) return false;
		AtofTerm.associate(d.u, prop, ws, ws.atof, s.rec, s.assoc);
		final int c = (int) s.assoc[3];
		if (!(s.assoc[2] > 0.5) || c < 0) return false;
		final double[] par = ws.atof.par;
		final double[] times = par[AtofTerm.P_TON] > 0.5 ? ev.aTime : (bankTimes ? ev.tBank : ev.tFinal);
		final double tres = times[ws.atof.crow[c]] - ev.startTime - (par[AtofTerm.P_TMODE] > 0.5 ? s.rec[11] : s.rec[10]) - par[AtofTerm.P_TOFF];
		return Math.abs(tres) < 1.0;      // NaN: false (finite(tres) of DR:38)
	}

	/** the profile around the sequence's final fit out; results into out. */
	private static void profile(GFitEvent ev, SwitchConstants.DedxConstants dx, Scratch s, GFitResult out, boolean prompt) {
		final double p = out.p();
		final double sig = p * Math.sqrt(Math.abs(out.cov[2 * GFitConfig.NPAR + 2]));      // est_sig_p of pack_est
		final DedxProfile f = s.prof;
		f.evaluate(s.x, s.lnu, s.use, s.nKept, p, sig, out.status, prompt, ev.half, dx);
		out.pDedx = f.pDedx; out.sigPDedx = f.sigPDedx; out.dedxLnE = f.lnE; out.dedxLnp = f.dedxLnp; out.dedxSigLnp = f.dedxSigLnp;
		out.dedxN = f.n; out.prompt = prompt;
		out.flags |= f.flags;
	}
}
