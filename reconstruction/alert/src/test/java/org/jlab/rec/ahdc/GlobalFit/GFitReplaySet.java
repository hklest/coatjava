package org.jlab.rec.ahdc.GlobalFit;

import java.io.File;
import java.util.Arrays;

import org.jlab.rec.ahdc.GlobalFit.GFitFixtures.Fx;

/**
 * Test-scope loader of a track fixture set (t500/t60/replay layout of export_fixtures.py): tracks.txt, atof.txt and
 * expected_before.txt, optionally pool.txt, expected_after.txt and expected_recover.txt, into one {@link GFitTrack}
 * and one {@link GFitEvent} per track plus the raw fixture arrays. Files may be .txt or .txt.gz.
 */
public final class GFitReplaySet {
	public int n;
	public GFitTrack[] tracks;
	public GFitEvent[] events;
	/** raw fixtures: after, recover and pool are null when not loaded. */
	public Fx tracksFx, atofFx, before, after, recover, poolFx;
	public int[] tier, tagEvent;

	private GFitReplaySet() {}

	/** all files of a test resource directory, e.g. {@code "t60"} (resource /gfit/t60/*.txt.gz). */
	public static GFitReplaySet resource(String prefix) {
		return build(name -> GFitFixtures.resource(prefix + "/" + name + ".txt.gz"), true);
	}

	/** tracks, atof and expected_before of a directory on disk. */
	public static GFitReplaySet file(String dir) { return file(dir, false); }

	/** with {@code extras}: also pool, expected_after and expected_recover (large for the replay sets). */
	public static GFitReplaySet file(String dir, boolean extras) {
		return build(name -> {
			File f = new File(dir, name + ".txt");
			return GFitFixtures.file(f.exists() ? f.getPath() : f.getPath() + ".gz");
		}, extras);
	}

	private interface Src { Fx get(String name); }

	private static GFitReplaySet build(Src src, boolean extras) {
		GFitReplaySet s = new GFitReplaySet();
		Fx t = s.tracksFx = src.get("tracks");
		Fx a = s.atofFx = src.get("atof");
		s.before = src.get("expected_before");
		if (extras) {
			s.poolFx = src.get("pool");
			s.after = src.get("expected_after");
			s.recover = src.get("expected_recover");
		}
		int n = s.n = t.shape("sel")[0];
		s.tier = t.i("tier");
		s.tagEvent = t.i("tag_event");
		int[] trackid = t.i("trackid"), hoff = t.i("hoff"), eoff = t.i("eoff"), hg = t.i("hg"), hism = t.i("hism");
		int[] hcls = t.i("hcls"), hev = t.i("h_ev"), eg = t.i("eg");
		double[] tc = t.d("TC"), start = t.d("start_time"), hm = t.d("hm"), hsig = t.d("hsig"), hh = t.d("hh");
		int[] aoff = a.i("aoff"), coff = a.i("coff");
		double[] ax = a.d("a_x"), ay = a.d("a_y"), az = a.d("a_z"), atb = a.d("a_t_bank"), atf = a.d("a_t_final");
		long[] abar = a.l("a_bar");
		int[] asec = a.i("a_sector"), alay = a.i("a_layer"), acomp = a.i("a_comp");
		int[] crow = a.i("c_row"), cbar = a.i("cbar"), csec = a.i("c_sec"), ccomp = a.i("c_comp"), fnone = a.i("fnone");
		double[] cz = a.d("cz"), cphi = a.d("cphi"), ct = a.d("ct"), ctime = a.d("c_time"), ctb = a.d("c_time_bank"), seed = a.d("atof_seed");
		int[] poff = null, pev = null, pg = null, pli = null, psrc = null, padc = null, pwf = null;
		double[] ptime = null, pdoca = null, plet = null, ptot = null, pped = null;
		if (s.poolFx != null) {
			Fx p = s.poolFx;
			poff = p.i("poff"); pev = p.i("p_ev"); pg = p.i("p_g"); pli = p.i("p_li"); psrc = p.i("p_src"); padc = p.i("p_adc"); pwf = p.i("p_wf");
			ptime = p.d("p_time"); pdoca = p.d("p_doca"); plet = p.d("p_let"); ptot = p.d("p_tot"); pped = p.d("p_ped");
		}
		s.tracks = new GFitTrack[n];
		s.events = new GFitEvent[n];
		for (int k = 0; k < n; k++) {
			int h0 = hoff[k], h1 = hoff[k + 1];
			s.tracks[k] = new GFitTrack(trackid[k], Arrays.copyOfRange(hg, h0, h1), Arrays.copyOfRange(hm, h0, h1),
					Arrays.copyOfRange(hsig, h0, h1), Arrays.copyOfRange(hh, h0, h1), Arrays.copyOfRange(hism, h0, h1),
					Arrays.copyOfRange(hcls, h0, h1), Arrays.copyOfRange(hev, h0, h1),
					tc[k * 7 + 1], tc[k * 7 + 2], tc[k * 7 + 3], tc[k * 7 + 4]);
			GFitEvent e = new GFitEvent();
			e.ze = tc[k * 7];
			e.half = (int) tc[k * 7 + 5];
			e.startTime = start[k];
			int a0 = aoff[k], a1 = aoff[k + 1];
			e.ax = Arrays.copyOfRange(ax, a0, a1); e.ay = Arrays.copyOfRange(ay, a0, a1); e.az = Arrays.copyOfRange(az, a0, a1);
			e.tBank = Arrays.copyOfRange(atb, a0, a1); e.tFinal = Arrays.copyOfRange(atf, a0, a1);
			e.isBar = new boolean[a1 - a0];
			for (int j = a0; j < a1; j++) e.isBar[j - a0] = abar[j] != 0;
			e.sector = Arrays.copyOfRange(asec, a0, a1); e.layer = Arrays.copyOfRange(alay, a0, a1);
			e.comp = Arrays.copyOfRange(acomp, a0, a1);
			int c0 = coff[k], c1 = coff[k + 1];
			e.cRow = Arrays.copyOfRange(crow, c0, c1); e.cBar = Arrays.copyOfRange(cbar, c0, c1);
			e.cSec = Arrays.copyOfRange(csec, c0, c1); e.cComp = Arrays.copyOfRange(ccomp, c0, c1);
			e.cz = Arrays.copyOfRange(cz, c0, c1); e.cphi = Arrays.copyOfRange(cphi, c0, c1); e.ct = Arrays.copyOfRange(ct, c0, c1);
			e.cTime = Arrays.copyOfRange(ctime, c0, c1); e.cTimeBank = Arrays.copyOfRange(ctb, c0, c1);
			e.fnone = fnone[k];
			e.atofSeed = Arrays.copyOfRange(seed, k * 20, k * 20 + 20);
			e.eg = Arrays.copyOfRange(eg, eoff[k], eoff[k + 1]);
			if (poff != null) {
				int p0 = poff[k], p1 = poff[k + 1];
				e.pEv = Arrays.copyOfRange(pev, p0, p1); e.pG = Arrays.copyOfRange(pg, p0, p1);
				e.pLi = Arrays.copyOfRange(pli, p0, p1); e.pSrc = Arrays.copyOfRange(psrc, p0, p1);
				e.pAdc = Arrays.copyOfRange(padc, p0, p1); e.pWf = Arrays.copyOfRange(pwf, p0, p1);
				e.pTime = Arrays.copyOfRange(ptime, p0, p1); e.pDoca = Arrays.copyOfRange(pdoca, p0, p1);
				e.pLet = Arrays.copyOfRange(plet, p0, p1); e.pTot = Arrays.copyOfRange(ptot, p0, p1);
				e.pPed = Arrays.copyOfRange(pped, p0, p1);
			}
			s.events[k] = e;
		}
		return s;
	}
}
