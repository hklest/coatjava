package org.jlab.rec.ahdc.GlobalFit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.jlab.rec.ahdc.GlobalFit.GFitFixtures.Fx;
import org.junit.Test;

/** Checks of the compressed test resources (S1d): stream hashes, array counts, size bound, config and replay-set loader. */
public class GFitFixturesTest {

	private static final String[] T60 = {"tracks", "atof", "pool", "expected_before", "expected_after", "expected_recover", "MANIFEST"};

	private static String sha256(InputStream in) throws Exception {
		MessageDigest md = MessageDigest.getInstance("SHA-256");
		byte[] buf = new byte[1 << 16];
		for (int n; (n = in.read(buf)) > 0; ) md.update(buf, 0, n);
		StringBuilder sb = new StringBuilder();
		for (byte b : md.digest()) sb.append(String.format("%02x", b));
		return sb.toString();
	}

	private static InputStream open(String res) throws Exception {
		InputStream in = GFitFixturesTest.class.getResourceAsStream("/gfit/" + res);
		assertTrue("resource " + res, in != null);
		return res.endsWith(".gz") ? new GZIPInputStream(in) : in;
	}

	@Test
	public void fixtureRoundTrip() throws Exception {
		// SHA256SUMS: sha256 of the decompressed text of every unit/ and t60/ resource, written by subset_fixtures.py
		List<String[]> sums = new ArrayList<>();
		try (BufferedReader r = new BufferedReader(new InputStreamReader(open("SHA256SUMS"), StandardCharsets.UTF_8))) {
			for (String l; (l = r.readLine()) != null; ) if (!l.isEmpty()) sums.add(l.split("\\s+"));
		}
		assertTrue(sums.size() >= 14);
		int nT60 = 0, nUnit = 0;
		for (String[] s : sums) {
			try (InputStream in = open(s[1] + ".gz")) {
				assertEquals(s[1], s[0], sha256(in));
			}
			if (s[1].startsWith("t60/")) nT60++;
			if (s[1].startsWith("unit/")) nUnit++;
			if (s[1].endsWith("MANIFEST.txt")) continue;
			// array counts: header shapes against the parsed arrays (the reader fails on missing or extra tokens)
			Fx fx = GFitFixtures.resource(s[1] + ".gz");
			int na = 0;
			for (String name : fx.names()) {
				int n = 1;
				for (int d : fx.shape(name)) n *= d;
				assertEquals(s[1] + ":" + name, n, fx.isFloat(name) ? fx.d(name).length : fx.l(name).length);
				na++;
			}
			assertTrue(s[1], na > 0);
		}
		assertEquals(7, nT60);
		assertEquals(7, nUnit);
		// total size of the resource directory (test classpath)
		File dir;
		try {
			dir = new File(GFitFixturesTest.class.getResource("/gfit/SHA256SUMS").toURI()).getParentFile();
		} catch (URISyntaxException e) {
			throw new AssertionError(e);
		}
		long total = sizeOf(dir);
		assertTrue("resources " + total + " bytes", total <= 1_500_000L);
	}

	private static long sizeOf(File f) {
		if (f.isFile()) return f.length();
		long t = 0;
		File[] l = f.listFiles();
		if (l != null) for (File c : l) t += sizeOf(c);
		return t;
	}

	@Test
	public void configMatchesPython() {
		Fx c = GFitFixtures.resource("unit/config.txt.gz");
		GFitConfig g = GFitConfig.pythonDefault();
		assertArrayEquals(c.d("cfg_cf"), g.cf(), 0.0);
		assertArrayEquals(c.i("cfg_ci"), g.ci());
		double[] su = c.d("cfg_steps_u");
		for (int k = 0; k < 5; k++) assertEquals(su[k], GFitConfig.stepU(k), 0.0);
		assertEquals(su[4], GFitConfig.STEP_KINK, 0.0);
		double[] fms = c.d("cfg_f_ms");
		for (int k = 0; k < 6; k++) assertEquals(fms[k], g.fMs(k), 0.0);
		double[] gr = c.d("cfg_grid");
		assertEquals(gr[0], GFitConfig.NG, 0.0);
		assertEquals(gr[1], GFitConfig.SMAX, 0.0);
		assertEquals(gr[2], GFitConfig.DS, 0.0);
		assertEquals(gr[3], GFitConfig.NPL, 0.0);
		assertEquals(gr[4], GFitConfig.NK, 0.0);
		assertEquals(gr[5], GFitConfig.NPAR, 0.0);
		assertEquals(gr[6], GFitConfig.NADD, 0.0);
		assertEquals(gr[7], GFitConfig.NMISS, 0.0);
		assertEquals(gr[8], GFitConfig.SOL, 0.0);
		double[] rw = c.d("rec_windows");
		assertEquals(rw[0], g.recTLo(), 0.0);
		assertEquals(rw[1], g.recTHi(), 0.0);
		assertEquals(rw[2], g.recTotMin(), 0.0);
		assertEquals(rw[3], g.recAdcMin(), 0.0);
		assertEquals(rw[4], g.recWfMax(), 0.0);
		assertEquals(rw[5], g.recKRes(), 0.0);
		assertEquals(rw[6], g.recQ(), 0.0);
		assertEquals(rw[7], g.recRounds(), 0.0);
		assertEquals(c.i0("rec_pool_w"), g.recPoolW());
		double[] ts = c.d("rec_t_stock");
		assertEquals(ts[0], g.recTStockLo(), 0.0);
		assertEquals(ts[1], g.recTStockHi(), 0.0);
		// copies
		assertEquals(1, g.withNFlip(1).nFlip());
		assertEquals(2, g.nFlip());
		assertEquals(2, g.withStarts(GFitConfig.Starts.B0).ci()[12]);
		assertEquals(3, g.ci()[12]);
		assertEquals(GFitConfig.Search.FAST, g.withSearch(GFitConfig.Search.FAST).search());
		assertEquals(GFitConfig.Search.FULL, g.search());
	}

	@Test
	public void replaySetT60() {
		GFitReplaySet s = GFitReplaySet.resource("t60");
		assertEquals(60, s.n);
		assertEquals(60, s.tracks.length);
		assertEquals(60, s.events.length);
		int[] hoff = s.tracksFx.i("hoff"), eoff = s.tracksFx.i("eoff"), aoff = s.atofFx.i("aoff"), poff = s.poolFx.i("poff");
		int[] cnt = new int[4];
		for (int k = 0; k < 60; k++) {
			assertEquals(hoff[k + 1] - hoff[k], s.tracks[k].nHits());
			assertEquals(eoff[k + 1] - eoff[k], s.events[k].eg.length);
			assertEquals(aoff[k + 1] - aoff[k], s.events[k].nAtof());
			assertEquals(poff[k + 1] - poff[k], s.events[k].nPool());
			assertTrue(s.events[k].half == 0 || s.events[k].half == 1);
			assertEquals(s.tagEvent[k] % 2, s.events[k].half);
			cnt[s.tier[k]]++;
		}
		assertEquals(2, cnt[1]);
		assertEquals(40, cnt[2]);
		assertEquals(18, cnt[3]);
		// before-fit hit lists have the same track count, and the recovery table indexes tracks 0..59 only
		assertEquals(60, s.before.shape("est_status")[0]);
		long[] pt = s.recover.l("pt_t");
		for (long t : pt) assertTrue(t >= 0 && t < 60);
		long[] acc = s.after.l("acc_t");
		for (long t : acc) assertTrue(t >= 0 && t < 60);
	}

	@Test
	public void resultAccessors() {
		GFitResult r = new GFitResult(10);
		r.u[1] = Math.toRadians(190.0);
		r.u[2] = Math.log(0.25);
		r.u[3] = 3.0; r.u[4] = Math.toRadians(80.0);
		r.cov[2 * 15 + 2] = 4.0;
		r.cov[0] = 1.0;
		assertEquals(0.25, r.p(), 1e-15);
		assertEquals(-170.0, r.phiDeg(), 1e-12);
		assertEquals(80.0, r.thetaDeg(), 1e-12);
		assertEquals(3.0, r.zv(), 0.0);
		double[] c = new double[25];
		r.covX(c);
		assertEquals(1.0, c[0], 0.0);
		assertEquals(4.0 / (0.25 * 0.25), c[24], 1e-12);
		r.clear();
		assertEquals(0.0, r.u[2], 0.0);
		assertEquals(18, r.hitD.length);
	}
}
