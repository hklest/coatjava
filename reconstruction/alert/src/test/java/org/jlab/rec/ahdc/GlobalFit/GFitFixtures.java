package org.jlab.rec.ahdc.GlobalFit;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Reader of the fixture text format written by the pinned Python exporter (tools/gr_java/export_fixtures.py).
 * Line 1: {@code # name:shape:dtype ...}; then all numbers of all arrays in header order, row-major, whitespace
 * separated. dtype f8 (Python repr, with nan, inf, -inf), i8, b1 (0/1). Files ending in .gz are decompressed.
 */
public final class GFitFixtures {

	/** One fixture file: arrays by name, in header order. */
	public static final class Fx {
		final Map<String, int[]> shapes = new LinkedHashMap<>();
		final Map<String, double[]> f8 = new LinkedHashMap<>();
		final Map<String, long[]> i8 = new LinkedHashMap<>();

		public boolean has(String name) { return shapes.containsKey(name); }

		public int[] shape(String name) { return need(shapes, name); }

		public Iterable<String> names() { return shapes.keySet(); }

		public boolean isFloat(String name) { need(shapes, name); return f8.containsKey(name); }

		/** f8 array; an i8 or b1 array is converted. */
		public double[] d(String name) {
			double[] a = f8.get(name);
			if (a != null) return a;
			long[] l = need(i8, name);
			double[] o = new double[l.length];
			for (int k = 0; k < l.length; k++) o[k] = l[k];
			return o;
		}

		public long[] l(String name) { return need(i8, name); }

		public int[] i(String name) {
			long[] l = need(i8, name);
			int[] o = new int[l.length];
			for (int k = 0; k < l.length; k++) o[k] = Math.toIntExact(l[k]);
			return o;
		}

		public double d0(String name) { return d(name)[0]; }

		public int i0(String name) { return Math.toIntExact(l(name)[0]); }

		/** row {@code r} of a 2-d (or higher, flattened after the first axis) array. */
		public double[] row(String name, int r) {
			int[] s = shape(name);
			int w = 1;
			for (int k = 1; k < s.length; k++) w *= s[k];
			double[] a = d(name);
			double[] o = new double[w];
			System.arraycopy(a, r * w, o, 0, w);
			return o;
		}

		private static <T> T need(Map<String, T> m, String name) {
			T v = m.get(name);
			if (v == null) throw new IllegalArgumentException("fixture has no array " + name);
			return v;
		}
	}

	private GFitFixtures() {}

	/** test resource {@code /gfit/<name>} (for example {@code unit/prop_states.txt.gz}). */
	public static Fx resource(String name) {
		String path = "/gfit/" + name;
		try (InputStream in = GFitFixtures.class.getResourceAsStream(path)) {
			if (in == null) throw new IllegalArgumentException("no test resource " + path);
			return read(name.endsWith(".gz") ? new GZIPInputStream(in, 1 << 16) : in);
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}

	/** fixture file on disk (harness). */
	public static Fx file(String path) {
		try (InputStream in = new FileInputStream(path)) {
			return read(path.endsWith(".gz") ? new GZIPInputStream(in, 1 << 16) : in);
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}

	public static Fx read(InputStream raw) throws IOException {
		Tok t = new Tok(new BufferedInputStream(raw, 1 << 16));
		String header = t.line();
		if (header == null || !header.startsWith("#")) throw new IOException("fixture header missing");
		Fx fx = new Fx();
		for (String e : header.substring(1).trim().split("\\s+")) {
			String[] p = e.split(":");
			if (p.length != 3) throw new IOException("bad header entry " + e);
			String[] dims = p[1].split("x");
			int[] shape = new int[dims.length];
			int n = 1;
			for (int k = 0; k < dims.length; k++) {
				shape[k] = Integer.parseInt(dims[k]);
				n *= shape[k];
			}
			fx.shapes.put(p[0], shape);
			if (p[2].equals("f8")) {
				double[] a = new double[n];
				for (int k = 0; k < n; k++) a[k] = parseF8(t.need());
				fx.f8.put(p[0], a);
			} else if (p[2].equals("i8") || p[2].equals("b1")) {
				long[] a = new long[n];
				for (int k = 0; k < n; k++) a[k] = Long.parseLong(t.need());
				fx.i8.put(p[0], a);
			} else {
				throw new IOException("unknown dtype " + p[2]);
			}
		}
		if (t.next() != null) throw new IOException("fixture has trailing tokens");
		return fx;
	}

	static double parseF8(String s) {
		switch (s) {
			case "nan": return Double.NaN;
			case "inf": return Double.POSITIVE_INFINITY;
			case "-inf": return Double.NEGATIVE_INFINITY;
			default: return Double.parseDouble(s);
		}
	}

	/** whitespace tokenizer over bytes (ASCII). */
	private static final class Tok {
		private final InputStream in;
		private final StringBuilder sb = new StringBuilder(32);

		Tok(InputStream in) { this.in = in; }

		String line() throws IOException {
			sb.setLength(0);
			int c;
			while ((c = in.read()) >= 0 && c != '\n') sb.append((char) c);
			return (c < 0 && sb.length() == 0) ? null : sb.toString();
		}

		String need() throws IOException {
			String s = next();
			if (s == null) throw new IOException("fixture ends early");
			return s;
		}

		String next() throws IOException {
			int c;
			do {
				c = in.read();
			} while (c == ' ' || c == '\n' || c == '\r' || c == '\t');
			if (c < 0) return null;
			sb.setLength(0);
			while (c >= 0 && c != ' ' && c != '\n' && c != '\r' && c != '\t') {
				sb.append((char) c);
				c = in.read();
			}
			return sb.toString();
		}
	}
}
