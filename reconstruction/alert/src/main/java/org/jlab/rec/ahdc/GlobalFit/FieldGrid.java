package org.jlab.rec.ahdc.GlobalFit;

/**
 * Solenoid map in (rho, z_map), bilinear with edge clamping: port of prop._bfield, map mode (pinned prop.py:95-117), on
 * the grid of prop.load_field (prop.py:66-87). The grid is passed already scaled (the Python grid multiplies by the RCDB
 * scale); z_map = z + zsh. A point outside the grid takes the edge cell with the fractions clamped to [0, 1]; nothing is
 * rejected per step. Arrays are flat row-major [i * nz + j] with i the rho node, j the z node. Immutable and thread safe.
 */
public final class FieldGrid {

	public final double drho, dz, zmin, zsh;
	public final int nr, nz;
	private final double[] br, bz;

	private FieldGrid(double drho, double dz, double zmin, double zsh, int nr, int nz, double[] br, double[] bz) {
		this.drho = drho;
		this.dz = dz;
		this.zmin = zmin;
		this.zsh = zsh;
		this.nr = nr;
		this.nz = nz;
		this.br = br;
		this.bz = bz;
	}

	/** br, bz flat [i*nz + j], already scaled (copied). */
	public static FieldGrid fromArrays(double drho, double dz, double zmin, double zsh, int nr, int nz, double[] br, double[] bz) {
		if (br.length != nr * nz || bz.length != nr * nz) throw new IllegalArgumentException("field array size");
		return new FieldGrid(drho, dz, zmin, zsh, nr, nz, br.clone(), bz.clone());
	}

	/** field (T) at (x, y, z) mm of the AHDC frame into b3 = (Bx, By, Bz). */
	public void field(double x, double y, double z, double[] b3) {
		double rho = Math.sqrt(x * x + y * y);
		double fr = rho / drho;
		int i = (int) fr;
		if (i > nr - 2) i = nr - 2;
		fr = Math.min(fr - i, 1.0);
		double fz = (z + zsh - zmin) / dz;
		int j = (int) Math.floor(fz);
		if (j < 0) j = 0;
		if (j > nz - 2) j = nz - 2;
		fz = Math.min(Math.max(fz - j, 0.0), 1.0);
		int a = i * nz + j, b = a + nz;
		double brv = (br[a] * (1 - fr) + br[b] * fr) * (1 - fz) + (br[a + 1] * (1 - fr) + br[b + 1] * fr) * fz;
		double bzv = (bz[a] * (1 - fr) + bz[b] * fr) * (1 - fz) + (bz[a + 1] * (1 - fr) + bz[b + 1] * fr) * fz;
		if (rho > 1e-9) {
			b3[0] = brv * x / rho;
			b3[1] = brv * y / rho;
			b3[2] = bzv;
		} else {
			b3[0] = 0.0;
			b3[1] = 0.0;
			b3[2] = bzv;
		}
	}
}
