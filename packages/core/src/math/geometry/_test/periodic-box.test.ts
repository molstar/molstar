/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '../../linear-algebra.js';
import { PeriodicBox } from '../periodic-box.js';
import { Cell } from '../spacegroup/cell.js';

/** The shortest of `d` moved by every combination of −n..n box vectors. */
function bruteForce(box: PeriodicBox, d: Vec3, n = 1) {
  const [a, b, c] = [0, 3, 6].map((o) => Vec3.create(box.basis[o], box.basis[o + 1], box.basis[o + 2]));
  let best = Vec3.clone(d);
  for (let i = -n; i <= n; i++) {
    for (let j = -n; j <= n; j++) {
      for (let k = -n; k <= n; k++) {
        const v = Vec3.clone(d);
        Vec3.scaleAndAdd(v, v, a, i);
        Vec3.scaleAndAdd(v, v, b, j);
        Vec3.scaleAndAdd(v, v, c, k);
        if (Vec3.magnitude(v) < Vec3.magnitude(best)) best = v;
      }
    }
  }
  return best;
}

// A rectangular box, and GROMACS's rhombic dodecahedron (xy-square, the MDAnalysis AdK box)
const Rectangular = PeriodicBox.create(Vec3.create(30, 0, 0), Vec3.create(0, 40, 0), Vec3.create(0, 0, 50));
const Dodecahedron = PeriodicBox.create(Vec3.create(80, 0, 0), Vec3.create(0, 80, 0), Vec3.create(40, 40, 56.57));

/** A deterministic sequence in [-1, 1) */
function* noise(n: number) {
  let s = 12345;
  for (let i = 0; i < n; i++) {
    s = (s * 1103515245 + 12345) % 2147483648;
    yield (s / 2147483648) * 2 - 1;
  }
}

describe('PeriodicBox', () => {
  it('is the nearest image in a rectangular box', () => {
    const r = [...noise(300)];
    for (let i = 0; i < r.length; i += 3) {
      const d = Vec3.create(r[i] * 45, r[i + 1] * 60, r[i + 2] * 75);
      const m = PeriodicBox.minimumImage(Vec3(), Rectangular, d);
      expect(Vec3.distance(m, bruteForce(Rectangular, d))).toBeLessThan(1e-9);
    }
  });

  it('is the nearest image in a dodecahedron for displacements under half a box', () => {
    const r = [...noise(300)];
    for (let i = 0; i < r.length; i += 3) {
      // a short displacement (a bond) moved by a few box vectors
      const bond = Vec3.create(r[i] * 2, r[i + 1] * 2, r[i + 2] * 2);
      const d = Vec3.clone(bond);
      Vec3.add(d, d, Vec3.create(80 * Math.round(r[i] * 2), 0, 0));
      Vec3.add(
        d,
        d,
        Vec3.create(40 * Math.round(r[i + 2] * 2), 40 * Math.round(r[i + 2] * 2), 56.57 * Math.round(r[i + 2] * 2)),
      );
      const m = PeriodicBox.minimumImage(Vec3(), Dodecahedron, d);
      expect(Vec3.distance(m, bond)).toBeLessThan(1e-9);
    }
  });

  it('is the nearest image in a dodecahedron and a truncated octahedron for any displacement (perf-audit-3 02/B7)', () => {
    const L = 10;
    const dodecahedron = PeriodicBox.create(
      Vec3.create(L, 0, 0),
      Vec3.create(0, L, 0),
      Vec3.create(L / 2, L / 2, (L * Math.SQRT2) / 2),
    );
    // GROMACS -bt octahedron
    const octahedron = PeriodicBox.create(
      Vec3.create(L, 0, 0),
      Vec3.create(L / 3, (2 * Math.SQRT2 * L) / 3, 0),
      Vec3.create(-L / 3, (Math.SQRT2 * L) / 3, (Math.sqrt(6) * L) / 3),
    );
    // rounding alone gives (−5, −5, −2.57), |m| 7.52
    expect(Vec3.magnitude(PeriodicBox.minimumImage(Vec3(), dodecahedron, Vec3.create(0, 0, 4.5)))).toBeCloseTo(4.5, 9);

    const r = [...noise(3000)];
    for (const box of [dodecahedron, octahedron]) {
      for (let i = 0; i < r.length; i += 3) {
        // anywhere within a box length on each axis, against every shift by up to 3 box vectors
        const d = Vec3.create(r[i] * L, r[i + 1] * L, r[i + 2] * L);
        const m = PeriodicBox.minimumImage(Vec3(), box, d);
        expect(Vec3.magnitude(m)).toBeCloseTo(Vec3.magnitude(bruteForce(box, d, 3)), 9);
        // and it is an image of d
        const f = Vec3.transformMat3(Vec3(), Vec3.sub(Vec3(), m, d), box.inverse);
        for (let k = 0; k < 3; k++) expect(Math.abs(f[k] - Math.round(f[k]))).toBeLessThan(1e-9);
      }
    }
  });

  it('is the nearest image in thin, tilted boxes (perf-audit-3 02/B7)', () => {
    // GROMACS-valid (|bx|, |cx| ≤ ax/2, |cy| ≤ by/2) but flat: the shortest image is more than one shift from the
    // rounded one (27 shifts gave 5.38 here)
    const slab = PeriodicBox.create(Vec3.create(10, 0, 0), Vec3.create(3.53, 10, 0), Vec3.create(-4.11, 2.12, 3));
    const m = PeriodicBox.minimumImage(Vec3(), slab, Vec3.create(15.3, 4.42, -13.5));
    expect(Vec3.magnitude(m)).toBeCloseTo(Vec3.magnitude(bruteForce(slab, Vec3.create(15.3, 4.42, -13.5), 8)), 9);

    const r = [...noise(1200)];
    for (let b = 0; b < 4; b++) {
      const tilt = r.slice(b * 4, b * 4 + 4);
      const box = PeriodicBox.create(
        Vec3.create(10, 0, 0),
        Vec3.create(tilt[0] * 5, 10, 0),
        Vec3.create(tilt[1] * 5, tilt[2] * 5, 1.5 + tilt[3]),
      );
      for (let i = 300 * b + 16; i + 2 < 300 * (b + 1); i += 6) {
        // within 5 Å of the origin: every image as short is within 10 shifts
        const d = Vec3.create(r[i] * 5, r[i + 1] * 5, r[i + 2] * 5);
        const m = PeriodicBox.minimumImage(Vec3(), box, d);
        expect(Vec3.magnitude(m)).toBeCloseTo(Vec3.magnitude(bruteForce(box, d, 10)), 9);
      }
    }
  });

  it('takes a cell in its standard frame', () => {
    const cell = Cell.create(Vec3.create(80, 80, 80), Vec3.create(Math.PI / 3, Math.PI / 3, Math.PI / 2));
    const box = PeriodicBox.fromCell(cell);
    // a along x, b along y (gamma 90°), c at 60° to both
    expect(box.basis[0]).toBeCloseTo(80);
    expect(box.basis[4]).toBeCloseTo(80);
    expect(box.basis[6]).toBeCloseTo(40);
    expect(box.basis[7]).toBeCloseTo(40);
    expect(PeriodicBox.isValid(box)).toBe(true);
    expect(PeriodicBox.isValid(PeriodicBox.fromCell(Cell.empty()))).toBe(false);
  });
});
