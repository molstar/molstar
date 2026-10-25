/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { RuntimeContext } from '@molstar/core/task';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { getBoundary } from '../boundary.js';
import { calcMolecularSurface, type MolecularSurfaceCalculationProps } from '../molecular-surface.js';

/** Deterministic LCG in [0, 1). */
function rng(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

const VdW = [1.52, 1.55, 1.7, 1.8];

/**
 * A compact blob of `count` atoms (min distance 1.3 Å, inside a sphere of `extent` Å) with gaps in the
 * index set: every fourth slot holds a decoy atom outside `indices`, placed on top of a real one, so a
 * mix-up of ordinal and resolved indices changes the field.
 */
function blob(count: number, extent: number, seed: number, probeRadius: number) {
  const next = rng(seed);
  const pts: number[] = [];
  let tries = 0;
  while (pts.length / 3 < count && tries++ < 200000) {
    const x = (next() * 2 - 1) * extent,
      y = (next() * 2 - 1) * extent,
      z = (next() * 2 - 1) * extent;
    if (x * x + y * y + z * z > extent * extent) continue;
    let ok = true;
    for (let i = 0; i < pts.length; i += 3) {
      const dx = pts[i] - x,
        dy = pts[i + 1] - y,
        dz = pts[i + 2] - z;
      if (dx * dx + dy * dy + dz * dz < 1.3 * 1.3) {
        ok = false;
        break;
      }
    }
    if (ok) pts.push(x, y, z);
  }
  const n = pts.length / 3;
  const slots = n + Math.floor(n / 3) + 2;
  const x = new Float32Array(slots),
    y = new Float32Array(slots),
    z = new Float32Array(slots);
  const radius = new Float32Array(slots);
  const indices: number[] = [];
  let maxRadius = 0;
  for (let s = 1, a = 0; a < n; s++) {
    if (s % 4 === 0) {
      // decoy: big sphere on the previous real atom
      x[s] = x[s - 1];
      y[s] = y[s - 1];
      z[s] = z[s - 1];
      radius[s] = 5;
      continue;
    }
    x[s] = pts[3 * a];
    y[s] = pts[3 * a + 1];
    z[s] = pts[3 * a + 2];
    const r = VdW[Math.floor(next() * VdW.length)];
    if (r > maxRadius) maxRadius = r;
    radius[s] = r + probeRadius;
    indices.push(s);
    a++;
  }
  const id = new Int32Array(n);
  for (let i = 0; i < n; i++) id[i] = 1000 + 3 * i;
  const position = { indices: OrderedSet.ofSortedArray(indices), x, y, z, id, radius };
  return { position, boundary: getBoundary(position), maxRadius };
}

/**
 * Atoms that clash on purpose, over an `Interval` index set: random centres with no minimum distance,
 * coincident pairs (same and different radii), a small atom nested in a big one, and a pair close to
 * external tangency.
 */
function clashes(probeRadius: number) {
  const next = rng(42);
  const xs: number[] = [],
    ys: number[] = [],
    zs: number[] = [],
    rs: number[] = [];
  const add = (x: number, y: number, z: number, r: number) => {
    xs.push(x);
    ys.push(y);
    zs.push(z);
    rs.push(r);
  };
  for (let i = 0; i < 40; i++) {
    add((next() * 2 - 1) * 6, (next() * 2 - 1) * 6, (next() * 2 - 1) * 6, VdW[Math.floor(next() * VdW.length)]);
  }
  for (let i = 0; i < 5; i++) add(xs[i], ys[i], zs[i], i < 2 ? rs[i] : 1.2);
  add(xs[5] + 0.2, ys[5], zs[5], 0.5);
  add(8, 0, 0, 1.5);
  add(8 + 2 * (1.5 + probeRadius), 0, 0, 1.5);
  const n = xs.length;
  const radius = new Float32Array(n);
  let maxRadius = 0;
  for (let i = 0; i < n; i++) {
    if (rs[i] > maxRadius) maxRadius = rs[i];
    radius[i] = rs[i] + probeRadius;
  }
  const id = new Int32Array(n);
  for (let i = 0; i < n; i++) id[i] = Math.floor(next() * 100000);
  const position = {
    indices: OrderedSet.ofBounds(0, n),
    x: Float32Array.from(xs),
    y: Float32Array.from(ys),
    z: Float32Array.from(zs),
    id,
    radius,
  };
  return { position, boundary: getBoundary(position), maxRadius };
}

/** FNV-1a over the bytes (not the values) of a float array, so any bit change shows. */
function hash(a: Float32Array) {
  return hashFnv32a(new Uint8Array(a.buffer, a.byteOffset, a.byteLength))
    .toString(16)
    .padStart(8, '0');
}

type Case = {
  name: string;
  /** atoms of a blob, or the `clashes` input */
  atoms: number | 'clashes';
  props: MolecularSurfaceCalculationProps;
  /** overrides the max VdW radius (also sets the pad, the lookup cell and the update chunk) */
  maxRadius?: number;
  nullBox?: boolean;
  expected: { dim: number[]; field: string; idField: string; inside: number; ids: number };
};

// Pinned on the implementation before the occlusion rewrite (compact neighbour lists, per-torus pruning):
// the field and the id field must stay bit-identical. `ngPoints > 0` (probeRadius < 2 * resolution) is
// the extended mode.
const Cases: Case[] = [
  {
    name: 'default, res 0.5',
    atoms: 220,
    props: { probeRadius: 1.4, resolution: 0.5, probePositions: 36 },
    expected: { dim: [56, 55, 55], field: '51088ec4', idField: 'b85beb94', inside: 44011, ids: 125 },
  },
  {
    name: 'extended, res 1.3',
    atoms: 220,
    props: { probeRadius: 1.4, resolution: 1.3, probePositions: 24 },
    expected: { dim: [23, 23, 23], field: 'cc5e5ddc', idField: '7c1cd090', inside: 2501, ids: 123 },
  },
  {
    name: 'extended, res 0.7, probe 1.2, null box',
    atoms: 220,
    props: { probeRadius: 1.2, resolution: 0.7, probePositions: 30 },
    nullBox: true,
    expected: { dim: [40, 40, 40], field: 'ab3a38ea', idField: 'b8e58b8e', inside: 15005, ids: 127 },
  },
  {
    name: 'probe 0, res 0.4',
    atoms: 120,
    props: { probeRadius: 0, resolution: 0.4, probePositions: 12 },
    expected: { dim: [63, 63, 61], field: 'f2bc7a6c', idField: 'd4564a5d', inside: 28529, ids: 120 },
  },
  {
    name: 'large probe, res 0.6, max radius 3',
    atoms: 160,
    props: { probeRadius: 3, resolution: 0.6, probePositions: 48 },
    maxRadius: 3,
    expected: { dim: [56, 56, 55], field: '634c0741', idField: 'bbc20a6e', inside: 38308, ids: 83 },
  },
  {
    name: 'extended, res 1.3, max radius 3.2',
    atoms: 220,
    props: { probeRadius: 1.4, resolution: 1.3, probePositions: 24 },
    maxRadius: 3.2,
    expected: { dim: [25, 25, 25], field: 'ac830a06', idField: '0943936a', inside: 2499, ids: 123 },
  },
  {
    name: 'clashes, res 0.5',
    atoms: 'clashes',
    props: { probeRadius: 1.4, resolution: 0.5, probePositions: 36 },
    expected: { dim: [61, 45, 46], field: '61cb0056', idField: '40bc68da', inside: 21158, ids: 44 },
  },
  {
    name: 'clashes, extended, res 1, probe 0.6',
    atoms: 'clashes',
    props: { probeRadius: 0.6, resolution: 1, probePositions: 20 },
    expected: { dim: [29, 22, 22], field: 'e3fba7fa', idField: '5a27b030', inside: 1398, ids: 44 },
  },
  {
    name: 'two atoms',
    atoms: 2,
    props: { probeRadius: 1.4, resolution: 0.5, probePositions: 36 },
    expected: { dim: [23, 23, 23], field: '013eae2c', idField: '20fb81ab', inside: 1251, ids: 2 },
  },
  {
    name: 'one atom, extended',
    atoms: 1,
    props: { probeRadius: 1.4, resolution: 1, probePositions: 36 },
    expected: { dim: [11, 11, 11], field: '504fed0c', idField: '302f6478', inside: 100, ids: 1 },
  },
];

function input(c: Case) {
  if (c.atoms === 'clashes') return clashes(c.props.probeRadius);
  const extent = c.atoms <= 2 ? 1.2 : 9;
  const b = blob(c.atoms, extent, 7 + c.atoms, c.props.probeRadius);
  return { ...b, maxRadius: c.maxRadius ?? b.maxRadius };
}

async function compute(ctx: RuntimeContext, c: Case) {
  const { position, boundary, maxRadius } = input(c);
  const { field, idField } = await calcMolecularSurface(
    ctx,
    position,
    boundary,
    maxRadius,
    c.nullBox ? null : boundary.box,
    c.props,
  );
  const data = field.data as Float32Array,
    ids = idField.data as Float32Array;
  let inside = 0;
  for (let i = 0; i < data.length; i++) if (data[i] > 0) inside++;
  const seen = new Set<number>();
  for (let i = 0; i < ids.length; i++) if (ids[i] >= 0) seen.add(ids[i]);
  return {
    dim: Array.from(field.space.dimensions),
    field: hash(data),
    idField: hash(ids),
    inside,
    ids: seen.size,
  };
}

/** Yields on every check: concurrent calculations interleave between chunks and passes. */
function yieldingContext(): RuntimeContext {
  return {
    shouldUpdate: true,
    isSynchronous: false,
    update: () => new Promise<void>((resolve) => setTimeout(resolve, 0)),
  };
}

describe('calcMolecularSurface', () => {
  for (const c of Cases) {
    it(c.name, async () => {
      expect(await compute(RuntimeContext.Synchronous, c)).toEqual(c.expected);
    });
  }

  it('projects the torii without the Vec3 helpers (perf-audit-2 02/P11)', async () => {
    const normalize = vi.spyOn(Vec3, 'normalize');
    const cross = vi.spyOn(Vec3, 'cross');
    try {
      expect(await compute(RuntimeContext.Synchronous, Cases[0])).toEqual(Cases[0].expected);
      expect(normalize).not.toHaveBeenCalled();
      expect(cross).not.toHaveBeenCalled();
    } finally {
      normalize.mockRestore();
      cross.mockRestore();
    }
  });

  it('concurrent calculations with yields give the same fields', async () => {
    const results = await Promise.all(Cases.map((c) => compute(yieldingContext(), c)));
    expect(results).toEqual(Cases.map((c) => c.expected));
  });
});

/** `n`³ atoms on a lattice with `spacing` Å, ids = ordinals; with a spacing above 2 · (r + probe) no torus is made */
function lattice(n: number, spacing: number, r: number, probeRadius: number) {
  const count = n * n * n;
  const x = new Float32Array(count),
    y = new Float32Array(count),
    z = new Float32Array(count);
  for (let i = 0; i < count; i++) {
    x[i] = (i % n) * spacing;
    y[i] = (Math.floor(i / n) % n) * spacing;
    z[i] = Math.floor(i / (n * n)) * spacing;
  }
  const radius = new Float32Array(count).fill(r + probeRadius);
  const id = new Int32Array(count);
  for (let i = 0; i < count; i++) id[i] = i;
  const position = { indices: OrderedSet.ofBounds(0, count), x, y, z, id, radius };
  return { position, boundary: getBoundary(position), count };
}

/** Counts the progress updates by message; asks for one at every check */
function countingContext(counts: Map<string, number>): RuntimeContext {
  return {
    shouldUpdate: true,
    isSynchronous: false,
    update: async (p) => {
      const message = typeof p === 'string' ? p : (p?.message ?? '');
      counts.set(message, (counts.get(message) ?? 0) + 1);
    },
  };
}

describe('calcMolecularSurface progress', () => {
  it('yields within a pass for small radii (perf-audit-2 02/B2)', async () => {
    // 15,625 atoms of radius 0.3: 797 atoms of work per chunk at resolution 1 (1e5 / r⁹ gave billions)
    const { position, boundary } = lattice(25, 4, 0.3, 1.4);
    const props = { probeRadius: 1.4, resolution: 1, probePositions: 12 };
    const counts = new Map<string, number>();
    await calcMolecularSurface(countingContext(counts), position, boundary, 0.3, boundary.box, props);
    expect(counts.get('projecting points')).toBe(20);
  });

  it('counts the per-atom fixed cost at coarse resolution (perf-audit-3 02/B6)', async () => {
    // 8,000 atoms of radius 1.7 with tori at resolution 8: 903 atoms per chunk in both passes (voxels alone: 93.5k,
    // one synchronous chunk of seconds on a 186k-atom structure surface)
    const { position, boundary } = lattice(20, 4, 1.7, 1.4);
    const props = { probeRadius: 1.4, resolution: 8, probePositions: 12 };
    const counts = new Map<string, number>();
    await calcMolecularSurface(countingContext(counts), position, boundary, 1.7, boundary.box, props);
    expect(counts.get('projecting points')).toBe(9);
    expect(counts.get('projecting torii')).toBe(9);
  });
});
