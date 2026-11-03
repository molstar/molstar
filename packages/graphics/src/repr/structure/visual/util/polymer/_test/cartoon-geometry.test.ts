/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ChunkedArray, hashFnv32a } from '@molstar/core/data/util';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { addTube } from '@molstar/graphics/geo/geometry/mesh/builder/tube';
import { addSheet } from '@molstar/graphics/geo/geometry/mesh/builder/sheet';
import { addRibbon } from '@molstar/graphics/geo/geometry/mesh/builder/ribbon';
import type { CurveSegmentControls } from '../curve-segment.js';

// `curve-segment.ts` keeps `Vec3.orthogonalize` in a module constant, so the spy goes in before it is loaded: it counts
// the calls, two per segment (first and last normal) plus the in-loop fallbacks for a zero or parallel tangent
const orthogonalize = vi.spyOn(Vec3, 'orthogonalize');
const { createCurveSegmentState, interpolateCurveSegment, interpolateSizes } = await import('../curve-segment.js');

// Pins the cartoon geometry bit for bit: the curve interpolation and the tube, sheet and ribbon builders, over a
// jittered helix with all the profile, cap and arrow variants the polymer trace visual uses. A change to the hashes
// is a change to the rendered cartoon.

function rng(seed: number) {
  return () => (seed = (seed * 16807) % 2147483647) / 2147483647;
}

/** a jittered alpha-helix-like trace (CA positions) and per-residue carbonyl-like directions */
function createTrace(n: number, degenerate: boolean) {
  const rnd = rng(11);
  const points: Vec3[] = [];
  const dirs: Vec3[] = [];
  for (let i = 0; i < n; ++i) {
    const a = (i * 100 * Math.PI) / 180;
    points.push(Vec3.create(2.3 * Math.cos(a) + rnd() * 0.3, 2.3 * Math.sin(a) + rnd() * 0.3, 1.5 * i + rnd() * 0.3));
    dirs.push(Vec3.normalize(Vec3(), Vec3.create(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5)));
  }
  // a zero and a tangent-parallel direction exercise the `orthogonalize` fallbacks
  dirs[5] = Vec3.create(0, 0, 0);
  dirs[9] = Vec3.sub(Vec3(), points[10], points[9]);
  if (degenerate) {
    // coincident points give zero tangents and the in-loop `orthogonalize` fallbacks
    for (let i = 13; i < 18; ++i) points[i] = points[12];
    dirs[14] = Vec3.create(0, 0, 0);
  }
  return { points, dirs };
}

function controlsAt(trace: ReturnType<typeof createTrace>, i: number, first: boolean, last: boolean) {
  const { points, dirs } = trace;
  const c: CurveSegmentControls = {
    secStrucFirst: first,
    secStrucLast: last,
    p0: points[i - 2],
    p1: points[i - 1],
    p2: points[i],
    p3: points[i + 1],
    p4: points[i + 2],
    d12: dirs[i - 1],
    d23: dirs[i],
  };
  return c;
}

const bits = (a: ArrayLike<number>) => hashFnv32a(new Uint32Array(Float32Array.from(a).buffer));

type Kind =
  | { kind: 'tube'; profile: 'elliptical' | 'rounded'; radialSegments: number; aspect: number; roundCap: boolean }
  | { kind: 'sheet'; arrow: boolean }
  | { kind: 'ribbon'; arrow: boolean };

function build(kind: Kind, linearSegments: number, tension: number, shift: number, degenerate = false) {
  const trace = createTrace(24, degenerate);
  const state = createCurveSegmentState(linearSegments);
  const builder = MeshBuilder.createState(64, 64);
  let curveHash = 0;
  for (let i = 2; i < 22; ++i) {
    // secondary structure elements of 5 residues; caps on their ends
    const first = i % 5 === 2;
    const last = i % 5 === 1;
    builder.currentGroup = i;
    interpolateCurveSegment(state, controlsAt(trace, i, first, last), tension, shift);
    const w = 0.2 + 0.01 * i;
    const h = kind.kind === 'tube' ? w * kind.aspect : w * 5;
    interpolateSizes(state, w * 0.9, w, w * 1.1, h * 0.9, h, h * 1.1, shift);
    const { curvePoints, tangentVectors, normalVectors, binormalVectors, widthValues, heightValues } = state;
    curveHash = hashFnv32a([
      curveHash,
      bits(curvePoints),
      bits(tangentVectors),
      bits(normalVectors),
      bits(binormalVectors),
    ]);
    if (kind.kind === 'tube') {
      const capped = kind.roundCap ? first || last || i % 7 === 0 : first;
      addTube(
        builder,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        kind.radialSegments,
        widthValues,
        heightValues,
        capped,
        kind.roundCap ? last || i % 7 === 0 : last,
        kind.profile,
        kind.roundCap,
      );
    } else if (kind.kind === 'sheet') {
      const arrowHeight = kind.arrow && last ? h * 1.5 : 0;
      addSheet(
        builder,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        widthValues,
        heightValues,
        arrowHeight,
        first,
        last,
      );
    } else {
      const arrowHeight = kind.arrow && last ? h * 1.5 : 0;
      addRibbon(
        builder,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        widthValues,
        heightValues,
        arrowHeight,
      );
    }
  }
  const { vertices, normals, indices, groups } = builder;
  return {
    curve: curveHash,
    vertexCount: vertices.elementCount,
    triangleCount: indices.elementCount,
    vertices: bits(ChunkedArray.compact(vertices, true)),
    normals: bits(ChunkedArray.compact(normals, true)),
    indices: hashFnv32a(ChunkedArray.compact(indices, true)),
    groups: bits(ChunkedArray.compact(groups, true)),
  };
}

const cases: [string, Kind, number, number, number][] = [
  [
    'helix elliptical',
    { kind: 'tube', profile: 'elliptical', radialSegments: 16, aspect: 5, roundCap: false },
    8,
    0.9,
    0.5,
  ],
  ['coil tube', { kind: 'tube', profile: 'elliptical', radialSegments: 16, aspect: 1, roundCap: false }, 8, 0.5, 0.5],
  ['helix rounded', { kind: 'tube', profile: 'rounded', radialSegments: 12, aspect: 5, roundCap: false }, 7, 0.9, 0.5],
  [
    'rounded, flat',
    { kind: 'tube', profile: 'rounded', radialSegments: 10, aspect: 0.5, roundCap: false },
    5,
    0.5,
    0.5,
  ],
  [
    'tubular helix, round caps',
    { kind: 'tube', profile: 'elliptical', radialSegments: 16, aspect: 1, roundCap: true },
    8,
    0.5,
    0.5,
  ],
  ['nucleic tube', { kind: 'tube', profile: 'elliptical', radialSegments: 6, aspect: 3, roundCap: false }, 6, 0.5, 0.3],
  ['sheet with arrows', { kind: 'sheet', arrow: true }, 8, 0.5, 0.5],
  ['sheet, square', { kind: 'sheet', arrow: false }, 3, 0.9, 0.3],
  ['ribbon with arrows', { kind: 'ribbon', arrow: true }, 4, 0.5, 0.5],
  [
    'one linear segment',
    { kind: 'tube', profile: 'elliptical', radialSegments: 4, aspect: 2, roundCap: false },
    1,
    0.5,
    0.5,
  ],
];

describe('cartoon geometry', () => {
  it('is pinned', () => {
    const out: Record<string, string> = {};
    for (const [name, kind, linearSegments, tension, shift] of cases) {
      const r = build(kind, linearSegments, tension, shift);
      out[name] =
        `${r.vertexCount}v ${r.triangleCount}t curve=${r.curve} v=${r.vertices} n=${r.normals} i=${r.indices} g=${r.groups}`;
    }
    // 20 segments: 40 calls without the fallback
    orthogonalize.mockClear();
    build(cases[0][1], 8, 0.9, 0.5);
    expect(orthogonalize).toHaveBeenCalledTimes(40);
    orthogonalize.mockClear();
    const d = build(cases[0][1], 8, 0.9, 0.5, true);
    expect(orthogonalize.mock.calls.length).toBeGreaterThan(40);
    out.degenerate = `${d.vertexCount}v ${d.triangleCount}t curve=${d.curve} v=${d.vertices} n=${d.normals} i=${d.indices} g=${d.groups}`;
    expect(out).toMatchInlineSnapshot(`
      {
        "coil tube": "3016v 5248t curve=2333589030 v=3862293639 n=1162244021 i=2528217029 g=3823400037",
        "degenerate": "3016v 5248t curve=3223941467 v=562689290 n=1333435979 i=2528217029 g=3823400037",
        "helix elliptical": "3016v 5248t curve=1144807048 v=863378944 n=3706819048 i=2528217029 g=3823400037",
        "helix rounded": "2024v 3456t curve=844614140 v=3736767465 n=3484826495 i=1912354837 g=2503719653",
        "nucleic tube": "896v 1488t curve=4158886629 v=514964450 n=4097688871 i=2968838797 g=2020398021",
        "one linear segment": "200v 192t curve=699304310 v=3443622594 n=3867938960 i=239117461 g=3150690405",
        "ribbon with arrows": "400v 320t curve=1420342767 v=562058045 n=545023109 i=1027626085 g=1653499653",
        "rounded, flat": "1288v 2080t curve=2069462369 v=3299690048 n=1572208531 i=526155157 g=2695361893",
        "sheet with arrows": "1488v 1304t curve=2333589030 v=1296133905 n=1943066435 i=910468645 g=3331605509",
        "sheet, square": "672v 496t curve=1618106690 v=3719697382 n=1909322117 i=3034569877 g=2020009541",
        "tubular helix, round caps": "3135v 5360t curve=2333589030 v=2844756929 n=2691703814 i=961893269 g=206096711",
      }
    `);
  });
});

describe('cartoon normals on tight turns of 4V6X (#649)', () => {
  const linearSegments = 8;
  const angle = (x: Vec3, y: Vec3) => (Math.acos(Math.min(1, Math.abs(Vec3.dot(x, y)))) * 180) / Math.PI;

  it.each([
    {
      name: 'RNA',
      shift: 0.3,
      p: [
        Vec3.create(93.775, -37.06, -14.431),
        Vec3.create(90.078, -32.656, -16.342),
        Vec3.create(86.01, -27.076, -14.678),
        Vec3.create(83.877, -32.335, -16.233),
        Vec3.create(82.599, -28.903, -10.817),
      ],
      d12: Vec3.create(-0.42075, 1.34875, 0.132),
      d23: Vec3.create(0.371, -1.19675, -0.15),
    },
    {
      name: 'protein',
      shift: 0.5,
      p: [
        Vec3.create(-97.194, 4.796, -58.202),
        Vec3.create(-100.441, 3.941, -59.373),
        Vec3.create(-100.373, 4.589, -63.03),
        Vec3.create(-101.99, 1.229, -62.481),
        Vec3.create(-105.258, 1.813, -64.168),
      ],
      d12: Vec3.create(-0.5275, -0.0985, -0.823),
      d23: Vec3.create(-0.39325, 0.025, 0.81325),
    },
  ])('$name: end normals on the direction lines, at most 90° of roll', ({ shift, p, d12, d23 }) => {
    const state = createCurveSegmentState(linearSegments);
    const controls: CurveSegmentControls = {
      secStrucFirst: false,
      secStrucLast: false,
      p0: p[0],
      p1: p[1],
      p2: p[2],
      p3: p[3],
      p4: p[4],
      d12,
      d23,
    };
    interpolateCurveSegment(state, controls, 0.5, shift);
    const { tangentVectors, normalVectors } = state;
    const t = (i: number) => Vec3.fromArray(Vec3(), tangentVectors, i * 3);
    const n = (i: number) => Vec3.fromArray(Vec3(), normalVectors, i * 3);

    expect(angle(n(0), Vec3.orthogonalize(Vec3(), t(0), d12))).toBeLessThan(0.1);
    expect(angle(n(linearSegments), Vec3.orthogonalize(Vec3(), t(linearSegments), d23))).toBeLessThan(0.1);
    for (let i = 0; i < linearSegments; ++i) {
      expect(angle(n(i), n(i + 1)) - angle(t(i), t(i + 1))).toBeLessThan(90 / linearSegments + 0.5);
    }
  });
});
