/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ChunkedArray, hashFnv32a } from '@molstar/core/data/util';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { addTube } from '@molstar/graphics/geo/geometry/mesh/builder/tube';
import { addSheet } from '@molstar/graphics/geo/geometry/mesh/builder/sheet';
import { addRibbon } from '@molstar/graphics/geo/geometry/mesh/builder/ribbon';
import {
  type CurveSegmentControls,
  createCurveSegmentState,
  interpolateCurveSegment,
  interpolateSizes,
} from '../curve-segment.js';

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
    const d = build(cases[0][1], 8, 0.9, 0.5, true);
    out.degenerate = `${d.vertexCount}v ${d.triangleCount}t curve=${d.curve} v=${d.vertices} n=${d.normals} i=${d.indices} g=${d.groups}`;
    expect(out).toMatchInlineSnapshot(`
      {
        "coil tube": "3016v 5248t curve=22430911 v=607549172 n=1795098989 i=2528217029 g=3823400037",
        "degenerate": "3016v 5248t curve=281995042 v=3482627721 n=3611377989 i=2528217029 g=3823400037",
        "helix elliptical": "3016v 5248t curve=3662355717 v=981845174 n=1525901127 i=2528217029 g=3823400037",
        "helix rounded": "2024v 3456t curve=3749712538 v=1711728876 n=1683080486 i=1912354837 g=2503719653",
        "nucleic tube": "896v 1488t curve=3214162229 v=3865116821 n=3567090349 i=2968838797 g=2020398021",
        "one linear segment": "200v 192t curve=4068756566 v=2239648245 n=1270629583 i=239117461 g=3150690405",
        "ribbon with arrows": "400v 320t curve=3092466637 v=2628284509 n=1848946621 i=1027626085 g=1653499653",
        "rounded, flat": "1288v 2080t curve=979255584 v=177480482 n=955150066 i=526155157 g=2695361893",
        "sheet with arrows": "1488v 1304t curve=22430911 v=810127339 n=2760085901 i=910468645 g=3331605509",
        "sheet, square": "672v 496t curve=3886727563 v=679132635 n=2103513057 i=3034569877 g=2020009541",
        "tubular helix, round caps": "3135v 5360t curve=22430911 v=4244870438 n=95562765 i=961893269 g=206096711",
      }
    `);
  });
});
