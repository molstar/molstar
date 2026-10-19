/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ChunkedArray } from '@molstar/core/data/util';
import { MeshBuilder } from '../mesh/mesh-builder.js';
import { addTube } from '../mesh/builder/tube.js';
import { addSheet } from '../mesh/builder/sheet.js';
import { addRibbon } from '../mesh/builder/ribbon.js';
import {
  type CurveSegmentControls,
  createCurveSegmentState,
  interpolateCurveSegment,
  interpolateSizes,
} from '../../../repr/structure/visual/util/polymer/curve-segment.js';

// Pins the cartoon geometry bit for bit: the hashes are of the output before the array reads were inlined.

/** FNV-1a over the Float32 bit patterns */
function hash(...arrays: ArrayLike<number>[]) {
  let h = 0x811c9dc5;
  for (const a of arrays) {
    const u = new Uint32Array(Float32Array.from(a).buffer);
    for (let i = 0; i < u.length; ++i) {
      h ^= u[i];
      h = Math.imul(h, 16777619) >>> 0;
    }
  }
  return h.toString(16);
}

function stateHash(state: MeshBuilder.State) {
  const used = (c: ChunkedArray<number, number>) =>
    Array.from(ChunkedArray.compact(c, true)).slice(0, c.elementCount * c.elementSize);
  return hash(used(state.vertices), used(state.normals), used(state.indices), used(state.groups));
}

/** A bent, twisting segment like a helix turn */
function controls(k: number): CurveSegmentControls {
  const p = (i: number) => Vec3.create(Math.cos(i * 0.9 + k) * 2.3, Math.sin(i * 0.9 + k) * 2.3, i * 1.5 + k * 0.1);
  const [p0, p1, p2, p3, p4] = [p(0), p(1), p(2), p(3), p(4)];
  return {
    secStrucFirst: k % 2 === 0,
    secStrucLast: k % 3 === 0,
    p0,
    p1,
    p2,
    p3,
    p4,
    d12: Vec3.normalize(Vec3(), Vec3.create(0.3, 1, -0.2 + k * 0.05)),
    d23: Vec3.normalize(Vec3(), Vec3.create(-0.4, 0.8, 0.5)),
  };
}

function curve(k: number, linearSegments = 12) {
  const state = createCurveSegmentState(linearSegments);
  interpolateCurveSegment(state, controls(k), 0.5, k % 2 ? 0.5 : 0.25);
  interpolateSizes(state, 0.2, 0.4, 0.3, 0.2, 0.9, 0.5, 0.5);
  return state;
}

describe('cartoon geometry', () => {
  it('interpolateCurveSegment', () => {
    const hashes = [0, 1, 2, 3].map((k) => {
      const s = curve(k);
      return hash(s.curvePoints, s.tangentVectors, s.normalVectors, s.binormalVectors, s.widthValues, s.heightValues);
    });
    expect(hashes).toEqual(['462fbb7', 'a65e2a19', 'd3140ee', '69a8bfde']);
  });

  it('addTube: elliptical and rounded, caps, round caps', () => {
    const state = MeshBuilder.createState(16, 16);
    for (let k = 0; k < 4; ++k) {
      const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues, linearSegments } = curve(k);
      state.currentGroup = k;
      addTube(
        state,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        16,
        widthValues,
        heightValues,
        true,
        true,
        'elliptical',
      );
      addTube(
        state,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        12,
        widthValues,
        heightValues,
        k % 2 === 0,
        true,
        'rounded',
      );
      addTube(
        state,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        8,
        widthValues,
        heightValues,
        true,
        k % 2 === 1,
        'elliptical',
        true,
      );
    }
    expect(stateHash(state)).toBe('27435d87');
  });

  it('addTube: odd segment counts, double round caps', () => {
    const state = MeshBuilder.createState(16, 16);
    for (const linearSegments of [1, 7, 9]) {
      const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues } = curve(
        linearSegments,
        linearSegments,
      );
      state.currentGroup = linearSegments;
      const args = [curvePoints, normalVectors, binormalVectors, linearSegments] as const;
      addTube(state, ...args, 10, widthValues, heightValues, true, true, 'elliptical', true);
      addTube(state, ...args, 6, widthValues, heightValues, true, true, 'rounded', true);
    }
    expect(stateHash(state)).toBe('b630c356');
  });

  it('addSheet: with and without arrow, caps', () => {
    const state = MeshBuilder.createState(16, 16);
    for (let k = 0; k < 4; ++k) {
      const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues, linearSegments } = curve(k);
      state.currentGroup = k;
      addSheet(
        state,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        widthValues,
        heightValues,
        0,
        true,
        true,
      );
      addSheet(
        state,
        curvePoints,
        normalVectors,
        binormalVectors,
        linearSegments,
        widthValues,
        heightValues,
        1.7,
        k % 2 === 0,
        true,
      );
    }
    expect(stateHash(state)).toBe('63f8e3fa');
  });

  it('addRibbon: with and without arrow', () => {
    const state = MeshBuilder.createState(16, 16);
    for (let k = 0; k < 4; ++k) {
      const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues, linearSegments } = curve(k);
      state.currentGroup = k;
      addRibbon(state, curvePoints, normalVectors, binormalVectors, linearSegments, widthValues, heightValues, 0);
      addRibbon(state, curvePoints, normalVectors, binormalVectors, linearSegments, widthValues, heightValues, 1.7);
    }
    expect(stateHash(state)).toBe('1ec27261');
  });
});
