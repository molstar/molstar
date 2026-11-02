/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CylindersBuilder } from '../cylinders/cylinders-builder.js';
import { CylinderStride } from '../cylinders/cylinders.js';
import { perVertexCylinders, perVertexLines } from './item-vertices.js';
import { LineStride } from '../lines/lines.js';
import { BoxItems, createItemIndexArray, LineItems } from '../../../gl/quad-indices.js';
import { LinesBuilder } from '../lines/lines-builder.js';
import { MeshBuilder } from '../mesh/mesh-builder.js';
import { addSphere } from '../mesh/builder/sphere.js';
import { addCylinder } from '../mesh/builder/cylinder.js';
import { PointsBuilder } from '../points/points-builder.js';
import { SpheresBuilder } from '../spheres/spheres-builder.js';

// Pins the geometry the builders give (counts and the used part of every buffer), whatever the size estimates the
// callers pass and however much of each buffer is allocated (perf-audit 08/P3, 08/P8).

function hash(array: ArrayLike<number> & { buffer?: ArrayBufferLike }, length: number) {
  const a = array as Float32Array | Uint32Array;
  const used = a.subarray(0, length);
  return (
    hashFnv32a(new Uint8Array(used.buffer, used.byteOffset, used.byteLength) as unknown as number[]) >>> 0
  ).toString(16);
}

/** Centers (xyz) and groups of `count` items of a position group array (xyz and group each), as separate buffers */
function deinterleave(positionGroup: Float32Array, count: number) {
  const centers = new Float32Array(count * 3);
  const groups = new Float32Array(count);
  for (let i = 0; i < count; i++) {
    centers.set(positionGroup.subarray(i * 4, i * 4 + 3), i * 3);
    groups[i] = positionGroup[i * 4 + 3];
  }
  return { centers, groups };
}

const point = (i: number) => Vec3.create(Math.sin(i) * 10, Math.cos(1.3 * i) * 10, (i % 17) - 8);

function cylinders(n: number, initialCount: number, chunkSize: number) {
  const builder = CylindersBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++)
    builder.add(...point(i), ...point(i + 1), 1 + (i % 3) * 0.1, i % 2 === 0, i % 3 === 0, 0, i);
  const c = builder.getCylinders();
  const vertexCount = c.cylinderCount * 6;
  // the per-vertex buffers cylinders had (6 vertices each), from the data per cylinder: the hashes are theirs
  const v = perVertexCylinders(c.cylinderData.ref.value.array, c.cylinderCount);
  return {
    counts: [c.cylinderCount],
    buffers: [
      hash(v.start, vertexCount * 3),
      hash(v.end, vertexCount * 3),
      hash(v.group, vertexCount),
      hash(v.scale, vertexCount),
      hash(v.cap, vertexCount),
      hash(v.mapping, vertexCount * 3),
      // the shared box indices are the indices cylinders had
      hash(createItemIndexArray(c.cylinderCount, BoxItems), c.cylinderCount * 12),
    ],
    capacity: c.cylinderData.ref.value.array.length / CylinderStride,
  };
}

function lines(n: number, initialCount: number, chunkSize: number) {
  const builder = LinesBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++) builder.add(...point(i), ...point(i + 1), i);
  const l = builder.getLines();
  const vertexCount = l.lineCount * 4;
  // the per-vertex buffers lines had (4 vertices each), from the data per line: the hashes are theirs
  const v = perVertexLines(l.lineData.ref.value.array, l.lineCount);
  return {
    counts: [l.lineCount],
    buffers: [
      hash(v.start, vertexCount * 3),
      hash(v.end, vertexCount * 3),
      hash(v.group, vertexCount),
      hash(v.mapping, vertexCount * 2),
      // the shared line indices are the indices lines had
      hash(createItemIndexArray(l.lineCount, LineItems), l.lineCount * 6),
    ],
    capacity: l.lineData.ref.value.array.length / LineStride,
  };
}

function spheres(n: number, initialCount: number, chunkSize: number) {
  const builder = SpheresBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++) builder.add(...point(i), i);
  const s = builder.getSpheres();
  // the hashes of the centers and groups as separate buffers: the contents of the position groups are theirs
  const { centers, groups } = deinterleave(s.positionGroup.ref.value.array, s.sphereCount);
  return {
    counts: [s.sphereCount],
    buffers: [hash(centers, s.sphereCount * 3), hash(groups, s.sphereCount)],
    capacity: s.positionGroup.ref.value.array.length / 4,
  };
}

function points(n: number, initialCount: number, chunkSize: number) {
  const builder = PointsBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++) builder.add(...point(i), i);
  const p = builder.getPoints();
  const { centers, groups } = deinterleave(p.positionGroup.ref.value.array, p.pointCount);
  return {
    counts: [p.pointCount],
    buffers: [hash(centers, p.pointCount * 3), hash(groups, p.pointCount)],
    capacity: p.positionGroup.ref.value.array.length / 4,
  };
}

function mesh(n: number, initialCount: number, chunkSize: number) {
  const state = MeshBuilder.createState(initialCount, chunkSize);
  for (let i = 0; i < n; i++) {
    state.currentGroup = i;
    if (i % 2) addSphere(state, point(i), 1 + (i % 4) * 0.2, 1);
    else addCylinder(state, point(i), point(i + 1), 1, { radiusTop: 0.3, radiusBottom: 0.4, radialSegments: 8 });
  }
  const m = MeshBuilder.getMesh(state);
  return {
    counts: [m.vertexCount, m.triangleCount],
    buffers: [
      hash(m.vertexBuffer.ref.value, m.vertexCount * 3),
      hash(m.normalBuffer.ref.value, m.vertexCount * 3),
      hash(m.groupBuffer.ref.value, m.vertexCount),
      hash(m.indexBuffer.ref.value, m.triangleCount * 3),
    ],
    vertexCapacity: m.vertexBuffer.ref.value.length / 3,
    triangleCapacity: m.indexBuffer.ref.value.length / 3,
  };
}

// item counts as the callers pass them: estimates that fit, overshoot 3-6x, or fall short
const sizes: [string, number, number][] = [
  ['exact', 1000, 250],
  ['overshoot', 6000, 1500],
  ['short', 100, 25],
];

describe('geometry builders', () => {
  it.each(sizes)('cylinders, %s estimate', (_, initialCount, chunkSize) => {
    const { counts, buffers } = cylinders(1000, initialCount, chunkSize);
    expect(counts).toEqual([1000]);
    expect(buffers).toEqual(['d2d2e5d', 'c0c0e6e9', 'aed38295', '90fbc9d9', 'c8acd315', '792a8d45', '5251c38d']);
  });

  it.each(sizes)('lines, %s estimate', (_, initialCount, chunkSize) => {
    const { counts, buffers } = lines(1000, initialCount, chunkSize);
    expect(counts).toEqual([1000]);
    expect(buffers).toEqual(['6666485', '63621c0d', 'a8e68ca5', 'e1d9cfc5', '2a9cd985']);
  });

  it.each(sizes)('spheres, %s estimate', (_, initialCount, chunkSize) => {
    const { counts, buffers } = spheres(1000, initialCount, chunkSize);
    expect(counts).toEqual([1000]);
    expect(buffers).toEqual(['91985a36', '9fa62ea8']);
  });

  it.each(sizes)('points, %s estimate', (_, initialCount, chunkSize) => {
    const { counts, buffers } = points(1000, initialCount, chunkSize);
    expect(counts).toEqual([1000]);
    expect(buffers).toEqual(['91985a36', '9fa62ea8']);
  });

  it.each(sizes)('mesh, %s estimate', (_, initialCount, chunkSize) => {
    const { counts, buffers } = mesh(200, initialCount * 20, chunkSize * 2);
    expect(counts).toEqual([5800, 9600]);
    expect(buffers).toEqual(['36ebb248', '66875299', 'fc9eaf45', '1f7a316f']);
  });

  it('trims buffers an estimate overshoots', () => {
    const m = mesh(200, 20000, 1000);
    expect([m.vertexCapacity, m.triangleCapacity]).toEqual([5800, 9600]);
    const fitting = mesh(200, 6000, 1000);
    // 97% and 80% used: kept as they are
    expect([fitting.vertexCapacity, fitting.triangleCapacity]).toEqual([6000, 12000]);
  });
});
