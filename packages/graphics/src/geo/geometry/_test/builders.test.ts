/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { CylindersBuilder } from '../cylinders/cylinders-builder.js';
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

const point = (i: number) => Vec3.create(Math.sin(i) * 10, Math.cos(1.3 * i) * 10, (i % 17) - 8);

function cylinders(n: number, initialCount: number, chunkSize: number) {
  const builder = CylindersBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++)
    builder.add(...point(i), ...point(i + 1), 1 + (i % 3) * 0.1, i % 2 === 0, i % 3 === 0, 0, i);
  const c = builder.getCylinders();
  const vertexCount = c.cylinderCount * 6;
  return {
    counts: [c.cylinderCount],
    buffers: [
      hash(c.startBuffer.ref.value, vertexCount * 3),
      hash(c.endBuffer.ref.value, vertexCount * 3),
      hash(c.groupBuffer.ref.value, vertexCount),
      hash(c.scaleBuffer.ref.value, vertexCount),
      hash(c.capBuffer.ref.value, vertexCount),
      hash(c.mappingBuffer.ref.value, vertexCount * 3),
      hash(c.indexBuffer.ref.value, c.cylinderCount * 12),
    ],
    capacity: c.startBuffer.ref.value.length / 3 / 6,
  };
}

function lines(n: number, initialCount: number, chunkSize: number) {
  const builder = LinesBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++) builder.add(...point(i), ...point(i + 1), i);
  const l = builder.getLines();
  const vertexCount = l.lineCount * 4;
  return {
    counts: [l.lineCount],
    buffers: [
      hash(l.startBuffer.ref.value, vertexCount * 3),
      hash(l.endBuffer.ref.value, vertexCount * 3),
      hash(l.groupBuffer.ref.value, vertexCount),
      hash(l.mappingBuffer.ref.value, vertexCount * 2),
      hash(l.indexBuffer.ref.value, l.lineCount * 6),
    ],
    capacity: l.startBuffer.ref.value.length / 3 / 4,
  };
}

function spheres(n: number, initialCount: number, chunkSize: number) {
  const builder = SpheresBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++) builder.add(...point(i), i);
  const s = builder.getSpheres();
  return {
    counts: [s.sphereCount],
    buffers: [hash(s.centerBuffer.ref.value, s.sphereCount * 3), hash(s.groupBuffer.ref.value, s.sphereCount)],
    capacity: s.centerBuffer.ref.value.length / 3,
  };
}

function points(n: number, initialCount: number, chunkSize: number) {
  const builder = PointsBuilder.create(initialCount, chunkSize);
  for (let i = 0; i < n; i++) builder.add(...point(i), i);
  const p = builder.getPoints();
  return {
    counts: [p.pointCount],
    buffers: [hash(p.centerBuffer.ref.value, p.pointCount * 3), hash(p.groupBuffer.ref.value, p.pointCount)],
    capacity: p.centerBuffer.ref.value.length / 3,
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
