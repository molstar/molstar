/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CylinderStride } from '../cylinders/cylinders.js';
import { LineStride } from '../lines/lines.js';

/** The corners of a cylinder's box as the per-vertex mapping had them */
const BoxMapping = [-1, 1, -1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1];

/**
 * The per-vertex buffers of `count` cylinders (6 vertices each) as cylinders had them, from their data (`CylinderStride`
 * floats each): tests written against vertices read them so (the items' data of cylinders and lines as vertices)
 */
export function perVertexCylinders(data: Float32Array, count: number) {
  const start = new Float32Array(count * 18);
  const end = new Float32Array(count * 18);
  const group = new Float32Array(count * 6);
  const scale = new Float32Array(count * 6);
  const cap = new Float32Array(count * 6);
  const colorMode = new Float32Array(count * 6);
  const mapping = new Float32Array(count * 18);
  for (let i = 0; i < count; i++) {
    const o = i * CylinderStride;
    for (let k = 0; k < 6; k++) {
      const v = i * 6 + k;
      start.set(data.subarray(o, o + 3), v * 3);
      end.set(data.subarray(o + 4, o + 7), v * 3);
      group[v] = data[o + 3];
      scale[v] = data[o + 7];
      cap[v] = data[o + 8];
      colorMode[v] = data[o + 9];
    }
    mapping.set(BoxMapping, i * 18);
  }
  return { start, end, group, scale, cap, colorMode, mapping };
}

/** `perVertexCylinders` of the values of a cylinders render object (12 indices per cylinder) */
export function cylinderVertices(values: any) {
  return perVertexCylinders(values.tCylinder.ref.value.array as Float32Array, values.drawCount.ref.value / 12);
}

/** The corners of a line as the per-vertex mapping had them (left and right of the start, then of the end) */
const LineMapping = [-1, -1, 1, -1, -1, 1, 1, 1];

/**
 * The per-vertex buffers of `count` lines (4 vertices each) as segments had them, from their data (`LineStride` floats
 * each): tests written against vertices read them so
 */
export function perVertexLines(data: Float32Array, count: number) {
  const start = new Float32Array(count * 12);
  const end = new Float32Array(count * 12);
  const group = new Float32Array(count * 4);
  const mapping = new Float32Array(count * 8);
  for (let i = 0; i < count; i++) {
    const o = i * LineStride;
    for (let k = 0; k < 4; k++) {
      const v = i * 4 + k;
      start.set(data.subarray(o, o + 3), v * 3);
      end.set(data.subarray(o + 4, o + 7), v * 3);
      group[v] = data[o + 3];
    }
    mapping.set(LineMapping, i * 8);
  }
  return { start, end, group, mapping };
}
