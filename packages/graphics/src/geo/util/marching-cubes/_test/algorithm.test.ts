/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Tensor } from '@molstar/core/math/linear-algebra';
import { computeMarchingCubesLines, computeMarchingCubesMesh, type MarchingCubesParams } from '../algorithm.js';
import { createWrappedTensor } from '../../../../repr/volume/util.js';
import type { Mesh } from '../../../geometry/mesh/mesh.js';
import type { Lines } from '../../../geometry/lines/lines.js';

// Pins the exact output (counts and bit-level checksums) so that changes to how the field is read cannot move a vertex.

const Dims = [17, 14, 11];
const Orders = [
  [0, 1, 2],
  [0, 2, 1],
  [1, 0, 2],
  [1, 2, 0],
  [2, 0, 1],
  [2, 1, 0],
];
const Iso = 0.3;

/** A shell around a cavity (for the floodfill), a blob cut by the grid faces (gradient clamping) and some ripple */
function value(i: number, j: number, k: number) {
  const r = Math.hypot(i - 7, j - 6, k - 5);
  const shell = 1 - Math.abs(r - 3.5) / 1.5;
  const blob = 1.2 * Math.exp(-((i - 16) ** 2 + (j - 1) ** 2 + (k - 0.5) ** 2) / 6);
  return Math.max(shell, blob) + 0.05 * Math.sin(1.3 * i + 0.7 * j + 2.1 * k);
}

function field(
  order: number[],
  f: (i: number, j: number, k: number) => number = value,
  ctor: Tensor.ArrayCtor = Float32Array,
) {
  const space = Tensor.Space(Dims, order, ctor);
  // these fields take the strided path
  expect(Tensor.strides(space)).toBeDefined();
  const data = space.create();
  for (let k = 0; k < Dims[2]; k++) {
    for (let j = 0; j < Dims[1]; j++) {
      for (let i = 0; i < Dims[0]; i++) space.set(data, i, j, k, f(i, j, k));
    }
  }
  return Tensor.create(space, data);
}

/** -2 (ignore) on the low x faces, -1 (no id) on a lattice, ids 0..4 elsewhere */
function idField(order: number[]) {
  return field(order, (i, j, k) => (i <= 2 ? -2 : (i + j + k) % 3 === 0 ? -1 : (7 * j + k) % 5), Int32Array);
}

/** FNV-1a over the bit patterns of the first `count` elements */
function hash(array: Float32Array | Uint32Array, count: number) {
  const bits = new Uint32Array(array.buffer, array.byteOffset, count);
  let h = 0x811c9dc5;
  for (let i = 0; i < count; i++) h = Math.imul(h ^ bits[i], 0x01000193) >>> 0;
  return h.toString(16).padStart(8, '0');
}

function meshSummary(mesh: Mesh) {
  const { vertexCount, triangleCount } = mesh;
  return {
    vertexCount,
    triangleCount,
    vertices: hash(mesh.vertexBuffer.ref.value, vertexCount * 3),
    normals: hash(mesh.normalBuffer.ref.value, vertexCount * 3),
    groups: hash(mesh.groupBuffer.ref.value, vertexCount),
    indices: hash(mesh.indexBuffer.ref.value, triangleCount * 3),
  };
}

function linesSummary(lines: Lines) {
  const { lineCount } = lines;
  return {
    lineCount,
    starts: hash(lines.startBuffer.ref.value, lineCount * 4 * 3),
    ends: hash(lines.endBuffer.ref.value, lineCount * 4 * 3),
    groups: hash(lines.groupBuffer.ref.value, lineCount * 4),
  };
}

async function mesh(params: MarchingCubesParams) {
  return meshSummary(await computeMarchingCubesMesh(params).run());
}

async function lines(params: MarchingCubesParams) {
  return linesSummary(await computeMarchingCubesLines(params).run());
}

describe('marching cubes', () => {
  it('gives the same mesh for every axis order', async () => {
    const reference = await mesh({ isoLevel: Iso, scalarField: field([0, 1, 2]) });
    expect(reference).toMatchInlineSnapshot(`
      {
        "groups": "843254cd",
        "indices": "b04fd405",
        "normals": "d6c184a8",
        "triangleCount": 1127,
        "vertexCount": 578,
        "vertices": "981276b5",
      }
    `);
    for (const order of Orders) {
      expect(await mesh({ isoLevel: Iso, scalarField: field(order) }), `order ${order}`).toEqual(reference);
    }
  });

  it('reads Float64 data', async () => {
    const reference = await mesh({ isoLevel: Iso, scalarField: field([0, 1, 2], value, Float64Array) });
    expect(reference).toMatchInlineSnapshot(`
      {
        "groups": "843254cd",
        "indices": "b04fd405",
        "normals": "4276cc80",
        "triangleCount": 1127,
        "vertexCount": 578,
        "vertices": "826f3e3b",
      }
    `);
    for (const order of Orders) {
      const scalarField = field(order, value, Float64Array);
      expect(await mesh({ isoLevel: Iso, scalarField }), `order ${order}`).toEqual(reference);
    }
  });

  it('flips normals and winding for a negative iso level', async () => {
    const negated = (i: number, j: number, k: number) => -value(i, j, k);
    const reference = await mesh({ isoLevel: -Iso, scalarField: field([0, 1, 2], negated) });
    expect(reference).toMatchInlineSnapshot(`
      {
        "groups": "843254cd",
        "indices": "57220765",
        "normals": "d6c184a8",
        "triangleCount": 1127,
        "vertexCount": 578,
        "vertices": "981276b5",
      }
    `);
    for (const order of Orders) {
      expect(await mesh({ isoLevel: -Iso, scalarField: field(order, negated) }), `order ${order}`).toEqual(reference);
    }
  });

  it('reads the floodfilled values, not the data', async () => {
    const plain = await mesh({ isoLevel: Iso, scalarField: field([0, 1, 2]) });
    for (const mode of ['inside', 'outside'] as const) {
      const reference = await mesh({
        isoLevel: Iso,
        scalarField: Tensor.createFloodfilled(field([0, 1, 2]), Iso, mode),
      });
      expect(reference, mode).not.toEqual(plain);
      for (const order of Orders) {
        const scalarField = Tensor.createFloodfilled(field(order), Iso, mode);
        expect(await mesh({ isoLevel: Iso, scalarField }), `${mode} ${order}`).toEqual(reference);
      }
    }
    expect(
      await mesh({ isoLevel: Iso, scalarField: Tensor.createFloodfilled(field([0, 1, 2]), Iso, 'inside') }),
    ).toMatchInlineSnapshot(`
      {
        "groups": "ce1a4015",
        "indices": "28d5a690",
        "normals": "a019c42b",
        "triangleCount": 879,
        "vertexCount": 452,
        "vertices": "3a466614",
      }
    `);
    expect(
      await mesh({ isoLevel: Iso, scalarField: Tensor.createFloodfilled(field([0, 1, 2]), Iso, 'outside') }),
    ).toMatchInlineSnapshot(`
      {
        "groups": "78fc8efd",
        "indices": "c3861887",
        "normals": "e57dcb6e",
        "triangleCount": 248,
        "vertexCount": 126,
        "vertices": "c4fc8764",
      }
    `);
  });

  it('reads a wrapped tensor through its get', async () => {
    expect(Tensor.strides(createWrappedTensor(field([0, 1, 2])).space)).toBeUndefined();
    const reference = await mesh({ isoLevel: Iso, scalarField: createWrappedTensor(field([0, 1, 2])) });
    expect(reference).toMatchInlineSnapshot(`
      {
        "groups": "b592799d",
        "indices": "3347ee72",
        "normals": "3e3b5e5c",
        "triangleCount": 1220,
        "vertexCount": 646,
        "vertices": "47ed6c0b",
      }
    `);
    for (const order of Orders) {
      const scalarField = createWrappedTensor(field(order));
      expect(await mesh({ isoLevel: Iso, scalarField }), `order ${order}`).toEqual(reference);
    }
  });

  it('assigns groups from the id field, skipping -2 and falling back on -1', async () => {
    const reference = await mesh({ isoLevel: Iso, scalarField: field([0, 1, 2]), idField: idField([0, 1, 2]) });
    expect(reference).toMatchInlineSnapshot(`
      {
        "groups": "29b8b82f",
        "indices": "0f128387",
        "normals": "c95610fe",
        "triangleCount": 1105,
        "vertexCount": 573,
        "vertices": "7b040041",
      }
    `);
    for (const order of Orders) {
      // the id field in its own layout, and in a different one from the scalar field
      const scalarField = field(order);
      expect(await mesh({ isoLevel: Iso, scalarField, idField: idField(order) }), `order ${order}`).toEqual(reference);
      const other = Orders[(Orders.indexOf(order) + 3) % Orders.length];
      expect(await mesh({ isoLevel: Iso, scalarField, idField: idField(other) }), `order ${order}/${other}`).toEqual(
        reference,
      );
    }
  });

  it('keeps to the bottomLeft/topRight range', async () => {
    const range = { bottomLeft: [2, 1, 1], topRight: [16, 13, 9] };
    const reference = await mesh({ isoLevel: Iso, scalarField: field([0, 1, 2]), ...range });
    expect(reference).toMatchInlineSnapshot(`
      {
        "groups": "7c2314a7",
        "indices": "00e4d65a",
        "normals": "6e1b1a69",
        "triangleCount": 907,
        "vertexCount": 487,
        "vertices": "7e1b2b2b",
      }
    `);
    for (const order of Orders) {
      expect(await mesh({ isoLevel: Iso, scalarField: field(order), ...range }), `order ${order}`).toEqual(reference);
    }
    const withIds = await mesh({ isoLevel: Iso, scalarField: field([1, 2, 0]), idField: idField([2, 0, 1]), ...range });
    expect(withIds).toMatchInlineSnapshot(`
      {
        "groups": "071a954d",
        "indices": "917d772a",
        "normals": "1d3d0c25",
        "triangleCount": 885,
        "vertexCount": 482,
        "vertices": "3d71b98f",
      }
    `);
  });

  it('builds lines with the edge filter', async () => {
    const reference = await lines({ isoLevel: Iso, scalarField: field([0, 1, 2]) });
    expect(reference).toMatchInlineSnapshot(`
      {
        "ends": "163f137d",
        "groups": "fb8f2e95",
        "lineCount": 1145,
        "starts": "6556d8b5",
      }
    `);
    const range = { bottomLeft: [2, 1, 1], topRight: [16, 13, 9] };
    const ranged = await lines({ isoLevel: Iso, scalarField: field([0, 1, 2]), idField: idField([0, 1, 2]), ...range });
    expect(ranged).toMatchInlineSnapshot(`
      {
        "ends": "7365e025",
        "groups": "f21ccef5",
        "lineCount": 919,
        "starts": "09dd4f59",
      }
    `);
    const negated = await lines({ isoLevel: -Iso, scalarField: field([0, 1, 2], (i, j, k) => -value(i, j, k)) });
    expect(negated).toMatchInlineSnapshot(`
      {
        "ends": "b1da13e1",
        "groups": "fb8f2e95",
        "lineCount": 1145,
        "starts": "69e95ad1",
      }
    `);
    for (const order of Orders) {
      expect(await lines({ isoLevel: Iso, scalarField: field(order) }), `order ${order}`).toEqual(reference);
      const scalarField = field(order);
      expect(await lines({ isoLevel: Iso, scalarField, idField: idField(order), ...range }), `order ${order}`).toEqual(
        ranged,
      );
    }
  });

  it('gives what get gives on random fields', async () => {
    // a copy of the space has no strides, so the copy is read through `get`, as before strides
    const viaGet = (t: Tensor) => Tensor.create({ ...t.space }, t.data);
    let seed = 7;
    const random = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
    const special = [Number.NaN, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, 0.5];
    for (let n = 0; n < 60; n++) {
      const dims = [1 + Math.floor(random() * 6), 1 + Math.floor(random() * 6), 1 + Math.floor(random() * 6)];
      const order = Orders[n % Orders.length];
      const space = Tensor.Space(dims, order, n % 2 ? Float64Array : Float32Array);
      const data = space.create();
      for (let i = 0; i < data.length; i++) {
        data[i] = random() < 0.05 ? special[Math.floor(random() * special.length)] : random();
      }
      const ids = Tensor.Space(dims, Orders[(n + 2) % Orders.length], Int32Array);
      const idData = ids.create();
      for (let i = 0; i < idData.length; i++) idData[i] = Math.floor(random() * 6) - 2;
      const scalarField = Tensor.create(space, data);
      const idField = n % 3 ? Tensor.create(ids, idData) : undefined;
      const isoLevel = n % 4 ? 0.5 : -0.5;
      const range =
        n % 5 === 0 ? { bottomLeft: [0, 0, 1], topRight: [dims[0], dims[1], Math.max(1, dims[2] - 1)] } : {};
      const params = { isoLevel, scalarField, idField, ...range };
      const label = `${n}: ${dims} ${order}`;
      const expected = await mesh({ ...params, scalarField: viaGet(scalarField), idField: idField && viaGet(idField) });
      expect(await mesh(params), label).toEqual(expected);
      const expectedLines = await lines({ ...params, scalarField: viaGet(scalarField) });
      expect(await lines(params), label).toEqual(expectedLines);
    }
  });

  it('reads through get outside the grid', async () => {
    // a range beyond the grid: `get`'s clamping of the central differences differs from the strided one there (with x
    // fastest, reads past the x face land in the next row instead of past the data)
    const scalarField = field([2, 1, 0]);
    const params = { isoLevel: Iso, scalarField, bottomLeft: [0, 0, 0], topRight: [Dims[0] + 1, Dims[1], Dims[2]] };
    const expected = await mesh({ ...params, scalarField: Tensor.create({ ...scalarField.space }, scalarField.data) });
    expect(await mesh(params)).toEqual(expected);
  });
});
