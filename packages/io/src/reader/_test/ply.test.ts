/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { parsePly } from '../ply/parser.js';
import { PlyTable, PlyList, type PlyFile } from '../ply/schema.js';
import { utf8Encode } from '@molstar/core/util/utf8';

const plyString = `ply
format ascii 1.0
comment file created by MegaMol
element vertex 6
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
property uchar alpha
property float nx
property float ny
property float nz
property int atomid
property uchar contactcount_r
property uchar contactcount_g
property uchar contactcount_b
property uchar contactsteps_r
property uchar contactsteps_g
property uchar contactsteps_b
property uchar hbonds_r
property uchar hbonds_g
property uchar hbonds_b
property uchar hbondsteps_r
property uchar hbondsteps_g
property uchar hbondsteps_b
property uchar molcount_r
property uchar molcount_g
property uchar molcount_b
property uchar spots_r
property uchar spots_g
property uchar spots_b
property uchar rmsf_r
property uchar rmsf_g
property uchar rmsf_b
element face 2
property list uchar int vertex_index
end_header
130.901 160.016 163.033 90 159 210 255 -0.382 -0.895 -0.231 181 21 100 150 24 102 151 20 100 150 20 100 150 30 106 154 20 100 150 171 196 212
131.372 159.778 162.83 90 159 210 255 -0.618 -0.776 -0.129 178 21 100 150 24 102 151 20 100 150 20 100 150 30 106 154 20 100 150 141 177 199
131.682 159.385 163.089 90 159 210 255 -0.773 -0.579 -0.259 180 21 100 150 24 102 151 20 100 150 20 100 150 30 106 154 20 100 150 172 196 212
131.233 160.386 162.11 90 159 210 255 -0.708 -0.383 -0.594 178 21 100 150 24 102 151 20 100 150 20 100 150 30 106 154 20 100 150 141 177 199
130.782 160.539 162.415 90 159 210 255 -0.482 -0.459 -0.746 181 21 100 150 24 102 151 20 100 150 20 100 150 30 106 154 20 100 150 171 196 212
131.482 160.483 161.621 90 159 210 255 -0.832 -0.431 -0.349 179 21 100 150 24 102 151 20 100 150 20 100 150 30 106 154 20 100 150 171 196 212
3 0 2 1
3 3 5 4
`;

const plyCubeString = `ply
format ascii 1.0
comment test cube
element vertex 24
property float32 x
property float32 y
property float32 z
property uint32 material_index
element face 6
property list uint8 int32 vertex_indices
element material 6
property uint8 red
property uint8 green
property uint8 blue
end_header
-1 -1 -1 0
1 -1 -1 0
1 1 -1 0
-1 1 -1 0
1 -1 1 1
-1 -1 1 1
-1 1 1 1
1 1 1 1
1 1 1 2
1 1 -1 2
1 -1 -1 2
1 -1 1 2
-1 1 -1 3
-1 1 1 3
-1 -1 1 3
-1 -1 -1 3
-1 1 1 4
-1 1 -1 4
1 1 -1 4
1 1 1 4
1 -1 1 5
1 -1 -1 5
-1 -1 -1 5
-1 -1 1 5
4 0 1 2 3
4 4 5 6 7
4 8 9 10 11
4 12 13 14 15
4 16 17 18 19
4 20 21 22 23
255 0 0
0 255 0
0 0 255
255 255 0
0 255 255
255 0 255
`;

describe('ply reader', () => {
  it('basic', async () => {
    const parsed = await parsePly(plyString).run();
    if (parsed.isError) return;
    const plyFile = parsed.result;

    const vertex = plyFile.getElement('vertex') as PlyTable;
    if (!vertex) return;
    const x = vertex.getProperty('x');
    if (!x) return;
    expect(x.value(0)).toEqual(130.901);

    const face = plyFile.getElement('face') as PlyList;
    if (!face) return;
    expect(face.value(0)).toEqual({ count: 3, entries: [0, 2, 1] });
    expect(face.value(1)).toEqual({ count: 3, entries: [3, 5, 4] });

    expect.assertions(3);
  });

  it('material', async () => {
    const parsed = await parsePly(plyCubeString).run();
    if (parsed.isError) return;
    const plyFile = parsed.result;

    const vertex = plyFile.getElement('vertex') as PlyTable;
    if (!vertex) return;
    expect(vertex.rowCount).toBe(24);

    const face = plyFile.getElement('face') as PlyList;
    if (!face) return;
    expect(face.rowCount).toBe(6);

    const material = plyFile.getElement('face') as PlyTable;
    if (!material) return;
    expect(face.rowCount).toBe(6);

    expect.assertions(3);
  });

  it('ascii as Uint8Array', async () => {
    const bytes = new TextEncoder().encode(plyString);
    const parsed = await parsePly(bytes).run();
    if (parsed.isError) throw new Error(parsed.message);
    const plyFile = parsed.result;

    const vertex = plyFile.getElement('vertex') as PlyTable;
    if (!vertex) return;
    const x = vertex.getProperty('x');
    if (!x) return;
    expect(x.value(0)).toEqual(130.901);

    const face = plyFile.getElement('face') as PlyList;
    if (!face) return;
    expect(face.value(0)).toEqual({ count: 3, entries: [0, 2, 1] });
    expect(face.value(1)).toEqual({ count: 3, entries: [3, 5, 4] });

    expect.assertions(3);
  });

  it('parses ASCII the same from a string and from its UTF-8 bytes', async () => {
    const dump = (ply: PlyFile) =>
      ply.elementNames.map((name) => {
        const e = ply.getElement(name)!;
        if (e.kind === 'list')
          return Array.from({ length: e.rowCount }, (_, i) => [...e.value(i).entries.slice(0, e.value(i).count)]);
        return e.propertyNames.map((p) => e.getProperty(p)!.toArray());
      });
    for (const text of [plyString, plyCubeString, plyString.replace(/\n/g, '\r\n')]) {
      const fromString = await parsePly(text).run();
      const fromBytes = await parsePly(utf8Encode(text)).run();
      if (fromString.isError || fromBytes.isError) throw new Error('parse failed');
      expect(fromBytes.result.comments).toEqual(fromString.result.comments);
      expect(fromBytes.result.elementNames).toEqual(fromString.result.elementNames);
      expect(dump(fromBytes.result)).toEqual(dump(fromString.result));
    }
  });

  it('reads a format line by its words (spaces, another version)', async () => {
    for (const line of ['format ascii 1.0 ', 'format  ascii 1.1']) {
      const parsed = await parsePly(utf8Encode(plyString.replace('format ascii 1.0', line))).run();
      if (parsed.isError) throw new Error(parsed.message);
      const x = (parsed.result.getElement('vertex') as PlyTable).getProperty('x')!;
      expect(x.value(0)).toEqual(130.901);
    }
  });

  it('reads a last list row without a line break', async () => {
    const text = plyString.trimEnd();
    const parsed = await parsePly(utf8Encode(text)).run();
    if (parsed.isError) throw new Error(parsed.message);
    const expected = await parsePly(utf8Encode(plyString)).run();
    if (expected.isError) throw new Error(expected.message);
    const face = parsed.result.getElement('face') as PlyList;
    const expectedFace = expected.result.getElement('face') as PlyList;
    expect(face.rowCount).toBe(expectedFace.rowCount);
    for (let i = 0; i < face.rowCount; i++) expect(face.value(i)).toEqual(expectedFace.value(i));
  });

  it('decodes a non-ASCII comment', async () => {
    const text = plyString.replace('format ascii 1.0\n', 'format ascii 1.0\ncomment Molécule Ångström\n');
    const parsed = await parsePly(utf8Encode(text)).run();
    if (parsed.isError) throw new Error(parsed.message);
    expect(parsed.result.comments[0]).toBe('Molécule Ångström');
    const x = (parsed.result.getElement('vertex') as PlyTable).getProperty('x')!;
    expect(x.value(0)).toEqual(130.901);
  });

  it('binary little-endian', async () => {
    // Build a minimal binary little-endian PLY with 3 vertices and 1 triangle face
    // Header (ASCII)
    const header =
      'ply\nformat binary_little_endian 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_index\nend_header\n';
    const headerBytes = new TextEncoder().encode(header);

    // Vertex data: 3 x (float32 x, float32 y, float32 z) = 36 bytes
    const vertexBuf = new ArrayBuffer(36);
    const vdv = new DataView(vertexBuf);
    vdv.setFloat32(0, 1.0, true);
    vdv.setFloat32(4, 0.0, true);
    vdv.setFloat32(8, 0.0, true);
    vdv.setFloat32(12, 0.0, true);
    vdv.setFloat32(16, 1.0, true);
    vdv.setFloat32(20, 0.0, true);
    vdv.setFloat32(24, 0.0, true);
    vdv.setFloat32(28, 0.0, true);
    vdv.setFloat32(32, 1.0, true);

    // Face data: uchar count (3), then 3 x int32 indices = 1 + 12 = 13 bytes
    const faceBuf = new ArrayBuffer(13);
    const fdv = new DataView(faceBuf);
    fdv.setUint8(0, 3);
    fdv.setInt32(1, 0, true);
    fdv.setInt32(5, 1, true);
    fdv.setInt32(9, 2, true);

    // Concatenate header + vertex + face bytes
    const totalLength = headerBytes.length + 36 + 13;
    const data = new Uint8Array(totalLength);
    data.set(headerBytes, 0);
    data.set(new Uint8Array(vertexBuf), headerBytes.length);
    data.set(new Uint8Array(faceBuf), headerBytes.length + 36);

    const parsed = await parsePly(data).run();
    if (parsed.isError) throw new Error(parsed.message);
    const plyFile = parsed.result;

    const vertex = plyFile.getElement('vertex') as PlyTable;
    if (!vertex) return;
    expect(vertex.rowCount).toBe(3);
    const x = vertex.getProperty('x');
    if (!x) return;
    expect(x.value(0)).toBeCloseTo(1.0);
    expect(x.value(1)).toBeCloseTo(0.0);

    const face = plyFile.getElement('face') as PlyList;
    if (!face) return;
    expect(face.rowCount).toBe(1);
    const faceVal = face.value(0);
    expect(faceVal.count).toBe(3);
    expect(faceVal.entries[0]).toBe(0);
    expect(faceVal.entries[1]).toBe(1);
    expect(faceVal.entries[2]).toBe(2);

    expect.assertions(8);
  });

  /** a binary little-endian PLY of `elements` (name, properties of a type each, rows), the header padded by `pad` */
  function binaryPly(
    elements: { name: string; props: [string, 'float' | 'uint' | 'int' | 'uchar'][]; rows: number[][] }[],
    pad = '',
  ) {
    let header = `ply\nformat binary_little_endian 1.0\ncomment ${pad}\n`;
    for (const e of elements) {
      header += `element ${e.name} ${e.rows.length}\n`;
      for (const [name, type] of e.props) header += `property ${type} ${name}\n`;
    }
    header += 'end_header\n';
    const headerBytes = new TextEncoder().encode(header);
    const size = { float: 4, uint: 4, int: 4, uchar: 1 };
    let length = headerBytes.length;
    for (const e of elements) for (const _ of e.rows) for (const [, t] of e.props) length += size[t];
    const data = new Uint8Array(length);
    data.set(headerBytes);
    const dv = new DataView(data.buffer);
    let o = headerBytes.length;
    for (const e of elements) {
      for (const row of e.rows) {
        e.props.forEach(([, t], j) => {
          if (t === 'float') dv.setFloat32(o, row[j], true);
          else if (t === 'uint') dv.setUint32(o, row[j], true);
          else if (t === 'int') dv.setInt32(o, row[j], true);
          else dv.setUint8(o, row[j]);
          o += size[t];
        });
      }
    }
    return data;
  }

  it('binary tables of one type: aligned or not, one after another, mixed ones value by value', async () => {
    const vertex = {
      name: 'vertex',
      props: [
        ['x', 'float'],
        ['y', 'float'],
        ['f_dc_0', 'float'],
      ] as [string, 'float'][],
      rows: [
        [1.5, -2, 0.25],
        [3, 4.5, -1],
      ],
    };
    const packed = {
      name: 'packed',
      props: [
        ['a', 'uint'],
        ['b', 'uint'],
      ] as [string, 'uint'][],
      rows: [
        [4294967295, 7],
        [0, 2147483648],
      ],
    };
    const mixed = {
      name: 'mixed',
      props: [
        ['u', 'uchar'],
        ['i', 'int'],
        ['v', 'uint'],
      ] as [string, 'uchar' | 'int' | 'uint'][],
      rows: [[200, -3, 9]],
    };
    // headers of every length modulo 4: the rows at every alignment
    for (const pad of ['', 'a', 'ab', 'abc']) {
      const parsed = await parsePly(binaryPly([vertex, packed, mixed], pad)).run();
      if (parsed.isError) throw new Error(parsed.message);
      const v = parsed.result.getElement('vertex') as PlyTable;
      expect(
        [0, 1].map((i) => [
          v.getProperty('x')!.value(i),
          v.getProperty('y')!.value(i),
          v.getProperty('f_dc_0')!.value(i),
        ]),
      ).toEqual(vertex.rows);
      const p = parsed.result.getElement('packed') as PlyTable;
      expect([0, 1].map((i) => [p.getProperty('a')!.value(i), p.getProperty('b')!.value(i)])).toEqual(packed.rows);
      const m = parsed.result.getElement('mixed') as PlyTable;
      expect([m.getProperty('u')!.value(0), m.getProperty('i')!.value(0), m.getProperty('v')!.value(0)]).toEqual(
        mixed.rows[0],
      );
    }
  });

  it("reads a binary table of a Node Buffer at an offset in a larger one (a pool's view)", async () => {
    const vertex = {
      name: 'vertex',
      props: [
        ['x', 'float'],
        ['y', 'float'],
      ] as [string, 'float'][],
      rows: [
        [1.5, -2],
        [3, 4.5],
      ],
    };
    for (const pad of ['', 'a', 'ab', 'abc']) {
      const bytes = binaryPly([vertex], pad);
      // the bytes from 1 in a buffer of other bytes before and after them: `slice` of a Buffer a view of it
      const pool = new Uint8Array(bytes.length + 16).fill(0xee);
      pool.set(bytes, 1);
      const data = Buffer.from(pool.buffer, 1, bytes.length);
      const parsed = await parsePly(data).run();
      if (parsed.isError) throw new Error(parsed.message);
      const v = parsed.result.getElement('vertex') as PlyTable;
      expect([0, 1].map((i) => [v.getProperty('x')!.value(i), v.getProperty('y')!.value(i)])).toEqual(vertex.rows);
    }
  });

  it('rejects a binary table past the end of the data', async () => {
    const data = binaryPly([{ name: 'vertex', props: [['x', 'float']], rows: [[1], [2]] }]);
    await expect(parsePly(data.subarray(0, data.length - 2)).run()).rejects.toThrow(/ends before its table/);
  });

  it('binary big-endian', async () => {
    const header =
      'ply\nformat binary_big_endian 1.0\nelement vertex 2\nproperty float x\nproperty float y\nend_header\n';
    const headerBytes = new TextEncoder().encode(header);

    // 2 vertices, each has float32 x + float32 y (big-endian)
    const buf = new ArrayBuffer(16);
    const dv = new DataView(buf);
    dv.setFloat32(0, 3.5, false);
    dv.setFloat32(4, -1.0, false);
    dv.setFloat32(8, 0.25, false);
    dv.setFloat32(12, 7.0, false);

    const data = new Uint8Array(headerBytes.length + 16);
    data.set(headerBytes, 0);
    data.set(new Uint8Array(buf), headerBytes.length);

    const parsed = await parsePly(data).run();
    if (parsed.isError) throw new Error(parsed.message);
    const plyFile = parsed.result;

    const vertex = plyFile.getElement('vertex') as PlyTable;
    if (!vertex) return;
    const x = vertex.getProperty('x');
    const y = vertex.getProperty('y');
    if (!x || !y) return;
    expect(x.value(0)).toBeCloseTo(3.5);
    expect(y.value(0)).toBeCloseTo(-1.0);
    expect(x.value(1)).toBeCloseTo(0.25);
    expect(y.value(1)).toBeCloseTo(7.0);

    expect.assertions(4);
  });
});
