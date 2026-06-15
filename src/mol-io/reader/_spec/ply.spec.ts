/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parsePly } from '../ply/parser';
import { PlyTable, PlyList } from '../ply/schema';

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

    it('binary little-endian', async () => {
        // Build a minimal binary little-endian PLY with 3 vertices and 1 triangle face
        // Header (ASCII)
        const header = 'ply\nformat binary_little_endian 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_index\nend_header\n';
        const headerBytes = new TextEncoder().encode(header);

        // Vertex data: 3 x (float32 x, float32 y, float32 z) = 36 bytes
        const vertexBuf = new ArrayBuffer(36);
        const vdv = new DataView(vertexBuf);
        vdv.setFloat32(0, 1.0, true); vdv.setFloat32(4, 0.0, true); vdv.setFloat32(8, 0.0, true);
        vdv.setFloat32(12, 0.0, true); vdv.setFloat32(16, 1.0, true); vdv.setFloat32(20, 0.0, true);
        vdv.setFloat32(24, 0.0, true); vdv.setFloat32(28, 0.0, true); vdv.setFloat32(32, 1.0, true);

        // Face data: uchar count (3), then 3 x int32 indices = 1 + 12 = 13 bytes
        const faceBuf = new ArrayBuffer(13);
        const fdv = new DataView(faceBuf);
        fdv.setUint8(0, 3);
        fdv.setInt32(1, 0, true); fdv.setInt32(5, 1, true); fdv.setInt32(9, 2, true);

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

    it('binary big-endian', async () => {
        const header = 'ply\nformat binary_big_endian 1.0\nelement vertex 2\nproperty float x\nproperty float y\nend_header\n';
        const headerBytes = new TextEncoder().encode(header);

        // 2 vertices, each has float32 x + float32 y (big-endian)
        const buf = new ArrayBuffer(16);
        const dv = new DataView(buf);
        dv.setFloat32(0, 3.5, false); dv.setFloat32(4, -1.0, false);
        dv.setFloat32(8, 0.25, false); dv.setFloat32(12, 7.0, false);

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
