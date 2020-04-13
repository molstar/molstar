/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Ply from '../ply/parser';
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
        const parsed = await Ply(plyString).run();
        if (parsed.isError) return;
        const plyFile = parsed.result;

        const vertex = plyFile.getElement('vertex') as PlyTable;
        if (!vertex) return;
        const x = vertex.getProperty('x');
        if (!x) return;
        expect(x.value(0)).toEqual(130.901);

        const face = plyFile.getElement('face') as PlyList;
        if (!face) return;
        expect(face.value(0)).toEqual({ count: 3, entries: [0, 2, 1]});
        expect(face.value(1)).toEqual({ count: 3, entries: [3, 5, 4]});

        expect.assertions(3);
    });

    it('material', async () => {
        const parsed = await Ply(plyCubeString).run();
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
});