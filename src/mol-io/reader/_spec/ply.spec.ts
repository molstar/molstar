/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Ply from '../ply/parser'
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
`

describe('ply reader', () => {
    it('basic', async () => {
        const parsed = await Ply(plyString).run();
        if (parsed.isError) return;
        const plyFile = parsed.result;

        const vertex = plyFile.getElement('vertex') as PlyTable
        if (!vertex) return
        const x = vertex.getProperty('x')
        if (!x) return
        console.log('x', x.toArray())

        const face = plyFile.getElement('face') as PlyList
        if (!face) return
        const f0 = face.value(0)
        console.log('f0', f0)
        const f1 = face.value(1)
        console.log('f1', f1)
    });
});