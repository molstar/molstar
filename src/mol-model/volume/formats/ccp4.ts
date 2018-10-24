/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData } from '../data'
import { Task } from 'mol-task';
import { SpacegroupCell, Box3D } from 'mol-math/geometry';
import { Tensor, Vec3 } from 'mol-math/linear-algebra';
import { Ccp4File } from 'mol-io/reader/ccp4/schema';
import { degToRad } from 'mol-math/misc';

// TODO implement voxelSize parameter

// TODO support for rescaling of values
// based on uglymol (https://github.com/uglymol/uglymol) by Marcin Wojdyr (wojdyr)
// if the file was converted by mapmode2to0 - scale the data
// if (intView[ 39 ] === -128 && intView[ 40 ] === 127) {
//     // scaling f(x)=b1*x+b0 such that f(-128)=min and f(127)=max
//     const b1 = (header.DMAX - header.DMIN) / 255.0
//     const b0 = 0.5 * (header.DMIN + header.DMAX + b1)
//     for (let j = 0, jl = data.length; j < jl; ++j) {
//         data[ j ] = b1 * data[ j ] + b0
//     }
// }

function volumeFromCcp4(source: Ccp4File): Task<VolumeData> {
    return Task.create<VolumeData>('Parse Volume Data', async ctx => {
        const { header, values } = source;

        const cell = SpacegroupCell.create(
            header.ISPG,
            Vec3.create(header.xLength, header.yLength, header.zLength),
            Vec3.create(degToRad(header.alpha), degToRad(header.beta), degToRad(header.gamma))
        );

        const axis_order_fast_to_slow = Vec3.create(header.MAPC - 1, header.MAPR - 1, header.MAPS - 1);
        const normalizeOrder = Tensor.convertToCanonicalAxisIndicesFastToSlow(axis_order_fast_to_slow);

        const grid = [header.NX, header.NY, header.NZ];
        const extent = normalizeOrder([header.NC, header.NR, header.NS]);

        const gridOrigin: number[] = normalizeOrder([header.NCSTART, header.NRSTART, header.NSSTART]);

        // TODO: this is code from LiteMol, but not sure about the part based on originXYZ. The 
        // if (header.originX === 0.0 && header.originY === 0.0 && header.originZ === 0.0) {
        //     gridOrigin = normalizeOrder([header.NCSTART, header.NRSTART, header.NSSTART]);
        // } else {
        //     // Use ORIGIN records rather than old n[xyz]start records
        //     //   http://www2.mrc-lmb.cam.ac.uk/image2000.html
        //     // XXX the ORIGIN field is only used by the EM community, and
        //     //     has undefined meaning for non-orthogonal maps and/or
        //     //     non-cubic voxels, etc.
        //     gridOrigin = [header.originX, header.originY, header.originZ];
        // }

        const origin_frac = Vec3.create(gridOrigin[0] / grid[0], gridOrigin[1] / grid[1], gridOrigin[2] / grid[2]);
        const dimensions_frac = Vec3.create(extent[0] / grid[0], extent[1] / grid[1], extent[2] / grid[2]);

        const space = Tensor.Space(extent, Tensor.invertAxisOrder(axis_order_fast_to_slow), header.MODE === 0 ? Int8Array : Float32Array);
        const data = Tensor.create(space, Tensor.Data1(values));

        // TODO Calculate stats? When to trust header data?
        // Min/max/mean are reliable (based on LiteMol/DensityServer usage)
        // These, however, calculate sigma, so no data on that.

        return {
            cell,
            fractionalBox: Box3D.create(origin_frac, Vec3.add(Vec3.zero(), origin_frac, dimensions_frac)),
            data,
            dataStats: {
                min: header.AMIN,
                max: header.AMAX,
                mean: header.AMEAN,
                sigma: header.ARMS
            }
        };
    });
}

export { volumeFromCcp4 }