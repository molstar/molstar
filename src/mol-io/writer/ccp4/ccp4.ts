/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Minimal MRC/CCP4 binary writer (MODE=2, float32).
 * Spec: https://www.ccpem.ac.uk/mrc_format/mrc2014.php
 */

import { Grid } from '../../../mol-model/volume';
import { Vec3 } from '../../../mol-math/linear-algebra';

const HEADER_BYTES = 1024;

export namespace CCP4Writer {
    export function writeMrc(grid: Grid, data: Uint8Array | Float32Array): ArrayBuffer {
        const [nx, ny, nz] = grid.cells.space.dimensions as [number, number, number];
        const voxelCount = nx * ny * nz;
        const buf = new ArrayBuffer(HEADER_BYTES + voxelCount * 4);
        const i32 = new Int32Array(buf);
        const f32 = new Float32Array(buf);

        i32[0] = nx; // NC
        i32[1] = ny; // NR
        i32[2] = nz; // NS
        i32[3] = 2;  // MODE = 2 (float32)

        i32[4] = 0; // NXSTART
        i32[5] = 0; // NYSTART
        i32[6] = 0; // NZSTART

        i32[7] = nx; // NX
        i32[8] = ny; // NY
        i32[9] = nz; // NZ

        const g2c = Grid.getGridToCartesianTransform(grid);
        const ax = Vec3(), ay = Vec3(), az = Vec3();
        Vec3.set(ax, g2c[0], g2c[1], g2c[2]);
        Vec3.set(ay, g2c[4], g2c[5], g2c[6]);
        Vec3.set(az, g2c[8], g2c[9], g2c[10]);

        f32[10] = Vec3.magnitude(ax) * nx; // XLEN
        f32[11] = Vec3.magnitude(ay) * ny; // YLEN
        f32[12] = Vec3.magnitude(az) * nz; // ZLEN

        const axN = Vec3.normalize(Vec3(), ax);
        const ayN = Vec3.normalize(Vec3(), ay);
        const azN = Vec3.normalize(Vec3(), az);
        const toDeg = 180 / Math.PI;
        const clamp = (v: number) => Math.min(1, Math.max(-1, v));
        f32[13] = Math.acos(clamp(Vec3.dot(ayN, azN))) * toDeg; // ALPHA
        f32[14] = Math.acos(clamp(Vec3.dot(axN, azN))) * toDeg; // BETA
        f32[15] = Math.acos(clamp(Vec3.dot(axN, ayN))) * toDeg; // GAMMA

        i32[16] = 1; // MAPC
        i32[17] = 2; // MAPR
        i32[18] = 3; // MAPS

        let min = Infinity, max = -Infinity, sum = 0;
        for (let i = 0; i < voxelCount; i++) {
            const v = data[i];
            if (v < min) min = v;
            if (v > max) max = v;
            sum += v;
        }
        f32[19] = min;              // AMIN
        f32[20] = max;              // AMAX
        f32[21] = sum / voxelCount; // AMEAN

        i32[22] = 0; // ISPG
        i32[23] = 0; // NSYMBT

        // MRC2014 origin (voxel 0,0,0 world position)
        f32[49] = g2c[12]; // XORIGIN
        f32[50] = g2c[13]; // YORIGIN
        f32[51] = g2c[14]; // ZORIGIN

        const bytes = new Uint8Array(buf);
        bytes[208] = 0x4D; // M
        bytes[209] = 0x41; // A
        bytes[210] = 0x50; // P
        bytes[211] = 0x20; // ' '
        bytes[212] = 0x44; // machine stamp: little-endian float
        bytes[213] = 0x44;
        bytes[214] = 0x00;
        bytes[215] = 0x00;

        const dataView = new Float32Array(buf, HEADER_BYTES, voxelCount);
        for (let i = 0; i < voxelCount; i++) dataView[i] = data[i];

        return buf;
    }
}
