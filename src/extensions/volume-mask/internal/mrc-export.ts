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

export function writeMrc(grid: Grid, maskData: Uint8Array | Float32Array): ArrayBuffer {
    const [nx, ny, nz] = grid.cells.space.dimensions as [number, number, number];
    const voxelCount = nx * ny * nz;
    const buf = new ArrayBuffer(HEADER_BYTES + voxelCount * 4);
    const i32 = new Int32Array(buf);
    const f32 = new Float32Array(buf);

    // Columns / rows / sections  (fast → slow)
    i32[0] = nx; // NC
    i32[1] = ny; // NR
    i32[2] = nz; // NS
    i32[3] = 2; // MODE = 2 (float32)

    // Grid start (0-based)
    i32[4] = 0; // NXSTART
    i32[5] = 0; // NYSTART
    i32[6] = 0; // NZSTART

    // Grid sampling
    i32[7] = nx; // NX
    i32[8] = ny; // NY
    i32[9] = nz; // NZ

    // Cell dimensions and angles in Ångströms / degrees
    const g2c = Grid.getGridToCartesianTransform(grid);
    const ax = Vec3(), ay = Vec3(), az = Vec3();
    Vec3.set(ax, g2c[0], g2c[1], g2c[2]);
    Vec3.set(ay, g2c[4], g2c[5], g2c[6]);
    Vec3.set(az, g2c[8], g2c[9], g2c[10]);
    const xLen = Vec3.magnitude(ax) * nx;
    const yLen = Vec3.magnitude(ay) * ny;
    const zLen = Vec3.magnitude(az) * nz;

    f32[10] = xLen; // XLEN
    f32[11] = yLen; // YLEN
    f32[12] = zLen; // ZLEN
    // Derive true cell angles from axis vectors (handles non-orthogonal volumes)
    const axN = Vec3.normalize(Vec3(), ax);
    const ayN = Vec3.normalize(Vec3(), ay);
    const azN = Vec3.normalize(Vec3(), az);
    const toDeg = 180 / Math.PI;
    const clamp = (v: number) => Math.min(1, Math.max(-1, v));
    f32[13] = Math.acos(clamp(Vec3.dot(ayN, azN))) * toDeg; // ALPHA (Y∧Z)
    f32[14] = Math.acos(clamp(Vec3.dot(axN, azN))) * toDeg; // BETA  (X∧Z)
    f32[15] = Math.acos(clamp(Vec3.dot(axN, ayN))) * toDeg; // GAMMA (X∧Y)

    // Axis mapping (standard: columns=X, rows=Y, sections=Z)
    i32[16] = 1; // MAPC
    i32[17] = 2; // MAPR
    i32[18] = 3; // MAPS

    // Statistics
    let min = Infinity, max = -Infinity, sum = 0;
    for (let i = 0; i < voxelCount; i++) { const v = maskData[i]; if (v < min) min = v; if (v > max) max = v; sum += v; }
    f32[19] = min; // AMIN
    f32[20] = max; // AMAX
    f32[21] = sum / voxelCount; // AMEAN

    i32[22] = 0; // ISPG (space group 0 = image stack)
    i32[23] = 0; // NSYMBT (no symmetry bytes)

    // MRC2014 origin in Ångströms — position of voxel (0,0,0) in world space.
    // g2c is column-major; translation is at indices 12–14.
    // ChimeraX and other MRC2014 readers use these fields to overlay maps correctly.
    f32[49] = g2c[12]; // XORIGIN  (byte 196)
    f32[50] = g2c[13]; // YORIGIN  (byte 200)
    f32[51] = g2c[14]; // ZORIGIN  (byte 204)

    // Machine stamp (0x44 0x44 = little-endian float)
    const bytes = new Uint8Array(buf);
    bytes[212] = 0x44;
    bytes[213] = 0x44;
    bytes[214] = 0x00;
    bytes[215] = 0x00;

    // "MAP " label at offset 208
    bytes[208] = 0x4D; // M
    bytes[209] = 0x41; // A
    bytes[210] = 0x50; // P
    bytes[211] = 0x20; // ' '

    // Voxel data (float32) starting at byte 1024
    const dataView = new Float32Array(buf, HEADER_BYTES, voxelCount);
    for (let i = 0; i < voxelCount; i++) dataView[i] = maskData[i];

    return buf;
}

export function downloadMrc(grid: Grid, maskData: Uint8Array | Float32Array, filename = 'mask.mrc') {
    const buf = writeMrc(grid, maskData);
    const blob = new Blob([buf], { type: 'application/octet-stream' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}
