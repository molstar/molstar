/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../../../mol-task';
import { Mutable } from '../../../mol-util/type-helpers';
import { NetcdfReader } from '../../common/netcdf/reader';
import { ReaderResult as Result } from '../result';

export interface NctrajFile {
    coordinates: number[][],
    velocities?: number[][],
    forces?: number[][],
    cell_lengths?: number[][],
    cell_angles?: number[][],
    time?: number[],
    timeOffset: number,
    deltaTime: number
}

async function parseInternal(data: Uint8Array) {
    // http://ambermd.org/netcdf/nctraj.xhtml

    const nc = new NetcdfReader(data);

    const f: Mutable<NctrajFile> = {
        coordinates: [],
        time: [],
        timeOffset: 0,
        deltaTime: 1
    };

    for (const c of nc.getDataVariable('coordinates')) f.coordinates.push(c);

    if (nc.hasDataVariable('velocities')) {
        const velocities: number[][] = [];
        for (const v of nc.getDataVariable('velocities')) velocities.push(v);
        f.velocities = velocities;
    }

    if (nc.hasDataVariable('forces')) {
        const forces: number[][] = [];
        for (const f of nc.getDataVariable('forces')) forces.push(f);
        f.forces = forces;
    }

    if (nc.hasDataVariable('cell_lengths')) {
        const cell_lengths: number[][] = [];
        for (const l of nc.getDataVariable('cell_lengths')) cell_lengths.push(l);
        f.cell_lengths = cell_lengths;
    }

    if (nc.hasDataVariable('cell_angles')) {
        const cell_angles: number[][] = [];
        for (const a of nc.getDataVariable('cell_angles')) cell_angles.push(a);
        f.cell_angles = cell_angles;
    }

    if (nc.hasDataVariable('time')) {
        const time: number[] = [];
        for (const t of nc.getDataVariable('time')) time.push(t);
        f.time = time;
    }

    if (f.time) {
        if (f.time.length >= 1) {
            f.timeOffset = f.time[0];
        }
        if (f.time.length >= 2) {
            f.deltaTime = f.time[1] - f.time[0];
        }
    }

    return f;
}

export function parseNctraj(data: Uint8Array) {
    return Task.create<Result<NctrajFile>>('Parse NCTRAJ', async ctx => {
        try {
            ctx.update({ canAbort: true, message: 'Parsing trajectory...' });
            const file = await parseInternal(data);
            return Result.success(file);
        } catch (e) {
            return Result.error('' + e);
        }
    });
}