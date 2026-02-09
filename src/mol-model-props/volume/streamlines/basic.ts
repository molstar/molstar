/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PCG } from '../../../mol-data/util/hash-functions';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Grid } from '../../../mol-model/volume/grid';
import { CustomProperty } from '../../common/custom-property';
import { Volume } from '../../../mol-model/volume';
import { StreamlinePoint, Streamlines } from './shared';
import { RuntimeContext } from '../../../mol-task/execution/runtime-context';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3fromArray = Vec3.fromArray;
const v3distance = Vec3.distance;
const v3setMagnitude = Vec3.setMagnitude;
const v3create = Vec3.create;
const v3copy = Vec3.copy;
const v3add = Vec3.add;

export const BasicStreamlineCalculationParams = {
    seedDensity: PD.Numeric(10, { min: 1, max: 30, step: 1 }, { description: 'Percentage of cells with seed.' }),
    stepSize: PD.Numeric(0.35, { min: 0.01, max: 1, step: 0.01 }, { description: 'Step size in grid space.' }),
};
export type BasicStreamlineCalculationParams = typeof BasicStreamlineCalculationParams
export type BasicStreamlineCalculationProps = PD.Values<BasicStreamlineCalculationParams>

type GetGradient = (gridCoords: Vec3, out: Vec3) => boolean

const d = Vec3();
const p = Vec3();
const g = Vec3();
const prev = Vec3();

/**
 * Basic tracer with fixed step size in grid space
 */
function traceOneDirection(out: StreamlinePoint[], grid: Grid, getGradient: GetGradient, seed: Vec3, stepSize: number, dir: 1 | -1): number {
    const { space } = grid.cells;
    const [nx, ny, nz] = space.dimensions;
    const o = space.dataOffset;

    const step = dir * stepSize;
    const maxSteps = Math.max(nx, ny, nz) * 5 / stepSize;
    const t = stepSize / 1.1; // tolerance for writing points

    const seenVoxel = new Map<number, number>();
    const maxVoxelVisits = Math.ceil(2 / stepSize);

    v3copy(p, seed);

    let written = 0;
    for (let i = 0; i < maxSteps; ++i) {
        // boundary check
        if (p[0] < 1 || p[0] > nx - 2 ||
            p[1] < 1 || p[1] > ny - 2 ||
            p[2] < 1 || p[2] > nz - 2) break;

        if (!getGradient(p, g)) break;

        const key = o(Math.round(p[0]), Math.round(p[1]), Math.round(p[2]));
        const c = (seenVoxel.get(key) || 0) + 1;
        seenVoxel.set(key, c);
        if (c > maxVoxelVisits) break;

        v3setMagnitude(d, g, step);
        // write midpoint (keeps lines smooth and avoids tiny box at start)
        const midX = p[0] + 0.5 * d[0];
        const midY = p[1] + 0.5 * d[1];
        const midZ = p[2] + 0.5 * d[2];
        const m = v3create(midX, midY, midZ);
        if (v3distance(prev, m) >= t) {
            out.push(m);
            v3copy(prev, m);
        }
        written++;

        v3add(p, p, d); // advance
    }

    return written;
}

function traceStreamlineBothDirs(grid: Grid, getGradient: GetGradient, seed: Vec3, stepSize: number): StreamlinePoint[] {
    const line: StreamlinePoint[] = [];
    const nBack = traceOneDirection(line, grid, getGradient, seed, stepSize, -1);
    if (nBack > 1) line.reverse();
    traceOneDirection(line, grid, getGradient, seed, stepSize, +1);
    return line;
}

async function computeBasicStreamlines(ctx: RuntimeContext, grid: Grid, props: BasicStreamlineCalculationProps): Promise<Streamlines> {
    const { space } = grid.cells;
    const [nx, ny, nz] = space.dimensions;

    const { seedDensity, stepSize } = props;
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));

    // bounds avoiding edges
    const xStart = 1, xEnd = nx - 2 - seedStep;
    const yStart = 1, yEnd = ny - 2 - seedStep;
    const zStart = 1, zEnd = nz - 2 - seedStep;

    const pcg = new PCG();

    const seeds: number[] = [];
    for (let z = zStart; z <= zEnd; z += seedStep) {
        for (let y = yStart; y <= yEnd; y += seedStep) {
            for (let x = xStart; x <= xEnd; x += seedStep) {
                seeds.push(
                    x + pcg.float() * seedStep,
                    y + pcg.float() * seedStep,
                    z + pcg.float() * seedStep
                );
            }
        }
    }

    // shuffle in-place by triplets [x,y,z]
    for (let i = (seeds.length / 3) - 1; i > 0; i--) {
        const j = Math.floor(pcg.float() * (i + 1));
        const ia = i * 3, ja = j * 3;
        // swap 3 numbers at once
        const tx = seeds[ia], ty = seeds[ia + 1], tz = seeds[ia + 2];
        seeds[ia] = seeds[ja]; seeds[ia + 1] = seeds[ja + 1]; seeds[ia + 2] = seeds[ja + 2];
        seeds[ja] = tx; seeds[ja + 1] = ty; seeds[ja + 2] = tz;
    }

    await ctx.update({ isIndeterminate: false, current: 0, max: seeds.length / 3 });

    const getGradient = Grid.makeGetInterpolatedGradient(grid);

    const lines: Streamlines = [];
    const pos = Vec3();
    for (let i = 0; i < seeds.length; i += 3) {
        if (ctx.shouldUpdate) await ctx.update({ current: i / 3 + 1 });

        v3fromArray(pos, seeds, i);
        const line = traceStreamlineBothDirs(grid, getGradient, pos, stepSize);
        if (line.length * stepSize >= 3) lines.push(line);
    }
    return lines;
}

//

export async function calculateBasicStreamlines(ctx: CustomProperty.Context, volume: Volume, props: BasicStreamlineCalculationProps): Promise<Streamlines> {
    return computeBasicStreamlines(ctx.runtime, volume.grid, props);
}
