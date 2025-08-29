/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 *
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { Cylinders } from '../../mol-geo/geometry/cylinders/cylinders';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
import { Volume, Grid } from '../../mol-model/volume';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider, VolumeKey } from './representation';
import { VisualUpdateState } from '../util';
import { VisualContext } from '../visual';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Theme, ThemeRegistryContext } from '../../mol-theme/theme';
import { CylindersBuilder } from '../../mol-geo/geometry/cylinders/cylinders-builder';
import { createVolumeCellLocationIterator, eachVolumeLoci } from './util';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../representation';
import Stream from 'stream';

// Constants
// const DEFAULT_FILTER_THRESHOLD = 0.2;
// const MIN_GRADIENT_MAGNITUDE = 1e-6;
const OBLATION_SPACING = 2;

export const VolumeGradientParams = {
    isoValue: Volume.IsoValueParam,
    seedDensity: PD.Numeric(10, { min: 1, max: 100, step: 1 }, { description: 'Seeds per dimension for streamlines.' }),  
    maxSteps: PD.Numeric(1000, { min: 1, max: 2000, step: 1 }, { description: 'Maximum number of steps for streamlines.' }),  
    stepSize: PD.Numeric(0.35, { min: 0.01, max: 10, step: 0.01 }, { description: 'Step size for streamlines.' }),  
    minSpeed: PD.Numeric(0.001, { min: 0, max: 1, step: 1e-6 }, { description: 'Minimum speed for streamlines.' }),  
    minLevel: PD.Numeric(-2.0, { min: -10, max: 10, step: 0.1 }, { description: 'Minimum level for streamline tracing.' }),  
    maxLevel: PD.Numeric(2.0, { min: -10, max: 10, step: 0.1 }, { description: 'Maximum level for streamline tracing.' }),  
    algorithm: PD.Select('simple', PD.arrayToOptions(['simple', 'advanced'] as const), { description: 'Streamline algorithm to use.' }),  
};

export type VolumeGradientParams = typeof VolumeGradientParams
export type VolumeGradientProps = PD.Values<VolumeGradientParams>

export const VolumeLinesParams = {
    ...Lines.Params,
    ...VolumeGradientParams,
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
};
export type VolumeLinesParams = typeof VolumeLinesParams
export type VolumeLinesProps = PD.Values<VolumeLinesParams>


export function VolumeLinesVisual(materialId: number): VolumeVisual<VolumeLinesParams> {
    return VolumeVisual<Lines, VolumeLinesParams>({
        defaultProps: PD.getDefaultValues(VolumeLinesParams),
        createGeometry: createVolumeLinesMesh,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<VolumeLinesParams>, currentProps: PD.Values<VolumeLinesParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel
            );
        },
        geometryUtils: Lines.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeLinesParams>, webgl?: WebGLContext) => {
            return !!webgl;
        }
    }, materialId);
}

export const VolumeCylindersParams = {
    ...Cylinders.Params,
    ...VolumeGradientParams,
    radius: PD.Numeric(1.0, { min: 0.1, max: 5, step: 0.1 }, { description: 'Radius scale of the cylinders.' }),
};
export type VolumeCylindersParams = typeof VolumeCylindersParams
export type VolumeCylindersProps = PD.Values<VolumeCylindersParams>

export function VolumeCylindersImpostorVisual(materialId: number): VolumeVisual<VolumeCylindersParams> {
    return VolumeVisual<Cylinders, VolumeCylindersParams>({
        defaultProps: PD.getDefaultValues(VolumeCylindersParams),
        createGeometry: createVolumeCylindersImpostor,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<VolumeCylindersParams>, currentProps: PD.Values<VolumeCylindersParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats)||
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.radius !== currentProps.radius ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel
            );
        },
        geometryUtils: Cylinders.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeCylindersParams>, webgl?: WebGLContext) => {
            return !webgl;
        }
    }, materialId);
}

export function createVolumeCylindersImpostor(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeCylindersProps, points?: Cylinders): Cylinders {
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const stream_lines : StreamlinePoint[][] = [];
    if (props.algorithm === 'advanced') {
        getStreamLineAdvanced(volume, props, stream_lines);
    } else {
        getStreamLineSimple(volume, props, stream_lines);
    }
    const count = stream_lines.length;
    const builder = CylindersBuilder.create(count, Math.ceil(count / 2), points);
    // --- iterate streamlines ---
    const start = Vec3();
    const start1 = Vec3();
    const end = Vec3();
    const end1 = Vec3();
    for (let s=0; s< stream_lines.length; s++) {
        for (let i=0; i< stream_lines[s].length -1; i++) {
            Vec3.set(start1, stream_lines[s][i].x, stream_lines[s][i].y, stream_lines[s][i].z);
            Vec3.set(end1, stream_lines[s][i + 1].x, stream_lines[s][i + 1].y, stream_lines[s][i + 1].z);
            Vec3.transformMat4(start, start1, gridToCartn);
            Vec3.transformMat4(end, end1, gridToCartn);
            builder.add(start[0], start[1], start[2],
                        end[0], end[1], end[2],
                        props.radius, true, true, 2, s);
        }
    }
    const pt = builder.getCylinders();
    pt.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return pt;
}

export function createVolumeLinesMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeLinesProps, points?: Lines): Lines {
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const stream_lines : StreamlinePoint[][] = [];
    if (props.algorithm === 'advanced') {
        getStreamLineAdvanced(volume, props, stream_lines);
    } else {
        getStreamLineSimple(volume, props, stream_lines);
    }
    const count = stream_lines.length;
    const builder = LinesBuilder.create(count, Math.ceil(count / 2), points);
    // --- iterate streamlines ---
    const start = Vec3();
    const start1 = Vec3();
    const end = Vec3();
    const end1 = Vec3();
    for (let s=0; s< stream_lines.length; s++) {
        for (let i=0; i< stream_lines[s].length -1; i++) {
            Vec3.set(start1, stream_lines[s][i].x, stream_lines[s][i].y, stream_lines[s][i].z);
            Vec3.set(end1, stream_lines[s][i + 1].x, stream_lines[s][i + 1].y, stream_lines[s][i + 1].z);
            Vec3.transformMat4(start, start1, gridToCartn);
            Vec3.transformMat4(end, end1, gridToCartn);
            builder.addVec(start, end, s);
        }
    }
    const pt = builder.getLines();
    pt.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return pt;
}

function getLoci(volume: Volume, props: VolumeGradientProps) {
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return Volume.Isosurface.Loci(volume, props.isoValue, instances);
}

function getGradientLoci(pickingId: PickingId, volume: Volume, key: number, props: VolumeGradientProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;

    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        if (granularity === 'volume') {
            return Volume.Loci(volume, instances);
        } else if (granularity === 'object') {
            return Volume.Isosurface.Loci(volume, props.isoValue, instances);
        } else {
            const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
            return Volume.Cell.Loci(volume, [{ indices, instances }]);
        }
    }
    return EmptyLoci;
}

function eachGradient(loci: Loci, volume: Volume, key: number, props: VolumeGradientProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, { isoValue: props.isoValue }, apply);
}

const GradientVisuals = {
    'lines': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumeLinesParams>) => VolumeRepresentation('Gradient lines', ctx, getParams, VolumeLinesVisual, getLoci),
    'cylinders': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumeCylindersParams>) => VolumeRepresentation('Gradient cylinders', ctx, getParams, VolumeCylindersImpostorVisual, getLoci),
};


export const GradientParams = {
    ...VolumeLinesParams,
    ...VolumeCylindersParams,
    visuals: PD.MultiSelect(['lines'], PD.objectToOptions(GradientVisuals)),
};
export type GradientParams = typeof GradientParams;

export function getGradientParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(GradientParams);
    p.isoValue = Volume.createIsoValueParam(Volume.IsoValue.relative(2), volume.grid.stats);
    return p;
}

export type GradientRepresentation = VolumeRepresentation<GradientParams>
export function GradientRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, GradientParams>): GradientRepresentation {
    return Representation.createMulti('Gradient', ctx, getParams, Representation.StateBuilder, GradientVisuals as unknown as Representation.Def<Volume, GradientParams>);
}

export const GradientRepresentationProvider = VolumeRepresentationProvider({
    name: 'gradient',
    label: 'Gradient',
    description: 'Displays gradient vector lines.',
    factory: GradientRepresentation,
    getParams: getGradientParams,
    defaultValues: PD.getDefaultValues(GradientParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume)
});

interface StreamlinePoint {
    x: number;
    y: number;
    z: number;
    value: number;
}

function traceOneDirection(
    out: StreamlinePoint[],
    space: Tensor.Space, data: ArrayLike<number>, seed: Vec3,
    maxSteps: number, ds: number, eps: number,
    hx: number, hy: number, hz: number, dirSign: 1 | -1, skipFirst: boolean,
    traceState: {
        p: Vec3,
        g: Vec3, g2: Vec3, g3: Vec3, g4: Vec3,
        v1: Vec3, v2: Vec3, v3: Vec3, v4: Vec3,
        k1: Vec3, k2: Vec3, k3: Vec3, k4: Vec3,
        dp: Vec3, t1: Vec3, t2: Vec3
    }
): number {

    // const line: StreamlinePoint[] = [];
    const [nx, ny, nz] = space.dimensions as Vec3;
    const p = traceState.p;
    Vec3.copy(p, seed);
    const inv6 = 1 / 6;
    let written = 0;

    for (let step = 0; step < maxSteps; step++) {
        if (p[0] < 1 || p[0] > (nx as number) - 2 ||
            p[1] < 1 || p[1] > (ny as number) - 2 ||
            p[2] < 1 || p[2] > (nz as number) - 2) break;

        getInterpolatedGradient(traceState.g, space, data, p, hx, hy, hz);
        const m = Vec3.magnitude(traceState.g);
        if (!(m > eps)) break;

        // direction only, signed
        Vec3.scale(traceState.v1, traceState.g, dirSign / m);

        // RK4 with constant arc-length ds
        Vec3.scale(traceState.k1, traceState.v1, ds);

        Vec3.add(traceState.t1, p, Vec3.scale(traceState.t2, traceState.k1, 0.5));
        getInterpolatedGradient(traceState.g2, space, data, traceState.t1, hx, hy, hz);
        Vec3.scale(traceState.v2, traceState.g2, dirSign / Math.max(Vec3.magnitude(traceState.g2), eps));
        Vec3.scale(traceState.k2, traceState.v2, ds);

        Vec3.add(traceState.t1, p, Vec3.scale(traceState.t2, traceState.k2, 0.5));
        getInterpolatedGradient(traceState.g3, space, data, traceState.t1, hx, hy, hz);
        Vec3.scale(traceState.v3, traceState.g3, dirSign / Math.max(Vec3.magnitude(traceState.g3), eps));
        Vec3.scale(traceState.k3, traceState.v3, ds);

        Vec3.add(traceState.t1, p, traceState.k3);
        getInterpolatedGradient(traceState.g4, space, data, traceState.t1, hx, hy, hz);
        Vec3.scale(traceState.v4, traceState.g4, dirSign / Math.max(Vec3.magnitude(traceState.g4), eps));
        Vec3.scale(traceState.k4, traceState.v4, ds);

        // dp = (k1 + 2*k2 + 2*k3 + k4) / 6, with zero allocs
        Vec3.copy(traceState.dp, traceState.k1);

        Vec3.scale(traceState.t1, traceState.k2, 2);
        Vec3.add(traceState.dp, traceState.dp, traceState.t1);

        Vec3.scale(traceState.t2, traceState.k3, 2);
        Vec3.add(traceState.dp, traceState.dp, traceState.t2);

        Vec3.add(traceState.dp, traceState.dp, traceState.k4);
        Vec3.scale(traceState.dp, traceState.dp, inv6);

        if (!(skipFirst && step === 0)) {
            const value = getInterpolatedValue(space, data, p);
            out.push({ x: p[0], y: p[1], z: p[2], value });
            written++;
        }
        Vec3.add(p, p, traceState.dp);
    }
    return written;
}

function reverseSegmentInPlace<T>(arr: T[], start: number, end: number) {
    for (let i = start, j = end; i < j; i++, j--) {
        const tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    }
}

function traceStreamlineBothDirs(
    space: Tensor.Space, data: ArrayLike<number>, seed: Vec3,
    maxSteps: number, ds: number, eps: number,
    hx: number, hy: number, hz: number
): StreamlinePoint[] {
    const line: StreamlinePoint[] = [];
    // one shared scratch block to avoid per-call allocations
    const traceState = {
        p: Vec3(),
        g: Vec3(), g2: Vec3(), g3: Vec3(), g4: Vec3(),
        v1: Vec3(), v2: Vec3(), v3: Vec3(), v4: Vec3(),
        k1: Vec3(), k2: Vec3(), k3: Vec3(), k4: Vec3(),
        dp: Vec3(), t1: Vec3(), t2: Vec3()
    };
    // 1) trace BACKWARD, including the seed
    const nBack = traceOneDirection(line,
        space, data, seed, maxSteps, ds, eps, hx, hy, hz, -1,
        /* skipFirst */ false, traceState
    );

    // 2) reverse that segment IN PLACE so it's ordered from farthest->seed
    if (nBack > 1) reverseSegmentInPlace(line, 0, nBack - 1);

    // 3) trace FORWARD, skipping the seed to avoid duplication; append directly
    traceOneDirection(line,
        space, data, seed, maxSteps, ds, eps, hx, hy, hz, +1,
        /* skipFirst */ true, traceState
    );
    return line;
}

function getInterpolatedGradient(out: Vec3, space: Tensor.Space, data: ArrayLike<number>,
                                 pos: Vec3, hx: number, hy: number, hz: number): Vec3 {
    const x = Math.floor(pos[0]), y = Math.floor(pos[1]), z = Math.floor(pos[2]);
    const fx = pos[0] - x, fy = pos[1] - y, fz = pos[2] - z;

    // reset accumulator
    out[0] = 0; out[1] = 0; out[2] = 0;

    let i = 0;
    const grads: Vec3[] = [
        Vec3(), Vec3(), Vec3(), Vec3(),
        Vec3(), Vec3(), Vec3(), Vec3(),
    ];
    for (let dz = 0; dz <= 1; dz++) {
        for (let dy = 0; dy <= 1; dy++) {
            for (let dx = 0; dx <= 1; dx++) {
                calculateGradient(
                    grads[i], // write directly into existing Vec3
                    space, data, x + dx, y + dy, z + dz, hx, hy, hz
                );
                i++;
            }
        }
    }

    for (let i = 0; i < 8; i++) {
        const wx = (i & 1) ? fx : (1 - fx);
        const wy = (i & 2) ? fy : (1 - fy);
        const wz = (i & 4) ? fz : (1 - fz);
        const w = wx * wy * wz;
        Vec3.scaleAndAdd(out, out, grads[i], w);
    }
    return out;
}

function getInterpolatedValue(space: Tensor.Space, data: ArrayLike<number>, pos: Vec3): number {
    // Trilinear interpolation of scalar value
    const x = Math.floor(pos[0]);
    const y = Math.floor(pos[1]);
    const z = Math.floor(pos[2]);

    const fx = pos[0] - x;
    const fy = pos[1] - y;
    const fz = pos[2] - z;

    let result = 0;
    for (let dz = 0; dz <= 1; dz++) {
        for (let dy = 0; dy <= 1; dy++) {
            for (let dx = 0; dx <= 1; dx++) {
                const weight =
                    (dx ? fx : (1 - fx)) *
                    (dy ? fy : (1 - fy)) *
                    (dz ? fz : (1 - fz));
                result += data[space.dataOffset(x + dx, y + dy, z + dz)] * weight;
            }
        }
    }

    return result;
}

function calculateGradient(out: Vec3, space: Tensor.Space, data: ArrayLike<number>, x: number, y: number, z: number, hx: number, hy: number, hz: number): Vec3 {
    // clamp indices inside [1..N-2] to avoid border reads
    const xi = Math.max(1, Math.min(x, (space.dimensions[0] as number) - 2));
    const yi = Math.max(1, Math.min(y, (space.dimensions[1] as number) - 2));
    const zi = Math.max(1, Math.min(z, (space.dimensions[2] as number) - 2));

    const gx = (data[space.dataOffset(xi + 1, yi, zi)] - data[space.dataOffset(xi - 1, yi, zi)]) / (2 * hx);
    const gy = (data[space.dataOffset(xi, yi + 1, zi)] - data[space.dataOffset(xi, yi - 1, zi)]) / (2 * hy);
    const gz = (data[space.dataOffset(xi, yi, zi + 1)] - data[space.dataOffset(xi, yi, zi - 1)]) / (2 * hz);

    // Electric field E = -∇φ
    out[0] = -gx;
    out[1] = -gy;
    out[2] = -gz;
    return out;
}

// Simple RK4 streamlines implementation
function getStreamLineSimple(volume: Volume, props: VolumeLinesProps, out: StreamlinePoint[][]) {
    const { cells: { space, data } } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const seedDensity = props.seedDensity ?? 8;
    const [nx, ny, nz] = space.dimensions as Vec3;
    const { hx, hy, hz } = getCellSize(gridToCartn);
    // Generate seed points
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));
    const maxSteps = props.maxSteps ?? 2000; // more steps
    const ds = (props.stepSize ?? 0.35); // in index units; 0.2–0.5 works well
    const eps = props.minSpeed ?? 1e-6; // MUCH smaller for APBS

    const pos = Vec3();
    const g = Vec3();
    const xStart = 1, xEnd = nx - 2;
    const yStart = 1, yEnd = ny - 2;
    const zStart = 1, zEnd = nz - 2;
    // --- build seed list ---
    const seedPoints: number[] = [];
    for (let z = zStart; z <= zEnd; z += seedStep) {
        for (let y = yStart; y <= yEnd; y += seedStep) {
            for (let x = xStart; x <= xEnd; x += seedStep) {
                seedPoints.push(x);
                seedPoints.push(y);
                seedPoints.push(z);
            }
        }
    }
    // --- shuffle (Fisher–Yates) ---
    for (let i = seedPoints.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        const tmp = seedPoints[i];
        seedPoints[i] = seedPoints[j];
        seedPoints[j] = tmp;
    }
    for (let i = 0; i < seedPoints.length; i += 3) {
        Vec3.set(pos, seedPoints[i], seedPoints[i + 1], seedPoints[i + 2]);

        // filtering
        const phi = getInterpolatedValue(space, data, pos);
        getInterpolatedGradient(g, space, data, pos, hx, hy, hz);
        const mag = Vec3.magnitude(g);

        // thresholds to adjust according to your APBS maps
        if ((phi < props.minLevel || phi > props.maxLevel) || mag < 1e-6) continue;

        // if ok, we trace the streamline
        const streamline = traceStreamlineBothDirs(
            space, data, pos,
            maxSteps, ds, eps, hx, hy, hz
        );

        if (streamline.length < 2) continue;
        out.push(streamline);
    }
}

function getCellSize(gridToCartn: Mat4) {
    const b = Mat4.extractBasis(gridToCartn);
    return {
        hx: Vec3.magnitude(b.x),
        hy: Vec3.magnitude(b.y),
        hz: Vec3.magnitude(b.z)
    };
}

// tiny LCG for reproducible shuffles (like std::mt19937 seeded with range_size)
function seededRand(seed: number) {
  let s = seed >>> 0;
  return () => {
    // xorshift32
    s ^= s << 13; s ^= s >>> 17; s ^= s << 5;
    return ((s >>> 0) / 4294967296);
  };
}


// Interpolate gradient & position in *index space* (cell units)
function interpolateValue(space: Tensor.Space, data: ArrayLike<number>, i: number, j: number, k: number, fx: number, fy: number, fz: number) {
  let acc = 0;
  for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
    const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
    acc += data[space.dataOffset(i + dx, j + dy, k + dz)] * w;
  }
  return acc;
}

// precompute gradient field in *index units* (like IsofieldComputeGradients)
// NOTE: this is centered diff; PyMOL stores gradient in separate field
function buildGradientField(space: Tensor.Space, data: ArrayLike<number>) {
  const [nx, ny, nz] = space.dimensions as Vec3;
  const gx = new Float32Array(nx * ny * nz);
  const gy = new Float32Array(nx * ny * nz);
  const gz = new Float32Array(nx * ny * nz);
  const idx = (x: number, y: number, z: number) => space.dataOffset(x, y, z);

  for (let x = 0; x < nx; x++) {
    for (let y = 0; y < ny; y++) {
      for (let z = 0; z < nz; z++) {
        const xi0 = Math.max(0, x - 1), xi1 = Math.min(nx - 1, x + 1);
        const yi0 = Math.max(0, y - 1), yi1 = Math.min(ny - 1, y + 1);
        const zi0 = Math.max(0, z - 1), zi1 = Math.min(nz - 1, z + 1);
        gx[idx(x, y, z)] = (data[idx(xi1, y, z)] - data[idx(xi0, y, z)]) * 0.5;
        gy[idx(x, y, z)] = (data[idx(x, yi1, z)] - data[idx(x, yi0, z)]) * 0.5;
        gz[idx(x, y, z)] = (data[idx(x, y, zi1)] - data[idx(x, y, zi0)]) * 0.5;
      }
    }
  }
  return { gx, gy, gz };
}

function interpolateGradient(space: Tensor.Space, grad: {gx: Float32Array, gy: Float32Array, gz: Float32Array},
                             i: number, j: number, k: number, fx: number, fy: number, fz: number, out: Vec3) {
  const { gx, gy, gz } = grad;
  let vx=0, vy=0, vz=0;
  for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
    const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
    const off = space.dataOffset(i + dx, j + dy, k + dz);
    vx += gx[off] * w;
    vy += gy[off] * w;
    vz += gz[off] * w;
  }
  out[0] = vx; out[1] = vy; out[2] = vz;
  return out;
}

function getStreamLineAdvanced(volume: Volume, props: VolumeLinesProps, out: StreamlinePoint[][]) {
    // PYMOL algo
    const { cells: { space, data } } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);

    const [nx, ny, nz] = space.dimensions as Vec3;
    const range: [number,number,number,number,number,number] = [0,0,0, nx,ny,nz];
    const seed = (nx * ny * nz) | 0;
    const i0 = Math.max(0, range[0]|0), j0 = Math.max(0, range[1]|0), k0 = Math.max(0, range[2]|0);
    const i1 = Math.min(nx, range[3]|0), j1 = Math.min(ny, range[4]|0), k1 = Math.min(nz, range[5]|0);

    const dx = Vec3.magnitude(Mat4.extractBasis(gridToCartn).x);
    const dy = Vec3.magnitude(Mat4.extractBasis(gridToCartn).y);
    const dz = Vec3.magnitude(Mat4.extractBasis(gridToCartn).z);
    const avgCell = (dx + dy + dz) / 3;

    // scale params into *index (cell) units* like C++
    const stepSize = Math.max(props.stepSize / avgCell, 0.01);
    const maxWalk = 10; // props.maxSteps / avgCell; // 20
    const minWalk = 2;
    const minSlope = Math.max(props.minSpeed * avgCell, 1e-5);
    const minDot = 0.0;
    let symmetry = 0.0;
    const symmetryFlag = symmetry !== 0;
    if (symmetry > 1) symmetry = 1 / symmetry;

    // gradients (index units, like PyMOL's IsofieldComputeGradients)
    const grad = buildGradientField(space, data);
    const pointsToWorld = (i:number,j:number,k:number, fx:number,fy:number,fz:number, out:Vec3) => {
        // trilinear in index space → then transform to world
        Vec3.set(out, i + fx, j + fy, k + fz);
        // return Vec3.transformMat4(out, out, gridToCartn);
    };

    // randomized order of all cells in range
    const rx = i1 - i0, ry = j1 - j0, rz = k1 - k0;
    const rangeSize = rx * ry * rz;
    const order = new Int32Array(rangeSize * 3);
    {
        let t = 0;
        for (let k = k0; k < k1; k++) for (let j = j0; j < j1; j++) for (let i = i0; i < i1; i++) {
        order[t++] = i; order[t++] = j; order[t++] = k;
        }
        const rand = seededRand(seed || rangeSize);
        for (let a = 0; a < rangeSize; a++) {
        const p = (Math.floor(rand() * rangeSize) * 3) | 0;
        const q = (Math.floor(rand() * rangeSize) * 3) | 0;
        const t0 = order[p], t1 = order[p+1], t2 = order[p+2];
        order[p] = order[q]; order[p+1] = order[q+1]; order[p+2] = order[q+2];
        order[q] = t0; order[q+1] = t1; order[q+2] = t2;
        }
    }

    // flag array for oblation
    const flag = new Uint8Array(rx * ry * rz); // 0/1
    const stride0 = 1, stride1 = rx, stride2 = rx * ry;
    const flagIndex = (i:number,j:number,k:number) => ((i - i0) * stride0) + ((j - j0) * stride1) + ((k - k0) * stride2);

    // outputs like PyMOL
    // const positions: number[] = [];
    const num: number[] = [];
    // const v: number[] = [];
    let nLine = 0;
    let nSeg = 0;
    num[nSeg] = nLine;

    // scratch
    const g = Vec3();
    const prevG = Vec3();
    out.push([]); // first segment
    // walk both directions from each start cell
    for (let a = 0; a < rangeSize; a++) {
        let walkRemain = maxWalk;
        const abortNLine = nLine;
        const abortNSeg = nSeg;
        let symmetryMin = Number.POSITIVE_INFINITY;
        let symmetryMax = Number.NEGATIVE_INFINITY;

        // track active cells visited for oblation
        const activeCells: number[] = [];

        const doPass = (sign: 1|-1) => {
        let havePrev = false;
        let li = order[a*3], lj = order[a*3+1], lk = order[a*3+2];
        let fx = 0, fy = 0, fz = 0; // start at cell corner (like C++)
        let nVert = 0;

        for (;;) {
            // normalize fract and adjust locus
            // keep fx,fy,fz in [0,1), move cell indices accordingly
            // if out of range, break
            const normComp = (f:number, l:number, lo:number, hi:number) => {
            while (f < 0) { f += 1; l--; }
            while (f >= 1) { f -= 1; l++; }
            if (l < lo || l > hi) return { f, l, done: true };
            return { f, l, done: false };
            };
            let r = normComp(fx, li, i0, i1-2); fx = r.f; li = r.l; if (r.done) break;
            r = normComp(fy, lj, j0, j1-2); fy = r.f; lj = r.l; if (r.done) break;
            r = normComp(fz, lk, k0, k1-2); fz = r.f; lk = r.l; if (r.done) break;

            // stop if flagged cell
            if (flag[flagIndex(li, lj, lk)]) break;

            // level test
            const level = interpolateValue(space, data, li, lj, lk, fx, fy, fz);
            if (level < props.minLevel || level > props.maxLevel) break;
            if (symmetryFlag) {
                if (level < symmetryMin) symmetryMin = level;
                if (level > symmetryMax) symmetryMax = level;
            }

            // gradient
            interpolateGradient(space, grad, li, lj, lk, fx, fy, fz, g);
            const gm = Vec3.magnitude(g);
            if (gm < minSlope) break;

            // add vertex (world coords, like PyMOL)
            {
            const pw = Vec3();
            pointsToWorld(li, lj, lk, fx, fy, fz, pw);
            // positions.push(pw[0], pw[1], pw[2]);
            // v.push(level);
            out[nSeg].push({ x: pw[0], y: pw[1], z: pw[2], value: level });
            nLine++; nVert++;
            }

            // record active cell when it changes
            if (!havePrev || havePrev && activeCells.length >= 3) {
            const n = activeCells.length;
            if (n < 3 || activeCells[n-3] !== li || activeCells[n-2] !== lj || activeCells[n-1] !== lk) {
                activeCells.push(li, lj, lk);
            }
            }

            // normalize and turn into step
            Vec3.scale(g, g, 1 / gm);
            if (havePrev) {
            const dp = g[0]*prevG[0] + g[1]*prevG[1] + g[2]*prevG[2];
            if (dp < minDot) break;
            }
            prevG[0] = g[0]; prevG[1] = g[1]; prevG[2] = g[2];

            // forward/backward
            const s = sign * stepSize;
            fx += g[0] * s;
            fy += g[1] * s;
            fz += g[2] * s;

            walkRemain -= stepSize;
            if (walkRemain < 0) break;
            havePrev = true;
        }

        // segment bookkeeping
        if (nVert < 2) {
            if (nVert) nLine = num[nSeg]; // revert last lone point
        } else if (nLine !== num[nSeg]) {
            num[nSeg] = nLine - num[nSeg];
            nSeg++;
            num[nSeg] = nLine;
            out.push([]);
        }
        };

        // pass 0: down gradient (+)
        doPass(+1);
        // pass 1: up gradient (-)
        doPass(-1);

        // symmetry/min length checks
        let abort = false;
        if (symmetryFlag) {
        if (symmetryMax * symmetryMin >= 0) abort = true;
        else {
            let ratio = Math.abs(symmetryMax) / Math.abs(symmetryMin);
            if (ratio > 1) ratio = 1 / ratio;
            if (ratio < symmetry) abort = true;
        }
        }
        if ((maxWalk / avgCell - walkRemain) < minWalk) abort = true;

        if (abort) {
            // rollback
            nSeg = abortNSeg;
            nLine = abortNLine;
            num[nSeg] = nLine;
        } else {
            // oblation (flag spherical neighborhood with radius = spacing)
            const R = OBLATION_SPACING|0;
            const R2 = R*R;
            for (let t = 0; t < activeCells.length; t += 3) {
                const ii = activeCells[t], jj = activeCells[t+1], kk = activeCells[t+2];
                for (let k = Math.max(k0, kk-R); k <= Math.min(k1-1, kk+R); k++) {
                    const dz2 = (kk-k)*(kk-k);
                    for (let j = Math.max(j0, jj-R); j <= Math.min(j1-1, jj+R); j++) {
                        const dy2 = (jj-j)*(jj-j) + dz2;
                        if (dy2 > R2) continue;
                        for (let i = Math.max(i0, ii-R); i <= Math.min(i1-1, ii+R); i++) {
                            const d2 = (ii-i)*(ii-i) + dy2;
                            if (d2 <= R2) flag[flagIndex(i,j,k)] = 1;
                        }
                    }
                }
            }
        }
    }

    // terminate like PyMOL
    num[nSeg] = 0;
    // return { positions, num, v };
}
