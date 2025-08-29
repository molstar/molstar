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

export const VolumeGradientParams = {
    isoValue: Volume.IsoValueParam,
    seedDensity: PD.Numeric(20, { min: 1, max: 100, step: 1 }, { description: 'Seeds per dimension for gradient lines.' }),
    maxSteps: PD.Numeric(1000, { min: 1, max: 2000, step: 1 }, { description: 'Maximum number of steps for gradient lines.' }),
    stepSize: PD.Numeric(0.35, { min: 0.01, max: 10, step: 0.01 }, { description: 'Step size for gradient lines.' }),
    minSpeed: PD.Numeric(0.001, { min: 0, max: 1, step: 1e-6 }, { description: 'Minimum speed for gradient lines.' }),
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
                newProps.minSpeed !== currentProps.minSpeed
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
                newProps.radius !== currentProps.radius
            );
        },
        geometryUtils: Cylinders.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeCylindersParams>, webgl?: WebGLContext) => {
            return !webgl;
        }
    }, materialId);
}

export function createVolumeCylindersImpostor(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeCylindersProps, geometry?: Cylinders): Cylinders {
    const { cells: { space, data } } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);

    const [nx, ny, nz] = space.dimensions as Vec3;

    const { hx, hy, hz } = getCellSize(gridToCartn);

    // seeds
    const seedDensity = props.seedDensity ?? 8;
    const maxSteps = props.maxSteps ?? 2000; // more steps
    const ds = (props.stepSize ?? 0.35); // in index units; 0.2–0.5 works well
    const eps = props.minSpeed ?? 1e-6; // MUCH smaller for APBS


    // Generate seed points
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));

    // const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const builder = CylindersBuilder.create(count, Math.ceil(count / 2), geometry);
    // initialize the Vec3
    const pos = Vec3();
    const start = Vec3();
    const end = Vec3();
    const g = Vec3();
    for (let z = seedStep; z < nz - seedStep; z += seedStep) {
        for (let y = seedStep; y < ny - seedStep; y += seedStep) {
            for (let x = seedStep; x < nx - seedStep; x += seedStep) {

                Vec3.set(pos, x, y, z);

                // filtering
                const phi = getInterpolatedValue(space, data, pos);
                getInterpolatedGradient(g, space, data, pos, hx, hy, hz);
                const mag = Vec3.magnitude(g);

                // thresholds to adjust according to your APBS maps
                if (Math.abs(phi) > 0.2 || mag < 1e-6) continue;

                // if ok, we trace the streamline
                const streamline = traceStreamlineBothDirs(
                    space, data, pos,
                    maxSteps, ds, eps, hx, hy, hz
                );

                if (streamline.length < 2) continue;
                const cellIdx = space.dataOffset(x, y, z);

                for (let i = 0; i < streamline.length - 1; i++) {
                    Vec3.set(start, streamline[i].x, streamline[i].y, streamline[i].z);
                    Vec3.set(end, streamline[i + 1].x, streamline[i + 1].y, streamline[i + 1].z);
                    Vec3.transformMat4(start, start, gridToCartn);
                    Vec3.transformMat4(end, end, gridToCartn);
                    builder.add(start[0], start[1], start[2],
                                end[0], end[1], end[2],
                                props.radius, true, true, 2, cellIdx);
                }
            }
        }
    }
    const pt = builder.getCylinders();
    pt.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return pt;
}


export function createVolumeLinesMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeLinesProps, points?: Lines): Lines {
    const { cells: { space, data } } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);

    const [nx, ny, nz] = space.dimensions as Vec3;

    const { hx, hy, hz } = getCellSize(gridToCartn);

    // seeds
    const seedDensity = props.seedDensity ?? 8;
    const maxSteps = props.maxSteps ?? 2000; // more steps
    const ds = (props.stepSize ?? 0.35); // in index units; 0.2–0.5 works well
    const eps = props.minSpeed ?? 1e-6; // MUCH smaller for APBS

    // Generate seed points
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));

    // const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const builder = LinesBuilder.create(count, Math.ceil(count / 2), points);

    // Precompute basis vectors and largest cell axis length
    // const basis = getBasis(gridToCartn);
    const pos = Vec3();
    const start = Vec3();
    const end = Vec3();
    const start1 = Vec3();
    const end1 = Vec3();
    const g = Vec3();
    for (let z = seedStep; z < nz - seedStep; z += seedStep) {
        for (let y = seedStep; y < ny - seedStep; y += seedStep) {
            for (let x = seedStep; x < nx - seedStep; x += seedStep) {

                Vec3.set(pos, x, y, z);

                // filtering
                const phi = getInterpolatedValue(space, data, pos);
                getInterpolatedGradient(g, space, data, pos, hx, hy, hz);
                const mag = Vec3.magnitude(g);

                // thresholds to adjust according to your APBS maps
                if (Math.abs(phi) > 0.2 || mag < 1e-6) continue;

                // if ok, we trace the streamline
                const streamline = traceStreamlineBothDirs(
                    space, data, pos,
                    maxSteps, ds, eps, hx, hy, hz
                );

                if (streamline.length < 2) continue;
                const cellIdx = space.dataOffset(x, y, z);

                for (let i = 0; i < streamline.length - 1; i++) {
                    Vec3.set(start1, streamline[i].x, streamline[i].y, streamline[i].z);
                    Vec3.set(end1, streamline[i + 1].x, streamline[i + 1].y, streamline[i + 1].z);
                    Vec3.transformMat4(start, start1, gridToCartn);
                    Vec3.transformMat4(end, end1, gridToCartn);
                    builder.addVec(start, end, cellIdx);
                }
            }
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


function getCellSize(gridToCartn: Mat4) {
    const b = Mat4.extractBasis(gridToCartn);
    return {
        hx: Vec3.magnitude(b.x),
        hy: Vec3.magnitude(b.y),
        hz: Vec3.magnitude(b.z)
    };
}
