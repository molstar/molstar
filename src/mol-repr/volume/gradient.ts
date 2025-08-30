/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @maintainer ChatGPT5
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

// ---------------------------------------------------------------------------------------
// Shared params
// ---------------------------------------------------------------------------------------

export const VolumeGradientParams = {
    isoValue: Volume.IsoValueParam,
    seedDensity: PD.Numeric(10, { min: 1, max: 100, step: 1 }, { description: 'Seeds per dimension for streamlines.' }),
    maxSteps: PD.Numeric(1000, { min: 1, max: 2000, step: 1 }, { description: 'Maximum number of steps for streamlines.' }),
    stepSize: PD.Numeric(0.35, { min: 0.01, max: 10, step: 0.01 }, { description: 'Step size (Å) for streamlines.' }),
    minSpeed: PD.Numeric(0.001, { min: 0, max: 1, step: 1e-6 }, { description: 'Minimum |E| (in 1/Å) to continue tracing.' }),
    minLevel: PD.Numeric(-5.0, { min: -10, max: 10, step: 0.1 }, { description: 'Minimum map value (φ) allowed when tracing.' }),
    maxLevel: PD.Numeric(5.0, { min: -10, max: 10, step: 0.1 }, { description: 'Maximum map value (φ) allowed when tracing.' }),
    algorithm: PD.Select('simple', PD.arrayToOptions(['simple', 'advanced'] as const), { description: 'Streamline algorithm to use.' }),
    mid_interpolation: PD.Boolean(true, { description: 'Use RK2 (midpoint) instead of RK4 for the simple tracer.' }),
    writeStride: PD.Numeric(4, { min: 1, max: 50, step: 1 }, { description: 'Sample every N integration steps.' }),
    geomStride: PD.Numeric(1, { min: 1, max: 10, step: 1 }, { description: 'Emit every Nth segment to geometry.' }),
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

export const VolumeCylindersParams = {
    ...Cylinders.Params,
    ...VolumeGradientParams,
    radius: PD.Numeric(1.0, { min: 0.1, max: 5, step: 0.1 }, { description: 'Cylinder radius (Å).' }),
};
export type VolumeCylindersParams = typeof VolumeCylindersParams
export type VolumeCylindersProps = PD.Values<VolumeCylindersParams>

// ---------------------------------------------------------------------------------------
// Visuals
// ---------------------------------------------------------------------------------------

export function VolumeLinesVisual(materialId: number): VolumeVisual<VolumeLinesParams> {
    return VolumeVisual<Lines, VolumeLinesParams>({
        defaultProps: PD.getDefaultValues(VolumeLinesParams),
        createGeometry: createVolumeLinesGeometry,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps, currentProps, _newTheme, _currentTheme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel ||
                newProps.writeStride !== currentProps.writeStride ||
                newProps.geomStride !== currentProps.geomStride
            );
        },
        geometryUtils: Lines.Utils,
        mustRecreate: (_: VolumeKey, __: PD.Values<VolumeLinesParams>, webgl?: WebGLContext) => !!webgl,
    }, materialId);
}

export function VolumeCylindersImpostorVisual(materialId: number): VolumeVisual<VolumeCylindersParams> {
    return VolumeVisual<Cylinders, VolumeCylindersParams>({
        defaultProps: PD.getDefaultValues(VolumeCylindersParams),
        createGeometry: createVolumeCylindersGeometry,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps, currentProps, _newTheme, _currentTheme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.radius !== currentProps.radius ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel ||
                newProps.writeStride !== currentProps.writeStride ||
                newProps.geomStride !== currentProps.geomStride
            );
        },
        geometryUtils: Cylinders.Utils,
        mustRecreate: (_: VolumeKey, __: PD.Values<VolumeCylindersParams>, webgl?: WebGLContext) => !webgl,
    }, materialId);
}

// ---------------------------------------------------------------------------------------
// Geometry builders (shared streamline collection)
// ---------------------------------------------------------------------------------------

function createVolumeLinesGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumeLinesProps, lines?: Lines): Lines {
    const streamLines = collectStreamlines(volume, props);
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);

    const segCount = countSegments(streamLines, props.geomStride);
    const builder = LinesBuilder.create(segCount, Math.ceil(segCount / 2), lines);

    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => builder.addVec(a, b, id));

    const g = builder.getLines();
    g.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return g;
}

function createVolumeCylindersGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumeCylindersProps, cylinders?: Cylinders): Cylinders {
    const streamLines = collectStreamlines(volume, props);
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);

    const segCount = countSegments(streamLines, props.geomStride);
    const builder = CylindersBuilder.create(segCount, Math.ceil(segCount / 2), cylinders);

    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => {
        builder.add(a[0], a[1], a[2], b[0], b[1], b[2], props.radius, true, true, 2, id);
    });

    const g = builder.getCylinders();
    g.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return g;
}

function countSegments(lines: StreamlinePoint[][], stride: number) {
    let c = 0;
    for (let s = 0; s < lines.length; s++) c += Math.max(0, Math.floor((lines[s].length - 1) / Math.max(1, stride)));
    return c;
}

function writeSegmentsToBuilder(
    streamLines: StreamlinePoint[][],
    gridToCartn: Mat4,
    geomStride: number,
    emit: (a: Vec3, b: Vec3, id: number) => void
) {
    const a = Vec3(), b = Vec3(), ai = Vec3(), bi = Vec3();
    for (let s = 0; s < streamLines.length; s++) {
        const L = streamLines[s];
        for (let i = 0, n = L.length - 1; i < n; i++) {
            if (i % Math.max(1, geomStride) !== 0) continue;
            ai[0] = L[i].x; ai[1] = L[i].y; ai[2] = L[i].z;
            bi[0] = L[i+1].x; bi[1] = L[i+1].y; bi[2] = L[i+1].z;
            Vec3.transformMat4(a, ai, gridToCartn);
            Vec3.transformMat4(b, bi, gridToCartn);
            emit(a, b, s);
        }
    }
}

// ---------------------------------------------------------------------------------------
// Loci plumbing
// ---------------------------------------------------------------------------------------

function getLoci(volume: Volume, props: VolumeGradientProps) {
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return Volume.Loci(volume, instances); //
    // return Volume.Isosurface.Loci(volume, props.isoValue, instances);
}

function getGradientLoci(pickingId: PickingId, volume: Volume, _key: number, props: VolumeGradientProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;
    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        if (granularity === 'volume') return Volume.Loci(volume, instances);
        // if (granularity === 'object') return Volume.Isosurface.Loci(volume, props.isoValue, instances);
        const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
        return Volume.Cell.Loci(volume, [{ indices, instances }]);
    }
    return EmptyLoci;
}

function eachGradient(loci: Loci, volume: Volume, _key: number, props: VolumeGradientProps, apply: (interval: Interval) => boolean) {
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

// ---------------------------------------------------------------------------------------
// Streamline integrators
// ---------------------------------------------------------------------------------------

interface StreamlinePoint { x: number; y: number; z: number; value: number; }

function collectStreamlines(volume: Volume, props: VolumeLinesProps | VolumeCylindersProps): StreamlinePoint[][] {
    const out: StreamlinePoint[][] = [];
    if (props.algorithm === 'advanced') getStreamLineAdvanced(volume, props, out);
    else getStreamLineSimple(volume, props, out);
    return out;
}

function getCellSize(gridToCartn: Mat4) {
    const b = Mat4.extractBasis(gridToCartn);
    return { hx: Vec3.magnitude(b.x), hy: Vec3.magnitude(b.y), hz: Vec3.magnitude(b.z) };
}

// ---------------- Simple tracer (RK2/RK4) in index space ----------------

function traceOneDirection(
    out: StreamlinePoint[],
    space: Tensor.Space, data: ArrayLike<number>, seed: Vec3,
    maxSteps: number, dsWorld: number, eps: number,
    hx: number, hy: number, hz: number, dirSign: 1 | -1, skipFirst: boolean,
    S: { p: Vec3, g: Vec3, g2: Vec3, g3: Vec3, g4: Vec3, v1: Vec3, v2: Vec3, v3: Vec3, v4: Vec3, k1: Vec3, k2: Vec3, k3: Vec3, k4: Vec3, dp: Vec3, t1: Vec3, t2: Vec3, t3: Vec3 },
    writeStride = 2,
    order: 2 | 4
): number {
    const [nx, ny, nz] = space.dimensions as Vec3;
    const p = S.p; Vec3.copy(p, seed);

    // convert a world step to index steps per-axis
    const sIx = dirSign * (dsWorld / hx), sIy = dirSign * (dsWorld / hy), sIz = dirSign * (dsWorld / hz);

    let written = 0;
    for (let step = 0; step < maxSteps; step++) {
        if (p[0] < 1 || p[0] > (nx as number) - 2 || p[1] < 1 || p[1] > (ny as number) - 2 || p[2] < 1 || p[2] > (nz as number) - 2) break;

        // gW = -∇φ in WORLD units (1/Å)
        gradientAtP_world(S.g, space, data, p, hx, hy, hz, S.t3);
        const m = Vec3.magnitude(S.g);
        if (!(m > eps)) break;
        Vec3.scale(S.v1, S.g, 1 / m); // unit world

        if (order === 2) {
            // RK2 midpoint in world, converted to index steps
            S.t1[0] = p[0] + 0.5 * sIx * S.v1[0];
            S.t1[1] = p[1] + 0.5 * sIy * S.v1[1];
            S.t1[2] = p[2] + 0.5 * sIz * S.v1[2];

            gradientAtP_world(S.g2, space, data, S.t1, hx, hy, hz, S.t3);
            const m2 = Math.max(Vec3.magnitude(S.g2), eps);
            S.v2[0] = S.g2[0] / m2; S.v2[1] = S.g2[1] / m2; S.v2[2] = S.g2[2] / m2;

            S.dp[0] = sIx * S.v2[0]; S.dp[1] = sIy * S.v2[1]; S.dp[2] = sIz * S.v2[2];
        } else {
            // RK4 in world, converted to index
            Vec3.set(S.k1, sIx * S.v1[0], sIy * S.v1[1], sIz * S.v1[2]);

            S.t1[0] = p[0] + 0.5 * S.k1[0]; S.t1[1] = p[1] + 0.5 * S.k1[1]; S.t1[2] = p[2] + 0.5 * S.k1[2];
            gradientAtP_world(S.g2, space, data, S.t1, hx, hy, hz, S.t3);
            const m2 = Math.max(Vec3.magnitude(S.g2), eps);
            S.v2[0] = S.g2[0] / m2; S.v2[1] = S.g2[1] / m2; S.v2[2] = S.g2[2] / m2;
            Vec3.set(S.k2, sIx * S.v2[0], sIy * S.v2[1], sIz * S.v2[2]);

            S.t1[0] = p[0] + 0.5 * S.k2[0]; S.t1[1] = p[1] + 0.5 * S.k2[1]; S.t1[2] = p[2] + 0.5 * S.k2[2];
            gradientAtP_world(S.g3, space, data, S.t1, hx, hy, hz, S.t3);
            const m3 = Math.max(Vec3.magnitude(S.g3), eps);
            S.v3[0] = S.g3[0] / m3; S.v3[1] = S.g3[1] / m3; S.v3[2] = S.g3[2] / m3;
            Vec3.set(S.k3, sIx * S.v3[0], sIy * S.v3[1], sIz * S.v3[2]);

            S.t1[0] = p[0] + S.k3[0]; S.t1[1] = p[1] + S.k3[1]; S.t1[2] = p[2] + S.k3[2];
            gradientAtP_world(S.g4, space, data, S.t1, hx, hy, hz, S.t3);
            const m4 = Math.max(Vec3.magnitude(S.g4), eps);
            S.v4[0] = S.g4[0] / m4; S.v4[1] = S.g4[1] / m4; S.v4[2] = S.g4[2] / m4;
            Vec3.set(S.k4, sIx * S.v4[0], sIy * S.v4[1], sIz * S.v4[2]);

            // dp = (k1 + 2*k2 + 2*k3 + k4) / 6
            Vec3.copy(S.dp, S.k1);
            Vec3.scaleAndAdd(S.dp, S.dp, S.k2, 2);
            Vec3.scaleAndAdd(S.dp, S.dp, S.k3, 2);
            Vec3.add(S.dp, S.dp, S.k4);
            Vec3.scale(S.dp, S.dp, 1 / 6);
        }

        if (!(skipFirst && step === 0) && (step % writeStride === 0)) {
            const value = getInterpolatedValue(space, data, p);
            out.push({ x: p[0], y: p[1], z: p[2], value });
            written++;
        }
        Vec3.add(p, p, S.dp);
    }
    return written;
}

function reverseSegmentInPlace<T>(arr: T[], start: number, end: number) {
    for (let i = start, j = end; i < j; i++, j--) { const tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp; }
}

function traceStreamlineBothDirs(space: Tensor.Space, data: ArrayLike<number>, seed: Vec3, maxSteps: number, dsWorld: number, eps: number, hx: number, hy: number, hz: number, writeStride: number, order: 2 | 4): StreamlinePoint[] {
    const line: StreamlinePoint[] = [];
    const S = { p: Vec3(), g: Vec3(), g2: Vec3(), g3: Vec3(), g4: Vec3(), v1: Vec3(), v2: Vec3(), v3: Vec3(), v4: Vec3(), k1: Vec3(), k2: Vec3(), k3: Vec3(), k4: Vec3(), dp: Vec3(), t1: Vec3(), t2: Vec3(), t3: Vec3() };

    const nBack = traceOneDirection(line, space, data, seed, maxSteps, dsWorld, eps, hx, hy, hz, -1, false, S, writeStride, order);
    if (nBack > 1) reverseSegmentInPlace(line, 0, nBack - 1);
    traceOneDirection(line, space, data, seed, maxSteps, dsWorld, eps, hx, hy, hz, +1, true, S, writeStride, order);
    return line;
}

function getStreamLineSimple(volume: Volume, props: VolumeLinesProps | VolumeCylindersProps, out: StreamlinePoint[][]) {
    const { cells: { space, data } } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const [nx, ny, nz] = space.dimensions as Vec3;
    const { hx, hy, hz } = getCellSize(gridToCartn);

    const seedDensity = props.seedDensity ?? 8;
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));

    const maxSteps = props.maxSteps ?? 2000;
    const dsWorld = props.stepSize ?? 0.35; // Å
    const eps = props.minSpeed ?? 1e-6; // |E| threshold in 1/Å
    const writeStride = Math.max(1, props.writeStride | 0);
    const order: 2 | 4 = props.mid_interpolation ? 2 : 4;

    const pos = Vec3(), g = Vec3();

    // bounds in index space avoiding edges
    const xStart = 1, xEnd = nx - 2;
    const yStart = 1, yEnd = ny - 2;
    const zStart = 1, zEnd = nz - 2;

    // seeds on a coarse lattice, shuffled
    const seeds: number[] = [];
    for (let z = zStart; z <= zEnd; z += seedStep)
        for (let y = yStart; y <= yEnd; y += seedStep)
            for (let x = xStart; x <= xEnd; x += seedStep) { seeds.push(x, y, z); }
    for (let i = seeds.length - 3; i > 0; i -= 3) {
        const j = (Math.floor(Math.random() * (i / 3 + 1)) * 3) | 0;
        const tx = seeds[i], ty = seeds[i+1], tz = seeds[i+2];
        seeds[i] = seeds[j]; seeds[i+1] = seeds[j+1]; seeds[i+2] = seeds[j+2];
        seeds[j] = tx; seeds[j+1] = ty; seeds[j+2] = tz;
    }

    for (let i = 0; i < seeds.length; i += 3) {
        Vec3.set(pos, seeds[i], seeds[i + 1], seeds[i + 2]);

        const phi = getInterpolatedValue(space, data, pos);
        gradientAtP_world(g, space, data, pos, hx, hy, hz);
        const mag = Vec3.magnitude(g);
        if ((phi < props.minLevel || phi > props.maxLevel) || mag < eps) continue;

        const line = traceStreamlineBothDirs(space, data, pos, maxSteps, dsWorld, eps, hx, hy, hz, writeStride, order);
        if (line.length >= 2) out.push(line);
    }
}

// ---------------- Advanced (PyMOL-style) tracer ----------------

const OBLATION_SPACING = 2; // cells

function getStreamLineAdvanced(volume: Volume, props: VolumeLinesProps | VolumeCylindersProps, out: StreamlinePoint[][]) {
    const { cells: { space, data } } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const [nx, ny, nz] = space.dimensions as Vec3;
    const { hx, hy, hz } = getCellSize(gridToCartn);

    // map range (full grid)
    const i0 = 0, j0 = 0, k0 = 0;
    const i1 = nx, j1 = ny, k1 = nz;

    // precompute WORLD gradients (units: 1/Å) on the whole grid
    const grad = buildGradientFieldWorld(space, data, hx, hy, hz);

    // index → world (we store index coords and transform at geometry stage)
    const toIndexPoint = (i:number, j:number, k:number, fx:number, fy:number, fz:number, outV:Vec3) => { Vec3.set(outV, i + fx, j + fy, k + fz); };

    // randomized order of cells
    const rx = i1 - i0, ry = j1 - j0, rz = k1 - k0;
    const rangeSize = rx * ry * rz;
    const order = new Int32Array(rangeSize * 3);
    {
        let t = 0;
        for (let k = k0; k < k1; k++) for (let j = j0; j < j1; j++) for (let i = i0; i < i1; i++) { order[t++] = i; order[t++] = j; order[t++] = k; }
        // deterministic but good shuffle
        let s = (nx * ny * nz) >>> 0; const rnd = () => { s ^= s << 13; s ^= s >>> 17; s ^= s << 5; return (s >>> 0) / 4294967296; };
        for (let a = 0; a < rangeSize; a++) {
            const p = (Math.floor(rnd() * rangeSize) * 3) | 0;
            const q = (Math.floor(rnd() * rangeSize) * 3) | 0;
            const t0 = order[p], t1 = order[p+1], t2 = order[p+2]; order[p] = order[q]; order[p+1] = order[q+1]; order[p+2] = order[q+2]; order[q] = t0; order[q+1] = t1; order[q+2] = t2;
        }
    }

    // oblation flags
    const flag = new Uint8Array(rx * ry * rz);
    const strideX = 1, strideY = rx, strideZ = rx * ry;
    const flagIndex = (i:number, j:number, k:number) => ((i - i0) * strideX) + ((j - j0) * strideY) + ((k - k0) * strideZ);

    const minDot = 0.0; // disallow sharp reversals
    const maxSteps = Math.max(1, props.maxSteps | 0);
    const dsWorld = props.stepSize; // Å
    const writeStride = Math.max(1, props.writeStride | 0);
    const minSlope = Math.max(props.minSpeed, 1e-6); // 1/Å

    // helpers to keep (fx,fy,fz) in [0,1) while moving voxel indices
    function normalizeFrac(f:number, l:number, lo:number, hi:number) {
        while (f < 0) { f += 1; l--; }
        while (f >= 1) { f -= 1; l++; }
        return (l < lo || l > hi) ? { f, l, done: true } : { f, l, done: false };
    }

    // emit segments similar to PyMOL bookkeeping but directly into out as separate polylines
    const pushNewPolyline = () => { out.push([]); return out.length - 1; };

    for (let a = 0; a < rangeSize; a++) {
        const activeCells: number[] = [];

        const walk = (sign: 1 | -1) => {
            let li = order[a*3], lj = order[a*3+1], lk = order[a*3+2];
            let fx = 0.5, fy = 0.5, fz = 0.5; // start at voxel center (stabler)
            let nVert = 0; let havePrev = false; const prev = Vec3();
            const segIdx = pushNewPolyline();

            for (let step = 0; step < maxSteps; step++) {
                // keep within valid trilinear range [0..N-2]
                let r = normalizeFrac(fx, li, i0, i1 - 2); fx = r.f; li = r.l; if (r.done) break;
                r = normalizeFrac(fy, lj, j0, j1 - 2); fy = r.f; lj = r.l; if (r.done) break;
                r = normalizeFrac(fz, lk, k0, k1 - 2); fz = r.f; lk = r.l; if (r.done) break;

                if (flag[flagIndex(li, lj, lk)]) break;

                const level = interpolateValueFrac(space, data, li, lj, lk, fx, fy, fz);
                if (level < props.minLevel || level > props.maxLevel) break;

                const gW = Vec3();
                interpolateGradientWorld(space, grad, li, lj, lk, fx, fy, fz, gW);
                const gm = Vec3.magnitude(gW);
                if (gm < minSlope) break;
                Vec3.scale(gW, gW, sign / gm); // unit world dir with sign

                // decimate writes in tracer: only record every Nth step
                if ((step % writeStride) === 0) {
                    // add point (index coords; geometry stage does gridToCartn)
                    const pw = Vec3();
                    toIndexPoint(li, lj, lk, fx, fy, fz, pw);
                    out[segIdx].push({ x: pw[0], y: pw[1], z: pw[2], value: level });
                    nVert++;
                }

                // record visited cell for oblation (store once per change)
                const n = activeCells.length;
                if (n < 3 || activeCells[n-3] !== li || activeCells[n-2] !== lj || activeCells[n-1] !== lk) {
                    activeCells.push(li, lj, lk);
                }

                // directional coherence (avoid flip-flop)
                if (havePrev) {
                    const dp = gW[0]*prev[0] + gW[1]*prev[1] + gW[2]*prev[2];
                    if (dp < minDot) break;
                }
                Vec3.copy(prev, gW); havePrev = true;

                // advance in INDEX coordinates using a world step
                fx += gW[0] * (dsWorld / hx);
                fy += gW[1] * (dsWorld / hy);
                fz += gW[2] * (dsWorld / hz);

                // decimate writes
                if ((step % writeStride) !== 0) {
                    // nothing; sampling throttled by outer segment build
                }
            }

            // discard degenerate polylines
            if (nVert < 2) {
                out[segIdx] = [] as StreamlinePoint[]; // keep slot empty; harmless to geometry stage
            }
        };

        walk(+1);
        walk(-1);

        // oblation (flag cells in a sphere of radius R around visited cells)
        const R = OBLATION_SPACING|0; const R2 = R*R;
        for (let t = 0; t < activeCells.length; t += 3) {
            const ii = activeCells[t], jj = activeCells[t+1], kk = activeCells[t+2];
            for (let k = Math.max(k0, kk-R); k <= Math.min(k1-1, kk+R); k++) {
                const dz2 = (kk-k)*(kk-k);
                for (let j = Math.max(j0, jj-R); j <= Math.min(j1-1, jj+R); j++) {
                    const dy2 = (jj-j)*(jj-j) + dz2; if (dy2 > R2) continue;
                    for (let i = Math.max(i0, ii-R); i <= Math.min(i1-1, ii+R); i++) {
                        const d2 = (ii-i)*(ii-i) + dy2; if (d2 <= R2) flag[flagIndex(i,j,k)] = 1;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------
// Field sampling (clamped + world-aware gradients)
// ---------------------------------------------------------------------------------------

function clampIndex(space: Tensor.Space, x: number, y: number, z: number) {
    const nx = space.dimensions[0] as number, ny = space.dimensions[1] as number, nz = space.dimensions[2] as number;
    return [Math.max(0, Math.min(nx - 1, x)), Math.max(0, Math.min(ny - 1, y)), Math.max(0, Math.min(nz - 1, z))] as [number, number, number];
}

// Trilinear interpolation from a Vec3 pos (index coords)
function getInterpolatedValue(space: Tensor.Space, data: ArrayLike<number>, pos: Vec3): number {
    const x = Math.floor(pos[0]), y = Math.floor(pos[1]), z = Math.floor(pos[2]);
    const fx = pos[0] - x, fy = pos[1] - y, fz = pos[2] - z;
    let acc = 0;
    for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
        const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
        acc += data[space.dataOffset(x + dx, y + dy, z + dz)] * w;
    }
    return acc;
}

// Same but with explicit (i,j,k) + fractional offsets
function interpolateValueFrac(space: Tensor.Space, data: ArrayLike<number>, i: number, j: number, k: number, fx: number, fy: number, fz: number) {
    let acc = 0;
    for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
        const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
        acc += data[space.dataOffset(i + dx, j + dy, k + dz)] * w;
    }
    return acc;
}

// Central-diff gradient of φ in WORLD units (1/Å), then E = -∇φ
function gradientAtP_world(out: Vec3, space: Tensor.Space, data: ArrayLike<number>, p: Vec3, hx: number, hy: number, hz: number, tmp?: Vec3) {
    const tmp0 = tmp ?? Vec3();
    // sample neighbors with clamping
    tmp0[0] = p[0] + 1; tmp0[1] = p[1]; tmp0[2] = p[2]; let [x1,y1,z1] = clampIndex(space, tmp0[0]|0, tmp0[1]|0, tmp0[2]|0); const φxp = data[space.dataOffset(x1,y1,z1)];
    tmp0[0] = p[0] - 1; tmp0[1] = p[1]; tmp0[2] = p[2]; [x1,y1,z1] = clampIndex(space, tmp0[0]|0, tmp0[1]|0, tmp0[2]|0); const φxm = data[space.dataOffset(x1,y1,z1)];
    tmp0[0] = p[0]; tmp0[1] = p[1] + 1; tmp0[2] = p[2]; [x1,y1,z1] = clampIndex(space, tmp0[0]|0, tmp0[1]|0, tmp0[2]|0); const φyp = data[space.dataOffset(x1,y1,z1)];
    tmp0[0] = p[0]; tmp0[1] = p[1] - 1; tmp0[2] = p[2]; [x1,y1,z1] = clampIndex(space, tmp0[0]|0, tmp0[1]|0, tmp0[2]|0); const φym = data[space.dataOffset(x1,y1,z1)];
    tmp0[0] = p[0]; tmp0[1] = p[1]; tmp0[2] = p[2] + 1; [x1,y1,z1] = clampIndex(space, tmp0[0]|0, tmp0[1]|0, tmp0[2]|0); const φzp = data[space.dataOffset(x1,y1,z1)];
    tmp0[0] = p[0]; tmp0[1] = p[1]; tmp0[2] = p[2] - 1; [x1,y1,z1] = clampIndex(space, tmp0[0]|0, tmp0[1]|0, tmp0[2]|0); const φzm = data[space.dataOffset(x1,y1,z1)];

    out[0] = - (φxp - φxm) / (2 * hx);
    out[1] = - (φyp - φym) / (2 * hy);
    out[2] = - (φzp - φzm) / (2 * hz);
    return out;
}

// Precompute WORLD gradient field (E = -∇φ) at voxel centers
function buildGradientFieldWorld(space: Tensor.Space, data: ArrayLike<number>, hx: number, hy: number, hz: number) {
    const nx = space.dimensions[0] as number, ny = space.dimensions[1] as number, nz = space.dimensions[2] as number;
    const gx = new Float32Array(nx * ny * nz);
    const gy = new Float32Array(nx * ny * nz);
    const gz = new Float32Array(nx * ny * nz);
    const idx = (x: number, y: number, z: number) => space.dataOffset(x, y, z);

    for (let x = 0; x < nx; x++) {
        const xm = Math.max(0, x - 1), xp = Math.min(nx - 1, x + 1);
        for (let y = 0; y < ny; y++) {
            const ym = Math.max(0, y - 1), yp = Math.min(ny - 1, y + 1);
            for (let z = 0; z < nz; z++) {
                const zm = Math.max(0, z - 1), zp = Math.min(nz - 1, z + 1);
                const off = idx(x,y,z);
                gx[off] = - (data[idx(xp, y, z)] - data[idx(xm, y, z)]) / (2 * hx);
                gy[off] = - (data[idx(x, yp, z)] - data[idx(x, ym, z)]) / (2 * hy);
                gz[off] = - (data[idx(x, y, zp)] - data[idx(x, y, zm)]) / (2 * hz);
            }
        }
    }
    return { gx, gy, gz };
}

// Trilinear interpolation of the WORLD gradient
function interpolateGradientWorld(space: Tensor.Space, grad: { gx: Float32Array, gy: Float32Array, gz: Float32Array }, i: number, j: number, k: number, fx: number, fy: number, fz: number, out: Vec3) {
    const { gx, gy, gz } = grad;
    let vx = 0, vy = 0, vz = 0;
    for (let dz = 0; dz <= 1; dz++) for (let dy = 0; dy <= 1; dy++) for (let dx = 0; dx <= 1; dx++) {
        const w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
        const off = space.dataOffset(i + dx, j + dy, k + dz);
        vx += gx[off] * w; vy += gy[off] * w; vz += gz[off] * w;
    }
    out[0] = vx; out[1] = vy; out[2] = vz; return out;
}
