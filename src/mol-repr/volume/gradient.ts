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
// import { Color } from '../../mol-util/color';
// import { ColorScale } from '../../mol-util/color';
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
    isoValue: Volume.IsoValueParam
};

export type VolumeGradientParams = typeof VolumeGradientParams
export type VolumeGradientProps = PD.Values<VolumeGradientParams>

export const VolumeLinesParams = {
    ...Lines.Params,
    ...Cylinders.Params,
    ...VolumeGradientParams,
    useCylinder: PD.Boolean(true),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    seedDensity: PD.Numeric(20, { min: 1, max: 100, step: 1 }, { description: 'Seeds per dimension for gradient lines.' }),
    maxSteps: PD.Numeric(1000, { min: 1, max: 2000, step: 1 }, { description: 'Maximum number of steps for gradient lines.' }),
    stepSize: PD.Numeric(0.35, { min: 0.01, max: 10, step: 0.01 }, { description: 'Step size for gradient lines.' }),
    minSpeed: PD.Numeric(0.001, { min: 0, max: 1, step: 1e-6 }, { description: 'Minimum speed for gradient lines.' }),
};
export type VolumeLinesParams = typeof VolumeLinesParams
export type VolumeLinesProps = PD.Values<VolumeLinesParams>

export function VolumeGradientVisual(materialId: number, volume: Volume, key: number, props: PD.Values<VolumeLinesParams>) {
    return props.useCylinder
        ? VolumeCylindersImpostorVisual(materialId)
        : VolumeLinesVisual(materialId);
}


export function VolumeCylindersImpostorVisual(materialId: number): VolumeVisual<VolumeLinesParams> {
    return VolumeVisual<Cylinders, VolumeLinesParams>({
        defaultProps: PD.getDefaultValues(VolumeLinesParams),
        createGeometry: createVolumeCylindersImpostor,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<VolumeLinesParams>, currentProps: PD.Values<VolumeLinesParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats)||
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed
            );
        },
        geometryUtils: Cylinders.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeLinesParams>, webgl?: WebGLContext) => {
            return !props.useCylinder || !webgl;
        }
    }, materialId);
}


export function VolumeLinesVisual(materialId: number): VolumeVisual<VolumeLinesParams> {
    return VolumeVisual<Lines, VolumeLinesParams>({
        defaultProps: PD.getDefaultValues(VolumeLinesParams),
        createGeometry: createVolumeLinesMesh,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<VolumeLinesParams>, currentProps: PD.Values<VolumeSphereParams>, newTheme: Theme, currentTheme: Theme) => {
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
            return props.useCylinder && !!webgl;
        }
    }, materialId);
}

type Basis = { x: Vec3, y: Vec3, z: Vec3, maxScale: number }
function getBasis(m: Mat4): Basis {
    return {
        ...Mat4.extractBasis(m),
        maxScale: Mat4.getMaxScaleOnAxis(m)
    };
}


const offset = Vec3();
function getRandomOffsetFromBasis({ x, y, z, maxScale }: Basis): Vec3 {
    const rx = (Math.random() - 0.5) * maxScale;
    const ry = (Math.random() - 0.5) * maxScale;
    const rz = (Math.random() - 0.5) * maxScale;

    Vec3.scale(offset, x, rx);
    Vec3.scaleAndAdd(offset, offset, y, ry);
    Vec3.scaleAndAdd(offset, offset, z, rz);

    return offset;
}


export function createVolumeCylindersImpostor(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeLinesProps, cylinders?: Cylinders): Cylinders {
    const { cells: { space, data }, stats } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    // const isoVal = Volume.IsoValue.toAbsolute(props.isoValue, stats).absoluteValue;

    // const { min, max } = stats;
    const [nx, ny, nz] = space.dimensions as Vec3;

    const { hx, hy, hz } = getCellSize(gridToCartn);

    // seeds
    const seedDensity = props.seedDensity ?? 8;
    const maxSteps = props.maxSteps ?? 2000; // more steps
    const ds = (props.stepSize ?? 0.35); // in index units; 0.2–0.5 works well
    const eps = props.minSpeed ?? 1e-6; // MUCH smaller for APBS

    // Create color scale
    // const colorScale = ColorScale.create({
    //     domain: [min, max],
    //     listOrName: [Color(0x0000ff), Color(0xff0000)] // blue to red
    // });

    // Generate seed points
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));

    // const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const builder = CylindersBuilder.create(count, Math.ceil(count / 2), points);

    // const invert = isoVal < 0;

    // Precompute basis vectors and largest cell axis length
    // const basis = getBasis(gridToCartn);

    for (let z = seedStep; z < nz - seedStep; z += seedStep) {
        for (let y = seedStep; y < ny - seedStep; y += seedStep) {
            for (let x = seedStep; x < nx - seedStep; x += seedStep) {

                const pos = Vec3.create(x, y, z);

                // Filtrage par potentiel et champ
                const phi = getInterpolatedValue(space, data, pos);
                const g = getInterpolatedGradient(space, data, pos, hx, hy, hz);
                const mag = Vec3.magnitude(g);

                // seuils à ajuster selon tes cartes APBS
                if (Math.abs(phi) > 0.2 || mag < 1e-6) continue;

                // si ok, on trace la streamline
                const streamline = traceStreamlineBothDirs(
                    space, data, pos,
                    maxSteps, ds, eps, hx, hy, hz
                );

                if (streamline.length < 2) continue;
                const cellIdx = space.dataOffset(x, y, z);

                for (let i = 0; i < streamline.length - 1; i++) {
                    const start = Vec3();
                    const end = Vec3();
                    Vec3.transformMat4(start, streamline[i].position, gridToCartn);
                    Vec3.transformMat4(end, streamline[i + 1].position, gridToCartn);
                    builder.add(start[0], start[1], start[2],
                                end[0], end[1], end[2],
                                2.0, true, true, 2, cellIdx);
                    // builder.addVec(start, end, cellIdx);
                }
            }
        }
    }


    const pt = builder.getCylinders();
    pt.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return pt;
}


export function createVolumeLinesMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeLinesProps, points?: Lines): Lines {
    const { cells: { space, data }, stats } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    // const isoVal = Volume.IsoValue.toAbsolute(props.isoValue, stats).absoluteValue;

    // const { min, max } = stats;
    const [nx, ny, nz] = space.dimensions as Vec3;

    const { hx, hy, hz } = getCellSize(gridToCartn);

    // seeds
    const seedDensity = props.seedDensity ?? 8;
    const maxSteps = props.maxSteps ?? 2000; // more steps
    const ds = (props.stepSize ?? 0.35); // in index units; 0.2–0.5 works well
    const eps = props.minSpeed ?? 1e-6; // MUCH smaller for APBS

    // Create color scale
    // const colorScale = ColorScale.create({
    //     domain: [min, max],
    //     listOrName: [Color(0x0000ff), Color(0xff0000)] // blue to red
    // });

    // Generate seed points
    const seedStep = Math.max(1, Math.floor(Math.min(nx, ny, nz) / seedDensity));

    // const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const builder = LinesBuilder.create(count, Math.ceil(count / 2), points);

    // const invert = isoVal < 0;

    // Precompute basis vectors and largest cell axis length
    // const basis = getBasis(gridToCartn);

    for (let z = seedStep; z < nz - seedStep; z += seedStep) {
        for (let y = seedStep; y < ny - seedStep; y += seedStep) {
            for (let x = seedStep; x < nx - seedStep; x += seedStep) {

                const pos = Vec3.create(x, y, z);

                // Filtrage par potentiel et champ
                const phi = getInterpolatedValue(space, data, pos);
                const g = getInterpolatedGradient(space, data, pos, hx, hy, hz);
                const mag = Vec3.magnitude(g);

                // seuils à ajuster selon tes cartes APBS
                if (Math.abs(phi) > 0.2 || mag < 1e-6) continue;

                // si ok, on trace la streamline
                const streamline = traceStreamlineBothDirs(
                    space, data, pos,
                    maxSteps, ds, eps, hx, hy, hz
                );

                if (streamline.length < 2) continue;
                const cellIdx = space.dataOffset(x, y, z);

                for (let i = 0; i < streamline.length - 1; i++) {
                    const start = Vec3();
                    const end = Vec3();
                    Vec3.transformMat4(start, streamline[i].position, gridToCartn);
                    Vec3.transformMat4(end, streamline[i + 1].position, gridToCartn);
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
};



export const GradientParams = {
    ...VolumeLinesParams,
    visuals: PD.MultiSelect(['gradient'], PD.objectToOptions(GradientVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
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
    position: Vec3;
    value: number;
}

function traceOneDirection(
    space: Tensor.Space, data: ArrayLike<number>, seed: Vec3,
    maxSteps: number, ds: number, eps: number,
    hx: number, hy: number, hz: number, dirSign: 1 | -1
): StreamlinePoint[] {

    const line: StreamlinePoint[] = [];
    const [nx, ny, nz] = space.dimensions as Vec3;
    const p = Vec3.clone(seed);

    for (let step = 0; step < maxSteps; step++) {
        if (p[0] < 1 || p[0] > (nx as number) - 2 ||
            p[1] < 1 || p[1] > (ny as number) - 2 ||
            p[2] < 1 || p[2] > (nz as number) - 2) break;

        const g = getInterpolatedGradient(space, data, p, hx, hy, hz);
        let m = Vec3.magnitude(g);
        if (!(m > eps)) break;

        // direction only, signed
        const v = Vec3.scale(Vec3(), g, dirSign / m);

        // RK4 with constant arc-length ds
        const k1 = Vec3.scale(Vec3(), v, ds);
        const g2 = getInterpolatedGradient(space, data, Vec3.add(Vec3(), p, Vec3.scale(Vec3(), k1, 0.5)), hx, hy, hz);
        const v2 = Vec3.scale(Vec3(), g2, dirSign / Math.max(Vec3.magnitude(g2), eps));
        const k2 = Vec3.scale(Vec3(), v2, ds);

        const g3 = getInterpolatedGradient(space, data, Vec3.add(Vec3(), p, Vec3.scale(Vec3(), k2, 0.5)), hx, hy, hz);
        const v3 = Vec3.scale(Vec3(), g3, dirSign / Math.max(Vec3.magnitude(g3), eps));
        const k3 = Vec3.scale(Vec3(), v3, ds);

        const g4 = getInterpolatedGradient(space, data, Vec3.add(Vec3(), p, k3), hx, hy, hz);
        const v4 = Vec3.scale(Vec3(), g4, dirSign / Math.max(Vec3.magnitude(g4), eps));
        const k4 = Vec3.scale(Vec3(), v4, ds);

        const dp = Vec3.scale(Vec3(),
            Vec3.add(Vec3(),
                Vec3.add(Vec3(), k1, Vec3.scale(Vec3(), k2, 2)),
                Vec3.add(Vec3(), Vec3.scale(Vec3(), k3, 2), k4)
            ), 1 / 6);

        const value = getInterpolatedValue(space, data, p);
        line.push({ position: Vec3.clone(p), value });

        Vec3.add(p, p, dp);
    }
    return line;
}

function traceStreamlineBothDirs(
    space: Tensor.Space, data: ArrayLike<number>, seed: Vec3,
    maxSteps: number, ds: number, eps: number,
    hx: number, hy: number, hz: number
): StreamlinePoint[] {
    const back = traceOneDirection(space, data, seed, maxSteps, ds, eps, hx, hy, hz, -1);
    const fwd  = traceOneDirection(space, data, seed, maxSteps, ds, eps, hx, hy, hz, +1);
    return back.reverse().concat(fwd.length ? fwd.slice(1) : []);
}

function getInterpolatedGradient(space: Tensor.Space, data: ArrayLike<number>, pos: Vec3, hx: number, hy: number, hz: number): Vec3 {
    const x = Math.floor(pos[0]), y = Math.floor(pos[1]), z = Math.floor(pos[2]);
    const fx = pos[0] - x, fy = pos[1] - y, fz = pos[2] - z;

    const grads: Vec3[] = [];
    for (let dz = 0; dz <= 1; dz++) {
        for (let dy = 0; dy <= 1; dy++) {
            for (let dx = 0; dx <= 1; dx++) {
                grads.push(calculateGradient(space, data, x + dx, y + dy, z + dz, hx, hy, hz));
            }
        }
    }

    const out = Vec3();
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

function calculateGradient(space: Tensor.Space, data: ArrayLike<number>, x: number, y: number, z: number, hx: number, hy: number, hz: number): Vec3 {
    // clamp indices inside [1..N-2] to avoid border reads
    const xi = Math.max(1, Math.min(x, (space.dimensions[0] as number) - 2));
    const yi = Math.max(1, Math.min(y, (space.dimensions[1] as number) - 2));
    const zi = Math.max(1, Math.min(z, (space.dimensions[2] as number) - 2));

    const gx = (data[space.dataOffset(xi + 1, yi, zi)] - data[space.dataOffset(xi - 1, yi, zi)]) / (2 * hx);
    const gy = (data[space.dataOffset(xi, yi + 1, zi)] - data[space.dataOffset(xi, yi - 1, zi)]) / (2 * hy);
    const gz = (data[space.dataOffset(xi, yi, zi + 1)] - data[space.dataOffset(xi, yi, zi - 1)]) / (2 * hz);

    // Electric field E = -∇φ
    return Vec3.create(-gx, -gy, -gz);
}


function getCellSize(gridToCartn: Mat4) {
    const b = Mat4.extractBasis(gridToCartn);
    return {
        hx: Vec3.magnitude(b.x),
        hy: Vec3.magnitude(b.y),
        hz: Vec3.magnitude(b.z)
    };
}
