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
    seedDensity: PD.Numeric(8, { min: 1, max: 100, step: 1 }, { description: 'Seeds per dimension for gradient lines.' }),
    maxSteps: PD.Numeric(100, { min: 1, max: 1000, step: 1 }, { description: 'Maximum number of steps for gradient lines.' }),
    stepSize: PD.Numeric(0.5, { min: 0.01, max: 10, step: 0.01 }, { description: 'Step size for gradient lines.' }),
    minSpeed: PD.Numeric(0.001, { min: 0, max: 1, step: 0.0001 }, { description: 'Minimum speed for gradient lines.' }),
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
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats)
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
                // newProps.perturbPositions !== currentProps.perturbPositions ||
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
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
    const isoVal = Volume.IsoValue.toAbsolute(props.isoValue, stats).absoluteValue;

    const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const builder = CylindersBuilder.create(count, Math.ceil(count / 2), cylinders);

    const invert = isoVal < 0;

    // Precompute basis vectors and largest cell axis length
    const basis = getBasis(gridToCartn);

    for (let z = 0; z < zn; ++z) {
        for (let y = 0; y < yn; ++y) {
            for (let x = 0; x < xn; ++x) {
                const value = space.get(data, x, y, z);
                if (!invert && value < isoVal || invert && value > isoVal) continue;

                const cellIdx = space.dataOffset(x, y, z);
                if (basis) {
                    Vec3.set(p, x, y, z);
                    Vec3.transformMat4(p, p, gridToCartn);
                    const offset = getRandomOffsetFromBasis(basis);
                    Vec3.add(p, p, offset);
                } else {
                    Vec3.set(p, x, y, z);
                    Vec3.transformMat4(p, p, gridToCartn);
                }
                // builder.add(p[0], p[1], p[2], cellIdx);
            }
        }
    }

    const s = builder.getCylinders();
    s.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return s;
}


export function createVolumeLinesMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeLinesProps, points?: Lines): Lines {
    const { cells: { space, data }, stats } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    // const isoVal = Volume.IsoValue.toAbsolute(props.isoValue, stats).absoluteValue;

    // const { min, max } = stats;
    const [nx, ny, nz] = space.dimensions as Vec3;

    // Streamline parameters
    const seedDensity = props.seedDensity || 8; // seeds per dimension
    const maxSteps = props.maxSteps || 100;
    const stepSize = props.stepSize || 0.5;
    const minSpeed = props.minSpeed || 0.001;

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
                // Trace streamline from this seed point
                const streamline = traceStreamline(
                    space, data, Vec3.create(x, y, z), 
                    maxSteps, stepSize, minSpeed
                );

                if (streamline.length < 2) continue;
                const cellIdx = space.dataOffset(x, y, z);

                // Convert to cartesian coordinates and add to builder
                for (let i = 0; i < streamline.length - 1; i++) {
                    const start = Vec3();
                    const end = Vec3();

                    Vec3.transformMat4(start, streamline[i].position, gridToCartn);
                    Vec3.transformMat4(end, streamline[i + 1].position, gridToCartn);

                    // Color based on volume value at start point
                    // const color = colorScale.color(streamline[i].value);
                    // const group = Color.toRgb(color)[0] * 255; // Use red component as group

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

function traceStreamline(
    space: Tensor.Space,
    data: ArrayLike<number>,
    seedPoint: Vec3,
    maxSteps: number,
    stepSize: number,
    minSpeed: number
): StreamlinePoint[] {
    const streamline: StreamlinePoint[] = [];
    const [nx, ny, nz] = space.dimensions as Vec3;

    const currentPos = Vec3.clone(seedPoint);

    for (let step = 0; step < maxSteps; step++) {
        // Check bounds
        if (currentPos[0] < 1 || currentPos[0] >= nx - 1 ||
            currentPos[1] < 1 || currentPos[1] >= ny - 1 ||
            currentPos[2] < 1 || currentPos[2] >= nz - 1) {
            break;
        }

        // Get interpolated gradient at current position
        const gradient = getInterpolatedGradient(space, data, currentPos);
        const speed = Vec3.magnitude(gradient);

        // Stop if speed is too low
        if (speed < minSpeed) break;

        // Add current point to streamline
        const value = getInterpolatedValue(space, data, currentPos);
        streamline.push({
            position: Vec3.clone(currentPos),
            value: value
        });

        // Runge-Kutta 4th order integration
        const k1 = Vec3.scale(Vec3(), gradient, stepSize / speed);

        const pos2 = Vec3.add(Vec3(), currentPos, Vec3.scale(Vec3(), k1, 0.5));
        const grad2 = getInterpolatedGradient(space, data, pos2);
        const k2 = Vec3.scale(Vec3(), grad2, stepSize / Vec3.magnitude(grad2));

        const pos3 = Vec3.add(Vec3(), currentPos, Vec3.scale(Vec3(), k2, 0.5));
        const grad3 = getInterpolatedGradient(space, data, pos3);
        const k3 = Vec3.scale(Vec3(), grad3, stepSize / Vec3.magnitude(grad3));

        const pos4 = Vec3.add(Vec3(), currentPos, k3);
        const grad4 = getInterpolatedGradient(space, data, pos4);
        const k4 = Vec3.scale(Vec3(), grad4, stepSize / Vec3.magnitude(grad4));

        // Update position using RK4
        const deltaPos = Vec3.scale(Vec3(), 
            Vec3.add(Vec3(), 
                Vec3.add(Vec3(), k1, Vec3.scale(Vec3(), k2, 2)),
                Vec3.add(Vec3(), Vec3.scale(Vec3(), k3, 2), k4)
            ), 
            1/6
        );

        Vec3.add(currentPos, currentPos, deltaPos);
    }

    return streamline;
}

function getInterpolatedGradient(space: Tensor.Space, data: ArrayLike<number>, pos: Vec3): Vec3 {
    // Trilinear interpolation of gradient
    const x = Math.floor(pos[0]);
    const y = Math.floor(pos[1]);
    const z = Math.floor(pos[2]);

    const fx = pos[0] - x;
    const fy = pos[1] - y;
    const fz = pos[2] - z;

    // Sample gradients at 8 corners of the cell
    const gradients = [];
    for (let dz = 0; dz <= 1; dz++) {
        for (let dy = 0; dy <= 1; dy++) {
            for (let dx = 0; dx <= 1; dx++) {
                gradients.push(calculateGradient(space, data, x + dx, y + dy, z + dz));
            }
        }
    }

    // Trilinear interpolation
    const result = Vec3();
    for (let i = 0; i < 8; i++) {
        const weight = 
            ((i & 1) ? fx : (1 - fx)) *
            ((i & 2) ? fy : (1 - fy)) *
            ((i & 4) ? fz : (1 - fz));
        Vec3.scaleAndAdd(result, result, gradients[i], weight);
    }

    return result;
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

function calculateGradient(space: Tensor.Space, data: ArrayLike<number>, x: number, y: number, z: number): Vec3 {
    const gradient = Vec3();

    // Central differences for gradient calculation
    gradient[0] = data[space.dataOffset(x + 1, y, z)] - data[space.dataOffset(x - 1, y, z)];
    gradient[1] = data[space.dataOffset(x, y + 1, z)] - data[space.dataOffset(x, y - 1, z)];
    gradient[2] = data[space.dataOffset(x, y, z + 1)] - data[space.dataOffset(x, y, z - 1)];

    return gradient;
}