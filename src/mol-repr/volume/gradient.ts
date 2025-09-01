/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { Cylinders } from '../../mol-geo/geometry/cylinders/cylinders';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
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
import { Points } from '../../mol-geo/geometry/points/points';
import { PointsBuilder } from '../../mol-geo/geometry/points/points-builder';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { sphereVertexCount } from '../../mol-geo/primitive/sphere';
import { addCylinder, BasicCylinderProps } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { collectStreamlines, StreamlineSet } from '../../mol-math/volume/streamlines';

// ---------------------------------------------------------------------------------------
// Shared params
// ---------------------------------------------------------------------------------------

export const VolumeGradientParams = {
    seedDensity: PD.Numeric(10, { min: 1, max: 100, step: 1 }, { description: 'Seeds per dimension for streamlines.' }),
    maxSteps: PD.Numeric(1000, { min: 1, max: 2000, step: 1 }, { description: 'Maximum number of steps for streamlines.' }),
    stepSize: PD.Numeric(0.35, { min: 0.01, max: 10, step: 0.01 }, { description: 'Step size (Å) for streamlines.' }),
    minSpeed: PD.Numeric(0.001, { min: 0, max: 1, step: 1e-6 }, { description: 'Minimum |E| (in 1/Å) to continue tracing.' }),
    minLevel: PD.Numeric(-5.0, { min: -10, max: 10, step: 0.1 }, { description: 'Minimum map value (φ) allowed when tracing.' }),
    maxLevel: PD.Numeric(5.0, { min: -10, max: 10, step: 0.1 }, { description: 'Maximum map value (φ) allowed when tracing.' }),
    algorithm: PD.Select('simple', PD.arrayToOptions(['simple', 'advanced'] as const), { description: 'Streamline algorithm to use.' }),
    mid_interpolation: PD.Boolean(true, { description: 'Use RK2 (midpoint) instead of RK4 for the simple tracer.', hideIf: (o) => o.algorithm !== 'simple' }),
    writeStride: PD.Numeric(4, { min: 1, max: 50, step: 1 }, { description: 'Sample every N integration steps.' }),
    geomStride: PD.Numeric(1, { min: 1, max: 10, step: 1 }, { description: 'Emit every Nth segment to geometry.' }),
};
export type VolumeGradientParams = typeof VolumeGradientParams
export type VolumeGradientProps = PD.Values<VolumeGradientParams>


// LINES
export const VolumeLinesParams = {
    ...Lines.Params,
    ...VolumeGradientParams,
};
export type VolumeLinesParams = typeof VolumeLinesParams
export type VolumeLinesProps = PD.Values<VolumeLinesParams>

// CYLINDERS
export const VolumeCylindersLinesParams = {
    ...Cylinders.Params,
    ...Mesh.Params,
    ...VolumeGradientParams,
    radius: PD.Numeric(1.0, { min: 0.1, max: 5, step: 0.1 }, { description: 'Cylinder radius (Å).' }),
    tryUseImpostor: PD.Boolean(true),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
};
export type VolumeCylindersLinesParams = typeof VolumeCylindersLinesParams
export type VolumeCylindersLinesProps = PD.Values<VolumeCylindersLinesParams>

// POINTS
export const VolumePointsLinesParams = {
    ...Points.Params,
    ...VolumeGradientParams,
};
export type VolumePointsLinesParams = typeof VolumePointsLinesParams
export type VolumePointsLinesProps = PD.Values<VolumePointsLinesParams>

// SPHERES
export const VolumeSpheresLinesParams = {
    ...Spheres.Params,
    ...Mesh.Params,
    ...VolumeGradientParams,
    tryUseImpostor: PD.Boolean(true),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
};
export type VolumeSpheresLinesParams = typeof VolumeSpheresLinesParams
export type VolumeSpheresLinesProps = PD.Values<VolumeSpheresLinesParams>

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
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel ||
                newProps.writeStride !== currentProps.writeStride ||
                newProps.geomStride !== currentProps.geomStride||
                newProps.mid_interpolation !== currentProps.mid_interpolation ||
                newProps.algorithm !== currentProps.algorithm
            );
        },
        geometryUtils: Lines.Utils,
    }, materialId);
}

export function VolumeCylindersLinesVisual(materialId: number, volume: Volume, key: number, props: PD.Values<VolumeCylindersLinesParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth && webgl.extensions.textureFloat
        ? VolumeCylindersImpostorVisual(materialId)
        : VolumeCylindersMeshLinesVisual(materialId);
}

export function VolumeCylindersImpostorVisual(materialId: number): VolumeVisual<VolumeCylindersLinesParams> {
    return VolumeVisual<Cylinders, VolumeCylindersLinesParams>({
        defaultProps: PD.getDefaultValues(VolumeCylindersLinesParams),
        createGeometry: createVolumeCylindersGeometry,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps, currentProps, _newTheme, _currentTheme) => {
            state.createGeometry = (
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.radius !== currentProps.radius ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel ||
                newProps.writeStride !== currentProps.writeStride ||
                newProps.geomStride !== currentProps.geomStride||
                newProps.mid_interpolation !== currentProps.mid_interpolation ||
                newProps.algorithm !== currentProps.algorithm
            );
        },
        geometryUtils: Cylinders.Utils,
        mustRecreate: (_: VolumeKey, __: PD.Values<VolumeCylindersLinesParams>, webgl?: WebGLContext) => !webgl,
    }, materialId);
}

export function VolumeCylindersMeshLinesVisual(materialId: number): VolumeVisual<VolumeCylindersLinesParams> {
    return VolumeVisual<Mesh, VolumeCylindersLinesParams>({
        defaultProps: PD.getDefaultValues(VolumeCylindersLinesParams),
        createGeometry: createVolumeCylindersMeshGeometry,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps, currentProps, _newTheme, _currentTheme) => {
            state.createGeometry = (
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.radius !== currentProps.radius ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel ||
                newProps.writeStride !== currentProps.writeStride ||
                newProps.geomStride !== currentProps.geomStride||
                newProps.mid_interpolation !== currentProps.mid_interpolation ||
                newProps.algorithm !== currentProps.algorithm ||
                newProps.detail !== currentProps.detail

            );
        },
        geometryUtils: Mesh.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeCylindersLinesParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}

export function VolumePointsLinesVisual(materialId: number): VolumeVisual<VolumePointsLinesParams> {
    return VolumeVisual<Points, VolumePointsLinesParams>({
        defaultProps: PD.getDefaultValues(VolumePointsLinesParams),
        createGeometry: createVolumePointsGeometry,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getGradientLoci,
        eachLocation: eachGradient,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps, currentProps, _newTheme, _currentTheme) => {
            state.createGeometry = (
                newProps.seedDensity !== currentProps.seedDensity ||
                newProps.maxSteps !== currentProps.maxSteps ||
                newProps.stepSize !== currentProps.stepSize ||
                newProps.minSpeed !== currentProps.minSpeed ||
                newProps.minLevel !== currentProps.minLevel ||
                newProps.maxLevel !== currentProps.maxLevel ||
                newProps.writeStride !== currentProps.writeStride ||
                newProps.geomStride !== currentProps.geomStride ||
                newProps.mid_interpolation !== currentProps.mid_interpolation ||
                newProps.algorithm !== currentProps.algorithm

            );
        },
        geometryUtils: Points.Utils,
        mustRecreate: (_: VolumeKey, __: PD.Values<VolumePointsLinesParams>, webgl?: WebGLContext) => !webgl,
    }, materialId);
}

export function VolumeSphereLinesVisual(materialId: number, volume: Volume, key: number, props: PD.Values<VolumeSpheresLinesParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth && webgl.extensions.textureFloat
        ? VolumeSpheresImpostorLinesVisual(materialId)
        : VolumeSphereMeshLinesVisual(materialId);
}

export function VolumeSpheresImpostorLinesVisual(materialId: number): VolumeVisual<VolumeSpheresLinesParams> {
    return VolumeVisual<Spheres, VolumeSpheresLinesParams>({
            defaultProps: PD.getDefaultValues(VolumeSpheresLinesParams),
            createGeometry: createVolumeSpheresLinesGeometry,
            createLocationIterator: createVolumeCellLocationIterator,
            getLoci: getGradientLoci,
            eachLocation: eachGradient,
            setUpdateState: (state: VisualUpdateState, volume: Volume, newProps, currentProps, _newTheme, _currentTheme) => {
                state.createGeometry = (
                    newProps.seedDensity !== currentProps.seedDensity ||
                    newProps.maxSteps !== currentProps.maxSteps ||
                    newProps.stepSize !== currentProps.stepSize ||
                    newProps.minSpeed !== currentProps.minSpeed ||
                    newProps.minLevel !== currentProps.minLevel ||
                    newProps.maxLevel !== currentProps.maxLevel ||
                    newProps.writeStride !== currentProps.writeStride ||
                    newProps.geomStride !== currentProps.geomStride ||
                    newProps.mid_interpolation !== currentProps.mid_interpolation ||
                    newProps.algorithm !== currentProps.algorithm

                );
            },
            geometryUtils: Spheres.Utils,
            mustRecreate: (_: VolumeKey, __: PD.Values<VolumeSpheresLinesParams>, webgl?: WebGLContext) => !webgl,
        }, materialId);
}

export function VolumeSphereMeshLinesVisual(materialId: number): VolumeVisual<VolumeSpheresLinesParams> {
    return VolumeVisual<Mesh, VolumeSpheresLinesParams>({
            defaultProps: PD.getDefaultValues(VolumeSpheresLinesParams),
            createGeometry: createVolumeSpheresMeshLinesGeometry,
            createLocationIterator: createVolumeCellLocationIterator,
            getLoci: getGradientLoci,
            eachLocation: eachGradient,
            setUpdateState: (state: VisualUpdateState, volume: Volume, newProps, currentProps, _newTheme, _currentTheme) => {
                state.createGeometry = (
                    newProps.seedDensity !== currentProps.seedDensity ||
                    newProps.maxSteps !== currentProps.maxSteps ||
                    newProps.stepSize !== currentProps.stepSize ||
                    newProps.minSpeed !== currentProps.minSpeed ||
                    newProps.minLevel !== currentProps.minLevel ||
                    newProps.maxLevel !== currentProps.maxLevel ||
                    newProps.writeStride !== currentProps.writeStride ||
                    newProps.geomStride !== currentProps.geomStride ||
                    newProps.mid_interpolation !== currentProps.mid_interpolation ||
                    newProps.algorithm !== currentProps.algorithm ||
                    newProps.detail !== currentProps.detail
                );
            },
            geometryUtils: Mesh.Utils,
            mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeSpheresLinesParams>, webgl?: WebGLContext) => {
                return props.tryUseImpostor && !!webgl;
            }
}, materialId);
}
// ---------------------------------------------------------------------------------------
// Geometry builders (shared streamline collection)
// ---------------------------------------------------------------------------------------

function createVolumeLinesGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumeLinesProps, lines?: Lines): Lines {
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const { cells: { space, data } } = volume.grid;
    const streamLines = collectStreamlines({
        space: space,
        data: data,
        gridToCartn: gridToCartn,
        dims: Vec3.zero(),
        cellSize: Vec3.zero()
    }, props, props.algorithm);

    const segCount = countSegments(streamLines, props.geomStride);
    const builder = LinesBuilder.create(segCount, Math.ceil(segCount / 2), lines);

    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => builder.addVec(a, b, id));

    const g = builder.getLines();
    return g;
}

function createVolumeCylindersGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumeCylindersLinesProps, cylinders?: Cylinders): Cylinders {
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const { cells: { space, data } } = volume.grid;
    const streamLines = collectStreamlines({
        space: space,
        data: data,
        gridToCartn: gridToCartn,
        dims: Vec3.zero(),
        cellSize: Vec3.zero()
    }, props, props.algorithm);

    const segCount = countSegments(streamLines, props.geomStride);
    const builder = CylindersBuilder.create(segCount, Math.ceil(segCount / 2), cylinders);

    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => {
        builder.add(a[0], a[1], a[2], b[0], b[1], b[2], props.radius, true, true, 2, id);
    });

    const g = builder.getCylinders();
    return g;
}


function createVolumeCylindersMeshGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumeCylindersLinesProps, mesh?: Mesh): Mesh {
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const { cells: { space, data } } = volume.grid;
    const streamLines = collectStreamlines({
        space: space,
        data: data,
        gridToCartn: gridToCartn,
        dims: Vec3.zero(),
        cellSize: Vec3.zero()
    }, props, props.algorithm);

    const segCount = countSegments(streamLines, props.geomStride);
    const builder = MeshBuilder.createState(segCount, Math.ceil(segCount / 2), mesh);
    const cylProps: BasicCylinderProps = {
        radiusBottom: props.radius,
        radiusTop: props.radius,
        topCap: true,
        bottomCap: true,
    };
    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => {
        addCylinder(builder, a, b, 1.0, cylProps);
    });

    const g = MeshBuilder.getMesh(builder);
    return g;
}

function createVolumePointsGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumePointsLinesProps, points?: Points): Points {
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const { cells: { space, data } } = volume.grid;
    const streamLines = collectStreamlines({
        space: space,
        data: data,
        gridToCartn: gridToCartn,
        dims: Vec3.zero(),
        cellSize: Vec3.zero()
    }, props, props.algorithm);

    const builder = PointsBuilder.create(streamLines.length, Math.ceil(streamLines.length / 2), points);

    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => {
        builder.add(a[0], a[1], a[2], id);
    });

    const g = builder.getPoints();
    return g;
}

function createVolumeSpheresLinesGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumeSpheresLinesProps, spheres?: Spheres): Spheres {
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const { cells: { space, data } } = volume.grid;
    const streamLines = collectStreamlines({
        space: space,
        data: data,
        gridToCartn: gridToCartn,
        dims: Vec3.zero(),
        cellSize: Vec3.zero()
    }, props, props.algorithm);
    const segCount = countSegments(streamLines, props.geomStride);
    const builder = SpheresBuilder.create(segCount, Math.ceil(segCount / 2), spheres);

    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => {
        builder.add(a[0], a[1], a[2], id);
    });

    const g = builder.getSpheres();
    return g;
}

function createVolumeSpheresMeshLinesGeometry(ctx: VisualContext, volume: Volume, _key: number, _theme: Theme, props: VolumeSpheresLinesProps, mesh?: Mesh): Mesh {
    const { detail, sizeFactor } = props;

    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const { cells: { space, data } } = volume.grid;
    const streamLines = collectStreamlines({
        space: space,
        data: data,
        gridToCartn: gridToCartn,
        dims: Vec3.zero(),
        cellSize: Vec3.zero()
    }, props, props.algorithm);

    const segCount = countSegments(streamLines, props.geomStride) * sphereVertexCount(detail);
    const builder = MeshBuilder.createState(segCount, Math.ceil(segCount / 2), mesh);

    writeSegmentsToBuilder(streamLines, gridToCartn, props.geomStride, (a, b, id) => {
        addSphere(builder, a, sizeFactor, id);
    });

    const g = MeshBuilder.getMesh(builder);
    return g;
}

function countSegments(lines: StreamlineSet, stride: number) {
    let c = 0;
    for (let s = 0; s < lines.length; s++) c += Math.max(0, Math.floor((lines[s].length - 1) / Math.max(1, stride)));
    return c;
}

function writeSegmentsToBuilder(
    streamLines: StreamlineSet,
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
    return Volume.Loci(volume, instances);
}

function getGradientLoci(pickingId: PickingId, volume: Volume, _key: number, props: VolumeGradientProps, id: number) {
    const { objectId, groupId, instanceId } = pickingId;
    if (id === objectId) {
        const granularity = Volume.PickingGranularity.get(volume);
        const instances = OrderedSet.ofSingleton(instanceId as Volume.InstanceIndex);
        if (granularity === 'volume') return Volume.Loci(volume, instances);
        const indices = Interval.ofSingleton(groupId as Volume.CellIndex);
        return Volume.Cell.Loci(volume, [{ indices, instances }]);
    }
    return EmptyLoci;
}

function eachGradient(loci: Loci, volume: Volume, _key: number, props: VolumeGradientProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, { }, apply);
}

const GradientVisuals = {
    'lines': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumeLinesParams>) => VolumeRepresentation('Gradient lines', ctx, getParams, VolumeLinesVisual, getLoci),
    'cylinders': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumeCylindersLinesParams>) => VolumeRepresentation('Gradient cylinders', ctx, getParams, VolumeCylindersLinesVisual, getLoci),
    'points': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumePointsLinesParams>) => VolumeRepresentation('Gradient points', ctx, getParams, VolumePointsLinesVisual, getLoci),
    'spheres': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumeSpheresLinesParams>) => VolumeRepresentation('Gradient spheres', ctx, getParams, VolumeSphereLinesVisual, getLoci),
};

export const GradientParams = {
    ...VolumeLinesParams,
    ...VolumeCylindersLinesParams,
    visuals: PD.MultiSelect(['lines'], PD.objectToOptions(GradientVisuals)),
};
export type GradientParams = typeof GradientParams;

export function getGradientParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(GradientParams);
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
