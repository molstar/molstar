/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Grid, Volume } from '../../mol-model/volume';
import { VisualContext } from '../visual';
import { Theme, ThemeRegistryContext } from '../../mol-theme/theme';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { VolumeVisual, VolumeRepresentation, VolumeRepresentationProvider, VolumeKey } from './representation';
import { VisualUpdateState } from '../util';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../representation';
import { PickingId } from '../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { createVolumeCellLocationIterator, eachVolumeLoci } from './util';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { sphereVertexCount } from '../../mol-geo/primitive/sphere';
import { Points } from '../../mol-geo/geometry/points/points';
import { PointsBuilder } from '../../mol-geo/geometry/points/points-builder';
import { Mat4 } from '../../mol-math/linear-algebra';

export const VolumeDotParams = {
    isoValue: Volume.IsoValueParam,
    perturbPositions: PD.Boolean(false)
};
export type VolumeDotParams = typeof VolumeDotParams
export type VolumeDotProps = PD.Values<VolumeDotParams>

//

export const VolumeSphereParams = {
    ...Spheres.Params,
    ...Mesh.Params,
    ...VolumeDotParams,
    tryUseImpostor: PD.Boolean(true),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
};
export type VolumeSphereParams = typeof VolumeSphereParams
export type VolumeSphereProps = PD.Values<VolumeSphereParams>

export function VolumeSphereVisual(materialId: number, volume: Volume, key: number, props: PD.Values<VolumeSphereParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth && webgl.extensions.textureFloat
        ? VolumeSphereImpostorVisual(materialId)
        : VolumeSphereMeshVisual(materialId);
}

export function VolumeSphereImpostorVisual(materialId: number): VolumeVisual<VolumeSphereParams> {
    return VolumeVisual<Spheres, VolumeSphereParams>({
        defaultProps: PD.getDefaultValues(VolumeSphereParams),
        createGeometry: createVolumeSphereImpostor,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getDotLoci,
        eachLocation: eachDot,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<VolumeSphereParams>, currentProps: PD.Values<VolumeSphereParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                newProps.perturbPositions !== currentProps.perturbPositions
            );
        },
        geometryUtils: Spheres.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeSphereParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

export function VolumeSphereMeshVisual(materialId: number): VolumeVisual<VolumeSphereParams> {
    return VolumeVisual<Mesh, VolumeSphereParams>({
        defaultProps: PD.getDefaultValues(VolumeSphereParams),
        createGeometry: createVolumeSphereMesh,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getDotLoci,
        eachLocation: eachDot,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<VolumeSphereParams>, currentProps: PD.Values<VolumeSphereParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                newProps.perturbPositions !== currentProps.perturbPositions ||
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            );
        },
        geometryUtils: Mesh.Utils,
        mustRecreate: (volumekey: VolumeKey, props: PD.Values<VolumeSphereParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
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

export function createVolumeSphereImpostor(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeSphereProps, spheres?: Spheres): Spheres {
    const { cells: { space, data }, stats } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const isoVal = Volume.IsoValue.toAbsolute(props.isoValue, stats).absoluteValue;

    const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const builder = SpheresBuilder.create(count, Math.ceil(count / 2), spheres);

    const invert = isoVal < 0;

    // Precompute basis vectors and largest cell axis length
    const basis = props.perturbPositions ? getBasis(gridToCartn) : undefined;

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
                builder.add(p[0], p[1], p[2], cellIdx);
            }
        }
    }

    const s = builder.getSpheres();
    s.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return s;
}

export function createVolumeSphereMesh(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumeSphereProps, mesh?: Mesh): Mesh {
    const { detail, sizeFactor } = props;

    const { cells: { space, data }, stats } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const isoVal = Volume.IsoValue.toAbsolute(props.isoValue, stats).absoluteValue;

    const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const vertexCount = count * sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, Math.ceil(vertexCount / 2), mesh);

    const l = Volume.Cell.Location(volume);
    const themeSize = theme.size.size;
    const invert = isoVal < 0;

    // Precompute basis vectors and largest cell axis length
    const basis = props.perturbPositions ? getBasis(gridToCartn) : undefined;

    for (let z = 0; z < zn; ++z) {
        for (let y = 0; y < yn; ++y) {
            for (let x = 0; x < xn; ++x) {
                const value = space.get(data, x, y, z);
                if (!invert && value < isoVal || invert && value > isoVal) continue;

                const cellIdx = space.dataOffset(x, y, z);
                l.cell = cellIdx as Volume.CellIndex;
                const size = themeSize(l) * sizeFactor;
                if (basis) {
                    Vec3.set(p, x, y, z);
                    Vec3.transformMat4(p, p, gridToCartn);
                    const offset = getRandomOffsetFromBasis(basis);
                    Vec3.add(p, p, offset);
                } else {
                    Vec3.set(p, x, y, z);
                    Vec3.transformMat4(p, p, gridToCartn);
                }
                builderState.currentGroup = cellIdx;
                addSphere(builderState, p, size, detail);
            }
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    m.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return m;
}

//

export const VolumePointParams = {
    ...Points.Params,
    ...VolumeDotParams,
};
export type VolumePointParams = typeof VolumePointParams
export type VolumePointProps = PD.Values<VolumePointParams>

export function VolumePointVisual(materialId: number): VolumeVisual<VolumePointParams> {
    return VolumeVisual<Points, VolumePointParams>({
        defaultProps: PD.getDefaultValues(VolumePointParams),
        createGeometry: createVolumePoint,
        createLocationIterator: createVolumeCellLocationIterator,
        getLoci: getDotLoci,
        eachLocation: eachDot,
        setUpdateState: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<VolumePointParams>, currentProps: PD.Values<VolumePointParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                !Volume.IsoValue.areSame(newProps.isoValue, currentProps.isoValue, volume.grid.stats) ||
                newProps.perturbPositions !== currentProps.perturbPositions
            );
        },
        geometryUtils: Points.Utils,
    }, materialId);
}

export function createVolumePoint(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: VolumePointProps, points?: Points): Points {
    const { cells: { space, data }, stats } = volume.grid;
    const gridToCartn = Grid.getGridToCartesianTransform(volume.grid);
    const isoVal = Volume.IsoValue.toAbsolute(props.isoValue, stats).absoluteValue;

    const p = Vec3();
    const [xn, yn, zn] = space.dimensions;

    const count = Math.ceil((xn * yn * zn) / 10);
    const builder = PointsBuilder.create(count, Math.ceil(count / 2), points);

    const invert = isoVal < 0;

    // Precompute basis vectors and largest cell axis length
    const basis = props.perturbPositions ? getBasis(gridToCartn) : undefined;

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
                builder.add(p[0], p[1], p[2], cellIdx);
            }
        }
    }

    const pt = builder.getPoints();
    pt.setBoundingSphere(Volume.Isosurface.getBoundingSphere(volume, props.isoValue));
    return pt;
}

//

function getLoci(volume: Volume, props: VolumeDotProps) {
    const instances = Interval.ofLength(volume.instances.length as Volume.InstanceIndex);
    return Volume.Isosurface.Loci(volume, props.isoValue, instances);
}

function getDotLoci(pickingId: PickingId, volume: Volume, key: number, props: VolumeDotProps, id: number) {
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

function eachDot(loci: Loci, volume: Volume, key: number, props: VolumeDotProps, apply: (interval: Interval) => boolean) {
    return eachVolumeLoci(loci, volume, { isoValue: props.isoValue }, apply);
}

//

const DotVisuals = {
    'sphere': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumeSphereParams>) => VolumeRepresentation('Dot sphere', ctx, getParams, VolumeSphereVisual, getLoci),
    'point': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, VolumePointParams>) => VolumeRepresentation('Dot point', ctx, getParams, VolumePointVisual, getLoci),
};

export const DotParams = {
    ...VolumeSphereParams,
    ...VolumePointParams,
    visuals: PD.MultiSelect(['sphere'], PD.objectToOptions(DotVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
};
export type DotParams = typeof DotParams
export function getDotParams(ctx: ThemeRegistryContext, volume: Volume) {
    const p = PD.clone(DotParams);
    p.isoValue = Volume.createIsoValueParam(Volume.IsoValue.relative(2), volume.grid.stats);
    return p;
}

export type DotRepresentation = VolumeRepresentation<DotParams>
export function DotRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Volume, DotParams>): DotRepresentation {
    return Representation.createMulti('Dot', ctx, getParams, Representation.StateBuilder, DotVisuals as unknown as Representation.Def<Volume, DotParams>);
}

export const DotRepresentationProvider = VolumeRepresentationProvider({
    name: 'dot',
    label: 'Dot',
    description: 'Displays dots of volumetric data.',
    factory: DotRepresentation,
    getParams: getDotParams,
    defaultValues: PD.getDefaultValues(DotParams),
    defaultColorTheme: { name: 'uniform' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (volume: Volume) => !Volume.isEmpty(volume) && !Volume.Segmentation.get(volume)
});