/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { ColorTheme } from '../../../mol-theme/color';
import { SizeTheme } from '../../../mol-theme/size';
import { UniformColorTheme } from '../../../mol-theme/color/uniform';
import { UniformSizeTheme } from '../../../mol-theme/size/uniform';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../mol-geo/geometry/lines/lines-builder';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../../mol-geo/geometry/spheres/spheres-builder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { sphereVertexCount } from '../../../mol-geo/primitive/sphere';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Particle, ParticleList } from '../../../mol-model/particles/particle-list';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { Loci, EmptyLoci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { VisualUpdateState } from '../../util';
import { VisualContext } from '../../visual';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Location } from '../../../mol-model/location';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../representation';
import { ParticleRepresentation, ParticleRepresentationProvider } from '../representation';
import { ParticleVisual, ParticleKey } from '../visual';

// ---- Shared constants -------------------------------------------------------

const AxisLengthOptions = { min: 0.1, max: 1000, step: 0.1 } as const;
const PointSizeOptions = { min: 0.1, max: 100, step: 0.1 } as const;

// ---- Position visual (sphere at world origin, instanced per particle) --------

export const PositionParticlesParams = {
    ...Spheres.Params,
    ...Mesh.Params,
    tryUseImpostor: PD.Boolean(true),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    pointSize: PD.Numeric(1, PointSizeOptions, { description: 'Radius used for the particle position marker.' }),
    positionColor: PD.Color(ColorNames.white),
};
export type PositionParticlesParams = typeof PositionParticlesParams;
export type PositionParticlesProps = PD.Values<PositionParticlesParams>;

function positionOverrideTheme(theme: Theme, props: PositionParticlesProps): Theme {
    const colorFromTheme = theme.color.factory !== ColorTheme.EmptyFactory;
    const sizeFromTheme = theme.size.factory !== SizeTheme.EmptyFactory;
    return {
        color: colorFromTheme ? theme.color : UniformColorTheme({} as any, { value: props.positionColor, lightness: 0, saturation: 0 }),
        size: sizeFromTheme ? theme.size : UniformSizeTheme({} as any, { value: props.pointSize }),
    };
}

// Impostor (Spheres geometry – GPU billboard) ---------------------------------

function createPositionSphereImpostor(_ctx: VisualContext, _particles: ParticleList, _theme: Theme, _props: PositionParticlesProps, spheres?: Spheres): Spheres {
    const builder = SpheresBuilder.create(1, 1, spheres);
    builder.add(0, 0, 0, 0);
    return builder.getSpheres();
}

export function PositionParticlesImpostorVisual(materialId: number): ParticleVisual<PositionParticlesParams> {
    return ParticleVisual<Spheres, PositionParticlesParams>({
        defaultProps: PD.getDefaultValues(PositionParticlesParams),
        createGeometry: createPositionSphereImpostor,
        createLocationIterator: createPositionLocationIterator,
        getLoci: getPositionLoci,
        eachLocation: eachPositionLoci,
        setUpdateState: (state: VisualUpdateState, _np: ParticleList, _cp: ParticleList, newProps: PositionParticlesProps, currentProps: PositionParticlesProps) => {
            if (newProps.pointSize !== currentProps.pointSize) state.updateSize = true;
            if (newProps.positionColor !== currentProps.positionColor) state.updateColor = true;
        },
        overrideTheme: positionOverrideTheme,
        geometryUtils: Spheres.Utils,
        mustRecreate: (_key: ParticleKey, props: PositionParticlesProps, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        },
    }, materialId);
}

// Mesh (Sphere mesh – CPU tessellation) ---------------------------------------

function createPositionSphereMesh(_ctx: VisualContext, _particles: ParticleList, _theme: Theme, props: PositionParticlesProps, mesh?: Mesh): Mesh {
    const { detail, pointSize } = props;
    const vertexCount = sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, Math.ceil(vertexCount / 2), mesh);
    builderState.currentGroup = 0;
    addSphere(builderState, Vec3(), pointSize, detail);
    return MeshBuilder.getMesh(builderState);
}

export function PositionParticlesMeshVisual(materialId: number): ParticleVisual<PositionParticlesParams> {
    return ParticleVisual<Mesh, PositionParticlesParams>({
        defaultProps: PD.getDefaultValues(PositionParticlesParams),
        createGeometry: createPositionSphereMesh,
        createLocationIterator: createPositionLocationIterator,
        getLoci: getPositionLoci,
        eachLocation: eachPositionLoci,
        setUpdateState: (state: VisualUpdateState, _np: ParticleList, _cp: ParticleList, newProps: PositionParticlesProps, currentProps: PositionParticlesProps) => {
            state.createGeometry = newProps.detail !== currentProps.detail || newProps.pointSize !== currentProps.pointSize;
            if (newProps.positionColor !== currentProps.positionColor) state.updateColor = true;
        },
        overrideTheme: positionOverrideTheme,
        geometryUtils: Mesh.Utils,
        mustRecreate: (_key: ParticleKey, props: PositionParticlesProps, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        },
    }, materialId);
}

/** Dispatch: returns an impostor visual when WebGL + required extensions are available, mesh otherwise. */
export function PositionParticlesVisual(materialId: number, _particles: ParticleList, props: PD.Values<PositionParticlesParams>, webgl?: WebGLContext): ParticleVisual<PositionParticlesParams> {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth && webgl.extensions.textureFloat
        ? PositionParticlesImpostorVisual(materialId)
        : PositionParticlesMeshVisual(materialId);
}

// Shared helpers for position visuals -----------------------------------------

function createPositionLocationIterator(particles: ParticleList, _geometry: Spheres | Mesh): LocationIterator {
    const { count } = particles;
    const location = Particle.Location(particles, 0);
    return LocationIterator(1, count, 1, (_groupIndex, instanceIndex) => {
        location.index = instanceIndex;
        return location;
    });
}

function getPositionLoci(pickingId: PickingId, particles: ParticleList, _props: PositionParticlesProps, id: number, _geometry: Spheres | Mesh): Loci {
    const { objectId, instanceId } = pickingId;
    if (id !== objectId) return EmptyLoci;
    if (instanceId < 0 || instanceId >= particles.count) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofSingleton(instanceId as any));
}

function eachPositionLoci(loci: Loci, particles: ParticleList, _props: PositionParticlesProps, apply: (interval: Interval) => boolean, _geometry: Spheres | Mesh): boolean {
    if (!Particle.isLoci(loci)) return false;
    if (Particle.isLociEmpty(loci)) return false;
    if (loci.particles !== particles) return false;
    let changed = false;
    OrderedSet.forEach(loci.indices, idx => {
        if (apply(Interval.ofSingleton(idx))) changed = true;
    });
    return changed;
}

// ---- Axes visual (3 orientation lines, instanced per particle) ---------------

export const AxesParticlesParams = {
    ...Lines.Params,
    axisLength: PD.Numeric(10, AxisLengthOptions, { description: 'Length of the particle orientation axes.' }),
    xColor: PD.Color(ColorNames.red),
    yColor: PD.Color(ColorNames.green),
    zColor: PD.Color(ColorNames.blue),
};
export type AxesParticlesParams = typeof AxesParticlesParams;
export type AxesParticlesProps = PD.Values<AxesParticlesParams>;

type AxesLocation = Particle.Location & { groupIndex: number }

function axisColorTheme(props: AxesParticlesProps): ColorTheme<{}> {
    const colors = [props.xColor, props.yColor, props.zColor] as const;
    const theme: ColorTheme<{}> = {
        factory: (() => theme) as any,
        granularity: 'group',
        color: (location: Location) => {
            const g = (location as any).groupIndex as number | undefined;
            if (g !== undefined && g >= 0 && g < 3) return colors[g];
            return Color(0xCCCCCC);
        },
        props: {},
    };
    return theme;
}

function createAxesGeometry(_ctx: VisualContext, _particles: ParticleList, _theme: Theme, props: AxesParticlesProps, lines?: Lines): Lines {
    const { axisLength } = props;
    const builder = LinesBuilder.create(3, 3, lines);
    builder.add(0, 0, 0, axisLength, 0, 0, 0); // X axis, group 0
    builder.add(0, 0, 0, 0, axisLength, 0, 1); // Y axis, group 1
    builder.add(0, 0, 0, 0, 0, axisLength, 2); // Z axis, group 2
    return builder.getLines();
}

function createAxesLocationIterator(particles: ParticleList, _geometry: Lines): LocationIterator {
    const { count } = particles;
    const location = Object.assign(Particle.Location(particles, 0), { groupIndex: 0 }) as AxesLocation;
    return LocationIterator(3, count, 1, (groupIndex, instanceIndex) => {
        location.index = instanceIndex;
        location.groupIndex = groupIndex;
        return location;
    });
}

function getAxesLoci(pickingId: PickingId, particles: ParticleList, _props: AxesParticlesProps, id: number, _geometry: Lines): Loci {
    const { objectId, instanceId } = pickingId;
    if (id !== objectId) return EmptyLoci;
    if (instanceId < 0 || instanceId >= particles.count) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofSingleton(instanceId as any));
}

function eachAxesLoci(loci: Loci, particles: ParticleList, _props: AxesParticlesProps, apply: (interval: Interval) => boolean, _geometry: Lines): boolean {
    if (!Particle.isLoci(loci)) return false;
    if (Particle.isLociEmpty(loci)) return false;
    if (loci.particles !== particles) return false;
    let changed = false;
    OrderedSet.forEach(loci.indices, idx => {
        // Each particle maps to 3 consecutive group slots (one per axis).
        const start = idx * 3;
        if (apply(Interval.ofBounds(start, start + 3))) changed = true;
    });
    return changed;
}

export function AxesParticlesVisual(materialId: number, _particles?: ParticleList, _props?: PD.Values<AxesParticlesParams>, _webgl?: WebGLContext): ParticleVisual<AxesParticlesParams> {
    return ParticleVisual<Lines, AxesParticlesParams>({
        defaultProps: PD.getDefaultValues(AxesParticlesParams),
        createGeometry: createAxesGeometry,
        createLocationIterator: createAxesLocationIterator,
        getLoci: getAxesLoci,
        eachLocation: eachAxesLoci,
        setUpdateState: (state: VisualUpdateState, _np: ParticleList, _cp: ParticleList, newProps: AxesParticlesProps, currentProps: AxesParticlesProps) => {
            state.createGeometry = newProps.axisLength !== currentProps.axisLength;
            if (newProps.xColor !== currentProps.xColor ||
                newProps.yColor !== currentProps.yColor ||
                newProps.zColor !== currentProps.zColor) {
                state.updateColor = true;
            }
        },
        overrideTheme: (_theme: Theme, props: AxesParticlesProps): Theme => ({
            color: axisColorTheme(props),
            size: UniformSizeTheme({} as any, { value: 1 }),
        }),
        geometryUtils: Lines.Utils,
    }, materialId);
}

// ---- Combined OrientationParticles representation ---------------------------

const OrientationParticlesVisualKinds = ['position', 'axes'] as const;
type OrientationParticlesVisualKind = typeof OrientationParticlesVisualKinds[number];

export const OrientationParticlesParams = {
    ...PositionParticlesParams,
    ...AxesParticlesParams,
    visuals: PD.MultiSelect(['position', 'axes'] as OrientationParticlesVisualKind[], PD.arrayToOptions(OrientationParticlesVisualKinds as unknown as OrientationParticlesVisualKind[])),
};
export type OrientationParticlesParams = typeof OrientationParticlesParams;
export type OrientationParticlesProps = PD.Values<OrientationParticlesParams>;

export function getOrientationParticlesParams(_ctx: ThemeRegistryContext, data: ParticleList) {
    const hasRotations = !!data.rotations;
    const visualKinds: OrientationParticlesVisualKind[] = hasRotations
        ? ['position', 'axes']
        : ['position'];
    const defaultVisuals: OrientationParticlesVisualKind[] = hasRotations
        ? ['position', 'axes']
        : ['position'];
    return {
        ...OrientationParticlesParams,
        visuals: PD.MultiSelect(defaultVisuals, PD.arrayToOptions(visualKinds)),
    };
}

function getAllLoci(particles: ParticleList): Loci {
    const { count } = particles;
    if (count === 0) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofBounds(0, count) as any);
}

const OrientationVisuals = {
    'position': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, PositionParticlesParams>) =>
        ParticleRepresentation('Position', ctx, getParams, PositionParticlesVisual, getAllLoci),
    'axes': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, AxesParticlesParams>) =>
        ParticleRepresentation('Axes', ctx, getParams, AxesParticlesVisual, getAllLoci),
};

export type OrientationParticlesRepresentation = ParticleRepresentation<OrientationParticlesParams>

export function OrientationParticlesRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, OrientationParticlesParams>): OrientationParticlesRepresentation {
    return Representation.createMulti('Orientation Particles', ctx, getParams, Representation.StateBuilder, OrientationVisuals as unknown as Representation.Def<ParticleList, OrientationParticlesParams>);
}

export const OrientationParticlesRepresentationProvider = ParticleRepresentationProvider({
    name: 'orientation',
    label: 'Orientation',
    description: 'Displays particles as position markers and orientation axes.',
    factory: OrientationParticlesRepresentation,
    getParams: getOrientationParticlesParams,
    defaultValues: PD.getDefaultValues(OrientationParticlesParams),
    defaultColorTheme: { name: 'particle-index' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (data: ParticleList) => data.count > 0,
});

