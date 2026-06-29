/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { ColorTheme } from '../../../mol-theme/color';
import { SizeTheme } from '../../../mol-theme/size';
import { UniformColorTheme } from '../../../mol-theme/color/uniform';
import { UniformSizeTheme } from '../../../mol-theme/size/uniform';
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
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../representation';
import { ParticleRepresentation, ParticleRepresentationProvider } from '../representation';
import { ParticleVisual, ParticleKey } from '../visual';

// ---- Params -----------------------------------------------------------------

const PointSizeOptions = { min: 0.1, max: 100, step: 0.1 } as const;

export const SpacefillParticlesParams = {
    ...Spheres.Params,
    ...Mesh.Params,
    tryUseImpostor: PD.Boolean(true),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    pointSize: PD.Numeric(1, PointSizeOptions, { description: 'Radius used for the particle position marker.' }),
    positionColor: PD.Color(ColorNames.white),
};
export type SpacefillParticlesParams = typeof SpacefillParticlesParams;
export type SpacefillParticlesProps = PD.Values<SpacefillParticlesParams>;

function spacefillOverrideTheme(theme: Theme, props: SpacefillParticlesProps): Theme {
    const colorFromTheme = theme.color.factory !== ColorTheme.EmptyFactory;
    const sizeFromTheme = theme.size.factory !== SizeTheme.EmptyFactory;
    return {
        color: colorFromTheme ? theme.color : UniformColorTheme({} as any, { value: props.positionColor, lightness: 0, saturation: 0 }),
        size: sizeFromTheme ? theme.size : UniformSizeTheme({} as any, { value: props.pointSize }),
    };
}

// ---- Impostor visual --------------------------------------------------------

function createSpacefillSphereImpostor(_ctx: VisualContext, _particles: ParticleList, _theme: Theme, _props: SpacefillParticlesProps, spheres?: Spheres): Spheres {
    const builder = SpheresBuilder.create(1, 1, spheres);
    builder.add(0, 0, 0, 0);
    return builder.getSpheres();
}

export function SpacefillParticlesImpostorVisual(materialId: number): ParticleVisual<SpacefillParticlesParams> {
    return ParticleVisual<Spheres, SpacefillParticlesParams>({
        defaultProps: PD.getDefaultValues(SpacefillParticlesParams),
        createGeometry: createSpacefillSphereImpostor,
        createLocationIterator: createSpacefillLocationIterator,
        getLoci: getSpacefillLoci,
        eachLocation: eachSpacefillLoci,
        setUpdateState: (state: VisualUpdateState, _np: ParticleList, _cp: ParticleList, newProps: SpacefillParticlesProps, currentProps: SpacefillParticlesProps) => {
            if (newProps.pointSize !== currentProps.pointSize) state.updateSize = true;
            if (newProps.positionColor !== currentProps.positionColor) state.updateColor = true;
        },
        overrideTheme: spacefillOverrideTheme,
        geometryUtils: Spheres.Utils,
        mustRecreate: (_key: ParticleKey, props: SpacefillParticlesProps, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        },
    }, materialId);
}

// ---- Mesh visual ------------------------------------------------------------

function createSpacefillSphereMesh(_ctx: VisualContext, _particles: ParticleList, _theme: Theme, props: SpacefillParticlesProps, mesh?: Mesh): Mesh {
    const { detail, pointSize } = props;
    const vertexCount = sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, Math.ceil(vertexCount / 2), mesh);
    builderState.currentGroup = 0;
    addSphere(builderState, Vec3(), pointSize, detail);
    return MeshBuilder.getMesh(builderState);
}

export function SpacefillParticlesMeshVisual(materialId: number): ParticleVisual<SpacefillParticlesParams> {
    return ParticleVisual<Mesh, SpacefillParticlesParams>({
        defaultProps: PD.getDefaultValues(SpacefillParticlesParams),
        createGeometry: createSpacefillSphereMesh,
        createLocationIterator: createSpacefillLocationIterator,
        getLoci: getSpacefillLoci,
        eachLocation: eachSpacefillLoci,
        setUpdateState: (state: VisualUpdateState, _np: ParticleList, _cp: ParticleList, newProps: SpacefillParticlesProps, currentProps: SpacefillParticlesProps) => {
            state.createGeometry = newProps.detail !== currentProps.detail || newProps.pointSize !== currentProps.pointSize;
            if (newProps.positionColor !== currentProps.positionColor) state.updateColor = true;
        },
        overrideTheme: spacefillOverrideTheme,
        geometryUtils: Mesh.Utils,
        mustRecreate: (_key: ParticleKey, props: SpacefillParticlesProps, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        },
    }, materialId);
}

/** Dispatch: returns an impostor visual when WebGL + required extensions are available, mesh otherwise. */
export function SpacefillParticlesVisual(materialId: number, _particles: ParticleList, props: PD.Values<SpacefillParticlesParams>, webgl?: WebGLContext): ParticleVisual<SpacefillParticlesParams> {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth && webgl.extensions.textureFloat
        ? SpacefillParticlesImpostorVisual(materialId)
        : SpacefillParticlesMeshVisual(materialId);
}

// ---- Shared helpers ---------------------------------------------------------

function createSpacefillLocationIterator(particles: ParticleList, _geometry: Spheres | Mesh): LocationIterator {
    const { count } = particles;
    const location = Particle.Location(particles, 0);
    return LocationIterator(1, count, 1, (_groupIndex, instanceIndex) => {
        location.index = instanceIndex;
        return location;
    });
}

function getSpacefillLoci(pickingId: PickingId, particles: ParticleList, _props: SpacefillParticlesProps, id: number, _geometry: Spheres | Mesh): Loci {
    const { objectId, instanceId } = pickingId;
    if (id !== objectId) return EmptyLoci;
    if (instanceId < 0 || instanceId >= particles.count) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofSingleton(instanceId as any));
}

function eachSpacefillLoci(loci: Loci, particles: ParticleList, _props: SpacefillParticlesProps, apply: (interval: Interval) => boolean, _geometry: Spheres | Mesh): boolean {
    if (!Particle.isLoci(loci)) return false;
    if (Particle.isLociEmpty(loci)) return false;
    if (loci.particles !== particles) return false;
    let changed = false;
    OrderedSet.forEach(loci.indices, idx => {
        if (apply(Interval.ofSingleton(idx))) changed = true;
    });
    return changed;
}

// ---- Representation ---------------------------------------------------------

const SpacefillVisuals = {
    'spacefill': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, SpacefillParticlesParams>) =>
        ParticleRepresentation('Spacefill', ctx, getParams, SpacefillParticlesVisual, getAllLoci),
};

function getAllLoci(particles: ParticleList): Loci {
    const { count } = particles;
    if (count === 0) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofBounds(0, count) as any);
}

export const SpacefillParticlesParamsGet = (_ctx: ThemeRegistryContext, _data: ParticleList) => SpacefillParticlesParams;

export type SpacefillParticlesRepresentation = ParticleRepresentation<SpacefillParticlesParams>

export function SpacefillParticlesRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<ParticleList, SpacefillParticlesParams>): SpacefillParticlesRepresentation {
    return Representation.createMulti('Spacefill Particles', ctx, getParams, Representation.StateBuilder, SpacefillVisuals as unknown as Representation.Def<ParticleList, SpacefillParticlesParams>);
}

export const SpacefillParticlesRepresentationProvider = ParticleRepresentationProvider({
    name: 'spacefill',
    label: 'Spacefill',
    description: 'Displays particles as spacefill spheres.',
    factory: SpacefillParticlesRepresentation,
    getParams: SpacefillParticlesParamsGet,
    defaultValues: PD.getDefaultValues(SpacefillParticlesParams),
    defaultColorTheme: { name: 'particle-entity' },
    defaultSizeTheme: { name: 'particle-size' },
    isApplicable: (data: ParticleList) => data.count > 0,
});
