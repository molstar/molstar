/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { SizeTheme } from '../../../mol-theme/size';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { SpheresBuilder } from '../../../mol-geo/geometry/spheres/spheres-builder';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { sphereVertexCount } from '../../../mol-geo/primitive/sphere';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Particle, ParticleList, getFiberParticleMask } from '../../../mol-model/particles/particle-list';
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

export const SpacefillParticlesParams = {
    ...Spheres.Params,
    ...Mesh.Params,
    tryUseImpostor: PD.Boolean(true),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    excludeFibers: PD.Boolean(true, { description: 'Do not show fiber particles (e.g. shown separately by the Fibers representation).' }),
    excludeTargets: PD.Boolean(false, { description: 'Do not show particles that have a reference object (e.g. shown separately by the Target representation).' }),
};
export type SpacefillParticlesParams = typeof SpacefillParticlesParams;
export type SpacefillParticlesProps = PD.Values<SpacefillParticlesParams>;

/** Marks particles that are excluded from the spacefill spheres. */
function getExcludedParticleMask(particles: ParticleList, props: SpacefillParticlesProps): Uint8Array | undefined {
    const fiberMask = props.excludeFibers ? getFiberParticleMask(particles) : undefined;
    const targetMapping = props.excludeTargets ? particles.targetMapping : undefined;
    if (!targetMapping || targetMapping.size === 0) return fiberMask;

    const { count, targets } = particles;
    // note: `fiberMask` is cached on the particle list and must not be mutated
    const mask = fiberMask ? new Uint8Array(fiberMask) : new Uint8Array(count);
    for (let i = 0; i < count; ++i) {
        if (targetMapping.has(targets[i])) mask[i] = 1;
    }
    return mask;
}

// ---- Impostor visual --------------------------------------------------------

function createSpacefillSphereImpostor(_ctx: VisualContext, particles: ParticleList, _theme: Theme, props: SpacefillParticlesProps, spheres?: Spheres): Spheres {
    const { count, coordinates } = particles;
    const excludedMask = getExcludedParticleMask(particles, props);
    const builder = SpheresBuilder.create(Math.max(1, count), Math.max(1, Math.ceil(count / 10)), spheres);
    for (let i = 0; i < count; ++i) {
        if (excludedMask && excludedMask[i]) continue;
        const o = i * 3;
        builder.add(coordinates[o], coordinates[o + 1], coordinates[o + 2], i);
    }
    const result = builder.getSpheres();
    result.setBoundingSphere(Particle.getBoundary(particles).sphere);
    return result;
}

export function SpacefillParticlesImpostorVisual(materialId: number): ParticleVisual<SpacefillParticlesParams> {
    return ParticleVisual<Spheres, SpacefillParticlesParams>({
        defaultProps: PD.getDefaultValues(SpacefillParticlesParams),
        createGeometry: createSpacefillSphereImpostor,
        createLocationIterator: createSpacefillLocationIterator,
        getLoci: getSpacefillLoci,
        eachLocation: eachSpacefillLoci,
        setUpdateState: (state: VisualUpdateState, newParticles: ParticleList, currentParticles: ParticleList, newProps: SpacefillParticlesProps, currentProps: SpacefillParticlesProps) => {
            // sphere positions are baked from particle coordinates, so any new particle data requires rebuilding the geometry
            if (newParticles !== currentParticles) state.createGeometry = true;
            if (newProps.excludeFibers !== currentProps.excludeFibers) state.createGeometry = true;
            if (newProps.excludeTargets !== currentProps.excludeTargets) state.createGeometry = true;
        },
        geometryUtils: Spheres.Utils,
        mustRecreate: (_key: ParticleKey, props: SpacefillParticlesProps, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        },
    }, materialId);
}

// ---- Mesh visual ------------------------------------------------------------

function createSpacefillSphereMesh(_ctx: VisualContext, particles: ParticleList, theme: Theme, props: SpacefillParticlesProps, mesh?: Mesh): Mesh {
    const { count, coordinates } = particles;
    const { detail } = props;
    const excludedMask = getExcludedParticleMask(particles, props);
    const perSphereVertexCount = sphereVertexCount(detail);
    const vertexCount = Math.max(1, perSphereVertexCount * count);
    const builderState = MeshBuilder.createState(vertexCount, Math.ceil(vertexCount / 10), mesh);
    const location = Particle.Location(particles);
    const center = Vec3();
    for (let i = 0; i < count; ++i) {
        if (excludedMask && excludedMask[i]) continue;
        location.index = i;
        const radius = theme.size.size(location);
        Vec3.fromArray(center, coordinates, i * 3);
        builderState.currentGroup = i;
        addSphere(builderState, center, radius, detail);
    }
    const m = MeshBuilder.getMesh(builderState);
    m.setBoundingSphere(Particle.getBoundary(particles).sphere);
    return m;
}

export function SpacefillParticlesMeshVisual(materialId: number): ParticleVisual<SpacefillParticlesParams> {
    return ParticleVisual<Mesh, SpacefillParticlesParams>({
        defaultProps: PD.getDefaultValues(SpacefillParticlesParams),
        createGeometry: createSpacefillSphereMesh,
        createLocationIterator: createSpacefillLocationIterator,
        getLoci: getSpacefillLoci,
        eachLocation: eachSpacefillLoci,
        setUpdateState: (state: VisualUpdateState, newParticles: ParticleList, currentParticles: ParticleList, newProps: SpacefillParticlesProps, currentProps: SpacefillParticlesProps, newTheme: Theme, currentTheme: Theme) => {
            // sphere positions and radii are baked from particle coordinates/theme, so any new particle data or size theme change requires rebuilding the geometry
            if (newParticles !== currentParticles ||
                newProps.detail !== currentProps.detail ||
                newProps.excludeFibers !== currentProps.excludeFibers ||
                newProps.excludeTargets !== currentProps.excludeTargets ||
                !SizeTheme.areEqual(newTheme.size, currentTheme.size)) {
                state.createGeometry = true;
            }
        },
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
    return LocationIterator(count, 1, 1, (groupIndex: number) => {
        location.index = groupIndex;
        return location;
    }, true /* nonInstanceable */);
}

function getSpacefillLoci(pickingId: PickingId, particles: ParticleList, _props: SpacefillParticlesProps, id: number, _geometry: Spheres | Mesh): Loci {
    const { objectId, groupId } = pickingId;
    if (id !== objectId) return EmptyLoci;
    if (groupId < 0 || groupId >= particles.count) return EmptyLoci;
    return Particle.Loci(particles, OrderedSet.ofSingleton(groupId as any));
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

function getAllLoci(particles: ParticleList, props: SpacefillParticlesProps): Loci {
    const { count } = particles;
    if (count === 0) return EmptyLoci;
    const excludedMask = getExcludedParticleMask(particles, props);
    if (!excludedMask) return Particle.Loci(particles, OrderedSet.ofBounds(0, count) as any);
    const indices: number[] = [];
    for (let i = 0; i < count; ++i) {
        if (!excludedMask[i]) indices.push(i);
    }
    return Particle.Loci(particles, OrderedSet.ofSortedArray(indices));
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
    defaultColorTheme: { name: 'particle-index' },
    defaultSizeTheme: { name: 'particle-size' },
    isApplicable: (data: ParticleList) => data.count > 0,
});
