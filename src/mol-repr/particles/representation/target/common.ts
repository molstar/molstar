/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { ParticleList, Particle, ParticleTarget, getParticleTransforms } from '../../../../mol-model/particles/particle-list';
import { Theme } from '../../../../mol-theme/theme';
import { GraphicsRenderObject, createRenderObject } from '../../../../mol-gl/render-object';
import { Loci as ModelLoci, isEveryLoci } from '../../../../mol-model/loci';
import { MarkerAction } from '../../../../mol-util/marker-action';
import { VisualUpdateState } from '../../../util';
import { Visual, VisualContext } from '../../../visual';
import { WebGLContext } from '../../../../mol-gl/webgl/context';
import { Geometry } from '../../../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { createTransform, createIdentityTransform, TransformData } from '../../../../mol-geo/geometry/transform-data';
import { createColors } from '../../../../mol-geo/geometry/color-data';
import { createSizes } from '../../../../mol-geo/geometry/size-data';
import { Interval, OrderedSet } from '../../../../mol-data/int';
import { Spheres } from '../../../../mol-geo/geometry/spheres/spheres';
import { Mesh } from '../../../../mol-geo/geometry/mesh/mesh';
import { Sphere3D } from '../../../../mol-math/geometry';
import { Mat4 } from '../../../../mol-math/linear-algebra';
import { isPromiseLike } from '../../../../mol-util/type-helpers';
import { CustomProperties } from '../../../../mol-model/custom-property';
import { CommonElementSphereParams } from '../../../structure/visual/element-sphere';
import { SizeValues } from '../../../../mol-gl/renderable/schema';
import { BaseGeometry } from '../../../../mol-geo/geometry/base';
import { ColorTheme } from '../../../../mol-theme/color';
import { SizeTheme } from '../../../../mol-theme/size';
import { createStructureTargetGeometry, structureTargetGeometryPropsChanged, structureTargetMustRecreate, structureTargetSizeFactorChanged } from './structure';
import { createShapeTargetGeometry } from './shape';

// ---- Params ----------------------------------------------------------------

export const ParticleTargetRepresentationParams = {
    ...CommonElementSphereParams,
    tryUseImpostor: PD.Boolean(true),
    ...Spheres.Params,
    ...Mesh.Params,
    instanceGranularity: PD.Boolean(true),
    lodLevels: PD.ObjectList({
        minDistance: PD.Numeric(0),
        maxDistance: PD.Numeric(0),
        overlap: PD.Numeric(0),
        stride: PD.Numeric(0),
        scaleBias: PD.Numeric(3, { min: 0.1, max: 10, step: 0.1 }),
    }, o => `${o.stride}`, {
        ...BaseGeometry.CullingLodCategory,
        defaultValue: Spheres.LodLevelsPresets['quality'],
        presets: Object.entries(Spheres.LodLevelsPresets).map(([k, v]) => [v, k])
    })
};
export type ParticleTargetRepresentationParams = typeof ParticleTargetRepresentationParams;
export type ParticleTargetRepresentationProps = PD.Values<ParticleTargetRepresentationParams>;

// ---- Filtered particle list helper -----------------------------------------

/** A lightweight view over a `ParticleList` restricted to particles with a given target ID. */
export function makeFilteredParticleList(particles: ParticleList, targetId: number): { list: ParticleList, indices: Int32Array } {
    const indicesList: number[] = [];
    const { targets, count } = particles;
    for (let i = 0; i < count; i++) {
        if (targets[i] === targetId) indicesList.push(i);
    }
    const n = indicesList.length;
    const indices = new Int32Array(indicesList);
    const keys = new Int32Array(n);
    const targets2 = new Int32Array(n);
    const coordinates = new Float32Array(n * 3);
    const rotations = particles.rotations ? new Float32Array(n * 4) : undefined;
    const radii = particles.radii ? new Float32Array(n) : undefined;

    for (let i = 0; i < n; i++) {
        const src = indices[i];
        keys[i] = particles.keys[src];
        targets2[i] = targetId;
        coordinates[i * 3 + 0] = particles.coordinates[src * 3 + 0];
        coordinates[i * 3 + 1] = particles.coordinates[src * 3 + 1];
        coordinates[i * 3 + 2] = particles.coordinates[src * 3 + 2];
        if (rotations && particles.rotations) {
            rotations[i * 4 + 0] = particles.rotations[src * 4 + 0];
            rotations[i * 4 + 1] = particles.rotations[src * 4 + 1];
            rotations[i * 4 + 2] = particles.rotations[src * 4 + 2];
            rotations[i * 4 + 3] = particles.rotations[src * 4 + 3];
        }
        if (radii && particles.radii) {
            radii[i] = particles.radii[src];
        }
    }

    return {
        list: {
            entryId: particles.entryId,
            label: particles.label,
            count: n,
            keys,
            targets: targets2,
            coordinates,
            rotations,
            radii,
            getParticleLabel: (idx: number) => particles.getParticleLabel(indices[idx]),
            sourceData: particles.sourceData,
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
        },
        indices,
    };
}

// ---- Transform helpers -----------------------------------------------------

/**
 * Build per-instance transforms placing the target geometry at each particle.
 *
 * The geometry is built at its actual coordinates, centered at `invariantBoundingSphere.center`
 * rather than at the origin. Each particle transform is therefore composed with a pre-translation
 * by `-center` so that the geometry renders centred on the particle position. When `scaleByRadius`
 * is set (e.g. for shape targets), each instance is additionally scaled by the particle radius;
 * the scaling is applied before the center subtraction so the basis already accounts for it.
 */
export function createFilteredParticleTransform(
    filtered: ParticleList,
    invariantBoundingSphere: Sphere3D,
    cellSize: number,
    batchSize: number,
    scaleByRadius: boolean,
    transformData?: TransformData
): TransformData {
    const transforms = getParticleTransforms(filtered);
    const instanceCount = transforms.length;
    const transformArray = new Float32Array(instanceCount * 16);

    const { radii } = filtered;
    if (scaleByRadius && radii) {
        for (let i = 0; i < instanceCount; i++) {
            const r = radii[i] > 0 ? radii[i] : 1;
            Mat4.scaleUniformly(transforms[i], transforms[i], r);
        }
    }

    const { center } = invariantBoundingSphere;
    const cx = center[0], cy = center[1], cz = center[2];

    for (let i = 0; i < instanceCount; i++) {
        const t = transforms[i];
        t[12] -= t[0] * cx + t[4] * cy + t[8] * cz;
        t[13] -= t[1] * cx + t[5] * cy + t[9] * cz;
        t[14] -= t[2] * cx + t[6] * cy + t[10] * cz;
        const offset = i * 16;
        for (let j = 0; j < 16; j++) transformArray[offset + j] = t[j];
    }
    return createTransform(transformArray, instanceCount, invariantBoundingSphere, cellSize, batchSize, transformData);
}

// ---- Per-target geometry dispatch ------------------------------------------

/** The reference data backing a target, used for change detection. */
function targetData(target: ParticleTarget): unknown {
    return target.kind === 'structure' ? target.structure : target.shape;
}

/** Whether the target geometry is scaled by the per-particle radius (shapes), or not (structures). */
function targetScaleByRadius(target: ParticleTarget): boolean {
    return target.kind === 'shape';
}

function createTargetGeometry(ctx: VisualContext, target: ParticleTarget, theme: Theme, props: ParticleTargetRepresentationProps, webgl: WebGLContext | undefined, existing?: Geometry): Geometry | Promise<Geometry> {
    switch (target.kind) {
        case 'structure': return createStructureTargetGeometry(ctx, target.structure, theme, props, webgl, existing as Spheres | Mesh | undefined);
        case 'shape': return createShapeTargetGeometry(target.shape, existing);
    }
}

function targetGeometryPropsChanged(target: ParticleTarget, oldProps: ParticleTargetRepresentationProps, newProps: ParticleTargetRepresentationProps, webgl: WebGLContext | undefined): boolean {
    return target.kind === 'structure' ? structureTargetGeometryPropsChanged(oldProps, newProps, webgl) : false;
}

function targetMustRecreate(target: ParticleTarget, oldProps: ParticleTargetRepresentationProps, newProps: ParticleTargetRepresentationProps): boolean {
    return target.kind === 'structure' ? structureTargetMustRecreate(oldProps, newProps) : false;
}

function targetSizeFactorChanged(target: ParticleTarget, oldProps: ParticleTargetRepresentationProps, newProps: ParticleTargetRepresentationProps, webgl: WebGLContext | undefined): boolean {
    return target.kind === 'structure' ? structureTargetSizeFactorChanged(oldProps, newProps, webgl) : false;
}

// ---- Per-target visual -----------------------------------------------------

/**
 * A single visual for one target (structure or shape), managing one renderObject.
 * Geometry = the reference geometry in its invariant frame.
 * Instances = particles from the filtered particle list for this target.
 */
export interface TargetVisual {
    renderObject: GraphicsRenderObject | undefined
    geometryVersion: number
    /** Maps local instance index → global particle index in the unfiltered ParticleList. */
    particleIndices: OrderedSet<number> | undefined
    createOrUpdate(ctx: VisualContext, theme: Theme, props: ParticleTargetRepresentationProps, particles: ParticleList, particleIndices: OrderedSet<number>, target: ParticleTarget, filtered: ParticleList): Promise<void>
    mark(loci: ModelLoci, action: MarkerAction): boolean
    destroy(): void
}

export function createTargetVisual(_targetId: number, materialId: number, webgl?: WebGLContext): TargetVisual {
    const updateState = VisualUpdateState.create();

    let renderObject: GraphicsRenderObject | undefined;
    let geometry: Geometry | undefined;
    let geometryVersion = -1;
    let locationIt: LocationIterator;
    let positionIt: LocationIterator;

    let currentTargetData: unknown;
    let currentFiltered: ParticleList | undefined;
    let currentProps: ParticleTargetRepresentationProps | undefined;
    let currentTheme: Theme | undefined;
    let currentParticleIndices: OrderedSet<number> | undefined;

    /** Location iterator: groupCount = 1 (particle-granularity), instanceCount = particle count for this target.
     *  Uses the original (unfiltered) `ParticleList` so that themes receive the correct particle index. */
    function createLocIt(particles: ParticleList, indices: OrderedSet<number>): LocationIterator {
        const loc = Particle.Location(particles, 0);
        return LocationIterator(1, OrderedSet.size(indices), 1, (_groupIndex, instanceIndex) => {
            loc.index = OrderedSet.getAt(indices, instanceIndex);
            return loc;
        }, false);
    }

    /** Create a renderObject from scratch for a geometry + filtered particles. */
    function createRO(geom: Geometry, particles: ParticleList, particleIndices: OrderedSet<number>, filtered: ParticleList, theme: Theme, props: ParticleTargetRepresentationProps, scaleByRadius: boolean): GraphicsRenderObject {
        const { createValues, createRenderableState } = Geometry.getUtils(geom);
        const instanceCount = filtered.count;
        const transform = instanceCount > 0
            ? createFilteredParticleTransform(filtered, geom.boundingSphere, props.cellSize, props.batchSize, scaleByRadius)
            : createIdentityTransform();
        locationIt = createLocIt(particles, particleIndices);
        const values = createValues(geom, transform, locationIt, theme, props as any);
        const state = createRenderableState(props as any);
        positionIt = Geometry.getUtils(geom).createPositionIterator(geom, values as any);
        return createRenderObject(geom.kind, values, state, materialId);
    }

    async function createOrUpdate(ctx: VisualContext, theme: Theme, props: ParticleTargetRepresentationProps, particles: ParticleList, particleIndices: OrderedSet<number>, target: ParticleTarget, filtered: ParticleList): Promise<void> {
        VisualUpdateState.reset(updateState);

        const data = targetData(target);
        const targetChanged = data !== currentTargetData;
        const particlesChanged = filtered !== currentFiltered || (currentFiltered && filtered.count !== currentFiltered.count);
        const scaleByRadius = targetScaleByRadius(target);

        const geometryPropsChanged = currentProps ? targetGeometryPropsChanged(target, currentProps, props, webgl) : false;
        const colorThemeChanged = !currentTheme || !ColorTheme.areEqual(theme.color, currentTheme.color);
        const sizeThemeChanged = !currentTheme || !SizeTheme.areEqual(theme.size, currentTheme.size);
        const sizeFactorChanged = !!currentProps && targetSizeFactorChanged(target, currentProps, props, webgl);

        if (!renderObject || targetChanged || (currentProps && targetMustRecreate(target, currentProps, props))) {
            updateState.createNew = true;
        } else {
            updateState.createGeometry = targetChanged || geometryPropsChanged;
            updateState.updateMatrix = particlesChanged || targetChanged;
            updateState.updateColor = colorThemeChanged || updateState.createGeometry;
            updateState.updateSize = sizeThemeChanged || sizeFactorChanged || updateState.createGeometry;
        }

        if (updateState.createNew) {
            updateState.createGeometry = true;
            updateState.updateColor = true;
            updateState.updateSize = true;
            updateState.updateMatrix = true;
        }

        let newGeom: Geometry | undefined;
        if (updateState.createGeometry) {
            const built = createTargetGeometry(ctx, target, theme, props, webgl, geometry);
            newGeom = isPromiseLike(built) ? await built : built;
        }

        if (updateState.createNew || newGeom) {
            const g = (newGeom ?? geometry)!;
            geometry = g;
            geometryVersion++;
            renderObject = createRO(g, particles, particleIndices, filtered, theme, props, scaleByRadius);
        } else if (renderObject) {
            if (updateState.updateMatrix && filtered.count > 0) {
                createFilteredParticleTransform(filtered, geometry!.boundingSphere, props.cellSize, props.batchSize, scaleByRadius, renderObject.values as unknown as TransformData);
                const geomUtils = Geometry.getUtils(geometry!);
                geomUtils.updateBoundingSphere(renderObject.values as any, geometry!);
                positionIt = geomUtils.createPositionIterator(geometry!, renderObject.values as any);
            }

            if (updateState.updateColor) {
                locationIt = createLocIt(particles, particleIndices);
                createColors(locationIt, positionIt, theme.color, renderObject.values as any);
            }

            if (updateState.updateSize && 'uSize' in renderObject.values) {
                locationIt = createLocIt(particles, particleIndices);
                createSizes(locationIt, positionIt, theme.size, renderObject.values as SizeValues);
            }

            const geomUtils = Geometry.getUtils(geometry!);
            geomUtils.updateValues(renderObject.values as any, props as any);
            geomUtils.updateRenderableState(renderObject.state, props as any);
        }

        currentTargetData = data;
        currentFiltered = filtered;
        currentProps = { ...props };
        currentTheme = theme;
        currentParticleIndices = particleIndices;
    }

    function mark(loci: ModelLoci, action: MarkerAction): boolean {
        if (!renderObject || !currentParticleIndices) return false;
        const instanceCount = OrderedSet.size(currentParticleIndices);
        const particleIndices = currentParticleIndices;
        function lociApply(loci: ModelLoci, apply: (interval: Interval) => boolean, _isMarking: boolean): boolean {
            // Fast path: whole visual selected
            if (isEveryLoci(loci)) {
                return apply(Interval.ofBounds(0, instanceCount));
            }
            if (!Particle.isLoci(loci) || Particle.isLociEmpty(loci)) return false;
            // Intersect loci indices with this visual's particle indices
            const intersection = OrderedSet.intersect(loci.indices, particleIndices);
            const intersectionSize = OrderedSet.size(intersection);
            if (intersectionSize === 0) return false;
            // Fast path: all instances in this target are covered
            if (intersectionSize === instanceCount) {
                return apply(Interval.ofBounds(0, instanceCount));
            }
            // Iterate intersection in order; indexOf gives monotonically increasing instance indices → merge runs
            let changed = false;
            let start = -1, end = -1;
            OrderedSet.forEach(intersection, particleIdx => {
                const instanceIdx = OrderedSet.indexOf(particleIndices, particleIdx);
                if (start === -1) {
                    start = end = instanceIdx;
                } else if (instanceIdx === end + 1) {
                    end = instanceIdx;
                } else {
                    if (apply(Interval.ofBounds(start, end + 1))) changed = true;
                    start = end = instanceIdx;
                }
            });
            if (start !== -1 && apply(Interval.ofBounds(start, end + 1))) changed = true;
            return changed;
        }
        return Visual.mark(renderObject, loci, action, lociApply);
    }

    function destroy() {
        renderObject = undefined;
        geometry = undefined;
        currentTargetData = undefined;
        currentFiltered = undefined;
        currentProps = undefined;
        currentTheme = undefined;
        currentParticleIndices = undefined;
    }

    return {
        get renderObject() { return renderObject; },
        get geometryVersion() { return geometryVersion; },
        get particleIndices() { return currentParticleIndices; },
        createOrUpdate,
        mark,
        destroy,
    };
}
