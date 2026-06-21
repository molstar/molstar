/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Representation, RepresentationContext, RepresentationParamsGetter, RepresentationProvider } from '../../representation';
import { ParticleList, Particle } from '../../../mol-model/particles/particle-list';
import { Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { GraphicsRenderObject, getNextMaterialId } from '../../../mol-gl/render-object';
import { Subject } from 'rxjs';
import { Task } from '../../../mol-task';
import { Loci as ModelLoci, EmptyLoci, isEveryLoci, EveryLoci } from '../../../mol-model/loci';
import { MarkerAction, MarkerActions } from '../../../mol-util/marker-action';
import { VisualUpdateState, LocationCallback, getQualityProps } from '../../util';
import { Visual, VisualContext } from '../../visual';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Geometry } from '../../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { createTransform, createIdentityTransform, TransformData } from '../../../mol-geo/geometry/transform-data';
import { createRenderObject } from '../../../mol-gl/render-object';
import { createColors } from '../../../mol-geo/geometry/color-data';
import { createSizes } from '../../../mol-geo/geometry/size-data';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Sphere3D } from '../../../mol-math/geometry';
import { isPromiseLike } from '../../../mol-util/type-helpers';
import { getParticleTransforms } from '../../../mol-model/particles/particle-list';
import { CustomProperties } from '../../../mol-model/custom-property';
import { CommonElementSphereParams } from '../../structure/visual/element-sphere';
import { ElementSphereImpostorProps, ElementSphereMeshProps } from '../../structure/visual/util/element';
import { createStructureElementSphereImpostor, createStructureElementSphereMesh } from '../../structure/visual/util/element';
import { checkSphereImpostorSupport } from '../../structure/visual/util/common';
import { SizeValues } from '../../../mol-gl/renderable/schema';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { ColorTheme } from '../../../mol-theme/color';
import { SizeTheme } from '../../../mol-theme/size';

// ---- Params ----------------------------------------------------------------

export const ParticlesStructureRepresentationParams = {
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
export type ParticlesStructureRepresentationParams = typeof ParticlesStructureRepresentationParams;
export type ParticlesStructureRepresentationProps = PD.Values<ParticlesStructureRepresentationParams>;

// ---- Filtered particle list helper -----------------------------------------

/** A lightweight view over a `ParticleList` restricted to particles with a given target ID. */
function makeFilteredParticleList(particles: ParticleList, targetId: number): { list: ParticleList, indices: Int32Array } {
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
            getParticleLabel: (idx: number) => particles.getParticleLabel(indices[idx]),
            sourceData: particles.sourceData,
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
        },
        indices,
    };
}

// ---- Transform helpers -----------------------------------------------------

function createFilteredParticleTransform(
    filtered: ParticleList,
    invariantBoundingSphere: Sphere3D,
    cellSize: number,
    batchSize: number,
    transformData?: TransformData
): TransformData {
    const transforms = getParticleTransforms(filtered);
    const instanceCount = transforms.length;
    const transformArray = new Float32Array(instanceCount * 16);

    // The geometry is built at the structure's actual coordinates, centered at
    // `invariantBoundingSphere.center` rather than at the origin.  Each particle
    // transform must therefore be composed with a pre-translation by -center so
    // that the structure renders centred on the particle position:
    //   T_final = T_particle * T(-center)
    // Translation component: t_final = t_particle - R_particle * center
    const { center } = invariantBoundingSphere;
    const cx = center[0], cy = center[1], cz = center[2];

    for (let i = 0; i < instanceCount; i++) {
        const t = transforms[i];
        t[12] -= t[0] * cx + t[4] * cy + t[8] * cz;
        t[13] -= t[1] * cx + t[5] * cy + t[9] * cz;
        t[14] -= t[2] * cx + t[6] * cy + t[10] * cz;
        const offset = i * 16;
        for (let j = 0; j < 16; j++) transformArray[offset + j] = t[j];
        // transformArray.set(t, i * 16); // ASR: much slower (Jun 2026)
    }
    return createTransform(transformArray, instanceCount, invariantBoundingSphere, cellSize, batchSize, transformData);
}

// ---- Per-target visual -----------------------------------------------------

/**
 * A single visual for one target structure, managing one renderObject.
 * Geometry = atoms of the reference structure in their invariant frame.
 * Instances = particles from the filtered particle list for this target.
 */
interface TargetVisual {
    renderObject: GraphicsRenderObject | undefined
    geometryVersion: number
    /** Maps local instance index → global particle index in the unfiltered ParticleList. */
    particleIndices: OrderedSet<number> | undefined
    createOrUpdate(ctx: VisualContext, theme: Theme, props: ParticlesStructureRepresentationProps, particles: ParticleList, particleIndices: OrderedSet<number>, structure: Structure, filtered: ParticleList): Promise<void>
    mark(loci: ModelLoci, action: MarkerAction): boolean
    destroy(): void
}

function createTargetVisual(targetId: number, materialId: number, webgl?: WebGLContext): TargetVisual {
    const updateState = VisualUpdateState.create();

    type GeomType = Spheres | Mesh;
    let renderObject: GraphicsRenderObject | undefined;
    let geometry: GeomType | undefined;
    let geometryVersion = -1;
    let locationIt: LocationIterator;
    let positionIt: LocationIterator;

    let currentStructure: Structure | undefined;
    let currentFiltered: ParticleList | undefined;
    let currentProps: ParticlesStructureRepresentationProps | undefined;
    let currentTheme: Theme | undefined;
    let currentParticleIndices: OrderedSet<number> | undefined;

    function createGeometry(ctx: VisualContext, structure: Structure, theme: Theme, props: ParticlesStructureRepresentationProps, existing?: GeomType): GeomType | Promise<GeomType> {
        const useImpostor = props.tryUseImpostor && checkSphereImpostorSupport(webgl);
        if (useImpostor) {
            const p: ElementSphereImpostorProps = {
                sizeFactor: props.sizeFactor,
                ignoreHydrogens: props.ignoreHydrogens,
                ignoreHydrogensVariant: props.ignoreHydrogensVariant,
                traceOnly: props.traceOnly,
                stride: props.stride,
            };
            return createStructureElementSphereImpostor(ctx, structure, theme, p,
                existing?.kind === 'spheres' ? existing as Spheres : undefined);
        } else {
            const p: ElementSphereMeshProps = {
                detail: props.detail,
                sizeFactor: props.sizeFactor,
                ignoreHydrogens: props.ignoreHydrogens,
                ignoreHydrogensVariant: props.ignoreHydrogensVariant,
                traceOnly: props.traceOnly,
                stride: props.stride,
            };
            return createStructureElementSphereMesh(ctx, structure, theme, p,
                existing?.kind === 'mesh' ? existing as Mesh : undefined);
        }
    }

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
    function createRO(geom: GeomType, particles: ParticleList, particleIndices: OrderedSet<number>, filtered: ParticleList, theme: Theme, props: ParticlesStructureRepresentationProps): GraphicsRenderObject {
        const { createValues, createRenderableState } = Geometry.getUtils(geom);
        const instanceCount = filtered.count;
        const transform = instanceCount > 0
            ? createFilteredParticleTransform(filtered, geom.boundingSphere, props.cellSize, props.batchSize)
            : createIdentityTransform();
        locationIt = createLocIt(particles, particleIndices);
        const values = createValues(geom, transform, locationIt, theme, props as any);
        const state = createRenderableState(props as any);
        positionIt = Geometry.getUtils(geom).createPositionIterator(geom, values as any);
        return createRenderObject(geom.kind, values, state, materialId);
    }

    async function createOrUpdate(ctx: VisualContext, theme: Theme, props: ParticlesStructureRepresentationProps, particles: ParticleList, particleIndices: OrderedSet<number>, structure: Structure, filtered: ParticleList): Promise<void> {
        VisualUpdateState.reset(updateState);

        const structureChanged = structure !== currentStructure;
        const particlesChanged = filtered !== currentFiltered || (currentFiltered && filtered.count !== currentFiltered.count);

        const useImpostor = props.tryUseImpostor && checkSphereImpostorSupport(webgl);
        let geometryPropsChanged = false;
        if (currentProps) {
            if (useImpostor) {
                geometryPropsChanged = (
                    props.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                    props.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                    props.traceOnly !== currentProps.traceOnly ||
                    props.stride !== currentProps.stride
                );
            } else {
                geometryPropsChanged = (
                    props.sizeFactor !== currentProps.sizeFactor ||
                    props.detail !== currentProps.detail ||
                    props.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                    props.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                    props.traceOnly !== currentProps.traceOnly ||
                    props.stride !== currentProps.stride
                );
            }
        }
        const colorThemeChanged = !currentTheme || !ColorTheme.areEqual(theme.color, currentTheme.color);
        const sizeThemeChanged = !currentTheme || !SizeTheme.areEqual(theme.size, currentTheme.size);
        const sizeFactorChanged = !!currentProps && useImpostor && props.sizeFactor !== currentProps.sizeFactor;

        if (!renderObject || structureChanged || (currentProps && props.tryUseImpostor !== currentProps.tryUseImpostor)) {
            updateState.createNew = true;
        } else {
            updateState.createGeometry = structureChanged || geometryPropsChanged;
            updateState.updateMatrix = particlesChanged || structureChanged;
            updateState.updateColor = colorThemeChanged || updateState.createGeometry;
            updateState.updateSize = sizeThemeChanged || sizeFactorChanged || updateState.createGeometry;
        }

        if (updateState.createNew) {
            updateState.createGeometry = true;
            updateState.updateColor = true;
            updateState.updateSize = true;
            updateState.updateMatrix = true;
        }

        let newGeom: GeomType | undefined;
        if (updateState.createGeometry) {
            const built = createGeometry(ctx, structure, theme, props, geometry);
            newGeom = isPromiseLike(built) ? await built : built;
        }

        if (updateState.createNew || newGeom) {
            const g = (newGeom ?? geometry)!;
            geometry = g;
            geometryVersion++;
            renderObject = createRO(g, particles, particleIndices, filtered, theme, props);
        } else if (renderObject) {
            if (updateState.updateMatrix && filtered.count > 0) {
                createFilteredParticleTransform(filtered, geometry!.boundingSphere, props.cellSize, props.batchSize, renderObject.values as unknown as TransformData);
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

        currentStructure = structure;
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
        currentStructure = undefined;
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

// ---- Public representation -------------------------------------------------

export interface ParticlesStructureRepresentation extends Representation<ParticleList, ParticlesStructureRepresentationParams> { }

export type ParticlesStructureRepresentationProvider<Id extends string = string> = RepresentationProvider<ParticleList, ParticlesStructureRepresentationParams, Representation.State, Id>

export function ParticlesStructureRepresentation(
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<ParticleList, ParticlesStructureRepresentationParams>
): ParticlesStructureRepresentation {
    const { webgl } = ctx;
    const updated = new Subject<number>();
    const geometryState = new Representation.GeometryState();
    const renderObjects: GraphicsRenderObject[] = [];
    const _state = Representation.createState();
    let _theme = Theme.createEmpty();

    let version = 0;
    let _particles: ParticleList | undefined;
    let _params: ParticlesStructureRepresentationParams;
    let _props: ParticlesStructureRepresentationProps;
    // Map from targetId → TargetVisual; re-used across updates.
    const targetVisuals = new Map<number, TargetVisual>();

    function createOrUpdate(props: Partial<ParticlesStructureRepresentationProps> = {}, input?: ParticleList) {
        if (input) {
            _particles = input;
            _params = getParams(ctx, _particles);
            if (!_props) _props = PD.getDefaultValues(_params);
        }
        const qualityProps = getQualityProps(Object.assign({}, _props, props));
        Object.assign(_props, props, qualityProps);

        return Task.create('ParticlesStructureRepresentation', async runtime => {
            if (!_particles) return;

            const particles = _particles;
            const targetStructures = Particle.getTargetStructures(particles) ?? new Map<number, Structure>();
            const { targets } = particles;

            // Collect unique target IDs present in this particle list.
            const presentTargets = new Set<number>();
            for (let i = 0; i < particles.count; i++) {
                const tid = targets[i];
                if (targetStructures.has(tid)) presentTargets.add(tid);
            }

            renderObjects.length = 0;

            for (const targetId of presentTargets) {
                const structure = targetStructures.get(targetId)!;
                const { list: filtered, indices: filteredIndices } = makeFilteredParticleList(particles, targetId);
                if (filtered.count === 0) continue;

                let visual = targetVisuals.get(targetId);
                if (!visual) {
                    visual = createTargetVisual(targetId, getNextMaterialId(), webgl);
                    targetVisuals.set(targetId, visual);
                }

                await visual.createOrUpdate({ webgl, runtime }, _theme, _props, particles, OrderedSet.ofSortedArray(filteredIndices), structure, filtered);

                if (visual.renderObject) {
                    renderObjects.push(visual.renderObject);
                    geometryState.add(visual.renderObject.id, visual.geometryVersion);
                }
            }

            // Remove visuals for targets no longer present.
            for (const [tid, visual] of targetVisuals) {
                if (!presentTargets.has(tid)) {
                    visual.destroy();
                    targetVisuals.delete(tid);
                }
            }

            geometryState.snapshot();
            updated.next(version++);
        });
    }

    function setState(state: Partial<Representation.State>) {
        Representation.updateState(_state, state);
        for (const visual of targetVisuals.values()) {
            if (!visual.renderObject) continue;
            if (state.visible !== undefined) visual.renderObject.state.visible = state.visible;
            if (state.alphaFactor !== undefined) visual.renderObject.state.alphaFactor = state.alphaFactor;
            if (state.pickable !== undefined) visual.renderObject.state.pickable = state.pickable;
        }
    }

    function setTheme(theme: Theme) {
        _theme = theme;
    }

    function getLoci(pickingId: PickingId): ModelLoci {
        if (!_particles) return EmptyLoci;
        const { objectId, instanceId } = pickingId;
        for (const visual of targetVisuals.values()) {
            if (visual.renderObject && visual.renderObject.id === objectId) {
                // instanceId is local to this target visual; map to global particle index
                const indices = visual.particleIndices;
                if (!indices || instanceId < 0 || instanceId >= OrderedSet.size(indices)) return EmptyLoci;
                return Particle.Loci(_particles, OrderedSet.ofSingleton(OrderedSet.getAt(indices, instanceId)));
            }
        }
        return EmptyLoci;
    }

    function getAllLoci(): ModelLoci[] {
        if (!_particles) return [];
        return [Particle.Loci(_particles, OrderedSet.ofRange(0, _particles.count - 1))];
    }

    function eachLocation(cb: LocationCallback) {
        if (!_particles) return;
        const particles = _particles;
        const loc = Particle.Location(particles, 0);
        for (let i = 0; i < particles.count; i++) {
            loc.index = i;
            cb(loc, false);
        }
    }

    function mark(loci: ModelLoci, action: MarkerAction): boolean {
        if (!_particles) return false;
        if (!MarkerActions.is(_state.markerActions, action)) return false;
        if (!isEveryLoci(loci) && (!Particle.isLoci(loci) || loci.particles !== _particles)) return false;
        if (Particle.isLoci(loci) && Particle.lociSize(loci) === _particles.count) {
            // Change to `EveryLoci` to allow for downstream optimizations
            loci = EveryLoci;
        }
        let changed = false;
        for (const visual of targetVisuals.values()) {
            changed = visual.mark(loci, action) || changed;
        }
        return changed;
    }

    function destroy() {
        for (const visual of targetVisuals.values()) {
            visual.destroy();
        }
        targetVisuals.clear();
        renderObjects.length = 0;
    }

    return {
        label: 'Particles Structure',
        get updated() { return updated; },
        get renderObjects(): ReadonlyArray<GraphicsRenderObject> { return renderObjects; },
        get props(): Readonly<ParticlesStructureRepresentationProps> { return _props; },
        get params(): Readonly<ParticlesStructureRepresentationParams> { return _params; },
        get state(): Readonly<Representation.State> { return _state; },
        get theme(): Readonly<Theme> { return _theme; },
        get groupCount() {
            return renderObjects.reduce((sum, ro) => sum + ((ro.values as any).uGroupCount?.ref?.value || 0), 0);
        },
        get geometryVersion() { return geometryState.version; },
        createOrUpdate,
        setState,
        setTheme,
        getLoci,
        getAllLoci,
        eachLocation,
        mark,
        destroy,
    };
}

export const ParticlesStructureRepresentationProvider: ParticlesStructureRepresentationProvider = {
    name: 'particles-structure',
    label: 'Particles Structure',
    description: 'Displays particles as instanced spacefill representations of a reference structure.',
    factory: ParticlesStructureRepresentation,
    getParams: (_ctx: RepresentationContext, _particles: ParticleList) => ParticlesStructureRepresentationParams,
    defaultValues: PD.getDefaultValues(ParticlesStructureRepresentationParams),
    defaultColorTheme: { name: 'particle-entity' },
    defaultSizeTheme: { name: 'uniform', props: { value: 1.6 } },
    isApplicable: (data: any) => (data as ParticleList).count > 0,
};


