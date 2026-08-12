/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { ValueCell } from '../../../../mol-util/value-cell';
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
import { isPromiseLike } from '../../../../mol-util/type-helpers';
import { SizeValues } from '../../../../mol-gl/renderable/schema';
import { BaseGeometry } from '../../../../mol-geo/geometry/base';
import { ColorTheme } from '../../../../mol-theme/color';
import { SizeTheme } from '../../../../mol-theme/size';
import { UniformSizeTheme } from '../../../../mol-theme/size/uniform';
import { createStructureTargetGeometry, createStructureTargetSizes, StructureTargetProps, structureTargetGeometryPropsChanged, structureTargetMustRecreate, structureTargetSizeFactorChanged, structureTargetUseImpostor } from './structure';
import { blobSurfaceTargetGeometryPropsChanged, BlobSurfaceTargetProps, createBlobSurfaceTargetGeometry } from './blob-surface';
import { createShapeTargetGeometry } from './shape';
import { createVolumeTargetGeometry, VolumeTargetProps, volumeTargetGeometryPropsChanged, volumeTargetMustRecreate, volumeTargetSizeFactorChanged } from './volume';
import { VolumeIsosurfaceParams } from '../../../volume/isosurface';
import { BlobSurfaceCoreParams } from '../../../structure/visual/util/blob-surface';
import { getQualityProps } from '../../../util';

// ---- Params ----------------------------------------------------------------

/** Removes the params every target kind shares (see `ParticleTargetCommonParams`). */
function withoutCommon<P extends PD.Params>(params: P): P {
    const out: PD.Params = { ...params };
    for (const k of Object.keys(ParticleTargetCommonParams)) delete out[k];
    return out as P;
}

/** Copies `params`, hiding each of them unless the group's `type` is `type`. */
function forType<P extends PD.Params>(type: string, params: P): P {
    const out: PD.Params = {};
    for (const k of Object.keys(params)) {
        out[k] = { ...params[k], hideIf: (p: { type: string }) => p.type !== type };
    }
    return out as P;
}

/** Params shared by all target kinds; the per-kind groups below are merged on top of these. */
export const ParticleTargetCommonParams = {
    ...BaseGeometry.Params,
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
export type ParticleTargetCommonParams = typeof ParticleTargetCommonParams;

/** Params for targets that are structures, rendered as sphere impostors, a sphere mesh or a blob surface. */
export const ParticleTargetStructureParams = {
    ...withoutCommon(Spheres.Params),
    ...withoutCommon(Mesh.Params),
    type: PD.Select<string>('spacefill', PD.arrayToOptions(['spacefill', 'blob-surface']), { isEssential: true, description: 'How the elements of the target are rendered.' }),
    ...forType('spacefill', {
        sizeFactor: PD.Numeric(1, { min: 0.01, max: 10, step: 0.01 }),
        detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
        tryUseImpostor: PD.Boolean(true),
    }),
    ...forType('blob-surface', {
        blobSize: BlobSurfaceCoreParams.blobSize,
        resolution: BlobSurfaceCoreParams.resolution,
        radiusOffset: BlobSurfaceCoreParams.radiusOffset,
        smoothness: BlobSurfaceCoreParams.smoothness,
    }),
};
export type ParticleTargetStructureParams = typeof ParticleTargetStructureParams;

/** Params for targets that are shapes, rendered using the shape's own geometry. */
export const ParticleTargetShapeParams = {
    scaleByRadius: PD.Boolean(true, { description: 'Scale each instance of the shape by the radius of the particle it is placed at.' }),
    ...withoutCommon(Mesh.Params),
};
export type ParticleTargetShapeParams = typeof ParticleTargetShapeParams;

/** Params for targets that are volumes, rendered as an isosurface mesh or as dots. */
export const ParticleTargetVolumeParams = {
    ...withoutCommon(Spheres.Params),
    ...withoutCommon(Mesh.Params),
    type: PD.Select<string>('isosurface', PD.arrayToOptions(['isosurface', 'dot']), { isEssential: true, description: 'How the volume of the target is rendered.' }),
    isoValue: VolumeIsosurfaceParams.isoValue,
    ...forType('isosurface', {
        wrap: VolumeIsosurfaceParams.wrap,
        floodfill: VolumeIsosurfaceParams.floodfill,
    }),
    ...forType('dot', {
        perturbPositions: PD.Boolean(false),
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, { description: 'Radius of a dot. Dots are sized uniformly, the representation size theme does not apply to them.' }),
        detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
        tryUseImpostor: PD.Boolean(true),
    }),
    scaleByRadius: PD.Boolean(false, { description: 'Scale each instance of the volume geometry by the radius of the particle it is placed at.' }),
};
export type ParticleTargetVolumeParams = typeof ParticleTargetVolumeParams;

export const ParticleTargetRepresentationParams = {
    ...ParticleTargetCommonParams,
    structure: PD.Group(ParticleTargetStructureParams, { isExpanded: false, description: 'Parameters applied to targets that are structures.' }),
    shape: PD.Group(ParticleTargetShapeParams, { isExpanded: false, description: 'Parameters applied to targets that are shapes.' }),
    volume: PD.Group(ParticleTargetVolumeParams, { isExpanded: false, description: 'Parameters applied to targets that are volumes.' }),
};
export type ParticleTargetRepresentationParams = typeof ParticleTargetRepresentationParams;
export type ParticleTargetRepresentationProps = PD.Values<ParticleTargetRepresentationParams>;

/**
 * The flat props used for a single target: the common props with the group of the target's
 * kind merged on top. Only the params of the matching kind are present, so kind-specific
 * props are optional here.
 */
export type ParticleTargetProps =
    & PD.Values<ParticleTargetCommonParams>
    & Partial<PD.Values<ParticleTargetStructureParams>>
    & Partial<PD.Values<ParticleTargetShapeParams>>
    & Partial<PD.Values<ParticleTargetVolumeParams>>;

/** Resolves the flat props for a single target, see `createParticleTargetPropsProvider`. */
export type ParticleTargetPropsProvider = (target: ParticleTarget) => ParticleTargetProps

/**
 * Create a provider that flattens the representation props for a single target: common props
 * with the group for the target's kind merged on top, with quality-derived props resolved
 * against the target's own data (e.g. sphere detail from the number of spheres).
 *
 * The per-kind flattening is done once up front and the quality resolution is memoized per
 * target data object, since targets are commonly shared between many particles and the work
 * would otherwise be repeated for every target on every update.
 *
 * Create one provider per update; it holds the props by value.
 */
export function createParticleTargetPropsProvider(props: ParticleTargetRepresentationProps): ParticleTargetPropsProvider {
    const { structure, shape, volume, ...common } = props;
    const bases = {
        structure: { ...common, ...structure } as ParticleTargetProps,
        shape: { ...common, ...shape } as ParticleTargetProps,
        volume: { ...common, ...volume } as ParticleTargetProps,
    };
    const cache = new Map<unknown, ParticleTargetProps>();

    return target => {
        const data = ParticleTarget.data(target);
        let resolved = cache.get(data);
        if (!resolved) {
            const flat = bases[target.kind];
            resolved = { ...flat, ...getQualityProps(flat, data) };
            cache.set(data, resolved);
        }
        return resolved;
    };
}


// ---- Transform helpers -----------------------------------------------------

/**
 * Build per-instance transforms placing the target geometry at each of `indices`.
 *
 * The geometry is built at its actual coordinates, centered at `invariantBoundingSphere.center`
 * rather than at the origin. Each particle transform is therefore composed with a pre-translation
 * by `-center` so that the geometry renders centred on the particle position. When `scaleByRadius`
 * is set (e.g. for shape targets), each instance is additionally scaled by the particle radius;
 * the scaling is applied before the center subtraction so the basis already accounts for it.
 * The bounding sphere used for grid/culling purposes is likewise inflated by the largest radius
 * in play, since it is otherwise unaware of the per-instance scale (see `calcInstanceGrid`).
 */
export function createTargetParticleTransform(
    particles: ParticleList,
    indices: OrderedSet<number>,
    invariantBoundingSphere: Sphere3D,
    cellSize: number,
    batchSize: number,
    scaleByRadius: boolean,
    transformData?: TransformData
): TransformData {
    // Shared across every target visual, so it must not be mutated here.
    const src = getParticleTransforms(particles);
    const instanceCount = OrderedSet.size(indices);
    const transformArray = new Float32Array(instanceCount * 16);

    const { radii } = particles;
    const { center } = invariantBoundingSphere;
    const cx = center[0], cy = center[1], cz = center[2];

    let maxRadius = 1;
    for (let i = 0; i < instanceCount; i++) {
        const so = OrderedSet.getAt(indices, i) * 16;
        const o = i * 16;

        let s = 1;
        if (scaleByRadius && radii) {
            const r = radii[OrderedSet.getAt(indices, i)];
            if (r > 0) s = r;
            if (s > maxRadius) maxRadius = s;
        }

        for (let j = 0; j < 12; j++) transformArray[o + j] = src[so + j] * s;
        transformArray[o + 15] = src[so + 15];

        transformArray[o + 12] = src[so + 12] - (transformArray[o + 0] * cx + transformArray[o + 4] * cy + transformArray[o + 8] * cz);
        transformArray[o + 13] = src[so + 13] - (transformArray[o + 1] * cx + transformArray[o + 5] * cy + transformArray[o + 9] * cz);
        transformArray[o + 14] = src[so + 14] - (transformArray[o + 2] * cx + transformArray[o + 6] * cy + transformArray[o + 10] * cz);
    }

    // `calcInstanceGrid` (used inside `createTransform`) builds its grid/culling bounds from
    // `invariantBoundingSphere.radius` alone and does not account for any scale baked into the
    // per-instance transform matrices. Since scaling-by-radius applies such a scale above, the
    // sphere passed on for grid/culling purposes must be inflated by the largest radius used,
    // otherwise larger particles get under-sized cells/LOD bounds and can be culled incorrectly.
    const gridBoundingSphere = scaleByRadius && maxRadius !== 1
        ? Sphere3D.create(center, invariantBoundingSphere.radius * maxRadius)
        : invariantBoundingSphere;

    return createTransform(transformArray, instanceCount, gridBoundingSphere, cellSize, batchSize, transformData);
}

// ---- Per-target geometry dispatch ------------------------------------------

/** The reference data backing a target, used for change detection. */
function targetData(target: ParticleTarget): unknown {
    return ParticleTarget.data(target);
}

/** Whether the target geometry is scaled by the per-particle radius. */
function targetScaleByRadius(target: ParticleTarget, props: ParticleTargetProps): boolean {
    return target.kind === 'structure' ? false : !!props.scaleByRadius;
}

/** Whether a volume target is rendered as dots rather than as an isosurface. */
function isVolumeDot(target: ParticleTarget, props: ParticleTargetProps): boolean {
    return target.kind === 'volume' && (props as VolumeTargetProps).type === 'dot';
}

let _dotSizeTheme: SizeTheme<any> | undefined;
function getDotSizeTheme(): SizeTheme<any> {
    if (!_dotSizeTheme) _dotSizeTheme = UniformSizeTheme({} as any, { value: 1 });
    return _dotSizeTheme;
}

export function createTargetGeometry(ctx: VisualContext, target: ParticleTarget, theme: Theme, props: ParticleTargetProps, webgl: WebGLContext | undefined, existing?: Geometry): Geometry | Promise<Geometry> {
    const spheresOrMesh = existing?.kind === 'spheres' || existing?.kind === 'mesh' ? existing as Spheres | Mesh : undefined;
    switch (target.kind) {
        case 'structure': return (props as StructureTargetProps).type === 'blob-surface'
            ? createBlobSurfaceTargetGeometry(ctx, target.structure, theme, props as BlobSurfaceTargetProps, existing?.kind === 'mesh' ? existing as Mesh : undefined)
            : createStructureTargetGeometry(ctx, target.structure, theme, props as StructureTargetProps, webgl, spheresOrMesh);
        case 'shape': return createShapeTargetGeometry(target.shape, existing);
        case 'volume': return createVolumeTargetGeometry(ctx, target.volume, theme, props as VolumeTargetProps, webgl, spheresOrMesh);
    }
}

export function targetGeometryPropsChanged(target: ParticleTarget, oldProps: ParticleTargetProps, newProps: ParticleTargetProps, webgl: WebGLContext | undefined): boolean {
    switch (target.kind) {
        case 'structure': return (newProps as StructureTargetProps).type === 'blob-surface'
            ? blobSurfaceTargetGeometryPropsChanged(oldProps as BlobSurfaceTargetProps, newProps as BlobSurfaceTargetProps)
            : structureTargetGeometryPropsChanged(oldProps as StructureTargetProps, newProps as StructureTargetProps, webgl);
        case 'shape': return false;
        case 'volume': return volumeTargetGeometryPropsChanged(target.volume, oldProps as VolumeTargetProps, newProps as VolumeTargetProps, webgl);
    }
}

export function targetMustRecreate(target: ParticleTarget, oldProps: ParticleTargetProps, newProps: ParticleTargetProps): boolean {
    switch (target.kind) {
        case 'structure': return structureTargetMustRecreate(oldProps as StructureTargetProps, newProps as StructureTargetProps);
        case 'shape': return false;
        case 'volume': return volumeTargetMustRecreate(oldProps as VolumeTargetProps, newProps as VolumeTargetProps);
    }
}

function targetSizeFactorChanged(target: ParticleTarget, oldProps: ParticleTargetProps, newProps: ParticleTargetProps, webgl: WebGLContext | undefined): boolean {
    switch (target.kind) {
        case 'structure': return structureTargetSizeFactorChanged(oldProps as StructureTargetProps, newProps as StructureTargetProps, webgl);
        case 'shape': return false;
        case 'volume': return volumeTargetSizeFactorChanged(oldProps as VolumeTargetProps, newProps as VolumeTargetProps, webgl);
    }
}

/** Whether a size theme change requires rebuilding the geometry (sphere meshes bake in their radii). */
function targetSizeThemeAffectsGeometry(target: ParticleTarget, props: ParticleTargetProps, webgl: WebGLContext | undefined): boolean {
    return target.kind === 'structure' && !structureTargetUseImpostor(props as StructureTargetProps, webgl);
}


// ---- Per-target visual -----------------------------------------------------

/**
 * A single visual for one target (structure or shape), managing one renderObject.
 * Geometry = the reference geometry in its invariant frame.
 * Instances = the particles of this target, given as indices into the unfiltered `ParticleList`.
 */
export interface TargetVisual {
    renderObject: GraphicsRenderObject | undefined
    geometryVersion: number
    /** Maps local instance index → global particle index in the unfiltered ParticleList. */
    particleIndices: OrderedSet<number> | undefined
    createOrUpdate(ctx: VisualContext, theme: Theme, props: ParticleTargetProps, particles: ParticleList, particleIndices: OrderedSet<number>, target: ParticleTarget): Promise<void>
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
    let currentProps: ParticleTargetProps | undefined;
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

    /** Create a renderObject from scratch for a geometry + the target's particles. */
    function createRO(geom: Geometry, particles: ParticleList, particleIndices: OrderedSet<number>, theme: Theme, props: ParticleTargetProps, target: ParticleTarget, scaleByRadius: boolean): GraphicsRenderObject {
        const { createValues, createRenderableState } = Geometry.getUtils(geom);
        const instanceCount = OrderedSet.size(particleIndices);
        const transform = instanceCount > 0
            ? createTargetParticleTransform(particles, particleIndices, geom.boundingSphere, props.cellSize, props.batchSize, scaleByRadius)
            : createIdentityTransform();
        locationIt = createLocIt(particles, particleIndices);
        const values = createValues(geom, transform, locationIt, theme, props as any);
        const state = createRenderableState(props as any);
        positionIt = Geometry.getUtils(geom).createPositionIterator(geom, values as any);
        // Sizes of a structure target are per-element, unlike the per-particle colors and picking.
        if (target.kind === 'structure' && 'uSize' in values) {
            createStructureTargetSizes(target.structure, theme, values as unknown as SizeValues);
        }
        return createRenderObject(geom.kind, values, state, materialId);
    }

    async function createOrUpdate(ctx: VisualContext, representationTheme: Theme, props: ParticleTargetProps, particles: ParticleList, particleIndices: OrderedSet<number>, target: ParticleTarget): Promise<void> {
        VisualUpdateState.reset(updateState);

        // Dots are sized uniformly and scaled by `sizeFactor`; the representation's per-particle
        // size theme would make every dot as large as the particle it belongs to.
        const theme = isVolumeDot(target, props)
            ? { color: representationTheme.color, size: getDotSizeTheme() }
            : representationTheme;

        const data = targetData(target);
        const targetChanged = data !== currentTargetData;
        const particlesChanged = particleIndices !== currentParticleIndices;
        const scaleByRadius = targetScaleByRadius(target, props);

        const geometryPropsChanged = currentProps ? targetGeometryPropsChanged(target, currentProps, props, webgl) : false;
        const colorThemeChanged = !currentTheme || !ColorTheme.areEqual(theme.color, currentTheme.color);
        const sizeThemeChanged = !currentTheme || !SizeTheme.areEqual(theme.size, currentTheme.size);
        const sizeFactorChanged = !!currentProps && targetSizeFactorChanged(target, currentProps, props, webgl);

        if (!renderObject || targetChanged || (currentProps && targetMustRecreate(target, currentProps, props))) {
            updateState.createNew = true;
        } else {
            updateState.createGeometry = targetChanged || geometryPropsChanged || (sizeThemeChanged && targetSizeThemeAffectsGeometry(target, props, webgl));
            updateState.updateMatrix = particlesChanged || targetChanged || targetScaleByRadius(target, currentProps!) !== scaleByRadius;
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
            renderObject = createRO(g, particles, particleIndices, theme, props, target, scaleByRadius);
        } else if (renderObject) {
            if (updateState.updateMatrix && OrderedSet.size(particleIndices) > 0) {
                createTargetParticleTransform(particles, particleIndices, geometry!.boundingSphere, props.cellSize, props.batchSize, scaleByRadius, renderObject.values as unknown as TransformData);
                const geomUtils = Geometry.getUtils(geometry!);
                // TODO: needs update? how to do efficeintly? geomUtils is too slow
                // geomUtils.updateBoundingSphere(renderObject.values as any, geometry!);
                positionIt = geomUtils.createPositionIterator(geometry!, renderObject.values as any);
                if ('lodLevels' in renderObject.values) {
                    // instanceGrid (and its cellCount) was just rebuilt above; bump the `lodLevels`
                    // version to force `renderable.cull` to resize/refresh its per-level multidraw
                    // data (`mdbDataList`) instead of reusing stale, potentially undersized arrays -
                    // otherwise cells beyond the old cellCount are silently dropped from culling,
                    // which can leave the wrong (smaller) LOD level's spheres rendered.
                    ValueCell.update(renderObject.values.lodLevels, renderObject.values.lodLevels.ref.value);
                }
            }

            if (updateState.updateColor) {
                locationIt = createLocIt(particles, particleIndices);
                createColors(locationIt, positionIt, theme.color, renderObject.values as any);
            }

            if (updateState.updateSize && 'uSize' in renderObject.values) {
                if (target.kind === 'structure') {
                    createStructureTargetSizes(target.structure, theme, renderObject.values as SizeValues);
                } else {
                    locationIt = createLocIt(particles, particleIndices);
                    createSizes(locationIt, positionIt, theme.size, renderObject.values as SizeValues);
                }
            }
            const geomUtils = Geometry.getUtils(geometry!);
            geomUtils.updateValues(renderObject.values as any, props as any);
            geomUtils.updateRenderableState(renderObject.state, props as any);
        }

        currentTargetData = data;
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
