/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Spherical-harmonic surface visual: reuses Mol*'s molecular-surface mesh, fits
 * a real spherical-harmonic expansion of the radial function r(theta, phi) about
 * the surface centroid, and reconstructs a smooth closed mesh from the fitted
 * coefficients. The parameter `sphericalHarmonicL` controls the level of detail.
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { CommonMolecularSurfaceCalculationParams, computeStructureMolecularSurface, computeUnitMolecularSurface } from './util/molecular-surface';
import { computeMarchingCubesMesh } from '../../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from './util/element';
import { VisualUpdateState } from '../../util';
import { CommonSurfaceParams, ensureReasonableResolution, getConformation, isHydrogen, isTrace } from './util/common';
import { Sphere3D } from '../../../mol-math/geometry';
import { MeshValues } from '../../../mol-gl/renderable/mesh';
import { Texture } from '../../../mol-gl/webgl/texture';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { applyMeshColorSmoothing } from '../../../mol-geo/geometry/mesh/color-smoothing';
import { ColorSmoothingParams, getColorSmoothingProps } from '../../../mol-geo/geometry/base';
import { ValueCell } from '../../../mol-util';
import { ComplexMeshVisual, ComplexVisual } from '../complex-visual';
import { Tensor } from '../../../mol-math/linear-algebra/tensor';
import { Vec3, Mat4 } from '../../../mol-math/linear-algebra';
import { Sphere } from '../../../mol-geo/primitive/sphere';
import { GridLookup3D } from '../../../mol-math/geometry/lookup3d/grid';
import { getBoundary } from '../../../mol-math/geometry/boundary';
import { OrderedSet, Interval, Segmentation } from '../../../mol-data/int';
import { fitSphericalHarmonics, reconstructRadius, shTermCount } from '../../../mol-math/geometry/spherical-harmonics';
import { Task } from '../../../mol-task';
import { calcMolecularSurface } from '../../../mol-math/geometry/molecular-surface';
import { createTransform, TransformData } from '../../../mol-geo/geometry/transform-data';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { SizeTheme } from '../../../mol-theme/size';

export const SphericalHarmonicSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...CommonMolecularSurfaceCalculationParams,
    ...CommonSurfaceParams,
    ...ColorSmoothingParams,
    sphericalHarmonicL: PD.Numeric(8, { min: 2, max: 20, step: 1 }, { description: 'Maximum spherical harmonic degree. Higher values follow the molecular surface more closely.' }),
    reconstructionDetail: PD.Numeric(3, { min: 1, max: 5, step: 1 }, { description: 'Triangulation detail of the reconstructed sphere mesh.' }),
    fitSource: PD.Select('surface', [['surface', 'Molecular surface'], ['atoms', 'Atom positions']] as const, { description: 'Fit the spherical harmonics to the molecular surface, or directly to atom positions (faster, blobbier, no surface step).' }),
    residueL: PD.Numeric(2, { min: 1, max: 4, step: 1 }, { description: 'Spherical harmonic degree for the per-residue blobs (kept low: a residue has too few atoms for high degrees).' }),
    residueDetail: PD.Numeric(2, { min: 1, max: 3, step: 1 }, { description: 'Triangulation detail of each per-residue blob.' }),
};
export type SphericalHarmonicSurfaceMeshParams = typeof SphericalHarmonicSurfaceMeshParams
export type SphericalHarmonicSurfaceMeshProps = PD.Values<SphericalHarmonicSurfaceMeshParams>

type SphericalHarmonicSurfaceMeta = {
    resolution?: number
    colorTexture?: Texture
    // cache of the L-independent base surface, keyed on surface-affecting props
    shCacheKey?: string
    shCenter?: Vec3
    shVertices?: Float32Array
    shGroups?: Float32Array
    shLookup?: GridLookup3D
    shBoundingSphere?: Sphere3D
}

// memoized sphere triangulations per detail level
const sphereCache = new Map<number, ReturnType<typeof Sphere>>();
function getSphere(detail: number) {
    const d = Math.round(detail);
    let s = sphereCache.get(d);
    if (!s) { s = Sphere(d); sphereCache.set(d, s); }
    return s;
}

/** Props that affect the base molecular surface (i.e. everything except L and reconstruction detail). */
function surfaceCacheKey(id: string, props: SphericalHarmonicSurfaceMeshProps) {
    return [
        id, props.fitSource, props.probeRadius, props.resolution, props.probePositions,
        props.ignoreHydrogens, props.ignoreHydrogensVariant, props.traceOnly,
        props.includeParent, props.floodfill,
    ].join('|');
}

const MinRadius = 1e-3;

/**
 * Reconstruct the SH surface from the cached base surface stored in `meta`.
 * Always refits at the current `sphericalHarmonicL` (cheap relative to the base
 * surface computation), so the result is independent of slider history.
 */
function reconstructSphericalHarmonicMesh(meta: SphericalHarmonicSurfaceMeta, props: SphericalHarmonicSurfaceMeshProps, ctx: VisualContext, mesh?: Mesh): Mesh {
    const { shVertices, shGroups, shLookup, shCenter, shBoundingSphere } = meta;
    if (!shVertices || !shGroups || !shLookup || !shCenter || !shBoundingSphere || shVertices.length === 0) {
        return Mesh.createEmpty(mesh);
    }

    const L = Math.round(props.sphericalHarmonicL);
    const K = shTermCount(L);
    const { coeffs } = fitSphericalHarmonics(shVertices, shCenter, L);

    const sphere = getSphere(props.reconstructionDetail);
    const sv = sphere.vertices;
    const vertexCount = sv.length / 3;
    const indices = sphere.indices instanceof Uint32Array ? sphere.indices : Uint32Array.from(sphere.indices);
    const triangleCount = indices.length / 3;

    const vertices = new Float32Array(vertexCount * 3);
    const groups = new Float32Array(vertexCount);
    const normals = new Float32Array(vertexCount * 3);

    const basisScratch = new Float64Array(K);
    const legendreScratch = new Float64Array(((L + 1) * (L + 2)) / 2);
    const cx = shCenter[0], cy = shCenter[1], cz = shCenter[2];

    for (let i = 0; i < vertexCount; ++i) {
        // sphere vertices are unit length (radius 1), so they are direction vectors
        const dx = sv[i * 3], dy = sv[i * 3 + 1], dz = sv[i * 3 + 2];
        const theta = Math.acos(Math.min(1, Math.max(-1, dz)));
        const phi = Math.atan2(dy, dx);
        let r = reconstructRadius(coeffs, L, theta, phi, basisScratch, legendreScratch);
        if (r < MinRadius) r = MinRadius;

        const px = cx + r * dx, py = cy + r * dy, pz = cz + r * dz;
        vertices[i * 3] = px; vertices[i * 3 + 1] = py; vertices[i * 3 + 2] = pz;

        // transfer the group (element id) of the nearest original surface vertex
        const res = shLookup.nearest(px, py, pz, 1);
        groups[i] = res.count > 0 ? shGroups[res.indices[0]] : 0;
    }

    const surface = Mesh.create(vertices, indices, normals, groups, vertexCount, triangleCount, mesh);
    Mesh.computeNormals(surface);

    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    surface.setBoundingSphere(shBoundingSphere);
    const newMeta = surface.meta as SphericalHarmonicSurfaceMeta;
    // preserve cache fields across the rebuild (Mesh.create reuses meta when reusing the mesh)
    if (newMeta !== meta) Object.assign(newMeta, meta);
    return surface;
}

/**
 * Populate the base cache in `meta` from an explicit point cloud (interleaved xyz) and per-point
 * group ids. Used by both the surface path (marching-cubes vertices) and the atom-cloud path.
 */
function populateCacheFromPoints(meta: SphericalHarmonicSurfaceMeta, points: Float32Array, groups: Float32Array, center: Vec3, boundingSphere: Sphere3D, resolution: number, cacheKey: string) {
    const n = points.length / 3;
    const x = new Float32Array(n);
    const y = new Float32Array(n);
    const z = new Float32Array(n);
    for (let i = 0; i < n; ++i) {
        x[i] = points[i * 3];
        y[i] = points[i * 3 + 1];
        z[i] = points[i * 3 + 2];
    }
    const position = { x, y, z, indices: OrderedSet.ofRange(0, n) };
    const shLookup = n > 0 ? GridLookup3D(position, getBoundary(position)) : undefined;

    meta.shCacheKey = cacheKey;
    meta.shCenter = Vec3.clone(center);
    meta.shVertices = points;
    meta.shGroups = groups;
    meta.shLookup = shLookup;
    meta.shBoundingSphere = Sphere3D.clone(boundingSphere);
    meta.resolution = resolution;
}

/** Populate the base-surface cache in `surface.meta` from a freshly computed molecular surface. */
function populateCache(surface: Mesh, center: Vec3, boundingSphere: Sphere3D, resolution: number, cacheKey: string) {
    const vertexCount = surface.vertexCount;
    const points = surface.vertexBuffer.ref.value.slice(0, vertexCount * 3);
    const groups = surface.groupBuffer.ref.value.slice(0, vertexCount);
    populateCacheFromPoints(surface.meta as SphericalHarmonicSurfaceMeta, points, groups, center, boundingSphere, resolution, cacheKey);
}

/**
 * Build an SH-fitting point cloud from atom positions: each atom is pushed outward by its radius
 * (size + probe) along the centroid->atom direction so the fit tracks the surface envelope rather
 * than the (inward-biased) atom centers. `ids[i]` is the group id carried to the reconstructed mesh.
 */
function buildAtomCloud(xs: number[], ys: number[], zs: number[], rs: number[], ids: number[]) {
    const n = xs.length;
    const position = { x: xs, y: ys, z: zs, radius: rs, indices: OrderedSet.ofRange(0, n) };
    const boundary = getBoundary(position);
    const center = boundary.sphere.center;
    const cx = center[0], cy = center[1], cz = center[2];

    const points = new Float32Array(n * 3);
    const groups = new Float32Array(n);
    let maxRadius = 0;
    for (let i = 0; i < n; ++i) {
        const dx = xs[i] - cx, dy = ys[i] - cy, dz = zs[i] - cz;
        const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
        const r = rs[i];
        if (r > maxRadius) maxRadius = r;
        if (len > 1e-3) {
            const push = r / len; // push outward by the atom radius
            points[i * 3] = xs[i] + dx * push;
            points[i * 3 + 1] = ys[i] + dy * push;
            points[i * 3 + 2] = zs[i] + dz * push;
        } else {
            points[i * 3] = xs[i]; points[i * 3 + 1] = ys[i]; points[i * 3 + 2] = zs[i];
        }
        groups[i] = ids[i];
    }
    const boundingSphere = Sphere3D.expand(Sphere3D(), boundary.sphere, maxRadius);
    return { points, groups, center: Vec3.clone(center), boundingSphere, maxRadius };
}

/** Reconstruct an SH blob directly from an atom point cloud (no molecular surface). */
function reconstructFromAtomCloud(cloud: ReturnType<typeof buildAtomCloud>, props: SphericalHarmonicSurfaceMeshProps, resolution: number, cacheKey: string, ctx: VisualContext, mesh?: Mesh): Mesh {
    const target = mesh ?? Mesh.createEmpty();
    populateCacheFromPoints(target.meta as SphericalHarmonicSurfaceMeta, cloud.points, cloud.groups, cloud.center, cloud.boundingSphere, resolution, cacheKey);
    return reconstructSphericalHarmonicMesh(target.meta as SphericalHarmonicSurfaceMeta, props, ctx, target);
}

/** Atom arrays (model coords) for one unit, with group id = index into unit.elements. */
function getUnitAtomArrays(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: SphericalHarmonicSurfaceMeshProps) {
    const { probeRadius, ignoreHydrogens, ignoreHydrogensVariant, traceOnly } = props;
    const { x, y, z } = getConformation(unit);
    const { elements } = unit;
    const l = StructureElement.Location.create(structure, unit);
    const xs: number[] = [], ys: number[] = [], zs: number[] = [], rs: number[] = [], ids: number[] = [];
    for (let j = 0, jl = elements.length; j < jl; ++j) {
        const eI = elements[j];
        if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
        if (traceOnly && !isTrace(unit, eI)) continue;
        xs.push(x[eI]); ys.push(y[eI]); zs.push(z[eI]);
        l.element = eI;
        rs.push(sizeTheme.size(l) + probeRadius);
        ids.push(j);
    }
    return { xs, ys, zs, rs, ids };
}

/** Atom arrays (world coords) for the whole structure, with group id = serial element index. */
function getStructureAtomArrays(structure: Structure, sizeTheme: SizeTheme<any>, props: SphericalHarmonicSurfaceMeshProps) {
    const { probeRadius, ignoreHydrogens, ignoreHydrogensVariant, traceOnly } = props;
    const { getSerialIndex } = structure.serialMapping;
    const l = StructureElement.Location.create(structure);
    const xs: number[] = [], ys: number[] = [], zs: number[] = [], rs: number[] = [], ids: number[] = [];
    for (const unit of structure.units) {
        const { elements, conformation: c } = unit;
        l.unit = unit;
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const eI = elements[j];
            if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
            if (traceOnly && !isTrace(unit, eI)) continue;
            xs.push(c.x(eI)); ys.push(c.y(eI)); zs.push(c.z(eI));
            l.element = eI;
            rs.push(sizeTheme.size(l) + probeRadius);
            ids.push(getSerialIndex(unit, eI));
        }
    }
    return { xs, ys, zs, rs, ids };
}

//

async function createSphericalHarmonicSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const cacheKey = surfaceCacheKey(`unit-${unit.id}`, props);
    const meta = mesh?.meta as SphericalHarmonicSurfaceMeta | undefined;
    if (mesh && meta && meta.shCacheKey === cacheKey && meta.shVertices) {
        return reconstructSphericalHarmonicMesh(meta, props, ctx, mesh);
    }

    if (props.fitSource === 'atoms') {
        const { xs, ys, zs, rs, ids } = getUnitAtomArrays(structure, unit, theme.size, props);
        const cloud = buildAtomCloud(xs, ys, zs, rs, ids);
        return reconstructFromAtomCloud(cloud, props, props.resolution, cacheKey, ctx, mesh);
    }

    const { transform, field, idField, resolution, maxRadius } = await computeUnitMolecularSurface(structure, unit, theme.size, props).runInContext(ctx.runtime);
    const params = {
        isoLevel: props.probeRadius,
        scalarField: props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
        idField,
    };
    const base = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);
    if (props.includeParent) {
        const iterations = Math.ceil(2 / props.resolution);
        Mesh.smoothEdges(base, { iterations, maxNewEdgeLength: Math.sqrt(2) });
    }
    Mesh.transform(base, transform);

    const boundingSphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius);
    populateCache(base, unit.boundary.sphere.center as Vec3, boundingSphere, resolution, cacheKey);
    return reconstructSphericalHarmonicMesh(base.meta as SphericalHarmonicSurfaceMeta, props, ctx, base);
}

export function SphericalHarmonicSurfaceMeshVisual(materialId: number): UnitsVisual<SphericalHarmonicSurfaceMeshParams> {
    return UnitsMeshVisual<SphericalHarmonicSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(SphericalHarmonicSurfaceMeshParams),
        createGeometry: createSphericalHarmonicSurfaceMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) => {
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.probeRadius !== currentProps.probeRadius ||
                newProps.probePositions !== currentProps.probePositions ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill ||
                newProps.sphericalHarmonicL !== currentProps.sphericalHarmonicL ||
                newProps.reconstructionDetail !== currentProps.reconstructionDetail ||
                newProps.fitSource !== currentProps.fitSource
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        processValues: (values: MeshValues, geometry: Mesh, props: SphericalHarmonicSurfaceMeshProps, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as SphericalHarmonicSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture?.destroy();
        },
    }, materialId);
}

//

async function createStructureSphericalHarmonicSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const cacheKey = surfaceCacheKey(`structure-${structure.hashCode}`, props);
    const meta = mesh?.meta as SphericalHarmonicSurfaceMeta | undefined;
    if (mesh && meta && meta.shCacheKey === cacheKey && meta.shVertices) {
        return reconstructSphericalHarmonicMesh(meta, props, ctx, mesh);
    }

    if (props.fitSource === 'atoms') {
        const { xs, ys, zs, rs, ids } = getStructureAtomArrays(structure, theme.size, props);
        const cloud = buildAtomCloud(xs, ys, zs, rs, ids);
        return reconstructFromAtomCloud(cloud, props, props.resolution, cacheKey, ctx, mesh);
    }

    const { transform, field, idField, resolution, maxRadius } = await computeStructureMolecularSurface(structure, theme.size, props).runInContext(ctx.runtime);
    const params = {
        isoLevel: props.probeRadius,
        scalarField: props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
        idField,
    };
    const base = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);
    if (props.includeParent) {
        const iterations = Math.ceil(2 / props.resolution);
        Mesh.smoothEdges(base, { iterations, maxNewEdgeLength: Math.sqrt(2) });
    }
    Mesh.transform(base, transform);

    const boundingSphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius);
    populateCache(base, structure.boundary.sphere.center as Vec3, boundingSphere, resolution, cacheKey);
    return reconstructSphericalHarmonicMesh(base.meta as SphericalHarmonicSurfaceMeta, props, ctx, base);
}

export function StructureSphericalHarmonicSurfaceMeshVisual(materialId: number): ComplexVisual<SphericalHarmonicSurfaceMeshParams> {
    return ComplexMeshVisual<SphericalHarmonicSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(SphericalHarmonicSurfaceMeshParams),
        createGeometry: createStructureSphericalHarmonicSurfaceMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) => {
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.probeRadius !== currentProps.probeRadius ||
                newProps.probePositions !== currentProps.probePositions ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.includeParent !== currentProps.includeParent ||
                newProps.floodfill !== currentProps.floodfill ||
                newProps.sphericalHarmonicL !== currentProps.sphericalHarmonicL ||
                newProps.reconstructionDetail !== currentProps.reconstructionDetail ||
                newProps.fitSource !== currentProps.fitSource
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        processValues: (values: MeshValues, geometry: Mesh, props: SphericalHarmonicSurfaceMeshProps, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as SphericalHarmonicSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture?.destroy();
        },
    }, materialId);
}

//
// Protomer variant: merge the chains of one asymmetric unit (in model coordinates),
// fit a single SH blob, then GPU-instance it across the assembly symmetry operators.
//

interface ProtomerMapping {
    /** representative unit of each symmetry group (one per chain), defining the ASU */
    asuUnits: Unit[]
    /** cumulative element counts over asuUnits, length asuUnits.length + 1 */
    offsets: Int32Array
    /** total serial element count over the ASU */
    elementCount: number
    /** unit.id -> index into asuUnits */
    unitToAsu: Map<number, number>
    /** assembly operators as a flat instanceCount * 16 matrix array */
    operators: Float32Array
    instanceCount: number
}

const protomerMappingCache = new WeakMap<Structure, ProtomerMapping>();
function getProtomerMapping(structure: Structure): ProtomerMapping {
    const cached = protomerMappingCache.get(structure);
    if (cached) return cached;

    const groups = structure.unitSymmetryGroups;
    const asuUnits = groups.map(g => g.units[0]);
    const offsets = new Int32Array(asuUnits.length + 1);
    const unitToAsu = new Map<number, number>();
    for (let k = 0; k < asuUnits.length; ++k) {
        offsets[k + 1] = offsets[k] + asuUnits[k].elements.length;
        unitToAsu.set(asuUnits[k].id, k);
    }
    const elementCount = offsets[asuUnits.length];

    // assembly operators: take the symmetry group with the most instances
    // (assumes a uniform assembly where all chains share the same operator set)
    let maxGroup = groups[0];
    for (const g of groups) if (g.units.length > maxGroup.units.length) maxGroup = g;
    const instanceCount = maxGroup ? maxGroup.units.length : 1;
    const operators = new Float32Array(instanceCount * 16);
    for (let i = 0; i < instanceCount; ++i) {
        Mat4.toArray(maxGroup.units[i].conformation.operator.matrix, operators, i * 16);
    }

    const mapping: ProtomerMapping = { asuUnits, offsets, elementCount, unitToAsu, operators, instanceCount };
    protomerMappingCache.set(structure, mapping);
    return mapping;
}

/** largest k with offsets[k] <= s (which asuUnit the serial element s belongs to) */
function asuUnitIndex(offsets: Int32Array, s: number): number {
    let lo = 0, hi = offsets.length - 1;
    while (lo + 1 < hi) {
        const mid = (lo + hi) >> 1;
        if (offsets[mid] <= s) lo = mid; else hi = mid;
    }
    return lo;
}

/** Gather model-space (operator-free) positions/radii for one ASU copy, with serial ids for picking. */
function getProtomerPositionDataAndMaxRadius(structure: Structure, sizeTheme: SizeTheme<any>, props: SphericalHarmonicSurfaceMeshProps) {
    const { probeRadius, ignoreHydrogens, ignoreHydrogensVariant, traceOnly } = props;
    const m = getProtomerMapping(structure);

    const xs: number[] = [];
    const ys: number[] = [];
    const zs: number[] = [];
    const rs: number[] = [];
    const id: number[] = [];

    const l = StructureElement.Location.create(structure);
    let maxRadius = 0;
    for (let k = 0; k < m.asuUnits.length; ++k) {
        const unit = m.asuUnits[k];
        const { x, y, z } = getConformation(unit); // model conformation (no operator applied)
        const { elements } = unit;
        const base = m.offsets[k];
        l.unit = unit;
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const eI = elements[j];
            if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
            if (traceOnly && !isTrace(unit, eI)) continue;
            xs.push(x[eI]); ys.push(y[eI]); zs.push(z[eI]);
            l.element = eI;
            const r = sizeTheme.size(l);
            if (r > maxRadius) maxRadius = r;
            rs.push(r + probeRadius);
            id.push(base + j); // serial index over the ASU, stable under filtering
        }
    }

    const position = { x: xs, y: ys, z: zs, radius: rs, id, indices: OrderedSet.ofRange(0, id.length) };
    const boundary = getBoundary(position);
    return { position, boundary, maxRadius };
}

async function createProtomerSphericalHarmonicSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const cacheKey = surfaceCacheKey(`protomer-${structure.hashCode}`, props);
    const meta = mesh?.meta as SphericalHarmonicSurfaceMeta | undefined;
    if (mesh && meta && meta.shCacheKey === cacheKey && meta.shVertices) {
        return reconstructSphericalHarmonicMesh(meta, props, ctx, mesh);
    }

    if (props.fitSource === 'atoms') {
        const { position } = getProtomerPositionDataAndMaxRadius(structure, theme.size, props);
        const cloud = buildAtomCloud(position.x, position.y, position.z, position.radius, position.id);
        return reconstructFromAtomCloud(cloud, props, props.resolution, cacheKey, ctx, mesh);
    }

    const { position, boundary, maxRadius } = getProtomerPositionDataAndMaxRadius(structure, theme.size, props);
    const p = ensureReasonableResolution(boundary.box, props);
    const { transform, field, idField, resolution } = await Task.create('Protomer Molecular Surface', async runtime => {
        return calcMolecularSurface(runtime, position, boundary, maxRadius, boundary.box, p);
    }).runInContext(ctx.runtime);

    const params = {
        isoLevel: props.probeRadius,
        scalarField: props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
        idField,
    };
    const base = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);
    Mesh.transform(base, transform); // grid -> model space

    const boundingSphere = Sphere3D.expand(Sphere3D(), boundary.sphere, maxRadius);
    populateCache(base, boundary.sphere.center as Vec3, boundingSphere, resolution, cacheKey);
    return reconstructSphericalHarmonicMesh(base.meta as SphericalHarmonicSurfaceMeta, props, ctx, base);
}

function protomerLocationIterator(structure: Structure): LocationIterator {
    const m = getProtomerMapping(structure);
    const location = StructureElement.Location.create(structure);
    const getLocation = (groupIndex: number) => {
        const k = asuUnitIndex(m.offsets, groupIndex);
        location.unit = m.asuUnits[k];
        location.element = m.asuUnits[k].elements[groupIndex - m.offsets[k]];
        return location;
    };
    return LocationIterator(m.elementCount, m.instanceCount, 1, getLocation);
}

function getProtomerLoci(pickingId: PickingId, structure: Structure, id: number): Loci {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const m = getProtomerMapping(structure);
        if (groupId === PickingId.Null) return Structure.Loci(structure);
        const k = asuUnitIndex(m.offsets, groupId);
        const unit = m.asuUnits[k];
        const idx = (groupId - m.offsets[k]) as StructureElement.UnitIndex;
        return StructureElement.Loci(structure, [{ unit, indices: OrderedSet.ofSingleton(idx) }]);
    }
    return EmptyLoci;
}

function eachProtomerLocation(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean): boolean {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    const m = getProtomerMapping(structure);
    for (const e of loci.elements) {
        const k = m.unitToAsu.get(e.unit.id);
        if (k === undefined) continue;
        const base = m.offsets[k];
        if (Interval.is(e.indices)) {
            const start = base + Interval.start(e.indices);
            const end = base + Interval.end(e.indices);
            if (apply(Interval.ofBounds(start, end))) changed = true;
        } else {
            for (let i = 0, il = e.indices.length; i < il; ++i) {
                if (apply(Interval.ofSingleton(base + e.indices[i]))) changed = true;
            }
        }
    }
    return changed;
}

function createProtomerInstances(structure: Structure, props: SphericalHarmonicSurfaceMeshProps, geometry: Mesh, transformData?: TransformData): TransformData {
    const { operators, instanceCount } = getProtomerMapping(structure);
    return createTransform(Float32Array.from(operators), instanceCount, geometry.boundingSphere, props.cellSize, props.batchSize, transformData);
}

export function ProtomerSphericalHarmonicSurfaceMeshVisual(materialId: number): ComplexVisual<SphericalHarmonicSurfaceMeshParams> {
    return ComplexMeshVisual<SphericalHarmonicSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(SphericalHarmonicSurfaceMeshParams),
        createGeometry: createProtomerSphericalHarmonicSurfaceMesh,
        createLocationIterator: (structure: Structure) => protomerLocationIterator(structure),
        getLoci: getProtomerLoci,
        eachLocation: eachProtomerLocation,
        createInstances: createProtomerInstances,
        setUpdateState: (state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) => {
            state.createGeometry = (
                newProps.resolution !== currentProps.resolution ||
                newProps.probeRadius !== currentProps.probeRadius ||
                newProps.probePositions !== currentProps.probePositions ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.floodfill !== currentProps.floodfill ||
                newProps.sphericalHarmonicL !== currentProps.sphericalHarmonicL ||
                newProps.reconstructionDetail !== currentProps.reconstructionDetail ||
                newProps.fitSource !== currentProps.fitSource
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        processValues: (values: MeshValues, geometry: Mesh, props: SphericalHarmonicSurfaceMeshProps, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as SphericalHarmonicSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture?.destroy();
        },
    }, materialId);
}

//
// Per-residue blobs: one small low-degree SH blob per residue, baked into a single mesh.
// A residue has too few atoms for high L, so this is an ellipsoidal-blob mode (residueL clamped low).
//

function clampInt(v: number, min: number, max: number) {
    return Math.max(min, Math.min(max, Math.round(v)));
}

async function createResidueSphericalHarmonicSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const { probeRadius, ignoreHydrogens, ignoreHydrogensVariant, traceOnly } = props;
    const L = clampInt(props.residueL, 1, 4);
    const sphere = getSphere(clampInt(props.residueDetail, 1, 3));
    const sv = sphere.vertices;
    const svCount = sv.length / 3;
    const sIndices = sphere.indices;

    const { getSerialIndex } = structure.serialMapping;
    const l = StructureElement.Location.create(structure);

    const vertices: number[] = [];
    const groups: number[] = [];
    const indices: number[] = [];

    const basisScratch = new Float64Array(shTermCount(L));
    const legendreScratch = new Float64Array(((L + 1) * (L + 2)) / 2);

    for (let u = 0; u < structure.units.length; ++u) {
        const unit = structure.units[u];
        if (!Unit.isAtomic(unit)) continue;
        if (ctx.runtime.shouldUpdate) await ctx.runtime.update({ message: 'Residue spherical harmonic blobs', current: u, max: structure.units.length });

        const { elements, conformation: c } = unit;
        const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);
        l.unit = unit;
        while (residueIt.hasNext) {
            const seg = residueIt.move();
            const xs: number[] = [], ys: number[] = [], zs: number[] = [], rs: number[] = [];
            let repSerial = -1;
            for (let j = seg.start; j < seg.end; ++j) {
                const eI = elements[j];
                if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
                if (traceOnly && !isTrace(unit, eI)) continue;
                if (repSerial < 0) repSerial = getSerialIndex(unit, eI);
                xs.push(c.x(eI)); ys.push(c.y(eI)); zs.push(c.z(eI));
                l.element = eI;
                rs.push(theme.size.size(l) + probeRadius);
            }
            if (xs.length === 0 || repSerial < 0) continue;

            const cloud = buildAtomCloud(xs, ys, zs, rs, xs.map(() => 0));
            const cx = cloud.center[0], cy = cloud.center[1], cz = cloud.center[2];

            // angular spread of the (envelope-pushed) atoms about the centroid; monatomic /
            // centroid-coincident residues (waters, ions, single-atom ligands) have none, so an
            // SH fit collapses to ~0 radius. Fall back to a plain vdW-sized sphere in that case.
            let spread = 0;
            let rConst = 0;
            for (let p = 0; p < cloud.points.length; p += 3) {
                const ddx = cloud.points[p] - cx, ddy = cloud.points[p + 1] - cy, ddz = cloud.points[p + 2] - cz;
                const d = Math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                if (d > spread) spread = d;
            }
            for (let ri = 0; ri < rs.length; ++ri) if (rs[ri] > rConst) rConst = rs[ri];
            const useFit = spread > 1e-3;
            const coeffs = useFit ? fitSphericalHarmonics(cloud.points, cloud.center, L).coeffs : undefined;

            const base = vertices.length / 3;
            for (let k = 0; k < svCount; ++k) {
                const dx = sv[k * 3], dy = sv[k * 3 + 1], dz = sv[k * 3 + 2];
                const theta = Math.acos(Math.min(1, Math.max(-1, dz)));
                const phi = Math.atan2(dy, dx);
                let r = useFit ? reconstructRadius(coeffs!, L, theta, phi, basisScratch, legendreScratch) : rConst;
                if (r < MinRadius) r = MinRadius;
                vertices.push(cx + r * dx, cy + r * dy, cz + r * dz);
                groups.push(repSerial);
            }
            for (let k = 0; k < sIndices.length; ++k) indices.push(base + sIndices[k]);
        }
    }

    const vb = new Float32Array(vertices);
    const ib = new Uint32Array(indices);
    const nb = new Float32Array(vb.length);
    const gb = new Float32Array(groups);
    const surface = Mesh.create(vb, ib, nb, gb, vb.length / 3, ib.length / 3, mesh);
    Mesh.computeNormals(surface);

    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    surface.setBoundingSphere(Sphere3D.clone(structure.boundary.sphere));
    (surface.meta as SphericalHarmonicSurfaceMeta).resolution = props.resolution;
    return surface;
}

export function ResidueSphericalHarmonicSurfaceMeshVisual(materialId: number): ComplexVisual<SphericalHarmonicSurfaceMeshParams> {
    return ComplexMeshVisual<SphericalHarmonicSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(SphericalHarmonicSurfaceMeshParams),
        createGeometry: createResidueSphericalHarmonicSurfaceMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) => {
            state.createGeometry = (
                newProps.probeRadius !== currentProps.probeRadius ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.traceOnly !== currentProps.traceOnly ||
                newProps.residueL !== currentProps.residueL ||
                newProps.residueDetail !== currentProps.residueDetail
            );

            if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
                state.updateColor = true;
            } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
                if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
                if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
            }
        },
        processValues: (values: MeshValues, geometry: Mesh, props: SphericalHarmonicSurfaceMeshProps, theme: Theme, webgl?: WebGLContext) => {
            const { resolution, colorTexture } = geometry.meta as SphericalHarmonicSurfaceMeta;
            const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
            if (csp) {
                applyMeshColorSmoothing(values, csp, webgl, colorTexture);
                (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
            }
        },
        dispose: (geometry: Mesh) => {
            (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture?.destroy();
        },
    }, materialId);
}
