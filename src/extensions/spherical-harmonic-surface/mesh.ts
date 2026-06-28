/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 *
 * Spherical-harmonic surface visual: fits a real spherical-harmonic expansion of
 * the radial function r(theta, phi) about the centroid of an atom cloud (each atom
 * pushed outward by its size + radiusOffset), and reconstructs a smooth closed mesh
 * from the fitted coefficients. The parameter `sphericalHarmonicL` controls the
 * level of detail.
 *
 * The fit is to atom positions, not to a computed molecular surface: at the low
 * degrees used here (L ~ 6-15) the truncated single-center expansion is low-pass,
 * well below the ~probe-radius (1.4 A) detail that distinguishes the solvent-excluded
 * surface from a vdW+probe atom envelope, so fitting the molecular surface yields the
 * same low-order coefficients at far higher cost. Surface-based SH / 3D-Zernike is
 * only needed for shape-complementarity docking or property mapping onto the real
 * surface, not for a coarse shape envelope (Ritchie & Kemp 1999, J Comput Chem
 * 20:383; Duncan & Olson 1993, Biopolymers 33:219).
 *
 * Note: a single-center radial expansion can only represent star-convex surfaces
 * (those single-valued in r from the centroid). Pockets, concavities and
 * multi-domain shapes fold back on a radial ray and collapse to the outermost
 * crossing, so a single lobe is best understood as a smooth surface *envelope*
 * rather than a faithful molecular surface. To handle non-star-shaped inputs the
 * cloud can be decomposed into several star-shaped lobes (the `lobes` param: by
 * contiguous residue sequence per chain, by spatial k-means clusters, or auto - a target
 * lobe size that sets the per-chain count for uniform blobs). Each lobe is
 * reconstructed as an analytic icosphere blob and the lobes are combined either as
 * overlapping blobs (the default, no marching cubes) or, with `watertight` on,
 * blended into one watertight surface through a marching-cubes implicit field.
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../../mol-repr/structure/units-visual';
import { VisualContext } from '../../mol-repr/visual';
import { Unit, Structure, StructureElement, ElementIndex } from '../../mol-model/structure';
import { Theme } from '../../mol-theme/theme';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '../../mol-geo/util/marching-cubes/algorithm';
import { ElementIterator, getElementLoci, eachElement, getSerialElementLoci, eachSerialElement } from '../../mol-repr/structure/visual/util/element';
import { VisualUpdateState } from '../../mol-repr/util';
import { getConformation, isHydrogen, isTrace } from '../../mol-repr/structure/visual/util/common';
import { Sphere3D } from '../../mol-math/geometry';
import { MeshValues } from '../../mol-gl/renderable/mesh';
import { Texture } from '../../mol-gl/webgl/texture';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { applyMeshColorSmoothing } from '../../mol-geo/geometry/mesh/color-smoothing';
import { BaseGeometry, ColorSmoothingParams, getColorSmoothingProps } from '../../mol-geo/geometry/base';
import { ValueCell } from '../../mol-util';
import { ComplexMeshVisual, ComplexVisual } from '../../mol-repr/structure/complex-visual';
import { Tensor } from '../../mol-math/linear-algebra/tensor';
import { Mat3, Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Sphere } from '../../mol-geo/primitive/sphere';
import { getBoundary } from '../../mol-math/geometry/boundary';
import { OrderedSet } from '../../mol-data/int';
import { fitSphericalHarmonicLobesByLabel, kmeansLabels, reconstructRadius, shTermCount, SphericalHarmonicLobe } from '../../mol-math/geometry/spherical-harmonics';
import { SizeTheme } from '../../mol-theme/size';
import { isTimingMode } from '../../mol-util/debug';
import { now } from '../../mol-util/now';

/**
 * Aggregate timing for the build path. `UnitsVisual` builds one geometry per unique unit/group, so
 * per-unit logs would be unreadable on a large assembly; instead accumulate across the burst of
 * `createGeometry` calls and flush one summary line shortly after the build settles. The unit count
 * is the decisive number (few instanced groups vs. thousands of distinct chains). isTimingMode-gated,
 * so it is inert in normal use.
 */
const shBuildTiming = { units: 0, points: 0, ms: 0, timer: undefined as ReturnType<typeof setTimeout> | undefined };
function shAccumulateTiming(start: number, points: number) {
    if (!isTimingMode) return;
    shBuildTiming.units += 1;
    shBuildTiming.points += points;
    shBuildTiming.ms += now() - start;
    if (shBuildTiming.timer !== undefined) clearTimeout(shBuildTiming.timer);
    shBuildTiming.timer = setTimeout(() => {
        console.log(`[SH surface] built ${shBuildTiming.units} unit(s), ${shBuildTiming.points} fit pts total in ${shBuildTiming.ms.toFixed(1)}ms (${(shBuildTiming.ms / shBuildTiming.units).toFixed(2)}ms/unit)`);
        shBuildTiming.units = 0; shBuildTiming.points = 0; shBuildTiming.ms = 0; shBuildTiming.timer = undefined;
    }, 300);
}

/** Default target atoms per auto lobe: most chains fall under this (one lobe), large chains (e.g. rRNA) split. */
export const DefaultTargetAtoms = 3000;
/** Fallback slider max for `auto.targetAtoms` when no structure is available; the representation overrides it with the biggest chain's atom count. */
const DefaultMaxTargetAtoms = 50000;

/**
 * The `lobes` param: how to split a chain/structure into star-convex lobes before fitting. `auto`
 * (the default) sizes the per-chain lobe count from a target atoms-per-lobe, so small chains stay a
 * single envelope and only large chains split. `maxTargetAtoms` bounds the target slider - the
 * representation passes the biggest chain's atom count so the range matches the loaded data.
 */
export function LobesParam(maxTargetAtoms: number) {
    return PD.MappedStatic('auto', {
        'single': PD.EmptyGroup(),
        'sequence': PD.Group({
            divisions: PD.Numeric(4, { min: 1, max: 50, step: 1 }, { description: 'Divide each chain into this many contiguous residue-sequence segments, fitting one star-convex lobe per segment. Higher values follow a threaded chain (e.g. rRNA) more closely; 1 is a single envelope.' }),
        }, { isFlat: true }),
        'kmeans': PD.Group({
            clusters: PD.Numeric(4, { min: 1, max: 32, step: 1 }, { description: 'Cluster atoms into this many spatial groups (k-means), fitting one star-convex lobe per cluster. Compact lobes that follow a threaded chain (e.g. rRNA) at far fewer pieces than sequence division.' }),
        }, { isFlat: true }),
        'auto': PD.Group({
            targetAtoms: PD.Numeric(DefaultTargetAtoms, { min: DefaultTargetAtoms, max: Math.max(DefaultTargetAtoms, maxTargetAtoms), step: 50 }, { description: 'Target atoms per lobe (never below 3000). Each chain is split (by k-means) into round(chainAtoms / targetAtoms) lobes, so lobes are uniform across chains: chains up to ~1.5x this stay a single envelope, larger chains (e.g. rRNA) split into k-means lobes. Max is the biggest chain in the loaded structure.' }),
        }, { isFlat: true }),
    }, { description: 'Split a chain/structure into star-convex lobes before fitting: auto (default - uniform lobe size, per-chain count from a target; single envelope for small chains), a single radial envelope, by residue sequence (per chain), or by a fixed number of spatial k-means clusters.' });
}

export const SphericalHarmonicSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...ColorSmoothingParams,
    ignoreHydrogens: PD.Boolean(false, { description: 'Exclude hydrogen atoms from the shape envelope.' }),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    traceOnly: PD.Boolean(false, { description: 'Use only trace atoms (e.g. CA, P) for the shape envelope.' }),
    radiusOffset: PD.Numeric(5, { min: 0, max: 10, step: 0.1 }, { description: 'Extra radius added to each atom when building the shape envelope; larger gives a thicker, rounder surface (analogous to a probe radius).' }),
    sphericalHarmonicL: PD.Numeric(8, { min: 2, max: 20, step: 1 }, { description: 'Maximum spherical harmonic degree. Higher values follow the shape envelope more closely (but cannot recover concavities: the fit is single-valued in radius about the center).' }),
    reconstructionDetail: PD.Numeric(3, { min: 1, max: 5, step: 1 }, { description: 'Triangulation detail of the reconstructed sphere mesh.' }),
    lobes: LobesParam(DefaultMaxTargetAtoms),
    watertight: PD.Boolean(false, { description: 'Blend the lobes into one watertight surface via marching cubes. Off: overlapping star-convex blobs (no grid / no marching cubes).' }),
    regularization: PD.Numeric(0.01, { min: 0, max: 0.5, step: 0.005 }, { description: 'Tikhonov smoothness prior on the spherical-harmonic fit (per-band l(l+1) damping). Higher values give a rounder, more stable surface; lower values follow the samples more closely but can ring or blow up on sparse/clustered clouds (e.g. trace-only). 0 disables it.' }),
    resolution: PD.Numeric(0.5, { min: 0.1, max: 20, step: 0.1 }, { description: 'Grid spacing for the watertight multi-lobe blend field (only used when watertight blend is on with more than one lobe).', ...BaseGeometry.CustomQualityParamInfo }),
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
    // per-cloud-point chain ordinal and normalized residue position within that chain ([0,1)),
    // used by the `sequence` lobe mode to label points into contiguous residue-run segments.
    // Stored independently of the division count so changing it re-fits without rebuilding the cloud.
    shSeqChain?: Int32Array
    shSeqFrac?: Float32Array
    // per-cloud-point radius (atom size + radiusOffset); the lobe fit inflates each point outward by
    // this amount about its lobe centroid so the envelope wraps the atoms without drifting off-center.
    shRadii?: Float32Array
    shBoundingSphere?: Sphere3D
    // cache of the fitted lobes, keyed by `${L}|${mode}|${k}|${reg}` so a detail- or blend-only change re-tessellates without refitting
    shFitKey?: string
    shFitLobes?: SphericalHarmonicLobe[]
    // scratch mesh holding the single ASU envelope for the assembly flavour, before it is replicated
    // across the assembly operators into the bound mesh (kept off the render object, reused).
    shAsuMesh?: Mesh
}

// memoized sphere triangulations per detail level
const sphereCache = new Map<number, ReturnType<typeof Sphere>>();
function getSphere(detail: number) {
    const d = Math.round(detail);
    let s = sphereCache.get(d);
    if (!s) { s = Sphere(d); sphereCache.set(d, s); }
    return s;
}

/** Props that affect the base atom cloud (i.e. everything except L and reconstruction detail). */
function surfaceCacheKey(id: string, props: SphericalHarmonicSurfaceMeshProps) {
    return [
        id, props.radiusOffset, props.resolution,
        props.ignoreHydrogens, props.ignoreHydrogensVariant, props.traceOnly,
    ].join('|');
}

/**
 * Cap on the number of atom-cloud points gathered for the fit and the group transfer. A smooth
 * degree-L envelope is determined by far fewer points (the fit strides to 8192 internally) and the
 * per-vertex group transfer is over a coarse envelope anyway, so for large structures the cloud is
 * uniformly strided: building the point arrays and running the brute-force group transfer over every
 * atom (millions for e.g. 3J3Q) would be wasteful. Uniform stride preserves the overall shape.
 */
const MaxFitCloudPoints = 50000;

const MinRadius = 1e-3;
/** Reconstruction radius is clamped to `rMax * RadiusMargin` so an ill-conditioned fit cannot send vertices far outside the data (which would blow up the bounding sphere and cull the surface). */
const RadiusMargin = 1.15;

/** True if any prop affecting the geometry changed (everything but pure shading/color props). */
function surfaceGeometryChanged(newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps): boolean {
    return (
        newProps.resolution !== currentProps.resolution ||
        newProps.radiusOffset !== currentProps.radiusOffset ||
        newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
        newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
        newProps.traceOnly !== currentProps.traceOnly ||
        newProps.sphericalHarmonicL !== currentProps.sphericalHarmonicL ||
        newProps.reconstructionDetail !== currentProps.reconstructionDetail ||
        lobesChanged(newProps.lobes, currentProps.lobes) ||
        newProps.watertight !== currentProps.watertight ||
        newProps.regularization !== currentProps.regularization
    );
}

type LobesMode = 'single' | 'sequence' | 'kmeans' | 'auto';

/** Lobe config: the splitting mode, the fixed lobe count `k` (sequence/kmeans), and the auto target atoms-per-lobe. */
function lobesConfig(lobes: SphericalHarmonicSurfaceMeshProps['lobes']): { mode: LobesMode, k: number, target: number } {
    if (lobes.name === 'sequence') return { mode: 'sequence', k: Math.max(1, Math.round(lobes.params.divisions)), target: 0 };
    if (lobes.name === 'kmeans') return { mode: 'kmeans', k: Math.max(1, Math.round(lobes.params.clusters)), target: 0 };
    if (lobes.name === 'auto') return { mode: 'auto', k: 0, target: Math.max(1, Math.round(lobes.params.targetAtoms)) };
    return { mode: 'single', k: 1, target: 0 };
}

/** True if the lobe mode or its active count/target changed (MappedStatic identity is not stable). */
function lobesChanged(a: SphericalHarmonicSurfaceMeshProps['lobes'], b: SphericalHarmonicSurfaceMeshProps['lobes']): boolean {
    const ca = lobesConfig(a), cb = lobesConfig(b);
    return ca.mode !== cb.mode || ca.k !== cb.k || ca.target !== cb.target;
}

/** Shared color-smoothing branch of `setUpdateState` (identical for every flavour). */
function updateSmoothColorsState(state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) {
    if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
        state.updateColor = true;
    } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
        if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor) state.updateColor = true;
        if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride) state.updateColor = true;
    }
}

/** Shared `processValues`: applies volumetric color smoothing and caches the color texture on the mesh meta. */
function shProcessValues(values: MeshValues, geometry: Mesh, props: SphericalHarmonicSurfaceMeshProps, theme: Theme, webgl?: WebGLContext) {
    const meta = geometry.meta as SphericalHarmonicSurfaceMeta;
    const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, meta.resolution);
    if (csp) {
        applyMeshColorSmoothing(values, csp, webgl, meta.colorTexture);
        meta.colorTexture = values.tColorGrid.ref.value;
    }
}

/** Shared `dispose`: releases the cached color-smoothing texture. */
function shDispose(geometry: Mesh) {
    (geometry.meta as SphericalHarmonicSurfaceMeta).colorTexture?.destroy();
}

/**
 * Mean triangle edge length of the reconstructed mesh, used as the color-smoothing grid spacing.
 *
 * The reconstructed SH surface is far coarser than the molecular-surface grid it was fit from (a
 * detail-3 icosphere has ~640 vertices spread over the whole envelope, vs. one marching-cubes vertex
 * per `resolution` Angstrom). Reusing that fine base resolution as the smoothing grid spacing builds a
 * grid whose cells *between* the sparse reconstructed vertices are never splatted into, so they stay at
 * 0 (black) and trilinear interpolation bleeds that black across the surface. Matching the grid to the
 * actual vertex spacing keeps every cell near the surface covered (and, when the envelope is coarse
 * enough that the spacing exceeds the smoothing cutoff, disables smoothing entirely — also black-free).
 */
function estimateVertexSpacing(mesh: Mesh): number {
    const triangleCount = mesh.triangleCount;
    if (triangleCount === 0) return 1;
    const vb = mesh.vertexBuffer.ref.value;
    const ib = mesh.indexBuffer.ref.value;
    const sampleStride = Math.max(1, Math.floor(triangleCount / 1000));
    const a = Vec3(), b = Vec3(), c = Vec3();
    let sum = 0, n = 0;
    for (let t = 0; t < triangleCount; t += sampleStride) {
        Vec3.fromArray(a, vb, ib[t * 3] * 3);
        Vec3.fromArray(b, vb, ib[t * 3 + 1] * 3);
        Vec3.fromArray(c, vb, ib[t * 3 + 2] * 3);
        sum += Vec3.distance(a, b) + Vec3.distance(b, c) + Vec3.distance(c, a);
        n += 3;
    }
    return n > 0 ? sum / n : 1;
}

/**
 * Copy the group id of the nearest atom-cloud point onto each reconstructed vertex (the group drives
 * per-vertex coloring and picking), via an exact uniform-grid nearest search over the full cloud.
 *
 * A plain brute-force scan is O(vertices * cloud) and becomes the build floor once a chain is split
 * into many lobes (each lobe adds ~640 vertices): e.g. a 64-lobe rRNA was ~5 s of pure transfer.
 * This bins the cloud into a uniform grid (~1 atom/cell) and, per vertex, expands cell rings from the
 * vertex's cell until the nearest found atom is closer than the next ring can be - exact (returns the
 * true Euclidean nearest, 0% deviation vs brute force) but ~130x faster on the 64-lobe case. Note this
 * is a purpose-built grid with exact ring termination, not `GridLookup3D.nearest`, whose cell-order
 * heap is both slower per call and only approximate for queries starting outside the data box.
 */
function transferGroups(vertices: Float32Array, vertexCount: number, cloud: Float32Array, cloudGroups: Float32Array, out: Float32Array) {
    const pointCount = cloud.length / 3;
    if (pointCount === 0) { out.fill(0, 0, vertexCount); return; }

    // cloud bounding box
    let minX = Infinity, minY = Infinity, minZ = Infinity, maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
    for (let j = 0; j < pointCount; ++j) {
        const x = cloud[j * 3], y = cloud[j * 3 + 1], z = cloud[j * 3 + 2];
        if (x < minX) minX = x; if (y < minY) minY = y; if (z < minZ) minZ = z;
        if (x > maxX) maxX = x; if (y > maxY) maxY = y; if (z > maxZ) maxZ = z;
    }
    // cell size ~ one atom per cell; clamped so degenerate clouds don't make a huge grid
    const h = Math.max(2, Math.cbrt(((maxX - minX + 1) * (maxY - minY + 1) * (maxZ - minZ + 1)) / pointCount));
    const nx = Math.max(1, Math.floor((maxX - minX) / h) + 1);
    const ny = Math.max(1, Math.floor((maxY - minY) / h) + 1);
    const nz = Math.max(1, Math.floor((maxZ - minZ) / h) + 1);
    const cellIndex = (x: number, y: number, z: number) => {
        let cx = Math.floor((x - minX) / h); cx = cx < 0 ? 0 : cx >= nx ? nx - 1 : cx;
        let cy = Math.floor((y - minY) / h); cy = cy < 0 ? 0 : cy >= ny ? ny - 1 : cy;
        let cz = Math.floor((z - minZ) / h); cz = cz < 0 ? 0 : cz >= nz ? nz - 1 : cz;
        return (cx * ny + cy) * nz + cz;
    };

    // CSR bucketing of atoms into cells
    const cellOf = new Int32Array(pointCount);
    const start = new Int32Array(nx * ny * nz + 1);
    for (let j = 0; j < pointCount; ++j) { const c = cellIndex(cloud[j * 3], cloud[j * 3 + 1], cloud[j * 3 + 2]); cellOf[j] = c; start[c + 1]++; }
    for (let c = 0; c < nx * ny * nz; ++c) start[c + 1] += start[c];
    const items = new Int32Array(pointCount);
    const cursor = start.slice(0, nx * ny * nz);
    for (let j = 0; j < pointCount; ++j) { const c = cellOf[j]; items[cursor[c]++] = j; }

    const maxRing = Math.max(nx, ny, nz);
    for (let i = 0; i < vertexCount; ++i) {
        const vx = vertices[i * 3], vy = vertices[i * 3 + 1], vz = vertices[i * 3 + 2];
        let cx = Math.floor((vx - minX) / h); cx = cx < 0 ? 0 : cx >= nx ? nx - 1 : cx;
        let cy = Math.floor((vy - minY) / h); cy = cy < 0 ? 0 : cy >= ny ? ny - 1 : cy;
        let cz = Math.floor((vz - minZ) / h); cz = cz < 0 ? 0 : cz >= nz ? nz - 1 : cz;
        let best = Infinity, bg = 0;
        for (let ring = 0; ring <= maxRing; ++ring) {
            const xlo = Math.max(0, cx - ring), xhi = Math.min(nx - 1, cx + ring);
            const ylo = Math.max(0, cy - ring), yhi = Math.min(ny - 1, cy + ring);
            const zlo = Math.max(0, cz - ring), zhi = Math.min(nz - 1, cz + ring);
            for (let gx = xlo; gx <= xhi; ++gx) for (let gy = ylo; gy <= yhi; ++gy) for (let gz = zlo; gz <= zhi; ++gz) {
                // only the shell of Chebyshev distance == ring (inner rings already scanned)
                if (Math.abs(gx - cx) !== ring && Math.abs(gy - cy) !== ring && Math.abs(gz - cz) !== ring) continue;
                const c = (gx * ny + gy) * nz + gz;
                for (let t = start[c], e = start[c + 1]; t < e; ++t) {
                    const j = items[t];
                    const dx = cloud[j * 3] - vx, dy = cloud[j * 3 + 1] - vy, dz = cloud[j * 3 + 2] - vz;
                    const d = dx * dx + dy * dy + dz * dz;
                    if (d < best) { best = d; bg = cloudGroups[j]; }
                }
            }
            // any atom in ring+1 or beyond is at least (ring*h) away, so we can stop once that exceeds best
            const reach = ring * h;
            if (best < Infinity && reach * reach >= best) break;
        }
        out[i] = bg;
    }
}

/**
 * Star-shaped lobes → analytic icosphere blobs, one per lobe, concatenated into a single mesh
 * (no grid / no marching cubes). Each blob scales the unit-sphere directions by the lobe's
 * reconstructed radius about its center; neighbouring lobes simply overlap (not a watertight union),
 * which reads as a smooth lumpy surface following the chain - the cheap MC-free reconstruction path.
 * One lobe degrades to a single radial envelope.
 */
function reconstructBlobs(lobes: ReadonlyArray<SphericalHarmonicLobe>, L: number, props: SphericalHarmonicSurfaceMeshProps, cloud: Float32Array, shGroups: Float32Array, mesh?: Mesh): Mesh {
    const K = shTermCount(L);
    const sphere = getSphere(props.reconstructionDetail);
    const sv = sphere.vertices;
    const baseVertexCount = sv.length / 3;
    const baseIndices = sphere.indices instanceof Uint32Array ? sphere.indices : Uint32Array.from(sphere.indices);
    const baseTriangleCount = baseIndices.length / 3;
    const nL = lobes.length;

    const vertexCount = baseVertexCount * nL;
    const triangleCount = baseTriangleCount * nL;
    const vertices = new Float32Array(vertexCount * 3);
    const groups = new Float32Array(vertexCount);
    const normals = new Float32Array(vertexCount * 3);
    const indices = new Uint32Array(triangleCount * 3);

    const basisScratch = new Float64Array(K);
    const legendreScratch = new Float64Array(((L + 1) * (L + 2)) / 2);

    for (let li = 0; li < nL; ++li) {
        const lobe = lobes[li];
        const cx = lobe.center[0], cy = lobe.center[1], cz = lobe.center[2];
        const rClamp = lobe.rMax * RadiusMargin;
        const voff = li * baseVertexCount * 3;
        const ioff = li * baseTriangleCount * 3;
        const vBase = li * baseVertexCount;
        for (let i = 0; i < baseVertexCount; ++i) {
            // sphere vertices are unit length (radius 1), so they are direction vectors
            const dx = sv[i * 3], dy = sv[i * 3 + 1], dz = sv[i * 3 + 2];
            const theta = Math.acos(Math.min(1, Math.max(-1, dz)));
            const phi = Math.atan2(dy, dx);
            let r = reconstructRadius(lobe.coeffs, L, theta, phi, basisScratch, legendreScratch);
            if (!(r > MinRadius)) r = MinRadius; // also catches NaN
            else if (r > rClamp) r = rClamp;
            vertices[voff + i * 3] = cx + r * dx; vertices[voff + i * 3 + 1] = cy + r * dy; vertices[voff + i * 3 + 2] = cz + r * dz;
        }
        for (let i = 0, il = baseTriangleCount * 3; i < il; ++i) indices[ioff + i] = baseIndices[i] + vBase;
    }
    transferGroups(vertices, vertexCount, cloud, shGroups, groups);

    const surface = Mesh.create(vertices, indices, normals, groups, vertexCount, triangleCount, mesh);
    Mesh.computeNormals(surface);
    return surface;
}

const FieldPad = 1.5;
/** Cap on the multi-lobe field grid edge so cost (O(dim^3 * lobes)) stays bounded and independent of the quality/resolution slider. The smooth-max field gives smooth shading even at this coarse spacing, so it is kept low for speed (only multi-lobe chains reach this path). */
const MaxFieldDim = 64;

/**
 * Polynomial smooth-max (Inigo Quilez's smooth-min, negated): a rounded `max(a, b)` that blends only
 * within `k` of the crossover and is exact (== max) beyond it. Unlike a log-sum-exp soft-max it is
 * LOCAL - no `log(nLobes)/beta` bias, so empty regions stay outside (no spurious halo) - it just
 * rounds the seam where two lobes meet and keeps the field gradient continuous (smooth normals).
 */
function smoothMax(a: number, b: number, k: number): number {
    if (k <= 0) return a > b ? a : b;
    const h = Math.max(0, Math.min(1, 0.5 + 0.5 * (a - b) / k));
    return a * h + b * (1 - h) + k * h * (1 - h);
}

/**
 * Multiple lobes → one watertight surface. Each lobe contributes an inside-ness scalar
 * f = r_fit(dir) - |p - center| (positive inside its SH envelope); the lobes are combined by a
 * smooth-max (their rounded union) and the field is marching-cubed once at iso-level 0.
 *
 * The smooth-max rounds the seam creases where lobes meet and, by keeping the field gradient
 * continuous, gives smooth shading instead of the faceted ridges of a hard `max`. NB: a log-sum-exp
 * soft-max was tried and rejected - its `log(nLobes)/beta` bias inflates the field where lobes are
 * comparable, adding a halo of spurious zero-crossings (measured +57% vertices on an 8-lobe rRNA);
 * the polynomial smooth-max is local and halo-free.
 */
async function reconstructLobesField(lobes: SphericalHarmonicLobe[], L: number, props: SphericalHarmonicSurfaceMeshProps, cloud: Float32Array, shGroups: Float32Array, boundingSphere: Sphere3D, ctx: VisualContext, mesh?: Mesh): Promise<Mesh> {
    const R = boundingSphere.radius + FieldPad + props.resolution;
    const c = boundingSphere.center;
    const min = Vec3.create(c[0] - R, c[1] - R, c[2] - R);
    // grid spacing follows the resolution slider but the edge count is capped so cost stays bounded
    const dimN = Math.min(MaxFieldDim, Math.max(2, Math.ceil((2 * R) / props.resolution) + 1));
    const resolution = (2 * R) / (dimN - 1);
    const space = Tensor.Space([dimN, dimN, dimN], [0, 1, 2], Float32Array);
    const data = space.create();

    const K = shTermCount(L);
    const basisScratch = new Float64Array(K);
    const legendreScratch = new Float64Array(((L + 1) * (L + 2)) / 2);

    // Seam-rounding width for the smooth-max union: a few Angstrom, but at least ~2 grid cells so it
    // is actually resolved by the field (and capped so it never over-rounds small lobes).
    const blend = Math.max(2, Math.min(6, 2 * resolution));

    const nL = lobes.length;
    const rClamp = lobes.map(l => l.rMax * RadiusMargin);
    for (let i = 0; i < dimN; ++i) {
        const px = min[0] + i * resolution;
        for (let j = 0; j < dimN; ++j) {
            const py = min[1] + j * resolution;
            for (let k = 0; k < dimN; ++k) {
                const pz = min[2] + k * resolution;
                // rounded union inside-ness: F = smooth-max over lobes of (r_fit(dir) - |p - center|)
                let acc = 0; let has = false;
                for (let li = 0; li < nL; ++li) {
                    const ctr = lobes[li].center;
                    const dx = px - ctr[0], dy = py - ctr[1], dz = pz - ctr[2];
                    const rr = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    // Far outside this lobe's radius the cell is unambiguously outside it, so the
                    // exact (negative) inside-ness doesn't affect the union; skip the expensive
                    // SH evaluation and use the cheap monotonic estimate rMax - rr.
                    let f: number;
                    if (rr > rClamp[li] + blend) {
                        f = rClamp[li] - rr;
                    } else {
                        const theta = rr > 1e-9 ? Math.acos(Math.min(1, Math.max(-1, dz / rr))) : 0;
                        const phi = Math.atan2(dy, dx);
                        let rfit = reconstructRadius(lobes[li].coeffs, L, theta, phi, basisScratch, legendreScratch);
                        if (!(rfit > MinRadius)) rfit = MinRadius; // clamp ringing/NaN low so it can't punch holes
                        else if (rfit > rClamp[li]) rfit = rClamp[li];
                        f = rfit - rr;
                    }
                    acc = has ? smoothMax(acc, f, blend) : f;
                    has = true;
                }
                space.set(data, i, j, k, acc);
            }
        }
    }

    const field = Tensor.create(space, data);
    const surface = await computeMarchingCubesMesh({ isoLevel: 0, scalarField: field }, mesh).runAsChild(ctx.runtime);

    const transform = Mat4.fromScaling(Mat4(), Vec3.create(resolution, resolution, resolution));
    Mat4.setTranslation(transform, min);
    Mesh.transform(surface, transform);

    const vertexCount = surface.vertexCount;
    const vb = surface.vertexBuffer.ref.value;
    const gb = surface.groupBuffer.ref.value;
    const groups = gb.length >= vertexCount ? gb : new Float32Array(vertexCount);
    transferGroups(vb, vertexCount, cloud, shGroups, groups);
    ValueCell.update(surface.groupBuffer, groups);

    return surface;
}

/** Cap on auto lobes per chain so the per-vertex group transfer (O(lobes * verts * cloud)) stays bounded. */
const MaxAutoLobesPerChain = 64;

/**
 * Auto labels: split each chain (by `seqChain`) into round(chainAtoms / targetAtoms) k-means clusters,
 * so every lobe is about `targetAtoms` atoms - uniform blob size with a per-chain count (small chains
 * stay one lobe, large chains split into many). Labels are unique across chains.
 */
function autoLabels(verts: Float32Array, seqChain: Int32Array | undefined, targetAtoms: number): Int32Array {
    const n = verts.length / 3;
    const labels = new Int32Array(n);
    const byChain = new Map<number, number[]>();
    for (let i = 0; i < n; ++i) {
        const c = seqChain ? seqChain[i] : 0;
        let arr = byChain.get(c);
        if (!arr) { arr = []; byChain.set(c, arr); }
        arr.push(i);
    }

    let base = 0;
    for (const idxArr of byChain.values()) {
        const m = idxArr.length;
        const kChain = Math.max(1, Math.min(MaxAutoLobesPerChain, Math.round(m / targetAtoms)));
        if (kChain <= 1) {
            for (const i of idxArr) labels[i] = base;
            base += 1;
            continue;
        }
        const sub = new Float32Array(m * 3);
        for (let t = 0; t < m; ++t) {
            const i = idxArr[t];
            sub[t * 3] = verts[i * 3]; sub[t * 3 + 1] = verts[i * 3 + 1]; sub[t * 3 + 2] = verts[i * 3 + 2];
        }
        const subLabels = kmeansLabels(sub, kChain);
        for (let t = 0; t < m; ++t) labels[idxArr[t]] = base + subLabels[t];
        base += kChain;
    }
    return labels;
}

/**
 * Fit the cloud into one or more star-shaped lobes per the chosen mode:
 * - `single`: one radial envelope about the cloud centroid.
 * - `sequence`: one lobe per contiguous residue-sequence segment (k per chain), labelled from the
 *   per-point chain/residue-fraction stored on the cloud - locally star-convex pieces that follow a
 *   threaded chain (e.g. rRNA) better than a single shell, but a contiguous run is not spatially
 *   compact so it needs many divisions on a heavily threaded chain.
 * - `kmeans`: k spatial clusters - compact, star-convex lobes at far fewer pieces than sequence.
 * - `auto`: a per-chain k-means count derived from `target` atoms-per-lobe, for uniform lobe size.
 */
function fitLobes(meta: SphericalHarmonicSurfaceMeta, mode: LobesMode, k: number, target: number, L: number, maxFitPoints: number, regularization: number): SphericalHarmonicLobe[] {
    const { shVertices, shRadii, shSeqChain, shSeqFrac } = meta;
    const verts = shVertices!;
    const n = verts.length / 3;

    let labels: Int32Array;
    if (mode === 'auto') {
        labels = autoLabels(verts, shSeqChain, target);
    } else if (k > 1 && mode === 'sequence' && shSeqChain && shSeqFrac) {
        labels = new Int32Array(n);
        for (let i = 0; i < n; ++i) {
            const run = Math.min(k - 1, Math.floor(shSeqFrac[i] * k));
            labels[i] = shSeqChain[i] * k + run;
        }
    } else if (k > 1 && mode === 'kmeans') {
        labels = kmeansLabels(verts, k);
    } else {
        labels = new Int32Array(n); // single lobe
    }

    return fitSphericalHarmonicLobesByLabel(verts, labels, L, { radii: shRadii, maxPoints: maxFitPoints, regularization }).lobes;
}

/**
 * Reconstruct the SH surface from the cached base cloud stored in `meta`.
 *
 * Fits one or more star-shaped lobes (cached by mode+count, so a detail- or blend-only change
 * re-tessellates without refitting). Lobes are combined either as overlapping icosphere blobs
 * (the MC-free default) or, when `watertight` is on, blended into one watertight surface via a
 * marching-cubes implicit field.
 */
async function reconstructSphericalHarmonicMesh(meta: SphericalHarmonicSurfaceMeta, props: SphericalHarmonicSurfaceMeshProps, ctx: VisualContext, mesh?: Mesh): Promise<Mesh> {
    const { shVertices, shGroups, shCenter, shBoundingSphere } = meta;
    if (!shVertices || !shGroups || !shCenter || !shBoundingSphere || shVertices.length === 0) {
        return Mesh.createEmpty(mesh);
    }

    const L = Math.round(props.sphericalHarmonicL);
    const { mode, k, target } = lobesConfig(props.lobes);
    const regularization = props.regularization;

    // Cap the points used for the least-squares fit. The fit cost is O(points * K^2), K = (L+1)^2,
    // and it runs once per unit, so for large chains (e.g. ~1800-atom CA in 3J3Q x ~1356 chains) the
    // accumulation dominates build/update time. A smooth degree-L envelope is over-determined by far
    // fewer points (~10 * K), so striding the cloud to that cap is near-lossless (sub-0.5 A on a
    // ~20 A envelope) while cutting the fit several-fold. Scales with L so high degrees keep enough
    // samples. The full cloud is still used for the per-vertex group transfer (accurate coloring).
    const maxFitPoints = Math.min(8192, Math.max(512, shTermCount(L) * 10));

    // The lobe fit is cached by mode+count (watertight is NOT here: it only changes how the cached
    // lobes are combined, not the fit), so toggling the blend re-tessellates without refitting.
    const fitKey = `${L}|${mode}|${k}|${target}|${regularization}`;
    let lobes = meta.shFitKey === fitKey ? meta.shFitLobes : undefined;
    if (!lobes) {
        lobes = fitLobes(meta, mode, k, target, L, maxFitPoints, regularization);
    }

    const surface = lobes.length <= 1
        ? reconstructBlobs(lobes, L, props, shVertices, shGroups, mesh)
        : props.watertight
            ? await reconstructLobesField(lobes, L, props, shVertices, shGroups, shBoundingSphere, ctx, mesh)
            : reconstructBlobs(lobes, L, props, shVertices, shGroups, mesh);

    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(surface);
        ValueCell.updateIfChanged(surface.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(surface.varyingGroup, true);
    }

    // NB: do NOT set the bounding sphere from the atom boundary here. The reconstructed
    // radial surface can extend past it, and a too-small bounding sphere makes the renderer
    // cull the whole mesh (blank viewport). Let Mesh derive it from the actual vertices.
    const newMeta = surface.meta as SphericalHarmonicSurfaceMeta;
    // preserve cache fields across the rebuild (Mesh.create reuses meta when reusing the mesh)
    if (newMeta !== meta) Object.assign(newMeta, meta);
    // cache the fit for detail-only re-tessellation
    newMeta.shFitKey = fitKey;
    newMeta.shFitLobes = lobes;
    // color smoothing must use the reconstructed mesh's own vertex spacing, NOT the (fine) base
    // surface resolution carried over above - otherwise the smoothing grid is full of empty (black)
    // cells between the sparse reconstructed vertices. Set unconditionally: on the mesh-reuse path
    // newMeta === meta, so the Object.assign above is skipped and the stale base resolution remains.
    newMeta.resolution = estimateVertexSpacing(surface);
    return surface;
}

/**
 * Populate the base cache in `meta` from a built atom cloud (positions, group ids and per-point
 * sequence labels).
 */
function populateCacheFromPoints(meta: SphericalHarmonicSurfaceMeta, cloud: ReturnType<typeof buildAtomCloud>, resolution: number, cacheKey: string) {
    meta.shCacheKey = cacheKey;
    meta.shCenter = Vec3.clone(cloud.center);
    meta.shVertices = cloud.points;
    meta.shGroups = cloud.groups;
    meta.shSeqChain = cloud.chains;
    meta.shSeqFrac = cloud.fracs;
    meta.shRadii = cloud.radii;
    meta.shBoundingSphere = Sphere3D.clone(cloud.boundingSphere);
    meta.resolution = resolution;
    // the point cloud changed: invalidate the cached fit so it is recomputed
    meta.shFitKey = undefined;
    meta.shFitLobes = undefined;
}

/**
 * Build an SH-fitting point cloud from atom positions. The atoms are stored at their true positions
 * together with a per-atom radius (size + radiusOffset); the outward inflation that makes the fit
 * track the surface envelope (rather than the inward-biased atom centers) is applied per lobe at fit
 * time, *about that lobe's own centroid*, so a larger offset thickens the envelope without drifting
 * its center off the atoms. `ids[i]` is the group id carried to the reconstructed mesh.
 */
function buildAtomCloud(xs: number[], ys: number[], zs: number[], rs: number[], ids: ArrayLike<number>, seqChain: ArrayLike<number>, seqFrac: ArrayLike<number>) {
    const n = xs.length;
    const position = { x: xs, y: ys, z: zs, radius: rs, indices: OrderedSet.ofRange(0, n) };
    const boundary = getBoundary(position);

    const points = new Float32Array(n * 3);
    const groups = new Float32Array(n);
    const chains = new Int32Array(n);
    const fracs = new Float32Array(n);
    const radii = new Float32Array(n);
    let maxRadius = 0;
    for (let i = 0; i < n; ++i) {
        points[i * 3] = xs[i]; points[i * 3 + 1] = ys[i]; points[i * 3 + 2] = zs[i];
        groups[i] = ids[i];
        chains[i] = seqChain[i];
        fracs[i] = seqFrac[i];
        radii[i] = rs[i];
        if (rs[i] > maxRadius) maxRadius = rs[i];
    }
    const boundingSphere = Sphere3D.expand(Sphere3D(), boundary.sphere, maxRadius);
    return { points, groups, chains, fracs, radii, center: Vec3.clone(boundary.sphere.center), boundingSphere, maxRadius };
}

/** Reconstruct an SH blob directly from an atom point cloud (no molecular surface). */
function reconstructFromAtomCloud(cloud: ReturnType<typeof buildAtomCloud>, props: SphericalHarmonicSurfaceMeshProps, resolution: number, cacheKey: string, ctx: VisualContext, mesh?: Mesh): Promise<Mesh> {
    const target = mesh ?? Mesh.createEmpty();
    populateCacheFromPoints(target.meta as SphericalHarmonicSurfaceMeta, cloud, resolution, cacheKey);
    return reconstructSphericalHarmonicMesh(target.meta as SphericalHarmonicSurfaceMeta, props, ctx, target);
}

/** Normalized residue position ([0,1)) of element `eI` within the contiguous residue range of `unit`. */
function residueFraction(unit: Unit, eI: ElementIndex): number {
    if (!Unit.isAtomic(unit)) return 0;
    const ri = unit.residueIndex;
    const els = unit.elements;
    const first = ri[els[0]];
    const last = ri[els[els.length - 1]];
    const span = last - first;
    return span > 0 ? (ri[eI] - first) / (span + 1) : 0;
}

/** Atom arrays (model coords) for one unit, with group id = index into unit.elements. */
function getUnitAtomArrays(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: SphericalHarmonicSurfaceMeshProps) {
    const { radiusOffset, ignoreHydrogens, ignoreHydrogensVariant, traceOnly } = props;
    const { x, y, z } = getConformation(unit);
    const { elements } = unit;
    const l = StructureElement.Location.create(structure, unit);
    const xs: number[] = [], ys: number[] = [], zs: number[] = [], rs: number[] = [], ids: number[] = [];
    const chains: number[] = [], fracs: number[] = [];
    const stride = Math.max(1, Math.floor(elements.length / MaxFitCloudPoints));
    let kept = 0;
    for (let j = 0, jl = elements.length; j < jl; ++j) {
        const eI = elements[j];
        if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
        if (traceOnly && !isTrace(unit, eI)) continue;
        if ((kept++ % stride) !== 0) continue; // uniformly subsample to MaxFitCloudPoints
        xs.push(x[eI]); ys.push(y[eI]); zs.push(z[eI]);
        l.element = eI;
        rs.push(sizeTheme.size(l) + radiusOffset);
        ids.push(j);
        chains.push(0); // one unit = one chain
        fracs.push(residueFraction(unit, eI));
    }
    return { xs, ys, zs, rs, ids, chains, fracs };
}

/** Atom arrays (world coords) for the given `units`, with group id = serial element index in `structure`. */
function getUnitsAtomArrays(structure: Structure, units: ReadonlyArray<Unit>, sizeTheme: SizeTheme<any>, props: SphericalHarmonicSurfaceMeshProps) {
    const { radiusOffset, ignoreHydrogens, ignoreHydrogensVariant, traceOnly } = props;
    const { getSerialIndex } = structure.serialMapping;
    const l = StructureElement.Location.create(structure);
    const xs: number[] = [], ys: number[] = [], zs: number[] = [], rs: number[] = [], ids: number[] = [];
    const chains: number[] = [], fracs: number[] = [];
    let total = 0;
    for (const unit of units) total += unit.elements.length;
    const stride = Math.max(1, Math.floor(total / MaxFitCloudPoints));
    let kept = 0;
    let chainOrdinal = 0;
    for (const unit of units) {
        const { elements, conformation: c } = unit;
        l.unit = unit;
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const eI = elements[j];
            if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
            if (traceOnly && !isTrace(unit, eI)) continue;
            if ((kept++ % stride) !== 0) continue; // uniformly subsample to MaxFitCloudPoints
            xs.push(c.x(eI)); ys.push(c.y(eI)); zs.push(c.z(eI));
            l.element = eI;
            rs.push(sizeTheme.size(l) + radiusOffset);
            ids.push(getSerialIndex(unit, eI));
            chains.push(chainOrdinal);
            fracs.push(residueFraction(unit, eI));
        }
        chainOrdinal++;
    }
    return { xs, ys, zs, rs, ids, chains, fracs };
}

/** Atom arrays (world coords) for the whole structure, with group id = serial element index. */
function getStructureAtomArrays(structure: Structure, sizeTheme: SizeTheme<any>, props: SphericalHarmonicSurfaceMeshProps) {
    return getUnitsAtomArrays(structure, structure.units, sizeTheme, props);
}

/** Group a structure's units by assembly operator id (`operId`); units of one id are one ASU copy. */
function groupUnitsByAssemblyOper(structure: Structure): Map<number, Unit[]> {
    const groups = new Map<number, Unit[]>();
    for (const u of structure.units) {
        const operId = u.conformation.operator.assembly?.operId ?? 0;
        let arr = groups.get(operId);
        if (!arr) { arr = []; groups.set(operId, arr); }
        arr.push(u);
    }
    return groups;
}

/**
 * Rigid transform mapping the `base` ASU copy onto the `other` copy. Found from any chain present in
 * both copies as `M_other * M_base^-1`; the chain's own (pre-assembly) operator cancels, so the result
 * is the pure inter-copy assembly transform regardless of which chain is used.
 */
function relativeAssemblyTransform(baseUnits: ReadonlyArray<Unit>, otherUnits: ReadonlyArray<Unit>): Mat4 | undefined {
    const baseByInvariant = new Map<number, Unit>();
    for (const u of baseUnits) baseByInvariant.set(u.invariantId, u);
    for (const uk of otherUnits) {
        const ub = baseByInvariant.get(uk.invariantId);
        if (ub) {
            const inv = Mat4.invert(Mat4(), ub.conformation.operator.matrix);
            return inv ? Mat4.mul(Mat4(), uk.conformation.operator.matrix, inv) : undefined;
        }
    }
    return undefined;
}

/** Replicate a single `envelope` mesh by each transform into one combined mesh (per-vertex groups preserved). */
function replicateMesh(envelope: Mesh, transforms: ReadonlyArray<Mat4>, mesh?: Mesh): Mesh {
    const V = envelope.vertexCount, T = envelope.triangleCount, n = transforms.length;
    const ev = envelope.vertexBuffer.ref.value, en = envelope.normalBuffer.ref.value;
    const ei = envelope.indexBuffer.ref.value, eg = envelope.groupBuffer.ref.value;

    const vertices = new Float32Array(V * 3 * n);
    const normals = new Float32Array(V * 3 * n);
    const groups = new Float32Array(V * n);
    const indices = new Uint32Array(T * 3 * n);

    const p = Vec3(), nrm = Vec3(), m3 = Mat3();
    for (let c = 0; c < n; ++c) {
        const t = transforms[c];
        Mat3.directionTransform(m3, t);
        const voff = c * V * 3, goff = c * V, ioff = c * T * 3, vBase = c * V;
        for (let i = 0; i < V; ++i) {
            Vec3.transformMat4(p, Vec3.set(p, ev[i * 3], ev[i * 3 + 1], ev[i * 3 + 2]), t);
            vertices[voff + i * 3] = p[0]; vertices[voff + i * 3 + 1] = p[1]; vertices[voff + i * 3 + 2] = p[2];
            Vec3.normalize(nrm, Vec3.transformMat3(nrm, Vec3.set(nrm, en[i * 3], en[i * 3 + 1], en[i * 3 + 2]), m3));
            normals[voff + i * 3] = nrm[0]; normals[voff + i * 3 + 1] = nrm[1]; normals[voff + i * 3 + 2] = nrm[2];
            groups[goff + i] = eg[i];
        }
        for (let i = 0, il = T * 3; i < il; ++i) indices[ioff + i] = ei[i] + vBase;
    }
    return Mesh.create(vertices, indices, normals, groups, V * n, T * n, mesh);
}

//

async function createSphericalHarmonicSurfaceMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const start = isTimingMode ? now() : 0;
    const cacheKey = surfaceCacheKey(`unit-${unit.id}`, props);
    const meta = mesh?.meta as SphericalHarmonicSurfaceMeta | undefined;
    if (mesh && meta && meta.shCacheKey === cacheKey && meta.shVertices) {
        const result = await reconstructSphericalHarmonicMesh(meta, props, ctx, mesh);
        shAccumulateTiming(start, meta.shVertices.length / 3);
        return result;
    }

    const { xs, ys, zs, rs, ids, chains, fracs } = getUnitAtomArrays(structure, unit, theme.size, props);
    const cloud = buildAtomCloud(xs, ys, zs, rs, ids, chains, fracs);
    const result = await reconstructFromAtomCloud(cloud, props, props.resolution, cacheKey, ctx, mesh);
    shAccumulateTiming(start, xs.length);
    return result;
}

export function SphericalHarmonicSurfaceMeshVisual(materialId: number): UnitsVisual<SphericalHarmonicSurfaceMeshParams> {
    return UnitsMeshVisual<SphericalHarmonicSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(SphericalHarmonicSurfaceMeshParams),
        createGeometry: createSphericalHarmonicSurfaceMesh,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) => {
            state.createGeometry = surfaceGeometryChanged(newProps, currentProps);
            updateSmoothColorsState(state, newProps, currentProps);
        },
        processValues: shProcessValues,
        dispose: shDispose,
    }, materialId);
}

//

/**
 * Structure flavour: fit ONE envelope to every atom of the whole structure merged together, for a
 * deliberately coarse single-object outline of the entire complex (vs. the per-chain units flavour or
 * the symmetry-replicated assembly flavour). Star-convex, so it reads as a smooth bounding shape of
 * the input rather than a faithful multi-chain surface.
 */
async function createStructureSphericalHarmonicSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const start = isTimingMode ? now() : 0;
    const cacheKey = surfaceCacheKey(`structure-${structure.hashCode}`, props);
    const meta = mesh?.meta as SphericalHarmonicSurfaceMeta | undefined;
    if (mesh && meta && meta.shCacheKey === cacheKey && meta.shVertices) {
        const result = await reconstructSphericalHarmonicMesh(meta, props, ctx, mesh);
        shAccumulateTiming(start, meta.shVertices.length / 3);
        return result;
    }

    const { xs, ys, zs, rs, ids, chains, fracs } = getStructureAtomArrays(structure, theme.size, props);
    const cloud = buildAtomCloud(xs, ys, zs, rs, ids, chains, fracs);
    const result = await reconstructFromAtomCloud(cloud, props, props.resolution, cacheKey, ctx, mesh);
    shAccumulateTiming(start, xs.length);
    return result;
}

export function StructureSphericalHarmonicSurfaceMeshVisual(materialId: number): ComplexVisual<SphericalHarmonicSurfaceMeshParams> {
    return ComplexMeshVisual<SphericalHarmonicSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(SphericalHarmonicSurfaceMeshParams),
        createGeometry: createStructureSphericalHarmonicSurfaceMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) => {
            state.createGeometry = surfaceGeometryChanged(newProps, currentProps);
            updateSmoothColorsState(state, newProps, currentProps);
        },
        processValues: shProcessValues,
        dispose: shDispose,
    }, materialId);
}

//

/**
 * Assembly flavour: fit a single SH envelope to ONE asymmetric-unit copy (the chains of the lowest
 * assembly operator id, in their world frame), then replicate that one mesh across the assembly's
 * inter-copy transforms. This is the right model for symmetric assemblies (one capsid protomer fit
 * once, drawn N times) and avoids both the per-chain blobs of the units flavour and the single merged
 * blob of the structure flavour. Uses the atom-cloud fit (smooth, robust); the per-vertex groups are
 * shared across copies, so chain-id coloring is consistent and each copy looks identical.
 *
 * Without an assembly (a single operator id) this degrades to one envelope over the whole input.
 */
async function createAssemblySphericalHarmonicSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const target = mesh ?? Mesh.createEmpty();
    const tmeta = target.meta as SphericalHarmonicSurfaceMeta;

    const operGroups = groupUnitsByAssemblyOper(structure);
    const operIds = Array.from(operGroups.keys()).sort((a, b) => a - b);
    const baseUnits = operGroups.get(operIds[0])!;

    // fit one ASU envelope from its atom cloud into a scratch mesh (kept off the render object)
    const cacheKey = surfaceCacheKey(`assembly-${structure.hashCode}-${operIds[0]}`, props);
    const { xs, ys, zs, rs, ids, chains, fracs } = getUnitsAtomArrays(structure, baseUnits, theme.size, props);
    const cloud = buildAtomCloud(xs, ys, zs, rs, ids, chains, fracs);
    const scratch = tmeta.shAsuMesh ?? Mesh.createEmpty();
    populateCacheFromPoints(scratch.meta as SphericalHarmonicSurfaceMeta, cloud, props.resolution, cacheKey);
    const envelope = await reconstructSphericalHarmonicMesh(scratch.meta as SphericalHarmonicSurfaceMeta, props, ctx, scratch);
    tmeta.shAsuMesh = envelope;

    // one transform per assembly copy, relative to the fitted base copy (base = identity)
    const transforms: Mat4[] = [];
    for (const oid of operIds) {
        if (oid === operIds[0]) { transforms.push(Mat4.identity()); continue; }
        const t = relativeAssemblyTransform(baseUnits, operGroups.get(oid)!);
        if (t) transforms.push(t);
    }

    const combined = replicateMesh(envelope, transforms, target);
    if (ctx.webgl && !ctx.webgl.isWebGL2) {
        Mesh.uniformTriangleGroup(combined);
        ValueCell.updateIfChanged(combined.varyingGroup, false);
    } else {
        ValueCell.updateIfChanged(combined.varyingGroup, true);
    }
    // keep the envelope scratch and the color-smoothing resolution on the bound mesh's meta.
    // Use the envelope's per-copy vertex spacing (set by reconstructSphericalHarmonicMesh), not
    // props.resolution: replication doesn't change per-copy edge length, and the fine base
    // resolution would leave the smoothing grid full of black cells (see estimateVertexSpacing).
    const cmeta = combined.meta as SphericalHarmonicSurfaceMeta;
    cmeta.shAsuMesh = envelope;
    cmeta.resolution = (envelope.meta as SphericalHarmonicSurfaceMeta).resolution;
    return combined;
}

export function AssemblySphericalHarmonicSurfaceMeshVisual(materialId: number): ComplexVisual<SphericalHarmonicSurfaceMeshParams> {
    return ComplexMeshVisual<SphericalHarmonicSurfaceMeshParams>({
        defaultProps: PD.getDefaultValues(SphericalHarmonicSurfaceMeshParams),
        createGeometry: createAssemblySphericalHarmonicSurfaceMesh,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps) => {
            state.createGeometry = surfaceGeometryChanged(newProps, currentProps);
            updateSmoothColorsState(state, newProps, currentProps);
        },
        processValues: shProcessValues,
        dispose: shDispose,
    }, materialId);
}
