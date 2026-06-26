/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 *
 * Spherical-harmonic surface visual: reuses Mol*'s molecular-surface mesh, fits
 * a real spherical-harmonic expansion of the radial function r(theta, phi) about
 * the surface centroid, and reconstructs a smooth closed mesh from the fitted
 * coefficients. The parameter `sphericalHarmonicL` controls the level of detail.
 *
 * Note: a single-center radial expansion can only represent star-convex surfaces
 * (those single-valued in r from the centroid). Pockets, concavities and
 * multi-domain shapes fold back on a radial ray and collapse to the outermost
 * crossing, so a single lobe is best understood as a smooth surface *envelope*
 * rather than a faithful molecular surface. To handle non-star-shaped inputs the
 * cloud can be decomposed into several star-shaped lobes (`maxLobes`); a single
 * lobe is reconstructed as an analytic icosphere blob, while multiple lobes are
 * blended into one watertight surface through a marching-cubes implicit field.
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
import { CommonSurfaceParams, getConformation, isHydrogen, isTrace } from './util/common';
import { Sphere3D } from '../../../mol-math/geometry';
import { MeshValues } from '../../../mol-gl/renderable/mesh';
import { Texture } from '../../../mol-gl/webgl/texture';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { applyMeshColorSmoothing } from '../../../mol-geo/geometry/mesh/color-smoothing';
import { ColorSmoothingParams, getColorSmoothingProps } from '../../../mol-geo/geometry/base';
import { ValueCell } from '../../../mol-util';
import { ComplexMeshVisual, ComplexVisual } from '../complex-visual';
import { Tensor } from '../../../mol-math/linear-algebra/tensor';
import { Mat3, Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Sphere } from '../../../mol-geo/primitive/sphere';
import { GridLookup3D } from '../../../mol-math/geometry/lookup3d/grid';
import { getBoundary } from '../../../mol-math/geometry/boundary';
import { OrderedSet } from '../../../mol-data/int';
import { fitSphericalHarmonics, fitSphericalHarmonicLobes, reconstructRadius, shTermCount, SphericalHarmonicLobe } from '../../../mol-math/geometry/spherical-harmonics';
import { SizeTheme } from '../../../mol-theme/size';
import { isTimingMode } from '../../../mol-util/debug';

export const SphericalHarmonicSurfaceMeshParams = {
    ...UnitsMeshParams,
    ...CommonMolecularSurfaceCalculationParams,
    ...CommonSurfaceParams,
    ...ColorSmoothingParams,
    sphericalHarmonicL: PD.Numeric(8, { min: 2, max: 20, step: 1 }, { description: 'Maximum spherical harmonic degree. Higher values follow the surface envelope more closely (but cannot recover concavities: the fit is single-valued in radius about the center).' }),
    reconstructionDetail: PD.Numeric(3, { min: 1, max: 5, step: 1 }, { description: 'Triangulation detail of the reconstructed sphere mesh.' }),
    fitSource: PD.Select('surface', [['surface', 'Molecular surface'], ['atoms', 'Atom positions']] as const, { description: 'Fit the spherical harmonics to the molecular surface, or directly to atom positions (faster, blobbier, no surface step).' }),
    maxLobes: PD.Numeric(1, { min: 1, max: 8, step: 1 }, { description: 'Split non-star-shaped inputs (elongated or multi-domain) into up to this many star-shaped lobes, fit separately, and blend into one watertight surface. 1 keeps a single radial envelope.' }),
    regularization: PD.Numeric(0.01, { min: 0, max: 0.5, step: 0.005 }, { description: 'Tikhonov smoothness prior on the spherical-harmonic fit (per-band l(l+1) damping). Higher values give a rounder, more stable surface; lower values follow the samples more closely but can ring or blow up on sparse/clustered clouds (e.g. trace-only). 0 disables it.' }),
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
    // cache of the fitted lobes, keyed by `${L}|${maxLobes}` so a detail-only change re-tessellates without refitting
    shFitKey?: string
    shFitLobes?: SphericalHarmonicLobe[]
    // scratch mesh holding the single ASU envelope for the assembly flavour, before it is replicated
    // across the assembly operators into the bound mesh (kept off the render object, reused).
    shAsuMesh?: Mesh
    // throwaway mesh holding the molecular-surface base for the 'surface' fit source; kept OFF the
    // render object (the bound mesh must only ever hold the final reconstruction) and reused across
    // recomputes. Mixing the base into the bound mesh briefly desyncs its vertex/index buffers
    // mid-update -> out-of-range index fetch -> the surface vanishes.
    shBaseMesh?: Mesh
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

/**
 * Cap on the number of atom-cloud points gathered for the fit and the group transfer. A smooth
 * degree-L envelope is determined by far fewer points (the fit strides to 8192 internally) and the
 * per-vertex group transfer is approximate, so for large structures the cloud is uniformly strided:
 * building the point arrays and the group-transfer `GridLookup3D` over every atom (millions for e.g.
 * 3J3Q) is the dominant cost of the atom-fit path and is unnecessary. Uniform stride preserves shape.
 */
const MaxFitCloudPoints = 50000;

const MinRadius = 1e-3;
/** Reconstruction radius is clamped to `rMax * RadiusMargin` so an ill-conditioned fit cannot send vertices far outside the data (which would blow up the bounding sphere and cull the surface). */
const RadiusMargin = 1.15;

/** True if any base-surface prop shared by the unit/structure surface flavours changed. */
function surfaceGeometryChanged(newProps: SphericalHarmonicSurfaceMeshProps, currentProps: SphericalHarmonicSurfaceMeshProps): boolean {
    return (
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
        newProps.fitSource !== currentProps.fitSource ||
        newProps.maxLobes !== currentProps.maxLobes ||
        newProps.regularization !== currentProps.regularization
    );
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

/** Copy the group id of the nearest original surface vertex onto each reconstructed vertex. */
function transferGroups(vertices: Float32Array, vertexCount: number, shLookup: GridLookup3D, shGroups: Float32Array, out: Float32Array) {
    for (let i = 0; i < vertexCount; ++i) {
        const res = shLookup.nearest(vertices[i * 3], vertices[i * 3 + 1], vertices[i * 3 + 2], 1);
        out[i] = res.count > 0 ? shGroups[res.indices[0]] : 0;
    }
}

/**
 * Single star-shaped lobe → analytic icosphere blob (fast, instanceable): scale
 * each unit-sphere direction by the reconstructed radius about the lobe center.
 */
function reconstructIcosphereBlob(lobe: SphericalHarmonicLobe, L: number, props: SphericalHarmonicSurfaceMeshProps, shLookup: GridLookup3D, shGroups: Float32Array, mesh?: Mesh): Mesh {
    const K = shTermCount(L);
    const { center, coeffs } = lobe;

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
    const cx = center[0], cy = center[1], cz = center[2];
    const rClamp = lobe.rMax * RadiusMargin;

    for (let i = 0; i < vertexCount; ++i) {
        // sphere vertices are unit length (radius 1), so they are direction vectors
        const dx = sv[i * 3], dy = sv[i * 3 + 1], dz = sv[i * 3 + 2];
        const theta = Math.acos(Math.min(1, Math.max(-1, dz)));
        const phi = Math.atan2(dy, dx);
        let r = reconstructRadius(coeffs, L, theta, phi, basisScratch, legendreScratch);
        if (!(r > MinRadius)) r = MinRadius; // also catches NaN
        else if (r > rClamp) r = rClamp;
        vertices[i * 3] = cx + r * dx; vertices[i * 3 + 1] = cy + r * dy; vertices[i * 3 + 2] = cz + r * dz;
    }
    transferGroups(vertices, vertexCount, shLookup, shGroups, groups);

    const surface = Mesh.create(vertices, indices, normals, groups, vertexCount, triangleCount, mesh);
    Mesh.computeNormals(surface);
    return surface;
}

const FieldPad = 1.5;
/** Smooth-union blend width (length units) for the multi-lobe implicit field. */
const FieldBlend = 2.0;
/** Cap on the multi-lobe field grid edge so cost (O(dim^3 * lobes)) stays bounded and independent of the quality/resolution slider. */
const MaxFieldDim = 64;

/**
 * Multiple lobes → one watertight surface. Each lobe contributes an inside-ness
 * scalar f = r_fit(dir) - |p - center| (positive inside its SH envelope); the
 * lobes are combined with a soft-max (log-sum-exp) so neighbours blend smoothly
 * instead of producing the seam crease / internal walls of overlaid blobs. The
 * field is marching-cubed once at iso-level 0.
 */
async function reconstructLobesField(lobes: SphericalHarmonicLobe[], L: number, props: SphericalHarmonicSurfaceMeshProps, shLookup: GridLookup3D, shGroups: Float32Array, boundingSphere: Sphere3D, ctx: VisualContext, mesh?: Mesh): Promise<Mesh> {
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
    const beta = 1 / FieldBlend;

    const nL = lobes.length;
    const fScratch = new Float64Array(nL); // per-lobe inside-ness at the current cell
    const rClamp = lobes.map(l => l.rMax * RadiusMargin);
    for (let i = 0; i < dimN; ++i) {
        const px = min[0] + i * resolution;
        for (let j = 0; j < dimN; ++j) {
            const py = min[1] + j * resolution;
            for (let k = 0; k < dimN; ++k) {
                const pz = min[2] + k * resolution;
                // f per lobe (computed once), combined via soft-max:
                // F = fmax + log(sum exp(beta (f - fmax))) / beta
                let fmax = -Infinity;
                for (let li = 0; li < nL; ++li) {
                    const ctr = lobes[li].center;
                    const dx = px - ctr[0], dy = py - ctr[1], dz = pz - ctr[2];
                    const rr = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    // Far outside this lobe's radius the cell is unambiguously outside it, so the
                    // exact (negative) inside-ness doesn't affect the soft-max union; skip the
                    // expensive SH evaluation and use the cheap monotonic estimate rMax - rr.
                    let f: number;
                    if (rr > rClamp[li]) {
                        f = rClamp[li] - rr;
                    } else {
                        const theta = rr > 1e-9 ? Math.acos(Math.min(1, Math.max(-1, dz / rr))) : 0;
                        const phi = Math.atan2(dy, dx);
                        let rfit = reconstructRadius(lobes[li].coeffs, L, theta, phi, basisScratch, legendreScratch);
                        if (rfit > rClamp[li]) rfit = rClamp[li];
                        f = rfit - rr;
                    }
                    fScratch[li] = f;
                    if (f > fmax) fmax = f;
                }
                let sum = 0;
                for (let li = 0; li < nL; ++li) sum += Math.exp(beta * (fScratch[li] - fmax));
                space.set(data, i, j, k, fmax + Math.log(sum) / beta);
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
    transferGroups(vb, vertexCount, shLookup, shGroups, groups);
    ValueCell.update(surface.groupBuffer, groups);

    return surface;
}

/**
 * Reconstruct the SH surface from the cached base surface stored in `meta`.
 *
 * Fits one or more star-shaped lobes (cached by `${L}|${maxLobes}`, so a
 * detail-only change re-tessellates without refitting). A single lobe is
 * rendered as an analytic icosphere blob; multiple lobes are blended into one
 * watertight surface via a marching-cubes implicit field.
 */
async function reconstructSphericalHarmonicMesh(meta: SphericalHarmonicSurfaceMeta, props: SphericalHarmonicSurfaceMeshProps, ctx: VisualContext, mesh?: Mesh): Promise<Mesh> {
    const { shVertices, shGroups, shLookup, shCenter, shBoundingSphere } = meta;
    if (!shVertices || !shGroups || !shLookup || !shCenter || !shBoundingSphere || shVertices.length === 0) {
        return Mesh.createEmpty(mesh);
    }

    const L = Math.round(props.sphericalHarmonicL);
    const maxLobes = Math.max(1, Math.round(props.maxLobes));
    const regularization = props.regularization;

    const fitKey = `${L}|${maxLobes}|${regularization}`;
    let lobes = meta.shFitKey === fitKey ? meta.shFitLobes : undefined;
    if (!lobes) {
        if (isTimingMode) console.time(`SphericalHarmonicSurface fit (${shVertices.length / 3} pts, L=${L}, lobes<=${maxLobes})`);
        if (maxLobes <= 1) {
            const { coeffs, rMax } = fitSphericalHarmonics(shVertices, shCenter, L, undefined, regularization);
            lobes = [{ center: [shCenter[0], shCenter[1], shCenter[2]], coeffs, rMax }];
        } else {
            // Split only genuinely non-star-shaped clouds (crescents, rings, centroid-outside
            // shapes) - NOT ordinary compact or two-domain blobs, which a single radial shell
            // already represents. The gap threshold scales with the cloud radius because a solid
            // cloud's interior sampling leaves gaps ~proportional to its size; below this even a
            // star-shaped blob looks "re-entrant". Empirically ~R/3 with tolerance ~0.3 keeps a
            // solid ball at one lobe while a true gap (~R) still triggers.
            const thickness = Math.max(3, shBoundingSphere.radius * 0.33);
            lobes = fitSphericalHarmonicLobes(shVertices, L, { maxLobes, tolerance: 0.3, thickness, regularization }).lobes;
        }
        if (isTimingMode) console.timeEnd(`SphericalHarmonicSurface fit (${shVertices.length / 3} pts, L=${L}, lobes<=${maxLobes})`);
    }

    const surface = lobes.length <= 1
        ? reconstructIcosphereBlob(lobes[0], L, props, shLookup, shGroups, mesh)
        : await reconstructLobesField(lobes, L, props, shLookup, shGroups, shBoundingSphere, ctx, mesh);

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
    // the point cloud changed: invalidate the cached fit so it is recomputed
    meta.shFitKey = undefined;
    meta.shFitLobes = undefined;
}

/** Populate the base-surface cache in `targetMeta` from a freshly computed molecular `surface` mesh. */
function populateCacheFromSurface(targetMeta: SphericalHarmonicSurfaceMeta, surface: Mesh, center: Vec3, boundingSphere: Sphere3D, resolution: number, cacheKey: string) {
    const vertexCount = surface.vertexCount;
    const points = surface.vertexBuffer.ref.value.slice(0, vertexCount * 3);
    const groups = surface.groupBuffer.ref.value.slice(0, vertexCount);
    populateCacheFromPoints(targetMeta, points, groups, center, boundingSphere, resolution, cacheKey);
}

/**
 * Build an SH-fitting point cloud from atom positions: each atom is pushed outward by its radius
 * (size + probe) along the centroid->atom direction so the fit tracks the surface envelope rather
 * than the (inward-biased) atom centers. `ids[i]` is the group id carried to the reconstructed mesh.
 */
function buildAtomCloud(xs: number[], ys: number[], zs: number[], rs: number[], ids: ArrayLike<number>) {
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
function reconstructFromAtomCloud(cloud: ReturnType<typeof buildAtomCloud>, props: SphericalHarmonicSurfaceMeshProps, resolution: number, cacheKey: string, ctx: VisualContext, mesh?: Mesh): Promise<Mesh> {
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
    const stride = Math.max(1, Math.floor(elements.length / MaxFitCloudPoints));
    let kept = 0;
    for (let j = 0, jl = elements.length; j < jl; ++j) {
        const eI = elements[j];
        if (ignoreHydrogens && isHydrogen(structure, unit, eI, ignoreHydrogensVariant)) continue;
        if (traceOnly && !isTrace(unit, eI)) continue;
        if ((kept++ % stride) !== 0) continue; // uniformly subsample to MaxFitCloudPoints
        xs.push(x[eI]); ys.push(y[eI]); zs.push(z[eI]);
        l.element = eI;
        rs.push(sizeTheme.size(l) + probeRadius);
        ids.push(j);
    }
    return { xs, ys, zs, rs, ids };
}

/** Atom arrays (world coords) for the given `units`, with group id = serial element index in `structure`. */
function getUnitsAtomArrays(structure: Structure, units: ReadonlyArray<Unit>, sizeTheme: SizeTheme<any>, props: SphericalHarmonicSurfaceMeshProps) {
    const { probeRadius, ignoreHydrogens, ignoreHydrogensVariant, traceOnly } = props;
    const { getSerialIndex } = structure.serialMapping;
    const l = StructureElement.Location.create(structure);
    const xs: number[] = [], ys: number[] = [], zs: number[] = [], rs: number[] = [], ids: number[] = [];
    let total = 0;
    for (const unit of units) total += unit.elements.length;
    const stride = Math.max(1, Math.floor(total / MaxFitCloudPoints));
    let kept = 0;
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
            rs.push(sizeTheme.size(l) + probeRadius);
            ids.push(getSerialIndex(unit, eI));
        }
    }
    return { xs, ys, zs, rs, ids };
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
    // the molecular-surface base goes into a throwaway mesh (reused via meta), never the bound mesh
    const target = mesh ?? Mesh.createEmpty();
    const tmeta = target.meta as SphericalHarmonicSurfaceMeta;
    const base = await computeMarchingCubesMesh(params, tmeta.shBaseMesh).runAsChild(ctx.runtime);
    tmeta.shBaseMesh = base;
    if (props.includeParent) {
        const iterations = Math.ceil(2 / props.resolution);
        Mesh.smoothEdges(base, { iterations, maxNewEdgeLength: Math.sqrt(2) });
    }
    Mesh.transform(base, transform);

    const boundingSphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius);
    populateCacheFromSurface(tmeta, base, unit.boundary.sphere.center as Vec3, boundingSphere, resolution, cacheKey);
    return reconstructSphericalHarmonicMesh(tmeta, props, ctx, target);
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

async function createStructureSphericalHarmonicSurfaceMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: SphericalHarmonicSurfaceMeshProps, mesh?: Mesh): Promise<Mesh> {
    const cacheKey = surfaceCacheKey(`structure-${structure.hashCode}`, props);
    const meta = mesh?.meta as SphericalHarmonicSurfaceMeta | undefined;
    if (mesh && meta && meta.shCacheKey === cacheKey && meta.shVertices) {
        return reconstructSphericalHarmonicMesh(meta, props, ctx, mesh);
    }

    if (props.fitSource === 'atoms') {
        if (isTimingMode) console.time('SphericalHarmonicSurface atom cloud');
        const { xs, ys, zs, rs, ids } = getStructureAtomArrays(structure, theme.size, props);
        const cloud = buildAtomCloud(xs, ys, zs, rs, ids);
        if (isTimingMode) console.timeEnd('SphericalHarmonicSurface atom cloud');
        return reconstructFromAtomCloud(cloud, props, props.resolution, cacheKey, ctx, mesh);
    }

    const { transform, field, idField, resolution, maxRadius } = await computeStructureMolecularSurface(structure, theme.size, props).runInContext(ctx.runtime);
    const params = {
        isoLevel: props.probeRadius,
        scalarField: props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
        idField,
    };
    // the molecular-surface base goes into a throwaway mesh (reused via meta), never the bound mesh
    const target = mesh ?? Mesh.createEmpty();
    const tmeta = target.meta as SphericalHarmonicSurfaceMeta;
    const base = await computeMarchingCubesMesh(params, tmeta.shBaseMesh).runAsChild(ctx.runtime);
    tmeta.shBaseMesh = base;
    if (props.includeParent) {
        const iterations = Math.ceil(2 / props.resolution);
        Mesh.smoothEdges(base, { iterations, maxNewEdgeLength: Math.sqrt(2) });
    }
    Mesh.transform(base, transform);

    const boundingSphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius);
    populateCacheFromSurface(tmeta, base, structure.boundary.sphere.center as Vec3, boundingSphere, resolution, cacheKey);
    return reconstructSphericalHarmonicMesh(tmeta, props, ctx, target);
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
    const { xs, ys, zs, rs, ids } = getUnitsAtomArrays(structure, baseUnits, theme.size, props);
    const cloud = buildAtomCloud(xs, ys, zs, rs, ids);
    const scratch = tmeta.shAsuMesh ?? Mesh.createEmpty();
    populateCacheFromPoints(scratch.meta as SphericalHarmonicSurfaceMeta, cloud.points, cloud.groups, cloud.center, cloud.boundingSphere, props.resolution, cacheKey);
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
