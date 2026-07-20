/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PositionData, DensityData, fillGridDim } from './common';
import { Boundary, getFastBoundary } from './boundary';
import { Box3D } from './primitives/box3d';
import { GridLookup3D } from './lookup3d/grid';
import { Result } from './lookup3d/common';
import { OrderedSet } from '../../mol-data/int/ordered-set';
import { Vec3 } from '../linear-algebra/3d/vec3';
import { Mat4 } from '../linear-algebra/3d/mat4';
import { Tensor } from '../linear-algebra/tensor';
import { PrincipalAxes } from '../linear-algebra/matrix/principal-axes';
import { fasterExp } from '../approx';
import { fitSphericalHarmonics, reconstructRadius, toSpherical, SphericalCoord, shTermCount } from './spherical-harmonics';

export const DefaultBlobSurfaceProps = {
    blobSize: 30,
    method: 'grid' as 'grid' | 'clustering',
    clusterIterations: 2,
    shape: 'ellipsoid' as 'ellipsoid' | 'sh',
    shDegree: 2,
    shRegularization: 0.05,
    resolution: 2,
    radiusOffset: 0,
    smoothness: 1.5,
};
export type BlobSurfaceProps = typeof DefaultBlobSurfaceProps

export type BlobSurfaceData = { radiusFactor: number } & DensityData

interface EllipsoidBlob {
    kind: 'ellipsoid'
    center: Vec3
    /** orthonormal local axis directions */
    u0: Vec3, u1: Vec3, u2: Vec3
    /** semi-axis lengths along u0/u1/u2 (already inflated by atom radii + radiusOffset) */
    a0: number, a1: number, a2: number
    /** representative atom id (nearest to `center`), same convention as `PositionData.id` */
    atomId: number
}

interface SHBlob {
    kind: 'sh'
    center: Vec3
    /** max supported degree of the fitted real spherical-harmonics radial boundary `R(theta, phi)` */
    degree: number
    coeffs: Float64Array
    /** hard ceiling for `R(theta, phi)` (also used for padding): a small fixed margin over the
     * largest observed per-atom target radius, bounding any fit overshoot in unsampled directions */
    maxR: number
    /** representative atom id (nearest to `center`), same convention as `PositionData.id` */
    atomId: number
}

/** margin applied over the largest observed per-atom target radius when clamping/padding an SH blob's `R(theta, phi)` */
const SH_OVERSHOOT_CLAMP = 1.25;

type Blob = EllipsoidBlob | SHBlob

function blobMaxExtent(blob: Blob) {
    return blob.kind === 'ellipsoid' ? Math.max(blob.a0, blob.a1, blob.a2) : blob.maxR;
}

/**
 * A group of atoms belonging to one blob, given as a `[offset, offset + count)` range into
 * `array` rather than its own copy - `array` may be shared/reused across many groups (e.g. all
 * of `GridLookup3D`'s buckets live in one flat array), so building a group is a zero-copy O(1)
 * operation. Values are ordinal positions (`t`, indices into `position.indices`), same convention
 * throughout this module.
 */
interface AtomGroup {
    readonly array: ArrayLike<number>
    readonly offset: number
    readonly count: number
}

/**
 * Groups atoms into a small number of spatial clusters via a uniform spatial hash grid (reusing
 * `GridLookup3D`'s existing bucket-grid construction with a fixed `blobSize` cell size, rather
 * than hand-rolling one with a `Map`): atoms are binned by cell, one group per occupied cell.
 * This is a single O(atoms) pass, no iteration. The returned groups are plain views into
 * `GridLookup3D`'s own flat `buckets.array` (no copying).
 * `GridLookup3D` also guards against pathologically large grids (tiny `blobSize` on
 * a huge structure) by capping total cell count, silently coarsening cell size if needed.
 */
function groupAtomsByGrid(position: PositionData, boundary: Boundary, blobSize: number): AtomGroup[] {
    const { offset, count, array } = GridLookup3D(position, boundary, Vec3.create(blobSize, blobSize, blobSize)).buckets;

    const groups: AtomGroup[] = new Array(offset.length);
    for (let i = 0; i < offset.length; ++i) {
        groups[i] = { array, offset: offset[i], count: count[i] };
    }
    return groups;
}

/**
 * Groups atoms into a small number of spatial clusters via k-means (Lloyd's algorithm), seeded
 * from `groupAtomsByGrid`'s bins (both for the target cluster count `k` and for well-spread
 * initial centroids, for free). Each iteration reassigns every atom to its nearest centroid and
 * recomputes centroids as the mean of their assigned atoms (O(atoms * k) per iteration); empty
 * clusters are reseeded to the atom farthest from its currently assigned centroid, a standard
 * empty-cluster remedy. Softens/avoids the grid method's atoms-near-a-bin-boundary-split
 * artifact, at the cost of `clusterIterations` extra O(atoms * k) passes.
 */
function groupAtomsByClustering(position: PositionData, boundary: Boundary, blobSize: number, iterations: number): AtomGroup[] {
    const initialGroups = groupAtomsByGrid(position, boundary, blobSize);
    const k = initialGroups.length;
    if (k <= 1) return initialGroups;

    const { indices, x, y, z } = position;
    const n = OrderedSet.size(indices);

    const centroids: Vec3[] = initialGroups.map(group => {
        const c = Vec3();
        for (let k = 0; k < group.count; ++k) {
            const t = group.array[group.offset + k];
            const i = OrderedSet.getAt(indices, t);
            c[0] += x[i]; c[1] += y[i]; c[2] += z[i];
        }
        return Vec3.scale(c, c, 1 / group.count);
    });

    const assignment = new Int32Array(n);
    for (let gi = 0; gi < initialGroups.length; ++gi) {
        const group = initialGroups[gi];
        for (let k = 0; k < group.count; ++k) {
            const t = group.array[group.offset + k];
            assignment[t] = gi;
        }
    }

    const cx = new Float64Array(k), cy = new Float64Array(k), cz = new Float64Array(k);
    const centroidPositions: PositionData = { x: cx, y: cy, z: cz, indices: OrderedSet.ofBounds(0, k) };
    const nearestResult = Result.create<number>();
    const cellSize = Vec3.create(blobSize, blobSize, blobSize);

    for (let iter = 0; iter < iterations; ++iter) {
        // console.time('k-means iteration');
        let changed = false;

        // nearest-centroid assignment via a spatial index over the centroids, rebuilt once per
        // iteration (cheap: O(k)) - a brute-force scan of all `k` centroids per atom is
        // O(atoms * k), and since `k` itself grows with structure size for a fixed `blobSize`
        // (both atom count and bin count scale with volume), that is effectively O(atoms^2) and
        // becomes very slow well before 1M atoms. `GridLookup3D.nearest` turns the per-atom
        // lookup into an O(1)-amortized query, dropping the whole pass to ~O(atoms) regardless
        // of how large `k` gets.
        for (let ci = 0; ci < k; ++ci) {
            const c = centroids[ci];
            cx[ci] = c[0]; cy[ci] = c[1]; cz[ci] = c[2];
        }
        const centroidLookup = GridLookup3D(centroidPositions, getFastBoundary(centroidPositions), cellSize);

        for (let t = 0; t < n; ++t) {
            const i = OrderedSet.getAt(indices, t);
            // TODO: rationalize radius value: blobSize / 3
            centroidLookup.approxNearest(x[i], y[i], z[i], blobSize / 3, nearestResult);
            const best = nearestResult.count > 0 ? nearestResult.indices[0] : assignment[t];
            if (assignment[t] !== best) {
                assignment[t] = best;
                changed = true;
            }
        }

        const sums = new Float64Array(k * 3);
        const counts = new Int32Array(k);
        for (let t = 0; t < n; ++t) {
            const i = OrderedSet.getAt(indices, t);
            const ci = assignment[t];
            sums[ci * 3] += x[i]; sums[ci * 3 + 1] += y[i]; sums[ci * 3 + 2] += z[i];
            counts[ci]++;
        }

        for (let ci = 0; ci < k; ++ci) {
            if (counts[ci] > 0) {
                Vec3.set(centroids[ci], sums[ci * 3] / counts[ci], sums[ci * 3 + 1] / counts[ci], sums[ci * 3 + 2] / counts[ci]);
                continue;
            }
            // empty cluster: reseed to the atom farthest from its currently assigned centroid
            let farthestT = 0, farthestDistSq = -1;
            for (let t = 0; t < n; ++t) {
                const i = OrderedSet.getAt(indices, t);
                const c = centroids[assignment[t]];
                const dx = x[i] - c[0], dy = y[i] - c[1], dz = z[i] - c[2];
                const dSq = dx * dx + dy * dy + dz * dz;
                if (dSq > farthestDistSq) {
                    farthestDistSq = dSq;
                    farthestT = t;
                }
            }
            const i = OrderedSet.getAt(indices, farthestT);
            Vec3.set(centroids[ci], x[i], y[i], z[i]);
            changed = true;
        }
        // console.timeEnd('k-means iteration');
        if (!changed) break;
    }

    const groups: number[][] = Array.from({ length: k }, () => []);
    for (let t = 0; t < n; ++t) groups[assignment[t]].push(t);
    return groups.filter(g => g.length > 0).map(array => ({ array, offset: 0, count: array.length }));
}

/**
 * Fits a single ellipsoid blob (center + local covariance axes via `PrincipalAxes`, inflated by
 * the group's atom radii) to one group of atoms. Groups with a single atom fall back to an
 * isotropic sphere (no covariance to fit).
 */
function fitEllipsoidBlobFromGroup(position: PositionData, radius: (index: number) => number, radiusOffset: number, group: AtomGroup): Blob {
    const { indices, x, y, z, id } = position;
    const m = group.count;

    const posBuf: number[] = [];
    let avgRadius = 0;
    for (let k = 0; k < group.count; ++k) {
        const t = group.array[group.offset + k];
        const i = OrderedSet.getAt(indices, t);
        posBuf.push(x[i], y[i], z[i]);
        avgRadius += radius(i);
    }
    avgRadius = avgRadius / m + radiusOffset;

    const axes = PrincipalAxes.calculateMomentsAxes(posBuf);
    const center = Vec3.clone(axes.origin);

    // representative atom: nearest to the fitted center
    let reprT = group.array[group.offset];
    let bestDistSq = Infinity;
    for (let k = 0; k < group.count; ++k) {
        const t = group.array[group.offset + k];
        const i = OrderedSet.getAt(indices, t);
        const dx = x[i] - center[0], dy = y[i] - center[1], dz = z[i] - center[2];
        const dSq = dx * dx + dy * dy + dz * dz;
        if (dSq < bestDistSq) {
            bestDistSq = dSq;
            reprT = t;
        }
    }
    const atomId = id ? id[reprT] : reprT;

    const u0 = Vec3.clone(axes.dirA), u1 = Vec3.clone(axes.dirB), u2 = Vec3.clone(axes.dirC);
    let a0: number, a1: number, a2: number;

    if (m === 1) {
        Vec3.set(u0, 1, 0, 0);
        Vec3.set(u1, 0, 1, 0);
        Vec3.set(u2, 0, 0, 1);
        a0 = a1 = a2 = avgRadius;
    } else {
        const s0 = Vec3.magnitude(u0), s1 = Vec3.magnitude(u1), s2 = Vec3.magnitude(u2);
        if (s0 > 1e-6) Vec3.scale(u0, u0, 1 / s0); else Vec3.set(u0, 1, 0, 0);
        if (s1 > 1e-6) Vec3.scale(u1, u1, 1 / s1); else Vec3.set(u1, 0, 1, 0);
        if (s2 > 1e-6) Vec3.scale(u2, u2, 1 / s2); else Vec3.set(u2, 0, 0, 1);
        a0 = (s0 > 1e-6 ? s0 : 0) + avgRadius;
        a1 = (s1 > 1e-6 ? s1 : 0) + avgRadius;
        a2 = (s2 > 1e-6 ? s2 : 0) + avgRadius;
    }

    return { kind: 'ellipsoid', center, u0, u1, u2, a0, a1, a2, atomId };
}

/**
 * Fits a spherical-harmonics radial boundary `R(theta, phi)` (about the group's centroid) to one
 * group of atoms via `fitSphericalHarmonics` (`./spherical-harmonics`), generalizing the
 * ellipsoid's fixed quadratic boundary into an angularly-varying one that can (for `degree > 2`)
 * hug non-ellipsoidal, even mildly concave, local atom-cluster shapes. Each atom contributes one
 * sample point placed at `center + direction * targetRadius` (`targetRadius` = distance to the
 * atom + its own radius + `radiusOffset`), so the fit already sees the desired inflated boundary
 * directly as its sample radii, rather than the bare atom positions. `fitSphericalHarmonics`'s
 * Tikhonov regularization (`regularization`) keeps degenerate/sparse groups (fewer atoms than
 * SH coefficients, collinear atoms, single-atom groups, etc.) well-defined instead of
 * oscillating/overfitting.
 */
function fitSHBlobFromGroup(position: PositionData, radius: (index: number) => number, radiusOffset: number, group: AtomGroup, degree: number, regularization: number): Blob {
    const { indices, x, y, z, id } = position;
    const m = group.count;

    const center = Vec3();
    for (let k = 0; k < group.count; ++k) {
        const t = group.array[group.offset + k];
        const i = OrderedSet.getAt(indices, t);
        center[0] += x[i]; center[1] += y[i]; center[2] += z[i];
    }
    Vec3.scale(center, center, 1 / m);

    // inflated sample points `center + direction * targetRadius`, so `fitSphericalHarmonics`
    // (which samples r = |point - center|) sees the desired boundary radius directly
    const samples = new Float64Array(m * 3);
    let reprT = group.array[group.offset];
    let bestDistSq = Infinity;

    for (let k = 0; k < group.count; ++k) {
        const t = group.array[group.offset + k];
        const i = OrderedSet.getAt(indices, t);
        const dx = x[i] - center[0], dy = y[i] - center[1], dz = z[i] - center[2];
        const dSq = dx * dx + dy * dy + dz * dz;
        if (dSq < bestDistSq) {
            bestDistSq = dSq;
            reprT = t;
        }

        const dist = Math.sqrt(dSq);
        const targetR = dist + radius(i) + radiusOffset;
        const ux = dist > 1e-8 ? dx / dist : 1, uy = dist > 1e-8 ? dy / dist : 0, uz = dist > 1e-8 ? dz / dist : 0;
        samples[k * 3] = center[0] + ux * targetR;
        samples[k * 3 + 1] = center[1] + uy * targetR;
        samples[k * 3 + 2] = center[2] + uz * targetR;
    }

    const fit = fitSphericalHarmonics(samples, center, degree, undefined, regularization);
    const atomId = id ? id[reprT] : reprT;

    // hard ceiling (used both for grid padding and, in computeBlobSurface, to clamp R(theta, phi)
    // itself) over the max observed target radius - bounds any residual fit overshoot in
    // unsampled directions to a small, fixed margin regardless of degree/sample distribution
    return { kind: 'sh', center, degree, coeffs: fit.coeffs, maxR: fit.rMax * SH_OVERSHOOT_CLAMP, atomId };
}

/**
 * Coarsens atoms into a small number of "blobs": atoms are first grouped either by a uniform
 * spatial hash grid (`method: 'grid'`, single O(atoms) pass) or refined further via k-means
 * clustering seeded from that same grid (`method: 'clustering'`, softens the grid's
 * atoms-near-a-bin-boundary-split artifact at the cost of `clusterIterations` extra passes), then
 * each group's atoms are fit to either an ellipsoid (`shape: 'ellipsoid'`) or a
 * spherical-harmonics radial boundary (`shape: 'sh'`).
 */
function fitBlobs(position: PositionData, boundary: Boundary, radius: (index: number) => number, props: BlobSurfaceProps): { blobs: Blob[], maxSemiAxis: number } {
    const { blobSize, radiusOffset, method, clusterIterations, shape, shDegree, shRegularization } = props;

    const groups = method === 'clustering'
        ? groupAtomsByClustering(position, boundary, blobSize, clusterIterations)
        : groupAtomsByGrid(position, boundary, blobSize);

    const blobs: Blob[] = shape === 'sh'
        ? groups.map(group => fitSHBlobFromGroup(position, radius, radiusOffset, group, shDegree, shRegularization))
        : groups.map(group => fitEllipsoidBlobFromGroup(position, radius, radiusOffset, group));

    let maxSemiAxis = 0;
    for (const b of blobs) {
        const maxA = blobMaxExtent(b);
        if (maxA > maxSemiAxis) maxSemiAxis = maxA;
    }

    return { blobs, maxSemiAxis };
}


/** half-widths, along world x/y/z, of the axis-aligned bounding box of an ellipsoid scaled by `cutoff` */
function ellipsoidWorldPad(out: Vec3, u0: Vec3, u1: Vec3, u2: Vec3, a0: number, a1: number, a2: number, cutoff: number) {
    for (let axis = 0; axis < 3; ++axis) {
        const e0 = a0 * u0[axis], e1 = a1 * u1[axis], e2 = a2 * u2[axis];
        out[axis] = cutoff * Math.sqrt(e0 * e0 + e1 * e1 + e2 * e2);
    }
    return out;
}

/**
 * Estimates the density of the atom point cloud with a small number of ellipsoid "blobs"
 * (`fitBlobs`), then accumulates their union directly into a dense grid. The grid box is the
 * TIGHT union of each individual blob's own padded extent (`center ± pad`, not the atoms'
 * bounding box expanded uniformly by the single largest blob) - this avoids one big/outlier blob
 * forcing wasted padding across the entire box, which matters most for sparse/spread-out
 * structures (elongated chains, multi-domain assemblies) where most of the atoms' own bounding
 * box is otherwise empty space between blobs. `resolution` (grid cell size) is left entirely to
 * the caller - unlike `blobSize`, it's not auto-derived here, so callers that want a smooth
 * high-resolution mesh over big blobs can still ask for one; keeping the box tight is what keeps
 * that affordable instead of also depending on a coarse resolution to avoid blowing up cell count.
 */
export function computeBlobSurface(position: PositionData, boundary: Boundary, radius: (index: number) => number, props: BlobSurfaceProps): BlobSurfaceData {
    const { resolution, smoothness, shDegree } = props;
    const alpha = smoothness;
    const scaleFactor = 1 / resolution;

    const { blobs } = fitBlobs(position, boundary, radius, props);

    const cutoff = 2; // matches the existing "2x radius" density cutoff convention
    const cutoffSq = cutoff * cutoff;

    // per-blob world-space padding, computed once and reused both for the tight box union below
    // and for the splat loop's local cell range (avoids recomputing it twice per blob)
    const blobPads: Vec3[] = [];
    let maxSemiAxis = 0;
    const boxMin = Vec3.create(Infinity, Infinity, Infinity);
    const boxMax = Vec3.create(-Infinity, -Infinity, -Infinity);

    for (const blob of blobs) {
        const maxExtent = blobMaxExtent(blob);
        if (maxExtent > maxSemiAxis) maxSemiAxis = maxExtent;

        const padVec = Vec3();
        if (blob.kind === 'ellipsoid') {
            ellipsoidWorldPad(padVec, blob.u0, blob.u1, blob.u2, blob.a0, blob.a1, blob.a2, cutoff);
        } else {
            Vec3.set(padVec, maxExtent * cutoff, maxExtent * cutoff, maxExtent * cutoff);
        }
        blobPads.push(padVec);

        for (let axis = 0; axis < 3; ++axis) {
            const lo = blob.center[axis] - padVec[axis];
            const hi = blob.center[axis] + padVec[axis];
            if (lo < boxMin[axis]) boxMin[axis] = lo;
            if (hi > boxMax[axis]) boxMax[axis] = hi;
        }
    }

    const expandedBox = blobs.length > 0
        ? Box3D.create(
            Vec3.create(boxMin[0] - resolution, boxMin[1] - resolution, boxMin[2] - resolution),
            Vec3.create(boxMax[0] + resolution, boxMax[1] + resolution, boxMax[2] + resolution)
        )
        : Box3D.expand(Box3D(), boundary.box, Vec3.create(resolution, resolution, resolution));

    const min = expandedBox.min;
    const scaledBox = Box3D.scale(Box3D(), expandedBox, scaleFactor);
    const dim = Box3D.size(Vec3(), scaledBox);
    Vec3.ceil(dim, dim);
    // console.log(`computeBlobSurface: ${blobs.length} blobs, grid dim = ${dim[0]} x ${dim[1]} x ${dim[2]}`);

    const space = Tensor.Space(dim, [0, 1, 2], Float32Array);
    const data = space.create();
    const field = Tensor.create(space, data);

    const idData = space.create();
    idData.fill(-1);
    const idField = Tensor.create(space, idData);

    const [dimX, dimY, dimZ] = dim;
    const iu = dimZ, iv = dimY, iuv = iu * iv;

    const gridx = fillGridDim(dim[0], min[0], resolution);
    const gridy = fillGridDim(dim[1], min[1], resolution);
    const gridz = fillGridDim(dim[2], min[2], resolution);

    const densData = space.create();
    const sphericalScratch: SphericalCoord = { r: 0, theta: 0, phi: 0 };
    const shBasisScratch = new Float64Array(shTermCount(shDegree));
    const legendreScratch = new Float64Array(((shDegree + 1) * (shDegree + 2)) / 2);

    for (let bi = 0; bi < blobs.length; ++bi) {
        const blob = blobs[bi];
        const { center, atomId } = blob;
        const maxExtent = blobMaxExtent(blob);
        const padVec = blobPads[bi];

        const iax = Math.floor(scaleFactor * (center[0] - min[0]));
        const iay = Math.floor(scaleFactor * (center[1] - min[1]));
        const iaz = Math.floor(scaleFactor * (center[2] - min[2]));

        const ngx = Math.ceil(padVec[0] * scaleFactor);
        const ngy = Math.ceil(padVec[1] * scaleFactor);
        const ngz = Math.ceil(padVec[2] * scaleFactor);

        const begX = Math.max(0, iax - ngx);
        const begY = Math.max(0, iay - ngy);
        const begZ = Math.max(0, iaz - ngz);
        const endX = Math.min(dimX, iax + ngx + 2);
        const endY = Math.min(dimY, iay + ngy + 2);
        const endZ = Math.min(dimZ, iaz + ngz + 2);

        if (blob.kind === 'ellipsoid') {
            const { u0, u1, u2, a0, a1, a2 } = blob;
            const a0Sq = a0 * a0, a1Sq = a1 * a1, a2Sq = a2 * a2;

            for (let xi = begX; xi < endX; ++xi) {
                const dx = gridx[xi] - center[0];
                const xIdx = xi * iuv;
                for (let yi = begY; yi < endY; ++yi) {
                    const dy = gridy[yi] - center[1];
                    const xyIdx = yi * iu + xIdx;
                    for (let zi = begZ; zi < endZ; ++zi) {
                        const dz = gridz[zi] - center[2];

                        const l0 = dx * u0[0] + dy * u0[1] + dz * u0[2];
                        const l1 = dx * u1[0] + dy * u1[1] + dz * u1[2];
                        const l2 = dx * u2[0] + dy * u2[1] + dz * u2[2];
                        const q = (l0 * l0) / a0Sq + (l1 * l1) / a1Sq + (l2 * l2) / a2Sq;

                        if (q <= cutoffSq) {
                            const dens = fasterExp(-alpha * q);
                            const idx = zi + xyIdx;
                            data[idx] += dens;
                            if (dens > densData[idx]) {
                                densData[idx] = dens;
                                idData[idx] = atomId;
                            }
                        }
                    }
                }
            }
        } else {
            const { degree, coeffs } = blob;
            const minR = 1e-3;

            for (let xi = begX; xi < endX; ++xi) {
                const dx = gridx[xi] - center[0];
                const xIdx = xi * iuv;
                for (let yi = begY; yi < endY; ++yi) {
                    const dy = gridy[yi] - center[1];
                    const xyIdx = yi * iu + xIdx;
                    for (let zi = begZ; zi < endZ; ++zi) {
                        const dz = gridz[zi] - center[2];

                        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                        if (dist > maxExtent * cutoff) continue;

                        toSpherical(dx, dy, dz, sphericalScratch);
                        let r = reconstructRadius(coeffs, degree, sphericalScratch.theta, sphericalScratch.phi, shBasisScratch, legendreScratch);
                        if (r < minR) r = minR;
                        else if (r > maxExtent) r = maxExtent; // hard ceiling: no lobe beyond a small margin over the real data extent

                        const q = (dist * dist) / (r * r);
                        if (q <= cutoffSq) {
                            const dens = fasterExp(-alpha * q);
                            const idx = zi + xyIdx;
                            data[idx] += dens;
                            if (dens > densData[idx]) {
                                densData[idx] = dens;
                                idData[idx] = atomId;
                            }
                        }
                    }
                }
            }
        }
    }

    const transform = Mat4.identity();
    Mat4.fromScaling(transform, Vec3.create(resolution, resolution, resolution));
    Mat4.setTranslation(transform, expandedBox.min);

    return { field, idField, transform, radiusFactor: 1, resolution, maxRadius: maxSemiAxis };
}
