/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../linear-algebra/3d/vec3';

/**
 * Incremental 3D convex hull for small point sets (typically 4–12 points).
 *
 * Returns triangle indices with outward-facing normals.
 * Degenerate cases (coplanar/collinear points) are handled gracefully:
 *   - < 4 points or all coplanar → returns empty hull (no triangles)
 *   - Duplicate points are tolerated; the resulting hull is still valid
 */

const EPSILON = 1e-8;

const _ab = Vec3();
const _ac = Vec3();
const _normal = Vec3();
const _ap = Vec3();

interface ConvexHullResult {
    /** Triangle indices into the original positions array (length = numTriangles * 3) */
    indices: number[]
}

/**
 * Compute 3D convex hull of a set of points.
 * @param positions Array of Vec3 positions
 * @returns Triangle indices with outward-facing normals, or empty if degenerate
 */
export function convexHull(positions: ArrayLike<Vec3>): ConvexHullResult {
    const n = positions.length;

    if (n < 4) {
        return { indices: [] };
    }

    // Find initial tetrahedron
    const tet = findInitialTetrahedron(positions);
    if (!tet) {
        return { indices: [] };
    }

    const [i0, i1, i2, i3] = tet;

    // Faces are stored as triplets of indices into `positions`.
    // Each face's vertices must be ordered so the outward normal = cross(v1-v0, v2-v0).
    // We orient the initial tetrahedron so all faces point outward.

    // Check orientation of first face relative to 4th point
    Vec3.sub(_ab, positions[i1], positions[i0]);
    Vec3.sub(_ac, positions[i2], positions[i0]);
    Vec3.cross(_normal, _ab, _ac);
    Vec3.sub(_ap, positions[i3], positions[i0]);
    const dot = Vec3.dot(_normal, _ap);

    // Faces of tetrahedron. If the 4th point is on the positive side of face (i0,i1,i2),
    // we need to flip that face.
    let faces: number[][];
    if (dot > 0) {
        // i3 is on positive side of (i0,i1,i2) → flip first face
        faces = [
            [i0, i2, i1],
            [i0, i1, i3],
            [i1, i2, i3],
            [i0, i3, i2],
        ];
    } else {
        faces = [
            [i0, i1, i2],
            [i0, i3, i1],
            [i1, i3, i2],
            [i0, i2, i3],
        ];
    }

    // Incrementally add each remaining point
    for (let pi = 0; pi < n; ++pi) {
        if (pi === i0 || pi === i1 || pi === i2 || pi === i3) continue;

        const p = positions[pi];

        // Find visible faces (point is on positive side of face)
        const visible: boolean[] = [];
        let anyVisible = false;
        for (let fi = 0; fi < faces.length; fi++) {
            const face = faces[fi];
            if (isFaceVisibleFrom(positions, face, p)) {
                visible.push(true);
                anyVisible = true;
            } else {
                visible.push(false);
            }
        }
        if (!anyVisible) continue; // Point is inside current hull

        // Find horizon edges: edges shared between a visible and non-visible face.
        // An edge [a,b] in a visible face that also appears as [b,a] in a non-visible face is a horizon edge.
        const horizon: number[][] = [];
        for (let fi = 0; fi < faces.length; fi++) {
            if (!visible[fi]) continue;

            const face = faces[fi];
            for (let ei = 0; ei < 3; ei++) {
                const a = face[ei];
                const b = face[(ei + 1) % 3];
                // Check if the reverse edge (b,a) belongs to a non-visible face
                let isHorizon = false;
                for (let fj = 0; fj < faces.length; fj++) {
                    if (visible[fj]) continue;

                    const other = faces[fj];
                    for (let ej = 0; ej < 3; ej++) {
                        if (other[ej] === b && other[(ej + 1) % 3] === a) {
                            isHorizon = true;
                            break;
                        }
                    }
                    if (isHorizon) break;
                }
                if (isHorizon) {
                    horizon.push([a, b]);
                }
            }
        }

        // Remove visible faces
        const newFaces: number[][] = [];
        for (let fi = 0; fi < faces.length; fi++) {
            if (!visible[fi]) {
                newFaces.push(faces[fi]);
            }
        }

        // Create new faces from horizon edges to the new point
        for (const [a, b] of horizon) {
            newFaces.push([a, b, pi]);
        }

        faces = newFaces;
    }

    // Flatten faces into indices
    const indices: number[] = [];
    for (const face of faces) {
        indices.push(face[0], face[1], face[2]);
    }

    return { indices };
}

function isFaceVisibleFrom(positions: ArrayLike<Vec3>, face: number[], point: Vec3): boolean {
    const a = positions[face[0]];
    const b = positions[face[1]];
    const c = positions[face[2]];
    Vec3.sub(_ab, b, a);
    Vec3.sub(_ac, c, a);
    Vec3.cross(_normal, _ab, _ac);
    Vec3.sub(_ap, point, a);
    return Vec3.dot(_normal, _ap) > EPSILON;
}

/**
 * Find 4 non-coplanar points for the initial tetrahedron.
 * Returns indices [i0, i1, i2, i3] or null if all points are coplanar.
 */
function findInitialTetrahedron(positions: ArrayLike<Vec3>): [number, number, number, number] | null {
    const n = positions.length;

    // Find two distinct points
    const i0 = 0;
    let i1 = -1;
    for (let i = 1; i < n; i++) {
        if (Vec3.distance(positions[i0], positions[i]) > EPSILON) {
            i1 = i;
            break;
        }
    }
    if (i1 < 0) return null;

    // Find a third point not collinear with the first two
    let i2 = -1;
    let maxArea = 0;
    for (let i = 0; i < n; i++) {
        if (i === i0 || i === i1) continue;

        Vec3.sub(_ab, positions[i1], positions[i0]);
        Vec3.sub(_ac, positions[i], positions[i0]);
        Vec3.cross(_normal, _ab, _ac);
        const area = Vec3.magnitude(_normal);
        if (area > maxArea) {
            maxArea = area;
            i2 = i;
        }
    }
    if (i2 < 0 || maxArea < EPSILON) return null;

    // Find a fourth point not coplanar with the first three
    let i3 = -1;
    let maxVol = 0;
    Vec3.sub(_ab, positions[i1], positions[i0]);
    Vec3.sub(_ac, positions[i2], positions[i0]);
    Vec3.cross(_normal, _ab, _ac);
    Vec3.normalize(_normal, _normal);

    for (let i = 0; i < n; i++) {
        if (i === i0 || i === i1 || i === i2) continue;

        Vec3.sub(_ap, positions[i], positions[i0]);
        const vol = Math.abs(Vec3.dot(_normal, _ap));
        if (vol > maxVol) {
            maxVol = vol;
            i3 = i;
        }
    }
    if (i3 < 0 || maxVol < EPSILON) return null;

    return [i0, i1, i2, i3];
}
