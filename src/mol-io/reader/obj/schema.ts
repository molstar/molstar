/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Intermediate representation of a parsed OBJ file.
 *
 * Positions and normals are stored as raw arrays from the file.
 * Faces have been triangulated (fan-triangulation) and are stored as
 * separate flat arrays of per-triangle-vertex position and normal indices.
 * Normal indices of -1 indicate that the vertex has no explicit normal.
 */
export interface ObjFile {
    /** Raw position data from `v` lines, interleaved [x0,y0,z0, x1,y1,z1, ...] */
    readonly positions: Float32Array
    /** Raw normal data from `vn` lines, interleaved [nx0,ny0,nz0, ...]. Length 0 if no normals. */
    readonly normals: Float32Array

    /**
     * Per-face-vertex position index (0-based), length = triangleCount * 3.
     * Three consecutive values define one triangle: [p0, p1, p2, p3, p4, p5, ...]
     */
    readonly positionIndices: Int32Array

    /**
     * Per-face-vertex normal index (0-based), length = triangleCount * 3.
     * -1 means no explicit normal for that face-vertex.
     */
    readonly normalIndices: Int32Array

    readonly positionCount: number
    readonly normalCount: number
    readonly triangleCount: number

    /**
     * Unique material names encountered via `usemtl` directives, in encounter order.
     * Empty when no `usemtl` directives are present.
     */
    readonly materialNames: readonly string[]

    /**
     * Per-triangle material index (0-based, indexes into materialNames), length = triangleCount.
     * 0 for faces before the first `usemtl` or when no materials are defined.
     */
    readonly faceGroups: Int32Array
}
