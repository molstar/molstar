/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Build a mol* Mesh / Lines from a parsed SffData (HFF).
 *
 * Concatenates all meshes from all segments into a single geometry, with the
 * per-vertex `groups` channel encoding the segment index. Per-group colour
 * and label lookups index back into SffData.segments. The output geometry is
 * consumed by the multi-visual representation in `representation.ts`.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
import { Color } from '../../mol-util/color';
import { Mat4 } from '../../mol-math/linear-algebra';
import { SffData, SffSegment, SffTransform } from '../../mol-io/reader/hff/schema';

function transformMatrix(t: SffTransform | undefined): Mat4 | undefined {
    if (!t || t.data.length === 0) return undefined;
    const m = Mat4.identity();
    if (t.rows === 3 && t.cols === 4) {
        for (let r = 0; r < 3; r++) {
            for (let c = 0; c < 4; c++) {
                Mat4.setValue(m, r, c, t.data[r * 4 + c]);
            }
        }
        return m;
    }
    if (t.rows === 4 && t.cols === 4) {
        for (let r = 0; r < 4; r++) {
            for (let c = 0; c < 4; c++) {
                Mat4.setValue(m, r, c, t.data[r * 4 + c]);
            }
        }
        return m;
    }
    return undefined;
}

function applyMat4(out: Float32Array, off: number, x: number, y: number, z: number, m: Mat4 | undefined) {
    if (!m) {
        out[off] = x; out[off + 1] = y; out[off + 2] = z;
        return;
    }
    out[off] = m[0] * x + m[4] * y + m[8] * z + m[12];
    out[off + 1] = m[1] * x + m[5] * y + m[9] * z + m[13];
    out[off + 2] = m[2] * x + m[6] * y + m[10] * z + m[14];
}

function rotateOnly(m: Mat4, scratch: Mat4): Mat4 {
    Mat4.copy(scratch, m);
    scratch[12] = 0; scratch[13] = 0; scratch[14] = 0;
    return scratch;
}

function findTransformById(transforms: SffTransform[], id: number | undefined): SffTransform | undefined {
    if (id === undefined) return undefined;
    return transforms.find(t => t.id === id);
}

export interface BuiltMesh {
    mesh: Mesh;
    /** Per-group-id segment lookup (group id == segment index in SffData.segments). */
    segmentByGroup: SffSegment[];
}

export function buildMesh(data: SffData): BuiltMesh {
    let totalV = 0, totalT = 0;
    for (const seg of data.segments) {
        for (const m of seg.meshes) {
            totalV += m.vertices.count;
            totalT += m.triangles.count;
        }
    }

    const vertices = new Float32Array(totalV * 3);
    const indices = new Uint32Array(totalT * 3);
    const normals = new Float32Array(totalV * 3);
    const groups = new Float32Array(totalV);

    const rotScratch = Mat4();
    let vBase = 0;
    let iOff = 0;
    let anyMissingNormals = false;
    const segmentByGroup: SffSegment[] = [];

    for (let segIdx = 0; segIdx < data.segments.length; segIdx++) {
        const seg = data.segments[segIdx];
        for (const mesh of seg.meshes) {
            const tr = transformMatrix(findTransformById(data.transforms, mesh.transformId));
            const rot = tr ? rotateOnly(tr, rotScratch) : undefined;

            const vSrc = mesh.vertices.data as ArrayLike<number>;
            const vCount = mesh.vertices.count;
            for (let i = 0; i < vCount; i++) {
                applyMat4(vertices, (vBase + i) * 3, vSrc[i * 3], vSrc[i * 3 + 1], vSrc[i * 3 + 2], tr);
                groups[vBase + i] = segIdx;
            }

            if (mesh.normals && mesh.normals.count === vCount) {
                const nSrc = mesh.normals.data as ArrayLike<number>;
                for (let i = 0; i < vCount; i++) {
                    applyMat4(normals, (vBase + i) * 3, nSrc[i * 3], nSrc[i * 3 + 1], nSrc[i * 3 + 2], rot);
                }
            } else {
                anyMissingNormals = true;
            }

            const tSrc = mesh.triangles.data as ArrayLike<number>;
            const tElems = mesh.triangles.count * 3;
            for (let i = 0; i < tElems; i++) {
                indices[iOff + i] = vBase + tSrc[i];
            }

            vBase += vCount;
            iOff += tElems;
        }
        segmentByGroup.push(seg);
    }

    const mesh = Mesh.create(vertices, indices, normals, groups, totalV, totalT);
    if (anyMissingNormals) Mesh.computeNormals(mesh);

    return { mesh, segmentByGroup };
}

export function buildLinesFromMesh(built: BuiltMesh): Lines {
    const m = built.mesh;
    const verts = m.vertexBuffer.ref.value;
    const tris = m.indexBuffer.ref.value;
    const groups = m.groupBuffer.ref.value;
    const triCount = m.triangleCount;

    const builder = LinesBuilder.create(triCount * 3);
    for (let t = 0; t < triCount; t++) {
        const a = tris[t * 3], b = tris[t * 3 + 1], c = tris[t * 3 + 2];
        const ax = verts[a * 3], ay = verts[a * 3 + 1], az = verts[a * 3 + 2];
        const bx = verts[b * 3], by = verts[b * 3 + 1], bz = verts[b * 3 + 2];
        const cx = verts[c * 3], cy = verts[c * 3 + 1], cz = verts[c * 3 + 2];
        const g = groups[a];
        builder.add(ax, ay, az, bx, by, bz, g);
        builder.add(bx, by, bz, cx, cy, cz, g);
        builder.add(cx, cy, cz, ax, ay, az, g);
    }
    return builder.getLines();
}

export function colourToColor(c: [number, number, number, number]): Color {
    return Color.fromNormalizedRgb(c[0], c[1], c[2]);
}

export function segmentLabel(seg: SffSegment): string {
    return seg.biologicalAnnotation?.name?.trim() || `Segment ${seg.id}`;
}

// re-exported for unit tests
export const _internals = { buildMesh, buildLinesFromMesh };
