/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Build a mol* Mesh from parsed EMDB-SFF data and expose it as a `ShapeProvider`,
 * so that HFF follows the same path as the built-in shape formats:
 *
 *   Binary -> ParseHff -> ShapeFromSff -> Shape.Provider -> ShapeRepresentation3D
 *
 * All meshes of all segments are concatenated into a single geometry, with the
 * per-vertex `groups` channel encoding the segment index; per-group colour and
 * label lookups index back into `SffData.segments`.
 *
 * The equivalent for the built-in formats lives in `mol-model-formats/shape/`
 * (see `ply.ts`, which likewise pairs mesh building with its shape provider).
 * HFF keeps it here because the reader is an opt-in extension.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { Shape } from '../../mol-model/shape';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { Material } from '../../mol-util/material';
import { Mat4 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RuntimeContext, Task } from '../../mol-task';
import { SffData, SffSegment, SffTransform } from './schema';



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

export function colourToColor(c: [number, number, number, number]): Color {
    return Color.fromNormalizedRgb(c[0], c[1], c[2]);
}

export function segmentLabel(seg: SffSegment): string {
    return seg.biologicalAnnotation?.name?.trim() || `Segment ${seg.id}`;
}


export const SffShapeParams = {
    ...Mesh.Params,
    // SFF mesh segments are typically thin oriented surfaces (membranes, organelles),
    // so both faces need to be drawn ...
    doubleSided: PD.Boolean(true, BaseGeometry.CustomQualityParamInfo),
    // ... and the back face should keep the segment colour instead of the default
    // grey, otherwise a closed segment reads as two different objects. Mirrors
    // `getInteriorParam()` with `colorStrength` defaulted to 0 instead of 1.
    interior: PD.Group({
        color: PD.Color(Color.fromRgb(76, 76, 76)),
        colorStrength: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
        substance: Material.getParam(),
        substanceStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    }),
};
export type SffShapeParams = typeof SffShapeParams


const DefaultSegmentColor = Color.fromNormalizedRgb(0.7, 0.7, 0.7);

function makeShapeGetter(transforms?: Mat4[]) {
    let _data: SffData | undefined;
    let _shape: Shape<Mesh>;

    return async (_ctx: RuntimeContext, data: SffData, _props: PD.Values<SffShapeParams>, _shape0?: Shape<Mesh>) => {
        if (_data !== data) {
            const built = buildMesh(data);
            const colors = built.segmentByGroup.map(s => colourToColor(s.colour));
            const labels = built.segmentByGroup.map(s => segmentLabel(s));
            _shape = Shape.create(
                data.name?.trim() || 'EMDB-SFF',
                data,
                built.mesh,
                (g: number) => colors[g] ?? DefaultSegmentColor,
                () => 1,
                (g: number) => labels[g] ?? `Segment ${g}`,
                transforms,
            );
            _data = data;
        }
        return _shape;
    };
}

export function shapeFromSff(data: SffData, params?: { transforms?: Mat4[], label?: string }) {
    return Task.create<ShapeProvider<SffData, Mesh, SffShapeParams>>('Shape Provider', async _ctx => {
        return {
            label: params?.label || data.name?.trim() || 'EMDB-SFF',
            data,
            params: SffShapeParams,
            getShape: makeShapeGetter(params?.transforms),
            geometryUtils: Mesh.Utils,
        };
    });
}
