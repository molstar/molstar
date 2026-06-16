/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Reader for EMDB-SFF segmentation files in their HDF5 serialization (.hff).
 * Walks the HDF5 tree using the vendored jsfive reader and produces a typed
 * SffData structure that downstream transforms turn into mol* shapes.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Task, RuntimeContext } from '../../../mol-task';
import { ReaderResult as Result } from '../result';
import { File as H5File } from '../../common/h5';
import { SffData, SffMesh, SffSegment, SffEncodedSequence, SffMode, SffEndianness, SffTransform, SffColour, SffPrimaryDescriptor, SffBiologicalAnnotation } from './schema';

const MODE_BYTES: Record<SffMode, number> = {
    int8: 1, uint8: 1,
    int16: 2, uint16: 2,
    int32: 4, uint32: 4,
    int64: 8, uint64: 8,
    float32: 4, float64: 8,
};

function decodeBase64(b64: string): Uint8Array {
    if (typeof Buffer !== 'undefined') {
        const buf = Buffer.from(b64, 'base64');
        return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
    }
    const bin = atob(b64);
    const out = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
}

function byteSwap(buf: Uint8Array, itemSize: number) {
    if (itemSize <= 1) return;
    for (let i = 0; i < buf.length; i += itemSize) {
        for (let a = 0, b = itemSize - 1; a < b; a++, b--) {
            const t = buf[i + a];
            buf[i + a] = buf[i + b];
            buf[i + b] = t;
        }
    }
}

function makeTypedArray(mode: SffMode, bytes: Uint8Array): SffEncodedSequence['data'] {
    // bytes is in machine-native order at this point (caller swapped if needed).
    // Use a fresh ArrayBuffer to avoid alignment issues against the b64-decoded view.
    const ab = new ArrayBuffer(bytes.byteLength);
    new Uint8Array(ab).set(bytes);
    switch (mode) {
        case 'int8': return new Int8Array(ab);
        case 'uint8': return new Uint8Array(ab);
        case 'int16': return new Int16Array(ab);
        case 'uint16': return new Uint16Array(ab);
        case 'int32': return new Int32Array(ab);
        case 'uint32': return new Uint32Array(ab);
        case 'float32': return new Float32Array(ab);
        case 'float64': return new Float64Array(ab);
        case 'int64':
        case 'uint64':
            // Mol* mesh consumers only want 32-bit indices; downcast 64-bit modes.
            // SFF rarely uses these for vertices/triangles, but be defensive.
            throw new Error(`SFF mode '${mode}' not supported in v1 (use 32-bit modes)`);
    }
}

/** Safe get: jsfive's Group.get throws on missing keys; we want a soft lookup. */
function tryGet(group: any, key: string): any {
    if (!group) return undefined;
    try {
        return group.get(key);
    } catch {
        return undefined;
    }
}

function readScalar(node: any): any {
    if (node === undefined || node === null) return undefined;
    const v = node.value;
    if (Array.isArray(v) || ArrayBuffer.isView(v)) {
        return (v as ArrayLike<any>).length === 1 ? (v as ArrayLike<any>)[0] : v;
    }
    return v;
}

function readString(node: any): string | undefined {
    const v = readScalar(node);
    if (v === undefined) return undefined;
    if (typeof v === 'string') return v;
    // jsfive sometimes returns Uint8Array for strings; decode UTF-8.
    if (v instanceof Uint8Array) return new TextDecoder('utf-8').decode(v);
    return String(v);
}

function readNumber(node: any): number | undefined {
    const v = readScalar(node);
    return v === undefined ? undefined : Number(v);
}

function readEncodedSequence(group: any): SffEncodedSequence | undefined {
    if (!group) return undefined;
    const mode = readString(tryGet(group, 'mode')) as SffMode | undefined;
    const endianness = (readString(tryGet(group, 'endianness')) as SffEndianness | undefined) ?? 'little';
    const dataStr = readString(tryGet(group, 'data'));
    if (!mode || !dataStr) return undefined;

    let count = readNumber(tryGet(group, 'num_vertices'));
    if (count === undefined) count = readNumber(tryGet(group, 'num_triangles'));
    if (count === undefined) count = readNumber(tryGet(group, 'num_normals'));
    if (count === undefined) throw new Error('SFF encoded sequence missing num_* count');

    const raw = decodeBase64(dataStr);
    const itemSize = MODE_BYTES[mode];
    const expected = count * 3 * itemSize;
    if (raw.byteLength !== expected) {
        throw new Error(`SFF data size mismatch: expected ${expected} bytes (${count} × 3 × ${itemSize}), got ${raw.byteLength}`);
    }
    if (endianness === 'big') byteSwap(raw, itemSize);
    return { count, mode, endianness, data: makeTypedArray(mode, raw) };
}

function readColour(group: any): SffColour {
    const c = tryGet(group, 'colour');
    if (!c) return [0.7, 0.7, 0.7, 1.0];
    const v = c.value;
    return [Number(v[0]), Number(v[1]), Number(v[2]), Number(v[3] ?? 1.0)];
}

function readBiologicalAnnotation(group: any): SffBiologicalAnnotation | undefined {
    const ba = tryGet(group, 'biological_annotation');
    if (!ba) return undefined;
    return {
        name: readString(tryGet(ba, 'name')),
        description: readString(tryGet(ba, 'description')),
        numberOfInstances: readNumber(tryGet(ba, 'number_of_instances')),
    };
}

function readMeshes(segmentGroup: any): SffMesh[] {
    const ml = tryGet(segmentGroup, 'mesh_list');
    if (!ml) return [];
    const keys = (ml.keys as string[]).slice().sort((a, b) => Number(a) - Number(b));
    const out: SffMesh[] = [];
    for (const k of keys) {
        const mg = ml.get(k);
        const id = readNumber(tryGet(mg, 'id')) ?? Number(k);
        const transformId = readNumber(tryGet(mg, 'transform_id'));
        const vertices = readEncodedSequence(tryGet(mg, 'vertices'));
        const triangles = readEncodedSequence(tryGet(mg, 'triangles'));
        const normals = readEncodedSequence(tryGet(mg, 'normals'));
        if (!vertices || !triangles) continue;
        out.push({ id, vertices, triangles, normals, transformId });
    }
    return out;
}

function readTransforms(root: any): SffTransform[] {
    const tl = tryGet(root, 'transform_list');
    if (!tl) return [];
    const keys = (tl.keys as string[]).slice().sort((a, b) => Number(a) - Number(b));
    const out: SffTransform[] = [];
    for (const k of keys) {
        const tg = tl.get(k);
        const id = readNumber(tryGet(tg, 'id')) ?? Number(k);
        const rows = readNumber(tryGet(tg, 'rows')) ?? 3;
        const cols = readNumber(tryGet(tg, 'cols')) ?? 4;
        const dataStr = readString(tryGet(tg, 'data')) ?? '';
        const nums = dataStr.trim().length === 0 ? [] : dataStr.trim().split(/\s+/).map(Number);
        out.push({ id, rows, cols, data: new Float64Array(nums) });
    }
    return out;
}

function readSegments(root: any): SffSegment[] {
    const sl = tryGet(root, 'segment_list');
    if (!sl) return [];
    const keys = (sl.keys as string[]).slice().sort((a, b) => Number(a) - Number(b));
    const out: SffSegment[] = [];
    for (const k of keys) {
        const sg = sl.get(k);
        const id = readNumber(tryGet(sg, 'id')) ?? Number(k);
        const parentId = readNumber(tryGet(sg, 'parent_id')) ?? 0;
        const colour = readColour(sg);
        const biologicalAnnotation = readBiologicalAnnotation(sg);
        const meshes = readMeshes(sg);
        out.push({ id, parentId, colour, biologicalAnnotation, meshes });
    }
    return out;
}

export async function parseHffData(buffer: ArrayBuffer): Promise<SffData> {
    const file = new (H5File as any)(buffer, '');
    return {
        version: readString(tryGet(file, 'version')),
        name: readString(tryGet(file, 'name')),
        details: readString(tryGet(file, 'details')),
        primaryDescriptor: readString(tryGet(file, 'primary_descriptor')) as SffPrimaryDescriptor | undefined,
        transforms: readTransforms(file),
        segments: readSegments(file),
    };
}

export function parseHff(data: ArrayBuffer | Uint8Array) {
    return Task.create<Result<SffData>>('Parse HFF (EMDB-SFF)', async (ctx: RuntimeContext) => {
        try {
            await ctx.update('Reading HDF5...');
            const src = data instanceof Uint8Array ? data : new Uint8Array(data);
            const ab = new ArrayBuffer(src.byteLength);
            new Uint8Array(ab).set(src);
            const result = await parseHffData(ab);
            return Result.success(result);
        } catch (e) {
            return Result.error(`HFF parse error: ${(e as Error).message ?? e}`);
        }
    });
}
