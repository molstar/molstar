/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { WebGPUTextureMeshGeometry } from './texture-mesh';
import { WebGPUTextureData } from './texture-data';
import { unpackRGBToInt } from '../../mol-util/number-packing';
import { GraphicsRenderObject } from '../render-object';
import { RenderableValues } from '../renderable/schema';

export interface WebGPUGeometry {
    /** Five vec4s: position, normal, end/radius direction, mapping, group/source/uv. */
    vertices: Float32Array
    indices: Uint32Array
    vertexCount: number
    kind: number
    sphereLods?: { first: number, count: number, min: number, max: number }[]
    native?: { device: GPUDevice, buffer: GPUBuffer }
}

export function value<T>(values: RenderableValues, key: string, fallback: T): T {
    return values[key]?.ref.value ?? fallback;
}

const empty = new Float32Array(0);
/** Color, size/marker/emission/instance, and substance overlay. */
export const WebGPUThemeStride = 12;

/** Build native vertex data; no GL resources or shader translation are involved. */
export function createWebGPUGeometry(object: GraphicsRenderObject): WebGPUGeometry {
    const v: RenderableValues = object.values;
    if (object.type === 'direct-volume') {
        throw new Error(`WebGPU ${object.type} rendering has not been implemented.`);
    }
    if (object.type === 'texture-mesh') return textureMesh(v);
    if (object.type === 'spheres') return spheres(v);
    if (object.type === 'cylinders') return cylinders(v);
    const positions = value(v, 'aPosition', empty);
    const starts = value(v, 'aStart', empty);
    const ends = value(v, 'aEnd', empty);
    const normals = value(v, 'aNormal', empty);
    const groups = value(v, 'aGroup', empty);
    const mappings = value(v, 'aMapping', empty);
    const coords = value(v, object.type === 'image' ? 'aUv' : 'aTexCoord', empty);
    const depths = value(v, 'aDepth', empty);
    const elements = value(v, 'elements', new Uint32Array(0));
    const drawCount = value(v, 'drawCount', 0);
    // Point sprites are explicit triangles: WebGPU has no programmable point size.
    const point = object.type === 'points';
    const sourceCount = point ? drawCount : (elements.length ? maxIndex(elements, drawCount) + 1 : drawCount);
    const vertexCount = sourceCount * (point ? 4 : 1);
    const vertices = new Float32Array(vertexCount * 20);
    const kind = object.type === 'mesh' ? 0 : object.type === 'lines' ? 1 : point ? 2 : object.type === 'text' ? 3 : 4;
    const corners = [-1, -1, 1, -1, -1, 1, 1, 1];
    for (let i = 0; i < vertexCount; i++) {
        const source = point ? Math.floor(i / 4) : i;
        const o = i * 20, p = source * 3;
        const pos = kind === 1 ? starts : positions;
        vertices.set(pos.subarray(p, p + 3), o);
        vertices[o + 3] = 1;
        vertices.set(normals.subarray(p, p + 3), o + 4);
        vertices.set(ends.subarray(p, p + 3), o + 8);
        vertices[o + 11] = depths[source] || 0;
        vertices[o + 12] = point ? corners[(i % 4) * 2] : mappings[source * 2] || 0;
        vertices[o + 13] = point ? corners[(i % 4) * 2 + 1] : mappings[source * 2 + 1] || 0;
        vertices[o + 16] = groups[source] || 0;
        vertices[o + 17] = source;
        vertices[o + 18] = coords[source * 2] || 0;
        vertices[o + 19] = coords[source * 2 + 1] || 0;
    }
    let indices: Uint32Array;
    if (point) {
        indices = new Uint32Array(sourceCount * 6);
        for (let i = 0; i < sourceCount; i++) indices.set([i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 2, i * 4 + 1, i * 4 + 3], i * 6);
    } else {
        indices = elements.length ? elements.slice(0, drawCount) : Uint32Array.from({ length: drawCount }, (_, i) => i);
    }
    return { vertices, indices, vertexCount, kind };
}

function textureMesh(v: RenderableValues): WebGPUGeometry {
    const count = value(v, 'drawCount', 0), dimensions = value(v, 'uGeoTexDim', [0, 0]);
    if (!Number.isSafeInteger(count) || count < 0 || count % 3 !== 0) throw new Error('WebGPU texture mesh must contain complete triangles.');
    const native = value<Record<string, unknown>>(v, 'meta', {}).webgpuGeometry;
    if (native instanceof WebGPUTextureMeshGeometry && native.matches(value(v, 'tPosition', undefined), value(v, 'tNormal', undefined), value(v, 'tGroup', undefined), count)) return native.geometry;
    const textures = ['tPosition', 'tNormal', 'tGroup'].map(key => {
        const texture = value<WebGPUTextureData | undefined>(v, key, undefined);
        if (!(texture instanceof WebGPUTextureData)) throw new Error('WebGPU texture mesh requires native texture data.');
        const data = texture.data;
        if (data.depth !== 1 || data.width !== dimensions[0] || data.height !== dimensions[1] || count > data.width * data.height) throw new Error('WebGPU texture mesh dimensions or vertex count are invalid.');
        return data.array;
    });
    const [positions, normals, groups] = textures;
    if (!(positions instanceof Float32Array) || !(normals instanceof Float32Array) || groups instanceof Uint16Array) throw new Error('WebGPU texture mesh requires float positions/normals and byte or normalized float groups.');
    const vertices = new Float32Array(count * 20), indices = new Uint32Array(count);
    const groupScale = groups instanceof Uint8Array ? 1 : 255;
    for (let i = 0; i < count; i++) {
        const offset = i * 20, texel = i * 4;
        for (let a = 0; a < 3; a++) {
            if (!Number.isFinite(positions[texel + a]) || !Number.isFinite(normals[texel + a])) throw new Error('WebGPU texture mesh positions and normals must be finite.');
            vertices[offset + a] = positions[texel + a]; vertices[offset + a + 4] = normals[texel + a];
        }
        for (let a = 0; a < 3; a++) if (!Number.isFinite(groups[texel + a]) || groups[texel + a] < 0 || groups[texel + a] * groupScale > 255) throw new Error('WebGPU texture mesh group channels must encode normalized RGB values.');
        const group = unpackRGBToInt(Math.round(groups[texel] * groupScale), Math.round(groups[texel + 1] * groupScale), Math.round(groups[texel + 2] * groupScale));
        if (!Number.isSafeInteger(group) || group < 0) throw new Error('WebGPU texture mesh group IDs must be encoded RGB integers.');
        vertices[offset + 3] = 1; vertices[offset + 16] = group; vertices[offset + 17] = i;
        indices[i] = i;
    }
    return { vertices, indices, vertexCount: count, kind: 0 };
}

function maxIndex(indices: Uint32Array, count: number) {
    if (count > indices.length) throw new Error('WebGPU draw count exceeds the index buffer.');
    let max = 0;
    for (let i = 0; i < count; i++) max = Math.max(max, indices[i]);
    return max;
}

function spheres(v: RenderableValues): WebGPUGeometry {
    const centers = value(v, 'centerBuffer', empty);
    const groups = value(v, 'groupBuffer', empty);
    const count = value(v, 'drawCount', 0) / 6;
    const levels = value<number[][]>(v, 'lodLevels', []);
    const packed = value(v, 'tPositionGroup', { array: empty }).array;
    const descriptors = levels.length ? levels.map(l => {
        if (l.length < 5 || !l.slice(0, 5).every(Number.isFinite) || l[0] > l[1] || l[2] < 0 || l[4] < 0 || (l[4] === 0 && l[3] > 0) || l[3] % 6 !== 0 || l[3] < 0 || l[3] / 6 > count) throw new Error('Native sphere distance levels are invalid.');
        return { count: l[3] / 6, lod: [l[0], l[1], l[2], l[4]] };
    }) : [{ count, lod: [0, 0, 0, 0] }];
    if (levels.length && packed.length < count * 4) throw new Error('Native sphere distance levels require packed positions and groups.');
    // Each level uses exactly the prefix prepared by Spheres.shaderData.
    const total = descriptors.reduce((n, l) => n + l.count, 0), perSphere = 12;
    const vertices = new Float32Array(total * perSphere * 20);
    const faces = [0, 2, 1, 1, 2, 3, 4, 5, 6, 5, 7, 6,
        0, 4, 2, 2, 4, 6, 1, 3, 5, 3, 7, 5,
        0, 1, 4, 1, 5, 4, 2, 6, 3, 3, 6, 7,
        8, 9, 10, 9, 11, 10];
    const indices = new Uint32Array(total * faces.length), sphereLods: NonNullable<WebGPUGeometry['sphereLods']> = [];
    let offset = 0;
    for (const level of descriptors) {
        sphereLods.push({ first: offset * faces.length, count: level.count * faces.length, min: level.lod[0], max: level.lod[1] });
        for (let s = 0; s < level.count; s++) {
            for (let c = 0; c < perSphere; c++) {
                const o = ((offset + s) * perSphere + c) * 20;
                vertices.set(levels.length ? packed.subarray(s * 4, s * 4 + 3) : centers.subarray(s * 3, s * 3 + 3), o); vertices[o + 3] = 1;
                vertices.set(level.lod, o + 4);
                vertices[o + 8] = c % 2 === 0 ? -1 : 1;
                vertices[o + 9] = (c % 4) < 2 ? -1 : 1;
                vertices[o + 10] = c < 4 ? -1 : 1;
                vertices[o + 12] = 1; vertices[o + 15] = c >= 8 ? -1 : 0;
                vertices[o + 16] = levels.length ? packed[s * 4 + 3] : groups[s] || 0; vertices[o + 17] = s;
            }
            for (let i = 0; i < faces.length; i++) indices[(offset + s) * faces.length + i] = (offset + s) * perSphere + faces[i];
        }
        offset += level.count;
    }
    return { vertices, indices, vertexCount: total * perSphere, kind: 5, sphereLods: levels.length ? sphereLods : undefined };
}

function cylinders(v: RenderableValues): WebGPUGeometry {
    const starts = value(v, 'aStart', empty), ends = value(v, 'aEnd', empty);
    const groups = value(v, 'aGroup', empty), scales = value(v, 'aScale', empty);
    const caps = value(v, 'aCap', empty);
    const colorModes = value(v, 'aColorMode', empty);
    // Mol* impostor attributes repeat six times per cylinder.
    const count = Math.floor(value(v, 'drawCount', 0) / 12);
    const vertices: number[] = [], indices: number[] = [];
    for (let c = 0; c < count; c++) {
        const source = c * 6, p = source * 3;
        const start = Array.from(starts.subarray(p, p + 3)), end = Array.from(ends.subarray(p, p + 3));
        const axis = end.map((a, i) => a - start[i]);
        const length = Math.hypot(...axis);
        if (length === 0) continue;
        // Keep small bonds inexpensive while increasing radial detail for the
        // larger cylinders used by surfaces and distance representations.
        // The result is deterministic and remains within a portable index size.
        const radius = Math.max(0.05, Math.abs(scales[source] ?? 1));
        const radialSegments = Math.max(8, Math.min(32, 8 + Math.ceil(radius * 8)));
        for (let i = 0; i < 3; i++) axis[i] /= length;
        const tangent = Math.abs(axis[1]) < 0.9 ? [axis[2], 0, -axis[0]] : [0, -axis[2], axis[1]];
        const tl = Math.hypot(...tangent);
        for (let i = 0; i < 3; i++) tangent[i] /= tl;
        const bitangent = [axis[1] * tangent[2] - axis[2] * tangent[1], axis[2] * tangent[0] - axis[0] * tangent[2], axis[0] * tangent[1] - axis[1] * tangent[0]];
        const base = vertices.length / 20;
        const cap = caps[source] || 0;
        for (let ring = 0; ring < 2; ring++) {
            for (let j = 0; j <= radialSegments; j++) {
                const a = j * 2 * Math.PI / radialSegments;
                const n = tangent.map((t, i) => t * Math.cos(a) + bitangent[i] * Math.sin(a));
                vertices.push(...(ring === 0 ? start : end), 1, ...n, cap, ...n, 0, scales[source] ?? 1, ring, colorModes[source] ?? 2, base + 1, groups[source] || 0, source, 0, 0);
                if (ring === 0 && j < radialSegments) {
                    const i = base + j, e = i + radialSegments + 1;
                    indices.push(i, i + 1, e, i + 1, e + 1, e);
                }
            }
        }
        for (let ring = 0; ring < 2; ring++) {
            const artificial = !(cap & (ring === 0 ? 1 : 2));
            if (artificial && !value(v, 'dSolidInterior', false)) continue;
            const center = vertices.length / 20;
            const pos = ring === 0 ? start : end;
            const normal = axis.map(a => a * (ring === 0 ? -1 : 1));
            vertices.push(...pos, 1, ...normal, artificial ? -1 : 0, 0, 0, 0, 0, scales[source] ?? 1, ring, colorModes[source] ?? 2, base + 1, groups[source] || 0, source, 0, 0);
            for (let j = 0; j <= radialSegments; j++) {
                const a = j * 2 * Math.PI / radialSegments;
                const n = tangent.map((t, i) => t * Math.cos(a) + bitangent[i] * Math.sin(a));
                vertices.push(...pos, 1, ...normal, artificial ? -1 : 0, ...n, 0, scales[source] ?? 1, ring, colorModes[source] ?? 2, base + 1, groups[source] || 0, source, 0, 0);
                if (j < radialSegments) {
                    if (ring === 0) indices.push(center, center + j + 2, center + j + 1);
                    else indices.push(center, center + j + 1, center + j + 2);
                }
            }
        }
        const patch = vertices.length / 20;
        for (let j = 0; j < 4; j++) {
            vertices.push(...start, 1, 0, 0, 0, base + 1, ...end, 0,
                j % 2 === 0 ? -1 : 1, j < 2 ? -1 : 1, scales[source] ?? 1, -1, groups[source] || 0, source, 0, 0);
        }
        indices.push(patch, patch + 1, patch + 2, patch + 1, patch + 3, patch + 2);
    }
    return { vertices: new Float32Array(vertices), indices: new Uint32Array(indices), vertexCount: vertices.length / 20, kind: 6 };
}

function locationIndex(type: string, group: number, instance: number, vertex: number, groups: number, vertices: number) {
    switch (type) {
        case 'uniform': return 0;
        case 'instance': return instance;
        case 'group': return group;
        case 'groupInstance': return instance * groups + group;
        case 'vertex': return vertex;
        case 'vertexInstance': return instance * vertices + vertex;
        default: throw new Error(`WebGPU theme granularity '${type}' has not been implemented.`);
    }
}

/** Per-vertex, per-instance color, size, marker, transparency and emissive values. */
export function createWebGPUThemes(object: GraphicsRenderObject, geometry: WebGPUGeometry): Float32Array {
    const v: RenderableValues = object.values;
    const instances = value(v, 'instanceCount', 1), groups = value(v, 'uGroupCount', 1), sources = value(v, 'uVertexCount', geometry.vertexCount);
    const colors = value(v, 'tColor', { array: new Uint8Array(0) }).array;
    const colorType = value<string>(v, 'dColorType', 'uniform');
    const spatialColor = colorType === 'volume' || colorType === 'volumeInstance';
    const uniformColor = value(v, 'uColor', [1, 1, 1]);
    const sizes = value(v, 'tSize', { array: new Uint8Array(0) }).array;
    const sizeType = value(v, 'dSizeType', 'uniform');
    const markers = value(v, 'tMarker', { array: new Uint8Array(0) }).array;
    const transparency = value(v, 'tTransparency', { array: new Uint8Array(0) }).array;
    const emissive = value(v, 'tEmissive', { array: new Uint8Array(0) }).array;
    const substance = value(v, 'tSubstance', { array: new Uint8Array(0) }).array;
    const material = [value(v, 'uMetalness', 0), value(v, 'uRoughness', 1), value(v, 'uBumpiness', 0)];
    const instanceIds = value(v, 'aInstance', empty);
    const result = new Float32Array(instances * geometry.vertexCount * WebGPUThemeStride);
    const alpha = Math.max(0, Math.min(1, value(v, 'alpha', value(v, 'uAlpha', 1)) * object.state.alphaFactor));
    for (let i = 0; i < instances; i++) {
        const instance = instanceIds[i] ?? i;
        for (let j = 0; j < geometry.vertexCount; j++) {
            const group = geometry.vertices[j * 20 + 16], source = geometry.vertices[j * 20 + 17];
            const o = (i * geometry.vertexCount + j) * WebGPUThemeStride;
            const colorMode = geometry.kind === 6 && geometry.vertices[j * 20 + 15] < 0 ? value(v, 'aColorMode', empty)[source] ?? 2 : geometry.vertices[j * 20 + 14];
            const dualColor = geometry.kind === 6 && value(v, 'dDualColor', false) && (colorType === 'group' || colorType === 'groupInstance') && colorMode !== 2;
            const ci = (spatialColor ? 0 : locationIndex(colorType, group, instance, source, groups, sources)) * (dualColor ? 6 : 3);
            for (let c = 0; c < 3; c++) result[o + c] = colorType === 'uniform' ? uniformColor[c] : (colors[ci + c] ?? 255) / 255;
            if (dualColor && (colorMode <= 1 || colorMode === 3)) {
                // Bonds contain two directed half-cylinders: each runs from its atom to the midpoint.
                const fraction = colorMode === 3 ? geometry.vertices[j * 20 + 13] * 0.5 : colorMode;
                for (let c = 0; c < 3; c++) result[o + c] = result[o + c] * (1 - fraction) + (colors[ci + 3 + c] ?? 255) / 255 * fraction;
            }
            const spatialTransparency = value<string>(v, 'dTransparencyType', 'groupInstance') === 'volumeInstance';
            const ti = spatialTransparency ? 0 : locationIndex(value(v, 'dTransparencyType', 'groupInstance'), group, instance, source, groups, sources);
            result[o + 3] = alpha * (1 - (value(v, 'dTransparency', false) && !spatialTransparency ? (transparency[ti] || 0) / 255 * value(v, 'uTransparencyStrength', 1) : 0));
            const si = locationIndex(sizeType, group, instance, source, groups, sources) * 3;
            result[o + 4] = (sizeType === 'uniform' ? value(v, 'uSize', 1) : ((sizes[si] || 0) * 65536 + (sizes[si + 1] || 0) * 256 + (sizes[si + 2] || 0) - 1) / 100) * value(v, 'uSizeFactor', 1);
            const mi = value<string>(v, 'dMarkerType', 'groupInstance') === 'instance' ? instance : instance * groups + group;
            result[o + 5] = value(v, 'uMarker', -1) === -1 ? markers[mi] || 0 : value(v, 'uMarker', 0);
            const sampledEmissive = value(v, 'dEmissive', false) && value<string>(v, 'dEmissiveType', 'groupInstance') !== 'volumeInstance';
            const ei = sampledEmissive ? locationIndex(value(v, 'dEmissiveType', 'groupInstance'), group, instance, source, groups, sources) : 0;
            result[o + 6] = value(v, 'uEmissive', 0) + (sampledEmissive ? (emissive[ei] || 0) / 255 * value(v, 'uEmissiveStrength', 1) : 0);
            result[o + 7] = instance;
            if (value(v, 'dSubstance', false) && value<string>(v, 'dSubstanceType', 'groupInstance') !== 'volumeInstance') {
                const index = locationIndex(value(v, 'dSubstanceType', 'groupInstance'), group, instance, source, groups, sources) * 4;
                const weight = (substance[index + 3] || 0) / 255;
                const strength = value(v, 'uSubstanceStrength', 1);
                // Pre-mix before interpolation, matching the existing material shader.
                for (let c = 0; c < 3; c++) result[o + 8 + c] = (material[c] * (1 - weight) + (substance[index + c] || 0) / 255 * weight) * strength;
                result[o + 11] = weight * strength;
            }
        }
    }
    return result;
}

/** Raw overpaint, pre-mixed after native color lookup and blended per fragment. */
export function createWebGPUOverpaintThemes(object: GraphicsRenderObject, geometry: WebGPUGeometry) {
    const v = object.values;
    if (!value(v, 'dOverpaint', false) || value<string>(v, 'dOverpaintType', '') === 'volumeInstance') return new Float32Array(4);
    const instances = value(v, 'instanceCount', 1), groups = value(v, 'uGroupCount', 1), sources = value(v, 'uVertexCount', geometry.vertexCount);
    const ids = value(v, 'aInstance', empty), data = value(v, 'tOverpaint', { array: new Uint8Array(0) }).array;
    const result = new Float32Array(instances * geometry.vertexCount * 4);
    for (let i = 0; i < instances; i++) for (let j = 0; j < geometry.vertexCount; j++) {
        const index = locationIndex(value(v, 'dOverpaintType', 'groupInstance'), geometry.vertices[j * 20 + 16], ids[i] ?? i, geometry.vertices[j * 20 + 17], groups, sources) * 4;
        for (let c = 0; c < 4; c++) result[(i * geometry.vertexCount + j) * 4 + c] = (data[index + c] || 0) / 255;
    }
    return result;
}
