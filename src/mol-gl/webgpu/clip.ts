/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { RenderableValues } from '../renderable/schema';
import { value } from './geometry';

export function createClipData(v: RenderableValues) {
    const count = value(v, 'dClipObjectCount', 0);
    const data = new Float32Array(Math.max(1, count) * 32);
    const types = value<number[]>(v, 'uClipObjectType', []), invert = value<boolean[]>(v, 'uClipObjectInvert', []);
    const positions = value<number[]>(v, 'uClipObjectPosition', []), rotations = value<number[]>(v, 'uClipObjectRotation', []);
    const scales = value<number[]>(v, 'uClipObjectScale', []), transforms = value<number[]>(v, 'uClipObjectTransform', []);
    for (let i = 0; i < count; i++) {
        const o = i * 32;
        data[o] = types[i]; data[o + 1] = invert[i] ? 1 : 0;
        data.set(positions.slice(i * 3, i * 3 + 3), o + 4);
        data.set(rotations.slice(i * 4, i * 4 + 4), o + 8);
        data.set(scales.slice(i * 3, i * 3 + 3), o + 12);
        data.set(transforms.slice(i * 16, i * 16 + 16), o + 16);
    }
    return data;
}

export const clipShader = /* wgsl */ `
struct ClipObject { info: vec4f, position: vec4f, rotation: vec4f, scale: vec4f, transform: mat4x4f };
fn quaternionRotate(q: vec4f, v: vec3f) -> vec3f { let t = 2.0 * cross(q.xyz, v); return v + q.w * t + cross(q.xyz, t); }
fn clipInside(point: vec3f, clip: ClipObject) -> bool {
    let transformed = (clip.transform * vec4f(point, 1.0)).xyz;
    let p = quaternionRotate(vec4f(-clip.rotation.xyz, clip.rotation.w), transformed - clip.position.xyz);
    let halfSize = clip.scale.xyz * 0.5;
    let kind = u32(clip.info.x);
    if (kind == 1u) { return dot(quaternionRotate(clip.rotation, vec3f(0.0, 1.0, 0.0)), transformed - clip.position.xyz) <= 0.0; }
    if (kind == 2u) { return length(p / halfSize) <= 1.0; }
    if (kind == 3u) { return all(abs(p) <= halfSize); }
    if (kind == 4u) { return length(p.xz) <= halfSize.x && abs(p.y) <= halfSize.y; }
    if (kind == 5u) { return dot(halfSize.xy, vec2f(length(p.xy), p.z)) <= 0.0; }
    return false;
}
`;
