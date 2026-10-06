/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { RenderableValues } from '../renderable/schema';
import { value } from './geometry';

/** Per-cell overlays; image colors themselves are already supplied by the representation. */
export function createImageStyles(v: RenderableValues) {
    const groups = value(v, 'uGroupCount', 1), instances = value(v, 'instanceCount', 1);
    const result = new Float32Array(Math.max(1, groups * instances) * 8);
    const marker = value(v, 'tMarker', { array: new Uint8Array(0) }).array;
    const transparency = value(v, 'tTransparency', { array: new Uint8Array(0) }).array;
    const emissive = value(v, 'tEmissive', { array: new Uint8Array(0) }).array;
    const overpaint = value(v, 'tOverpaint', { array: new Uint8Array(0) }).array;
    for (let i = 0; i < instances; i++) for (let g = 0; g < groups; g++) {
        const o = (i * groups + g) * 8;
        const index = (key: string) => value<string>(v, key, 'groupInstance') === 'instance' ? i : i * groups + g;
        const oi = index('dOverpaintType') * 4;
        if (value(v, 'dOverpaint', false)) {
            for (let c = 0; c < 3; c++) result[o + c] = (overpaint[oi + c] || 0) / 255;
            result[o + 3] = (overpaint[oi + 3] || 0) / 255 * value(v, 'uOverpaintStrength', 1);
        }
        result[o + 4] = value(v, 'uMarker', -1) === -1 ? marker[index('dMarkerType')] || 0 : value(v, 'uMarker', 0);
        result[o + 6] = value(v, 'uEmissive', 0) + (value(v, 'dEmissive', false) ? (emissive[index('dEmissiveType')] || 0) / 255 * value(v, 'uEmissiveStrength', 1) : 0);
        result[o + 5] = value(v, 'dTransparency', false) ? (transparency[index('dTransparencyType')] || 0) / 255 * value(v, 'uTransparencyStrength', 1) : 0;
    }
    return result;
}

export const imageShader = /* wgsl */ `
@group(1) @binding(6) var imageGroups: texture_2d<f32>;
@group(1) @binding(7) var imageValues: texture_2d<f32>;
@group(1) @binding(8) var imagePalette: texture_2d<f32>;
fn paletteColor(encoded: f32) -> vec3f {
    if (object.colorPalette.y > 0.0) { return textureSampleLevel(imagePalette, atlasSampler, vec2f(encoded, 0.5), 0.0).rgb; }
    let size = i32(textureDimensions(imagePalette).x);
    return textureLoad(imagePalette, vec2i(clamp(i32(floor(encoded * f32(size))), 0, size - 1), 0), 0).rgb;
}
struct ImageStyle { color: vec4f, style: vec4f };
@group(1) @binding(9) var<storage, read> imageStyles: array<ImageStyle>;
fn imageCoord(uv: vec2f) -> vec2i {
    let size = vec2i(textureDimensions(atlas));
    return clamp(vec2i(floor(uv * vec2f(size))), vec2i(0), size - 1);
}
fn cubic(x: f32, b: f32, c: f32) -> f32 {
    let a = abs(x);
    if (a < 1.0) { return ((12.0 - 9.0*b - 6.0*c)*a*a*a + (-18.0 + 12.0*b + 6.0*c)*a*a + 6.0 - 2.0*b) / 6.0; }
    if (a < 2.0) { return ((-b - 6.0*c)*a*a*a + (6.0*b + 30.0*c)*a*a + (-12.0*b - 48.0*c)*a + 8.0*b + 24.0*c) / 6.0; }
    return 0.0;
}
fn imageSample(uv: vec2f) -> vec4f {
    let nearest = textureLoad(atlas, imageCoord(uv), 0);
    if (object.image.y == 0.0 || (object.image.z > 0.0 && all(nearest.rgb == vec3f(1.0)))) { return nearest; }
    let size = vec2i(textureDimensions(atlas));
    let p = uv * vec2f(size) - 0.5;
    let cell = fract(p); let base = vec2i(floor(p));
    var b = 1.0; var c = 0.0;
    if (object.image.y == 1.0) { b = 0.0; c = 0.5; }
    if (object.image.y == 2.0) { b = 1.0/3.0; c = 1.0/3.0; }
    var color = vec4f(0.0); var weight = 0.0;
    for (var y = -1; y <= 2; y++) { for (var x = -1; x <= 2; x++) {
        let w = abs(cubic(f32(x) - cell.x, b, c) * cubic(f32(y) - cell.y, b, c));
        color += textureLoad(atlas, clamp(base + vec2i(x, y), vec2i(0), size - 1), 0) * w;
        weight += w;
    } }
    return color / max(weight, 0.00001);
}
fn imageOutside(position: vec3f) -> bool {
    let kind = u32(object.trimCenter.w);
    if (kind == 0u) { return false; }
    let center = (object.trimTransform * vec4f(position, 1.0)).xyz;
    let q = vec4f(-object.trimRotation.xyz, object.trimRotation.w);
    let p = quaternionRotate(q, center - object.trimCenter.xyz);
    let halfSize = object.trimScale.xyz * 0.5;
    if (kind == 1u) { return dot(quaternionRotate(object.trimRotation, vec3f(0.0, 1.0, 0.0)), center - object.trimCenter.xyz) > 0.0; }
    if (kind == 2u) { return length(p / halfSize) > 1.0; }
    if (kind == 3u) { return any(abs(p) > halfSize); }
    if (kind == 4u) { return length(p.xz) > halfSize.x || abs(p.y) > halfSize.y; }
    return dot(halfSize.xy, vec2f(length(p.xy), p.z)) > 0.0;
}
`;
