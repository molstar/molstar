/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './depth.wgsl.js';

/** Derived switches, constants and helpers used by most programs */
export const common_wgsl = /* wgsl */ `
#if defined(dRenderVariant_colorBlended) || defined(dRenderVariant_colorWboit) || defined(dRenderVariant_colorDpoit)
    #define dRenderVariant_color
#endif

#if defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance) || defined(dColorType_vertex) || defined(dColorType_vertexInstance)
    #define dColorType_texture
#endif

#if defined(dColorType_volume) || defined(dColorType_volumeInstance)
    #define dColorType_grid
#endif

#if defined(dColorType_texture) || defined(dColorType_grid)
    #define dColorType_varying
#endif

#if ((defined(dRenderVariant_color) || defined(dRenderVariant_tracing)) && defined(dColorMarker)) || defined(dRenderVariant_marking)
    #define dNeedsMarker
#endif

#if defined(dXrayShaded_on) || defined(dXrayShaded_inverted)
    #define dXrayShaded
#endif

const MaskAll = 0;
const MaskOpaque = 1;
const MaskTransparent = 2;

const MaxDpoitDepth = 99999.0; // NOTE constant also set in TypeScript
const PI = 3.14159265;
const RECIPROCAL_PI = 0.31830988618;
const EPSILON = 1e-6;

fn pow2(x: f32) -> f32 {
    return x * x;
}

fn luminance(c: vec3<f32>) -> f32 {
    // https://www.w3.org/TR/WCAG21/#dfn-relative-luminance
    return dot(c, vec3<f32>(0.2125, 0.7154, 0.0721));
}

${depth_wgsl}

fn packIntToRGB(v: f32) -> vec3<f32> {
    var value = clamp(round(v), 0.0, 16777216.0 - 1.0) + 1.0;
    var c = vec3<f32>(0.0);
    c.b = value % 256.0;
    value = floor(value / 256.0);
    c.g = value % 256.0;
    value = floor(value / 256.0);
    c.r = value % 256.0;
    return c / 255.0;
}

fn packMarkingMask(hidden: f32, highlight: bool, fogAlpha: f32, opacity: f32) -> vec4<f32> {
    let g = hidden + (1.0 - 2.0 * hidden) * (1.0 - opacity) * 0.25;
    return vec4<f32>(0.0, g, select(0.0, 1.0, highlight), fogAlpha);
}

fn unpackRGBToInt(rgb: vec3<f32>) -> f32 {
    return (floor(rgb.r * 255.0 + 0.5) * 256.0 * 256.0 + floor(rgb.g * 255.0 + 0.5) * 256.0 + floor(rgb.b * 255.0 + 0.5)) - 1.0;
}

/** Depth in [0, 1] of a view space position (the projection maps to [0, 1]) */
fn calcDepth(pos: vec3<f32>) -> f32 {
    let clipZW = pos.z * uGlobal.uProjection[2].zw + uGlobal.uProjection[3].zw;
    return clipZW.x / clipZW.y;
}

fn mat3From4(m: mat4x4<f32>) -> mat3x3<f32> {
    return mat3x3<f32>(m[0].xyz, m[1].xyz, m[2].xyz);
}

// see https://github.com/graphitemaster/normals_revisited and https://www.shadertoy.com/view/3s33zj
fn adjoint(m: mat4x4<f32>) -> mat3x3<f32> {
    return mat3x3<f32>(
        cross(m[1].xyz, m[2].xyz),
        cross(m[2].xyz, m[0].xyz),
        cross(m[0].xyz, m[1].xyz)
    );
}

fn hash(h: f32) -> f32 {
    return fract(sin(h) * 43758.5453123);
}

fn noise(x: vec3<f32>) -> f32 {
    let p = floor(x);
    var f = fract(x);
    f = f * f * (3.0 - 2.0 * f);

    let n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(
        mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
            mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
        mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
            mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}

fn fbm(x: vec3<f32>) -> f32 {
    var p = x;
    var f = 0.0;
    f += 0.5 * noise(p);
    p *= 2.01;
    f += 0.25 * noise(p);
    p *= 2.02;
    f += 0.125 * noise(p);
    return f;
}
`;
