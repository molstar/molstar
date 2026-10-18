/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Slightly adapted from https://github.com/mrdoob/three.js
 * MIT License Copyright (c) 2010-2020 three.js authors
 *
 * WebGL port of Subpixel Morphological Antialiasing (SMAA) v2.8
 * Preset: SMAA 1x Medium (with color edge detection)
 * https://github.com/iryoku/smaa/releases/tag/v2.8
 */

import { quad_vert_wgsl, texture_sample_wgsl } from './chunks/quad.wgsl.js';

// The quads cover the whole target, so the GLSL vertex stage uv (`vUv`) is `gl_FragCoord.xy * uTexSizeInv`. SMAA's
// searches and area lookups are oriented, so its programs keep GL's window space (bottom-left origin) of the GLSL port:
// uvs from `glCoords`, `read_<name>` flips at the reads of rendered textures, the area and search textures are read in
// upload order.

/** uv of the fragment with a bottom-left origin (`gl_FragCoord.xy / size` of WebGL) */
const gl_coords_wgsl = /* wgsl */ `
fn glCoords(position: vec4<f32>, size: vec2<f32>) -> vec2<f32> {
    return vec2<f32>(position.x / size.x, 1.0 - position.y / size.y);
}
`;

/** `read_<name>(uv)`: texel of a rendered texture at a uv with a bottom-left origin, the nearest or filtered */
function texture_read_wgsl(name: string) {
  return /* wgsl */ `
fn read_${name}(uv: vec2<f32>) -> vec4<f32> {
    let st = vec2<f32>(uv.x, 1.0 - uv.y);
    #if defined(dLinear_${name})
        return textureSampleLevel(${name}, sLinear, st, 0.0);
    #else
        let size = vec2<i32>(textureDimensions(${name}));
        return textureLoad(${name}, clamp(vec2<i32>(floor(st * vec2<f32>(size))), vec2<i32>(0), size - 1), 0);
    #endif
}
`;
}

/** Color edge detection */
export const smaa_edges_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${gl_coords_wgsl}
${texture_read_wgsl('tColor')}

override oEdgeThreshold: f32 = 0.1;

fn colorDelta(a: vec3<f32>, b: vec3<f32>) -> f32 {
    let t = abs(a - b);
    return max(max(t.r, t.g), t.b);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let texSizeInv = uObject.uTexSizeInv;
    let texCoord = glCoords(in.position, 1.0 / texSizeInv);
    let offset0 = texCoord.xyxy + texSizeInv.xyxy * vec4<f32>(-1.0, 0.0, 0.0, 1.0);
    let offset1 = texCoord.xyxy + texSizeInv.xyxy * vec4<f32>(1.0, 0.0, 0.0, -1.0);
    let offset2 = texCoord.xyxy + texSizeInv.xyxy * vec4<f32>(-2.0, 0.0, 0.0, 2.0);

    let threshold = vec2<f32>(oEdgeThreshold, oEdgeThreshold);

    // Calculate color deltas:
    var delta: vec4<f32>;
    let C = read_tColor(texCoord).rgb;
    delta.x = colorDelta(C, read_tColor(offset0.xy).rgb);
    delta.y = colorDelta(C, read_tColor(offset0.zw).rgb);

    // We do the usual threshold:
    var edges = step(threshold, delta.xy);

    // Then discard if there is no edge:
    if (dot(edges, vec2<f32>(1.0, 1.0)) == 0.0) {
        discard;
    }

    // Calculate right and bottom deltas:
    delta.z = colorDelta(C, read_tColor(offset1.xy).rgb);
    delta.w = colorDelta(C, read_tColor(offset1.zw).rgb);

    // Calculate the maximum delta in the direct neighborhood:
    var maxDelta = max(max(max(delta.x, delta.y), delta.z), delta.w);

    // Calculate left-left and top-top deltas:
    delta.z = colorDelta(C, read_tColor(offset2.xy).rgb);
    delta.w = colorDelta(C, read_tColor(offset2.zw).rgb);

    // Calculate the final maximum delta:
    maxDelta = max(max(maxDelta, delta.z), delta.w);

    // Local contrast adaptation in action:
    edges *= step(vec2<f32>(0.5 * maxDelta), delta.xy);

    return vec4<f32>(edges, 0.0, 0.0);
}
`;

/** Blending weight calculation */
export const smaa_weights_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${gl_coords_wgsl}
${texture_read_wgsl('tEdges')}
${texture_sample_wgsl('tArea')}
${texture_sample_wgsl('tSearch')}

override oMaxSearchSteps: i32 = 16;

const AreaTexMaxDistance = 16.0;
const AreaTexPixelSize = 1.0 / vec2<f32>(160.0, 560.0);
const AreaTexSubtexSize = 1.0 / 7.0;

fn searchLength(e: vec2<f32>, bias: f32, scale: f32) -> f32 {
    return 255.0 * sample_tSearch(vec2<f32>(bias + e.x * scale, e.y)).r;
}

fn searchXLeft(start: vec2<f32>, end: f32) -> f32 {
    let texSizeInv = uObject.uTexSizeInv;
    var texCoord = start;
    // @PSEUDO_GATHER4: the coordinate is offset by (-0.25, -0.125) to sample between edges, fetching four in a row
    var e = vec2<f32>(0.0, 1.0);
    for (var i = 0; i < oMaxSearchSteps; i++) {
        e = read_tEdges(texCoord).rg;
        texCoord -= vec2<f32>(2.0, 0.0) * texSizeInv;
        if (!(texCoord.x > end && e.g > 0.8281 && e.r == 0.0)) {
            break;
        }
    }
    // correct the offset, the bias by one, and the last step
    texCoord.x += 0.25 * texSizeInv.x;
    texCoord.x += texSizeInv.x;
    texCoord.x += 2.0 * texSizeInv.x;
    texCoord.x -= texSizeInv.x * searchLength(e, 0.0, 0.5);
    return texCoord.x;
}

fn searchXRight(start: vec2<f32>, end: f32) -> f32 {
    let texSizeInv = uObject.uTexSizeInv;
    var texCoord = start;
    var e = vec2<f32>(0.0, 1.0);
    for (var i = 0; i < oMaxSearchSteps; i++) {
        e = read_tEdges(texCoord).rg;
        texCoord += vec2<f32>(2.0, 0.0) * texSizeInv;
        if (!(texCoord.x < end && e.g > 0.8281 && e.r == 0.0)) {
            break;
        }
    }
    texCoord.x -= 0.25 * texSizeInv.x;
    texCoord.x -= texSizeInv.x;
    texCoord.x -= 2.0 * texSizeInv.x;
    texCoord.x += texSizeInv.x * searchLength(e, 0.5, 0.5);
    return texCoord.x;
}

fn searchYUp(start: vec2<f32>, end: f32) -> f32 {
    let texSizeInv = uObject.uTexSizeInv;
    var texCoord = start;
    var e = vec2<f32>(1.0, 0.0);
    for (var i = 0; i < oMaxSearchSteps; i++) {
        e = read_tEdges(texCoord).rg;
        texCoord += vec2<f32>(0.0, 2.0) * texSizeInv;
        if (!(texCoord.y > end && e.r > 0.8281 && e.g == 0.0)) {
            break;
        }
    }
    texCoord.y -= 0.25 * texSizeInv.y;
    texCoord.y -= texSizeInv.y;
    texCoord.y -= 2.0 * texSizeInv.y;
    texCoord.y += texSizeInv.y * searchLength(e.gr, 0.0, 0.5);
    return texCoord.y;
}

fn searchYDown(start: vec2<f32>, end: f32) -> f32 {
    let texSizeInv = uObject.uTexSizeInv;
    var texCoord = start;
    var e = vec2<f32>(1.0, 0.0);
    for (var i = 0; i < oMaxSearchSteps; i++) {
        e = read_tEdges(texCoord).rg;
        texCoord -= vec2<f32>(0.0, 2.0) * texSizeInv;
        if (!(texCoord.y < end && e.r > 0.8281 && e.g == 0.0)) {
            break;
        }
    }
    texCoord.y += 0.25 * texSizeInv.y;
    texCoord.y += texSizeInv.y;
    texCoord.y += 2.0 * texSizeInv.y;
    texCoord.y -= texSizeInv.y * searchLength(e.gr, 0.5, 0.5);
    return texCoord.y;
}

fn area(dist: vec2<f32>, e1: f32, e2: f32, offset: f32) -> vec2<f32> {
    // Rounding prevents precision errors of bilinear filtering:
    var texCoord = AreaTexMaxDistance * round(4.0 * vec2<f32>(e1, e2)) + dist;
    // We do a scale and bias for mapping to texel space:
    texCoord = AreaTexPixelSize * texCoord + (0.5 * AreaTexPixelSize);
    // Move to proper place, according to the subpixel offset:
    texCoord.y += AreaTexSubtexSize * offset;
    return sample_tArea(texCoord).rg;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let texSizeInv = uObject.uTexSizeInv;
    let texCoord = glCoords(in.position, 1.0 / texSizeInv);
    let pixCoord = texCoord / texSizeInv;

    // offsets for the searches (see @PSEUDO_GATHER4) and the ends of the loops
    let offset0 = texCoord.xyxy + texSizeInv.xyxy * vec4<f32>(-0.25, 0.125, 1.25, 0.125);
    let offset1 = texCoord.xyxy + texSizeInv.xyxy * vec4<f32>(-0.125, 0.25, -0.125, -1.25);
    let offset2 = vec4<f32>(offset0.xz, offset1.yw) + vec4<f32>(-2.0, 2.0, -2.0, 2.0) * texSizeInv.xxyy * f32(oMaxSearchSteps);

    var weights = vec4<f32>(0.0);
    let e = read_tEdges(texCoord).rg;

    if (e.g > 0.0) { // Edge at north
        // Find the distance to the left:
        var coords = vec2<f32>(searchXLeft(offset0.xy, offset2.x), offset1.y);
        var d = vec2<f32>(coords.x, 0.0);

        // Fetch the left crossing edges, two at a time using bilinear filtering:
        let e1 = read_tEdges(coords).r;

        // Find the distance to the right:
        coords.x = searchXRight(offset0.zw, offset2.y);
        d.y = coords.x;

        // We want the distances to be in pixel units:
        d = d / texSizeInv.x - pixCoord.x;

        // The areas texture is compressed quadratically:
        let sqrtD = sqrt(abs(d));

        // Fetch the right crossing edges:
        coords.y -= texSizeInv.y;
        // the GLSL offset macro converts ivec2(1, 0) with float(), taking the first component: (1, 1) texels
        let e2 = read_tEdges(coords + vec2<f32>(1.0, 1.0) * texSizeInv).r;

        let a = area(sqrtD, e1, e2, 0.0);
        weights = vec4<f32>(a, weights.ba);
    }

    if (e.r > 0.0) { // Edge at west
        // Find the distance to the top:
        var coords = vec2<f32>(offset0.x, searchYUp(offset1.xy, offset2.z));
        var d = vec2<f32>(coords.y, 0.0);

        // Fetch the top crossing edges:
        let e1 = read_tEdges(coords).g;

        // Find the distance to the bottom:
        coords.y = searchYDown(offset1.zw, offset2.w);
        d.y = coords.y;

        // We want the distances to be in pixel units:
        d = d / texSizeInv.y - pixCoord.y;

        let sqrtD = sqrt(abs(d));

        // Fetch the bottom crossing edges:
        coords.y -= texSizeInv.y;
        // the GLSL offset macro converts ivec2(0, 1) with float(), taking the first component: no offset
        let e2 = read_tEdges(coords).g;

        let a = area(sqrtD, e1, e2, 0.0);
        weights = vec4<f32>(weights.rg, a);
    }

    return weights;
}
`;

/** Neighborhood blending */
export const smaa_blend_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${gl_coords_wgsl}
${texture_read_wgsl('tColor')}
${texture_read_wgsl('tWeights')}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let texSizeInv = uObject.uTexSizeInv;
    var texCoord = glCoords(in.position, 1.0 / texSizeInv);
    let offset1 = texCoord.xyxy + texSizeInv.xyxy * vec4<f32>(1.0, 0.0, 0.0, -1.0);

    // Fetch the blending weights for current pixel:
    let w = read_tWeights(texCoord);
    let a = vec4<f32>(w.x, read_tWeights(offset1.zw).g, w.z, read_tWeights(offset1.xy).a);

    // Is there any blending weight with a value greater than 0.0?
    if (dot(a, vec4<f32>(1.0)) < 1e-5) {
        return read_tColor(texCoord);
    }

    // Up to 4 lines can be crossing a pixel (one through each edge). We favor blending by choosing the line with the
    // maximum weight for each direction:
    var offset = vec2<f32>(
        select(-a.b, a.a, a.a > a.b), // left vs. right
        select(a.r, -a.g, a.g > a.r) // top vs. bottom
    );

    // Then we go in the direction that has the maximum weight:
    if (abs(offset.x) > abs(offset.y)) { // horizontal vs. vertical
        offset.y = 0.0;
    } else {
        offset.x = 0.0;
    }

    // Fetch the opposite color and lerp by hand:
    var C = read_tColor(texCoord);
    texCoord += sign(offset) * texSizeInv;
    var Cop = read_tColor(texCoord);
    let s = select(abs(offset.y), abs(offset.x), abs(offset.x) > abs(offset.y));

    // gamma correction
    C = vec4<f32>(pow(C.xyz, vec3<f32>(2.2)), C.w);
    Cop = vec4<f32>(pow(Cop.xyz, vec3<f32>(2.2)), Cop.w);
    let mixed = mix(C, Cop, s);
    return vec4<f32>(pow(mixed.xyz, vec3<f32>(1.0 / 2.2)), mixed.w);
}
`;
