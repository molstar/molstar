/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

const DepthTextures = [
  'tDepth',
  'tDepthHalf',
  'tDepthQuarter',
  'tDepthTransparent',
  'tDepthHalfTransparent',
  'tDepthQuarterTransparent',
];

/**
 * Screen space ambient occlusion (StarCraft II AO, Filion and McNaughton 2008) into packed occlusion (rg) and depth
 * (ba). The math is in GL window space with the GL projection of the pass uniforms, as in GLSL.
 */
export const ssao_wgsl = /* wgsl */ `
${packing_wgsl}
${quad_vert_wgsl}
${DepthTextures.map(texture_read_wgsl).join('\n')}

override oNSamples: i32 = 32;
override oLevels: i32 = 3;

const PI = 3.14159265;
const QuarterThreshold = 0.1;
const HalfThreshold = 0.05;

fn smootherstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

fn noise(coords: vec2<f32>) -> f32 {
    let dt = dot(coords, vec2<f32>(12.9898, 78.233));
    let sn = dt - PI * floor(dt / PI); // GLSL mod
    return abs(fract(sin(sn) * 43758.5453));
}

fn getNoiseVec2(coords: vec2<f32>) -> vec2<f32> {
    return vec2<f32>(noise(coords), noise(coords + vec2<f32>(PI, 2.71828)));
}

fn isBackground(depth: f32) -> bool {
    return depth == 1.0;
}

// inset by half a texel to avoid sampling texels outside the viewport
fn clampToBounds(coords: vec2<f32>) -> vec2<f32> {
    let halfTexel = 0.5 / uObject.uTexSize;
    return clamp(coords, uObject.uBounds.xy + halfTexel, uObject.uBounds.zw - halfTexel);
}

fn getDepth(coords: vec2<f32>, transparentFlag: i32) -> f32 {
    if (transparentFlag == 1) {
        #if defined(dIncludeTransparent)
            return unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords)).x;
        #else
            return 1.0;
        #endif
    }
    #if defined(dPackedDepth)
        return unpackRGBAToDepth(read_tDepth(coords));
    #else
        return read_tDepth(coords).r;
    #endif
}

fn getLevelDepth(v: vec4<f32>) -> f32 {
    #if defined(dPackedDepthLevels)
        return unpackRGBAToDepth(v);
    #else
        return v.r;
    #endif
}

fn getMappedDepth(coords: vec2<f32>, selfCoords: vec2<f32>) -> f32 {
    let d = distance(coords, selfCoords);
    if (d > QuarterThreshold) {
        return getLevelDepth(read_tDepthQuarter(coords));
    } else if (d > HalfThreshold) {
        return getLevelDepth(read_tDepthHalf(coords));
    }
    return getDepth(coords, 0);
}

#if defined(dIncludeTransparent)
    fn getDepthTransparentWithAlpha(coords: vec2<f32>) -> vec2<f32> {
        return unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords));
    }

    fn getMappedDepthTransparentWithAlpha(coords: vec2<f32>, selfCoords: vec2<f32>) -> vec2<f32> {
        let d = distance(coords, selfCoords);
        if (d > QuarterThreshold) {
            return unpackRGBAToDepthWithAlpha(read_tDepthQuarterTransparent(coords));
        } else if (d > HalfThreshold) {
            return unpackRGBAToDepthWithAlpha(read_tDepthHalfTransparent(coords));
        }
        return unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords));
    }
#endif

fn viewPos(coords: vec2<f32>, depth: f32) -> vec3<f32> {
    return screenSpaceToViewSpace(vec3<f32>(coords, depth), uObject.uInvProjection);
}

// adapted from https://gist.github.com/bgolus/a07ed65602c009d5e2f753826e8078a0
fn viewNormalAtPixelPositionAccurate(vpos: vec2<f32>, transparentFlag: i32) -> vec3<f32> {
    let texel = 1.0 / uObject.uTexSize;
    let dx = vec2<f32>(texel.x, 0.0);
    let dy = vec2<f32>(0.0, texel.y);

    // current pixel's depth and view space position
    let c = getDepth(vpos, transparentFlag);
    let viewSpacePos_c = viewPos(vpos, c);

    // view space positions at 1 pixel offsets in each major direction
    let viewSpacePos_l = viewPos(vpos - dx, getDepth(vpos - dx, transparentFlag));
    let viewSpacePos_r = viewPos(vpos + dx, getDepth(vpos + dx, transparentFlag));
    let viewSpacePos_d = viewPos(vpos - dy, getDepth(vpos - dy, transparentFlag));
    let viewSpacePos_u = viewPos(vpos + dy, getDepth(vpos + dy, transparentFlag));

    let l = viewSpacePos_c - viewSpacePos_l;
    let r = viewSpacePos_r - viewSpacePos_c;
    let d = viewSpacePos_c - viewSpacePos_d;
    let u = viewSpacePos_u - viewSpacePos_c;

    // depth values at 1 & 2 pixels offsets from current along the horizontal and vertical axes
    let H = vec4<f32>(
        getDepth(vpos - dx, transparentFlag),
        getDepth(vpos + dx, transparentFlag),
        getDepth(vpos - 2.0 * dx, transparentFlag),
        getDepth(vpos + 2.0 * dx, transparentFlag)
    );
    let V = vec4<f32>(
        getDepth(vpos - dy, transparentFlag),
        getDepth(vpos + dy, transparentFlag),
        getDepth(vpos - 2.0 * dy, transparentFlag),
        getDepth(vpos + 2.0 * dy, transparentFlag)
    );

    // current pixel's depth difference from slope of offset depth samples (non-linear depth values)
    let he = abs((2.0 * H.xy - H.zw) - c);
    let ve = abs((2.0 * V.xy - V.zw) - c);

    // horizontal and vertical diff with the smallest depth difference from slopes
    let hDeriv = select(r, l, he.x < he.y);
    let vDeriv = select(u, d, ve.x < ve.y);

    return normalize(cross(hDeriv, vDeriv));
}

fn getPixelSize(coords: vec2<f32>, depth: f32) -> f32 {
    let viewPos0 = viewPos(coords, depth);
    let viewPos1 = viewPos(coords + vec2<f32>(1.0, 0.0) / uObject.uTexSize, depth);
    return distance(viewPos0, viewPos1);
}

/** GL window space of a view space position */
fn projectToWindow(p: vec3<f32>) -> vec2<f32> {
    let offset = uObject.uProjection * vec4<f32>(p, 1.0);
    return clampToBounds((offset.xy / offset.w) * 0.5 + 0.5);
}

/** occlusion of a sample, from opaque and (with dIncludeTransparent) transparent depth */
fn sampleOcclusion(sampleViewPos: vec3<f32>, coords: vec2<f32>, selfCoords: vec2<f32>, selfViewZ: f32, radius: f32, bias: f32, mapped: bool) -> f32 {
    var sampleOcc = 0.0;
    #if defined(dIllumination)
        if (uObject.uTransparencyFlag == 1) {
    #endif
            // NOTE: using the mapped depth in the single-scale case causes issues on some mobile devices
            let sampleDepth = select(getDepth(coords, 0), getMappedDepth(coords, selfCoords), mapped);
            if (!isBackground(sampleDepth)) {
                let sampleViewZ = viewPos(coords, sampleDepth).z;
                sampleOcc = step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, radius / abs(selfViewZ - sampleViewZ)) * bias;
            }
    #if defined(dIllumination)
        }
    #endif
    #if defined(dIncludeTransparent)
        let sampleDepthWithAlpha = select(getDepthTransparentWithAlpha(coords), getMappedDepthTransparentWithAlpha(coords, selfCoords), mapped);
        if (!isBackground(sampleDepthWithAlpha.x)) {
            let sampleViewZWithAlpha = viewPos(coords, sampleDepthWithAlpha.x).z;
            sampleOcc = max(sampleOcc, step(sampleViewPos.z + 0.025, sampleViewZWithAlpha) * smootherstep(0.0, 1.0, radius / abs(selfViewZ - sampleViewZWithAlpha)) * bias * sampleDepthWithAlpha.y);
        }
    #endif
    return sampleOcc;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let selfCoords = glCoords(in.position, uObject.uTexSize);
    let selfDepth = getDepth(selfCoords, uObject.uTransparencyFlag);
    let selfPackedDepth = packUnitIntervalToRG(selfDepth);

    if (isBackground(selfDepth)) {
        return vec4<f32>(packUnitIntervalToRG(1.0), selfPackedDepth);
    }

    let selfViewNormal = viewNormalAtPixelPositionAccurate(selfCoords, uObject.uTransparencyFlag);
    let selfViewPos = viewPos(selfCoords, selfDepth);

    let randomVec = normalize(vec3<f32>(getNoiseVec2(selfCoords) * 2.0 - 1.0, 0.0));
    let tangent = normalize(randomVec - selfViewNormal * dot(randomVec, selfViewNormal));
    let bitangent = cross(selfViewNormal, tangent);
    let TBN = mat3x3<f32>(tangent, bitangent, selfViewNormal);

    var occlusion = 0.0;
    #if defined(dMultiScale)
        let pixelSize = getPixelSize(selfCoords, selfDepth);
        for (var l = 0; l < oLevels; l++) {
            let radius = uObject.uLevelRadius[l].x;
            // TODO: smooth transition
            if (pixelSize * uObject.uNearThreshold > radius) {
                continue;
            }
            if (pixelSize * uObject.uFarThreshold < radius) {
                continue;
            }
            var levelOcclusion = 0.0;
            for (var i = 0; i < oNSamples; i++) {
                let sampleViewPos = selfViewPos + (TBN * uObject.uSamples[i]) * radius;
                let coords = projectToWindow(sampleViewPos);
                levelOcclusion += sampleOcclusion(sampleViewPos, coords, selfCoords, selfViewPos.z, radius, uObject.uLevelBias[l].x, true);
            }
            levelOcclusion /= f32(oNSamples);
            occlusion = max(occlusion, levelOcclusion);
        }
    #else
        for (var i = 0; i < oNSamples; i++) {
            let sampleViewPos = selfViewPos + (TBN * uObject.uSamples[i]) * uObject.uRadius;
            let coords = projectToWindow(sampleViewPos);
            occlusion += sampleOcclusion(sampleViewPos, coords, selfCoords, selfViewPos.z, uObject.uRadius, 1.0, false);
        }
        occlusion /= f32(oNSamples);
    #endif
    occlusion = 1.0 - (uObject.uBias * occlusion);

    return vec4<f32>(packUnitIntervalToRG(clamp(occlusion, 0.01, 1.0)), selfPackedDepth);
}
`;
