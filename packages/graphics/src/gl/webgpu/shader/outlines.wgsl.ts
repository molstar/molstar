/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/** 2 values in [0, 1] with 4 bits each, in one 8-bit channel */
export const pack_2x4_wgsl = /* wgsl */ `
fn pack2x4(v: vec2<f32>) -> f32 {
    let scaled = floor(clamp(v, vec2<f32>(0.0), vec2<f32>(1.0)) * 15.0 + 0.5); // round to 0-15
    return (scaled.x + scaled.y * 16.0) / 255.0;
}

fn unpack2x4(f: f32) -> vec2<f32> {
    let c = floor(f * 255.0 + 0.5);
    let lo = c - 16.0 * floor(c / 16.0);
    let hi = floor(c / 16.0);
    return vec2<f32>(lo, hi) / 15.0;
}
`;

/**
 * `getOutline(coords)`: the nearest opaque and transparent outlines within `oOutlineScale` texels, read from
 * `tOutlines` (`read_tOutlines`, with `pack_2x4_wgsl`); for the compose passes with outlines on
 */
export const outline_read_wgsl = /* wgsl */ `
override oOutlineScale: i32 = 1;

struct Outline {
    hasOpaque: bool,
    hasTransparent: bool,
    opaqueDepth: f32,
    transparentDepth: f32,
    alpha: f32,
};

fn getOutline(coords: vec2<f32>) -> Outline {
    let invTexSize = 1.0 / uObject.uTexSize;
    let squaredOutlineScale = oOutlineScale * oOutlineScale;

    var o = Outline(false, false, 0.0, 0.0, 0.0);
    for (var y = -oOutlineScale; y <= oOutlineScale; y++) {
        for (var x = -oOutlineScale; x <= oOutlineScale; x++) {
            if (x * x + y * y > squaredOutlineScale) {
                continue;
            }

            let sampleCoords = coords + vec2<f32>(f32(x), f32(y)) * invTexSize;

            let sampleOutlineCombined = read_tOutlines(sampleCoords);
            let sampleOpaqueDepth = unpackRGToUnitInterval(sampleOutlineCombined.gb);
            let sampleTransparentDepth = sampleOutlineCombined.a;
            let sampleFlagWithAlpha = unpack2x4(sampleOutlineCombined.r);

            let sampleFlag = sampleFlagWithAlpha.x;
            let sampleAlpha = clamp(sampleFlagWithAlpha.y * 0.5, 0.01, 1.0);

            let isOpaque = sampleFlag > 0.20 && sampleFlag < 0.30;
            let isTransparent = sampleFlag > 0.45 && sampleFlag < 0.55;
            let isBoth = sampleFlag > 0.70 && sampleFlag < 0.80;

            if ((isOpaque || isBoth) && sampleOpaqueDepth > o.opaqueDepth) {
                o.hasOpaque = true;
                o.opaqueDepth = sampleOpaqueDepth;
            }

            if ((isTransparent || isBoth) && sampleTransparentDepth > o.transparentDepth) {
                o.hasTransparent = true;
                o.transparentDepth = sampleTransparentDepth;
                o.alpha = sampleAlpha;
            }
        }
    }
    return o;
}
`;

/**
 * Outline edges of the opaque and transparent depth: packed outline type and transparent alpha
 * (r), the packed opaque or transparent edge depth (gb) and the transparent edge depth (a). Depth textures are always
 * available on WebGPU, so the opaque depth is never packed.
 */
export const outlines_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${pack_2x4_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tDepthOpaque')}
${texture_read_wgsl('tDepthTransparent')}

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

fn getDepthOpaque(coords: vec2<f32>) -> f32 {
    return read_tDepthOpaque(coords).r;
}

fn getDepthTransparentWithAlpha(coords: vec2<f32>) -> vec2<f32> {
    #if defined(dTransparentOutline)
        return unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords));
    #else
        return vec2<f32>(0.0, 0.0);
    #endif
}

// depth is reversed: the background (far plane) at 0, nearer is greater
fn isBackground(depth: f32) -> bool {
    return depth == 0.0;
}

fn getViewZOrBackground(depth: f32) -> f32 {
    return select(getViewZ(depth), 2.0 * uObject.uFar, isBackground(depth));
}

fn getPixelSize(coords: vec2<f32>, depth: f32) -> f32 {
    let viewPos0 = screenSpaceToViewSpace(vec3<f32>(coords, depth), uObject.uInvProjection);
    let viewPos1 = screenSpaceToViewSpace(vec3<f32>(coords + vec2<f32>(1.0, 0.0) / uObject.uTexSize, depth), uObject.uInvProjection);
    return distance(viewPos0, viewPos1);
}

/** Largest second difference of the view z around \`coords\` */
fn getCurvature(left: f32, right: f32, up: f32, down: f32, selfViewZ: f32) -> f32 {
    let ddx = abs(getViewZOrBackground(left) + getViewZOrBackground(right) - 2.0 * selfViewZ);
    let ddy = abs(getViewZOrBackground(up) + getViewZOrBackground(down) - 2.0 * selfViewZ);
    return max(ddx, ddy);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, uObject.uTexSize);
    let invTexSize = 1.0 / uObject.uTexSize;
    let threshold = uObject.uOutlineThreshold;

    let selfDepthOpaque = getDepthOpaque(coords);
    let selfViewZOpaque = getViewZOrBackground(selfDepthOpaque);
    let pixelSizeOpaque = getPixelSize(coords, selfDepthOpaque) * threshold;

    let selfDepthTransparentWithAlpha = getDepthTransparentWithAlpha(coords);
    let selfDepthTransparent = selfDepthTransparentWithAlpha.x;
    let selfViewZTransparent = getViewZOrBackground(selfDepthTransparent);
    let pixelSizeTransparent = getPixelSize(coords, selfDepthTransparent) * threshold;

    var bestOpaqueDepth = 0.0;
    var bestTransparentDepth = 0.0;
    var bestTransparentAlpha = 0.0;

    var opaqueOutlineFlag = 0.0;
    var transparentOutlineFlag = 0.0;

    for (var y = -1; y <= 1; y++) {
        for (var x = -1; x <= 1; x++) {
            let sampleCoords = coords + vec2<f32>(f32(x), f32(y)) * invTexSize;

            // opaque
            let sampleDepthOpaque = getDepthOpaque(sampleCoords);
            let sampleViewZOpaque = getViewZOrBackground(sampleDepthOpaque);
            if (abs(selfViewZOpaque - sampleViewZOpaque) > pixelSizeOpaque && selfDepthOpaque < sampleDepthOpaque && sampleDepthOpaque >= bestOpaqueDepth) {
                bestOpaqueDepth = sampleDepthOpaque;
                opaqueOutlineFlag = 1.0;
            }

            // transparent
            let sampleDepthTransparentWithAlpha = getDepthTransparentWithAlpha(sampleCoords);
            let sampleDepthTransparent = sampleDepthTransparentWithAlpha.x;
            let sampleViewZTransparent = getViewZOrBackground(sampleDepthTransparent);
            if (abs(selfViewZTransparent - sampleViewZTransparent) > pixelSizeTransparent && selfDepthTransparent < sampleDepthTransparent && sampleDepthTransparent >= bestTransparentDepth) {
                bestTransparentDepth = sampleDepthTransparent;
                bestTransparentAlpha = sampleDepthTransparentWithAlpha.y;
                transparentOutlineFlag = 1.0;
            }
        }
    }

    // disable transparent outline close to opaque elements
    if (transparentOutlineFlag > 0.0 && abs(selfViewZOpaque - getViewZ(bestTransparentDepth)) < pixelSizeOpaque) {
        transparentOutlineFlag = 0.0;
    }

    // lazy curvature veto: only computed if an outline was found and the center is not background
    let curvatureGate = 0.75;
    let dx = vec2<f32>(invTexSize.x, 0.0);
    let dy = vec2<f32>(0.0, invTexSize.y);

    if (opaqueOutlineFlag > 0.0 && !isBackground(selfDepthOpaque)) {
        let curvOpaque = getCurvature(
            getDepthOpaque(coords - dx), getDepthOpaque(coords + dx),
            getDepthOpaque(coords + dy), getDepthOpaque(coords - dy),
            selfViewZOpaque
        );
        if (curvOpaque < pixelSizeOpaque * curvatureGate) {
            opaqueOutlineFlag = 0.0;
            bestOpaqueDepth = 0.0;
        }
    }

    if (transparentOutlineFlag > 0.0 && !isBackground(selfDepthTransparent)) {
        let curvTransparent = getCurvature(
            getDepthTransparentWithAlpha(coords - dx).x, getDepthTransparentWithAlpha(coords + dx).x,
            getDepthTransparentWithAlpha(coords + dy).x, getDepthTransparentWithAlpha(coords - dy).x,
            selfViewZTransparent
        );
        if (curvTransparent < pixelSizeTransparent * curvatureGate) {
            transparentOutlineFlag = 0.0;
            bestTransparentDepth = 0.0;
            bestTransparentAlpha = 0.0;
        }
    }

    // pack depth in the g and b channels
    var depthPacked = vec2<f32>(0.0);
    var outlineTypeFlag = 0.0;
    if (opaqueOutlineFlag > 0.0 && transparentOutlineFlag > 0.0) {
        outlineTypeFlag = 0.75; // both
        depthPacked = packUnitIntervalToRG(bestOpaqueDepth);
    } else if (transparentOutlineFlag > 0.0) {
        outlineTypeFlag = 0.5; // transparent only
        depthPacked = packUnitIntervalToRG(bestTransparentDepth);
    } else if (opaqueOutlineFlag > 0.0) {
        outlineTypeFlag = 0.25; // opaque only
        depthPacked = packUnitIntervalToRG(bestOpaqueDepth);
    }

    // limited to [0.0, 0.5] to improve alpha precision since a wider range isn't needed
    let alpha = clamp(bestTransparentAlpha, 0.0, 0.5) * 2.0;
    let packedFlagWithAlpha = pack2x4(vec2<f32>(outlineTypeFlag, alpha));
    return vec4<f32>(packedFlagWithAlpha, depthPacked.x, depthPacked.y, bestTransparentDepth);
}
`;
