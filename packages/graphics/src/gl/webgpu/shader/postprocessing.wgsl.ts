/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';
import { pack_2x4_wgsl } from './outlines.wgsl.js';

const Textures = [
  'tColor',
  'tTransparentColor',
  'tDepthOpaque',
  'tDepthTransparent',
  'tSsaoDepth',
  'tSsaoDepthTransparent',
  'tShadows',
  'tOutlines',
  'tBloom',
];

/**
 * Composes occlusion, shadows, outlines, blended transparency and bloom onto the scene color. Depth textures are
 * always available on WebGPU, so the opaque depth is never packed.
 */
export const postprocessing_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${pack_2x4_wgsl}
${quad_vert_wgsl}
${Textures.map(texture_read_wgsl).join('\n')}

override oOutlineScale: i32 = 1;

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

fn isBackground(depth: f32) -> bool {
    // (2^24 - 1) / 2^24, max of 24-bit packed depth; also passes raw fp32
    return depth >= 0.99999994;
}

fn getFogFactor(depth: f32) -> f32 {
    return smoothstep(uObject.uFogNear, uObject.uFogFar, abs(getViewZ(depth)));
}

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

    var o = Outline(false, false, 1.0, 1.0, 0.0);
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

            if ((isOpaque || isBoth) && sampleOpaqueDepth < o.opaqueDepth) {
                o.hasOpaque = true;
                o.opaqueDepth = sampleOpaqueDepth;
            }

            if ((isTransparent || isBoth) && sampleTransparentDepth < o.transparentDepth) {
                o.hasTransparent = true;
                o.transparentDepth = sampleTransparentDepth;
                o.alpha = sampleAlpha;
            }
        }
    }
    return o;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = glCoords(in.position, uObject.uTexSize);
    var color = read_tColor(coords);
    let opaqueDepth = read_tDepthOpaque(coords).r;
    let transparentBackground = uObject.uTransparentBackground != 0u;

    #if defined(dBlendTransparency)
        var blendTransparency = true;
        var transparentColor = read_tTransparentColor(coords);
        #if defined(dOcclusionEnable) && defined(dOcclusionIncludeTransparency)
            let transparentDepth = unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords)).x;
        #endif
    #endif

    #if defined(dOcclusionEnable) || defined(dShadowEnable)
        let isOpaqueBackground = isBackground(opaqueDepth);
        let fogFactor = getFogFactor(opaqueDepth);
    #endif

    #if defined(dOcclusionEnable)
        if (!isOpaqueBackground) {
            let occlusionFactor = unpackSsao(read_tSsaoDepth(coords + uObject.uOcclusionOffset));
            if (!transparentBackground) {
                color = vec4<f32>(mix(mix(uObject.uOcclusionColor, uObject.uFogColor, fogFactor), color.rgb, occlusionFactor), color.a);
            } else {
                color = vec4<f32>(mix(uObject.uOcclusionColor * (1.0 - fogFactor), color.rgb, occlusionFactor), color.a);
            }
        }
        #if defined(dBlendTransparency) && defined(dOcclusionIncludeTransparency)
            if (!isBackground(transparentDepth)) {
                let transparentFogFactor = getFogFactor(transparentDepth);
                let occlusionFactor = unpackSsao(read_tSsaoDepthTransparent(coords + uObject.uOcclusionOffset));
                transparentColor = vec4<f32>(mix(uObject.uOcclusionColor * (1.0 - transparentFogFactor), transparentColor.rgb, occlusionFactor), transparentColor.a);
            }
        #endif
    #endif

    #if defined(dShadowEnable)
        if (!isOpaqueBackground) {
            let shadow = read_tShadows(coords);
            if (!transparentBackground) {
                color = vec4<f32>(mix(mix(vec3<f32>(0.0), uObject.uFogColor, fogFactor), color.rgb, shadow.a), color.a);
            } else {
                color = vec4<f32>(mix(vec3<f32>(0.0) * (1.0 - fogFactor), color.rgb, shadow.a), color.a);
            }
        }
    #endif

    // outline needs to be handled after occlusion and shadow to keep them clean
    #if defined(dOutlineEnable)
        let outline = getOutline(coords);

        if (outline.hasOpaque) {
            let outlineFogFactor = getFogFactor(outline.opaqueDepth);
            if (!transparentBackground) {
                color = vec4<f32>(mix(uObject.uOutlineColor, uObject.uFogColor, outlineFogFactor), color.a);
            } else {
                color = vec4<f32>(mix(uObject.uOutlineColor, vec3<f32>(0.0), outlineFogFactor), 1.0 - outlineFogFactor);
            }
        }

        #if defined(dBlendTransparency)
            if (outline.hasTransparent) {
                if (outline.hasOpaque && outline.opaqueDepth < outline.transparentDepth) {
                    blendTransparency = false;
                } else {
                    let finalOutlineAlpha = clamp(outline.alpha * 2.0, 0.0, 1.0);
                    let outlineFogFactor = getFogFactor(outline.transparentDepth);
                    let finalAlpha = max(transparentColor.a, finalOutlineAlpha * (1.0 - outlineFogFactor));
                    transparentColor = vec4<f32>(uObject.uOutlineColor * finalAlpha, finalAlpha);
                }
            }
        #endif
    #endif

    #if defined(dBlendTransparency)
        let alpha = transparentColor.a;
        if (blendTransparency && alpha != 0.0) {
            color = transparentColor + color * (1.0 - alpha);
        }
    #endif

    #if defined(dBloomEnable)
        let bloom = read_tBloom(coords);
        // additive glow over geometry and background alike; alpha just covers the premultiplied rgb
        color = vec4<f32>(color.rgb + bloom.rgb, min(max(color.a, max(color.r + bloom.r, max(color.g + bloom.g, color.b + bloom.b))), 1.0));
    #endif

    return color;
}
`;
