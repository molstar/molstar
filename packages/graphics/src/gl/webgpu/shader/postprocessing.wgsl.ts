/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';
import { outline_read_wgsl, pack_2x4_wgsl } from './outlines.wgsl.js';
import type { DefineLookup } from '../../renderable/schema.js';

/** The textures a variant reads through `read_<name>` (the upsampled occlusion samples its texture directly) */
function getTextures(d: DefineLookup) {
  const blend = d('dBlendTransparency');
  const occlusion = d('dOcclusionEnable');
  const transparentOcclusion = blend && occlusion && d('dOcclusionIncludeTransparency');
  const upsample = d('dOcclusionUpsample');
  return [
    'tColor',
    'tDepthOpaque',
    blend && 'tTransparentColor',
    transparentOcclusion && 'tDepthTransparent',
    occlusion && !upsample && 'tSsaoDepth',
    transparentOcclusion && !upsample && 'tSsaoDepthTransparent',
    d('dShadowEnable') && 'tShadows',
    d('dOutlineEnable') && 'tOutlines',
    d('dBloomEnable') && 'tBloom',
  ].filter((t): t is string => !!t);
}

/**
 * Composes occlusion, shadows, outlines, blended transparency and bloom onto the scene color. Depth textures are
 * always available on WebGPU, so the opaque depth is never packed. The texture reads and the outline helper are
 * included by the variant's defines (`WgslProgram.source`).
 */
export const postprocessing_wgsl = (d: DefineLookup) => /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${getTextures(d).map(texture_read_wgsl).join('\n')}
${d('dOutlineEnable') ? `${pack_2x4_wgsl}\n${outline_read_wgsl}` : ''}

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

fn isBackground(depth: f32) -> bool {
    // the far plane (depth is reversed), packed or raw
    return depth == 0.0;
}

fn getFogFactor(depth: f32) -> f32 {
    return smoothstep(uObject.uFogNear, uObject.uFogFar, abs(getViewZ(depth)));
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, uObject.uTexSize);
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
            #if defined(dOcclusionUpsample)
                let occlusionFactor = getSsaoUpsampled(tSsaoDepth, coords + uObject.uOcclusionOffset, getViewZ(opaqueDepth), uObject.uNear, uObject.uFar);
            #else
                let occlusionFactor = unpackSsao(read_tSsaoDepth(coords + uObject.uOcclusionOffset));
            #endif
            if (!transparentBackground) {
                color = vec4<f32>(mix(mix(uObject.uOcclusionColor, uObject.uFogColor, fogFactor), color.rgb, occlusionFactor), color.a);
            } else {
                color = vec4<f32>(mix(uObject.uOcclusionColor * (1.0 - fogFactor), color.rgb, occlusionFactor), color.a);
            }
        }
        #if defined(dBlendTransparency) && defined(dOcclusionIncludeTransparency)
            if (!isBackground(transparentDepth)) {
                let transparentFogFactor = getFogFactor(transparentDepth);
                #if defined(dOcclusionUpsample)
                    let occlusionFactor = getSsaoUpsampled(tSsaoDepthTransparent, coords + uObject.uOcclusionOffset, getViewZ(transparentDepth), uObject.uNear, uObject.uFar);
                #else
                    let occlusionFactor = unpackSsao(read_tSsaoDepthTransparent(coords + uObject.uOcclusionOffset));
                #endif
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
                if (outline.hasOpaque && outline.opaqueDepth > outline.transparentDepth) {
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
