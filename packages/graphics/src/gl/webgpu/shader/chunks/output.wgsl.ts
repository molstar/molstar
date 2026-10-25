/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Texel of `texture` at the fragment, as large as the target or a 1×1 placeholder (no texture to test against); the
 * coordinates are clamped, so that a placeholder reads its texel as on WebGL instead of an out-of-bounds value
 */
function fragment_texel_wgsl(texture: string, coords = 'vec2<i32>(in.position.xy)') {
  return `textureLoad(${texture}, min(${coords}, vec2<i32>(textureDimensions(${texture})) - 1), 0)`;
}

/**
 * `FragmentOutput` of the graphics programs: the color, or the object, instance and group ids and the packed depth
 * of the pick variant (one target each), or the shaded color, normal and material of the tracing variant, with the
 * program's own `extra` members
 */
export function fragment_output_wgsl(extra = '') {
  return /* wgsl */ `
struct FragmentOutput {
#if defined(dRenderVariant_pick)
    @location(0) object: vec4<f32>,
    @location(1) instance: vec4<f32>,
    @location(2) group: vec4<f32>,
    @location(3) depth: vec4<f32>,
#else
    @location(0) color: vec4<f32>,
    #if defined(dRenderVariant_colorWboit)
        @location(1) wboit: vec4<f32>,
    #elif defined(dRenderVariant_colorDpoit)
        @location(1) dpoitBack: vec4<f32>,
        @location(2) dpoitDepth: vec4<f32>,
    #elif defined(dRenderVariant_tracing)
        @location(1) normal: vec4<f32>,
        @location(2) material: vec4<f32>,
    #endif
#endif
${extra}
};
`;
}

/**
 * Inputs of the path tracer: the shaded `color`, the view space `normal` with `emissive` (8-bit, `packNormalToRGBA`),
 * and the unlit material color with the density (to estimate thickness); expects `color`, `material` and `emissive`
 */
export function tracing_output_wgsl(normal: string) {
  return /* wgsl */ `
    out.color = color;
    out.normal = packNormalToRGBA(${normal}, emissive);
    out.material = vec4<f32>(material.rgb, uMaterial.uDensity);
`;
}

/**
 * Weighted blended order-independent transparency: in the opaque pass only opaque fragments, in the transparent pass the
 * weighted color into `out.color` and the weight into `out.wboit` for transparent fragments in front of the opaque depth
 * (`tDepth`); expects `color`, `preFogAlpha`, `fragmentDepth` and `interior`
 */
export const wboit_output_wgsl = /* wgsl */ `
    #if defined(dRenderVariant_colorWboit)
        if (uGlobal.uRenderMask == MaskOpaque) {
            if (preFogAlpha < 1.0) {
                discard;
            }
        } else if (uGlobal.uRenderMask == MaskTransparent) {
            let opaqueDepth = ${fragment_texel_wgsl('tDepth')};
            // opaque fragments are drawn in the opaque pass; direct volumes have none: drawn here also where opaque
            #if defined(dGeometryType_directVolume)
                let wboitOpaque = false;
            #else
                let wboitOpaque = preFogAlpha == 1.0;
            #endif
            if (!wboitOpaque && fragmentDepth > opaqueDepth) {
                #if defined(dTransparentBackfaces_off)
                    if (interior) {
                        discard;
                    }
                #endif
                let alpha = color.a;
                // depth is reversed: larger in front
                let wboitWeight = alpha * clamp(pow(fragmentDepth, 2.0), 0.01, 1.0);
                color = vec4<f32>(color.rgb * alpha * wboitWeight, alpha);
                // extra alpha is to handle pre-multiplied alpha (direct volumes: always, their color is premultiplied)
                #if defined(dGeometryType_directVolume)
                    let backgroundAlpha = alpha;
                #else
                    let backgroundAlpha = select(1.0, alpha, uGlobal.uTransparentBackground != 0u);
                #endif
                out.wboit = vec4<f32>(backgroundAlpha * alpha * wboitWeight);
            } else {
                discard;
            }
        }
    #endif
`;

/**
 * Dual depth peeling: in the opaque pass only opaque fragments, in a transparent pass the front color (`color`, max
 * blended), the back color (`out.dpoitBack`) and the depth range still to peel (`out.dpoitDepth`, max blended) from the
 * last layer (`tDpoitDepth`, `tDpoitFrontColor`); expects `color`, `preFogAlpha`, `fragmentDepth` and `interior`
 *
 * adapted from https://github.com/tsherif/webgl2examples
 * The MIT License, Copyright 2017 Tarek Sherif, Shuai Shao
 */
export const dpoit_output_wgsl = /* wgsl */ `
    #if defined(dRenderVariant_colorDpoit)
        if (uGlobal.uRenderMask == MaskOpaque) {
            if (preFogAlpha < 1.0) {
                discard;
            }
        } else if (uGlobal.uRenderMask == MaskTransparent) {
            let dpoitCoord = vec2<i32>(in.position.xy);
            if (preFogAlpha != 1.0 && fragmentDepth > ${fragment_texel_wgsl('tDepth', 'dpoitCoord')}) {
                #if defined(dTransparentBackfaces_off)
                    if (interior) {
                        discard;
                    }
                #endif
                // peeling works with depth increasing away from the camera (depth is reversed)
                let peelDepth = 1.0 - fragmentDepth;
                let lastDepth = textureLoad(tDpoitDepth, dpoitCoord, 0).rg;
                let lastFrontColor = textureLoad(tDpoitFrontColor, dpoitCoord, 0);
                let fragColor = color;

                // depth value always increases, so we can use MAX blend equation
                out.dpoitDepth = vec4<f32>(vec2<f32>(-MaxDpoitDepth), 0.0, 0.0);
                // front color always increases, so we can use MAX blend equation
                color = lastFrontColor;
                // back color is separately blend afterwards each pass
                out.dpoitBack = vec4<f32>(0.0);

                let nearestDepth = -lastDepth.x;
                let furthestDepth = lastDepth.y;
                let alphaMultiplier = 1.0 - lastFrontColor.a;

                if (peelDepth < nearestDepth || peelDepth > furthestDepth) {
                    // skip this depth since it's been peeled
                } else if (peelDepth > nearestDepth && peelDepth < furthestDepth) {
                    // this needs to be peeled, the ones remaining after MAX blended for all need-to-peel will be
                    // peeled next pass
                    out.dpoitDepth = vec4<f32>(-peelDepth, peelDepth, 0.0, 0.0);
                } else if (peelDepth == nearestDepth) {
                    // write to front color buffer
                    color = vec4<f32>(color.rgb + fragColor.rgb * fragColor.a * alphaMultiplier, 1.0 - alphaMultiplier * (1.0 - fragColor.a));
                } else {
                    // write to back color buffer
                    out.dpoitBack += fragColor;
                }
            } else {
                discard;
            }
        }
    #endif
`;

/**
 * Depth passes: fragments behind the opaque depth (\`tDepth\`) are discarded; packed depth of opaque fragments for the
 * opaque mask, packed depth and alpha of transparent fragments for the transparent mask, packed depth of all fragments
 * otherwise, into \`out.color\`; expects \`fragmentDepth\`, and for x-ray shading, \`normal\`. \`materialOpacity\` is an
 * expression of the per-fragment opacity factor.
 */
export function depth_output_wgsl(materialOpacity = '1.0') {
  return /* wgsl */ `
    // behind the opaque depth (depth is reversed)
    if (fragmentDepth < ${fragment_texel_wgsl('tDepth')}) {
        discard;
    }
    let depthAlpha = uMaterial.uAlpha * ${materialOpacity};
    let depthRenderMask = uGlobal.uRenderMask;
    if (depthRenderMask == MaskOpaque) {
        #if defined(dXrayShaded)
            discard;
        #endif
        #if defined(dTransparency)
            var dta = 1.0 - in.vTransparency;
            #if defined(dVaryingGroup)
                if (in.vTransparency < 0.1) {
                    dta = 1.0; // hard cutoff to avoid artifacts
                }
            #endif
            if (depthAlpha * dta < 1.0) {
                discard;
            }
        #else
            if (depthAlpha < 1.0) {
                discard;
            }
        #endif
        out.color = packDepthToRGBA(fragmentDepth);
    } else if (depthRenderMask == MaskTransparent) {
        var alpha = depthAlpha;
        #if defined(dTransparency)
            alpha *= 1.0 - in.vTransparency;
        #endif
        #if defined(dXrayShaded)
            alpha = calcXrayShadedAlpha(alpha, normal);
        #else
            if (alpha == 1.0) {
                discard;
            }
        #endif
        out.color = packDepthWithAlphaToRGBA(fragmentDepth, alpha);
    } else {
        out.color = packDepthToRGBA(fragmentDepth);
    }
`;
}

/**
 * Emissive prepass for bloom: the premultiplied emissive glow in the object's own color, faded by fog and dimmed behind
 * a transparent blocker (`tDepthAlpha`); expects `fragmentDepth` and, for x-ray shading, `normal`
 */
export const emissive_output_wgsl = /* wgsl */ `
    var emissiveValue = uMaterial.uEmissive;
    #if defined(dEmissive)
        emissiveValue += in.vEmissive;
    #endif
    var emissiveAlpha = uMaterial.uAlpha;
    #if defined(dXrayShaded)
        emissiveAlpha = calcXrayShadedAlpha(emissiveAlpha, normal);
    #endif
    // fade emissive bloom with fog so the glow dims into the background like the geometry
    if (uGlobal.uFog != 0u) {
        let emissiveViewZ = depthToViewZ(uGlobal.uIsOrtho, fragmentDepth, uGlobal.uNear, uGlobal.uFar);
        emissiveAlpha *= 1.0 - smoothstep(uGlobal.uFogNear, uGlobal.uFogFar, abs(emissiveViewZ));
    }
    // dim emitters behind a transparent blocker by its coverage (tDepthAlpha packs front depth+alpha)
    let emissiveBlocker = unpackRGBAToDepthWithAlpha(${fragment_texel_wgsl('tDepthAlpha')});
    if (fragmentDepth < emissiveBlocker.x - 0.0001) {
        emissiveAlpha *= 1.0 - emissiveBlocker.y;
    }
    // glow with the object's own color so emissive isn't double-counted via the lit buffer
    #if defined(dUsePalette)
        var emissiveColor = sample_tPalette(vec2<f32>(in.vPaletteV, 0.5)).rgb;
    #elif defined(dColorType_uniform)
        var emissiveColor = uMaterial.uColor;
    #elif defined(dColorType_varying)
        var emissiveColor = in.vColor.rgb;
    #else
        var emissiveColor = vec3<f32>(1.0);
    #endif
    #if defined(dOverpaint)
        emissiveColor = mix(emissiveColor, in.vOverpaint.rgb, in.vOverpaint.a);
    #endif
    let emissiveFactor = emissiveValue * emissiveAlpha;
    var emissiveMaterial = vec4<f32>(emissiveColor * emissiveFactor, emissiveFactor);
    #if defined(dTransparency)
        var emissiveTa = 1.0 - in.vTransparency;
        if (in.vTransparency < 0.09) {
            emissiveTa = 1.0; // hard cutoff looks better
        }
        emissiveMaterial *= emissiveTa;
    #endif
    out.color = emissiveMaterial;
`;

/**
 * Marking depth (`uMarkingType` 1, packed depth and opacity of unmarked fragments) or mask (packed marking of marked
 * fragments, depth tested against `tDepthAlpha`) into `out.color`; expects `material`, `marker` and `fragmentDepth`
 */
export const marking_output_wgsl = /* wgsl */ `
    let markingOpacity = clamp(material.a, 0.0, 1.0);
    if (uGlobal.uMarkingType == 1) {
        if (marker > 0.0 || markingOpacity == 0.0) {
            discard;
        }
        out.color = packDepthWithAlphaToRGBA(fragmentDepth, markingOpacity);
    } else {
        if (marker == 0.0) {
            discard;
        }
        var depthTest = 1.0;
        if (uGlobal.uMarkingDepthTest != 0u) {
            let markingDepth = unpackRGBAToDepthWithAlpha(${fragment_texel_wgsl('tDepthAlpha')}).x;
            depthTest = select(0.0, 1.0, fragmentDepth <= markingDepth);
        }
        let priority = uGlobal.uMarkerPriority;
        let isHighlight = (priority == 1 && marker != 2.0) || (priority != 1 && marker == 1.0);
        let markingViewZ = depthToViewZ(uGlobal.uIsOrtho, fragmentDepth, uGlobal.uNear, uGlobal.uFar);
        let markingFogFactor = smoothstep(uGlobal.uFogNear, uGlobal.uFogFar, abs(markingViewZ));
        if (markingFogFactor == 1.0) {
            discard;
        }
        out.color = packMarkingMask(depthTest, isHighlight, 1.0 - markingFogFactor, markingOpacity);
    }
`;

/**
 * Pick ids and packed depth into `out`, discarding fragments that are not opaque enough to be picked; expects
 * `in.vPickIds`, `fragmentDepth`, `interior` and, for x-ray shading, `normal`
 */
export const pick_output_wgsl = /* wgsl */ `
    let pickViewZ = depthToViewZ(uGlobal.uIsOrtho, fragmentDepth, uGlobal.uNear, uGlobal.uFar);
    let pickFogFactor = smoothstep(uGlobal.uFogNear, uGlobal.uFogFar, abs(pickViewZ));
    let pickFogAlpha = (1.0 - pickFogFactor) * uMaterial.uAlpha;
    var pickAlpha = uMaterial.uAlpha;
    #if defined(dXrayShaded)
        // add bias to make picking xray shaded elements easier
        pickAlpha = calcXrayShadedAlpha(pickAlpha, normal) + (0.3 * uGlobal.uPickingAlphaThreshold);
    #endif
    // if not opaque enough ignore so the element below can be picked
    if (pickAlpha < uGlobal.uPickingAlphaThreshold || pickFogAlpha < 0.1) {
        #if defined(dTransparentBackfaces_opaque)
            if (!interior) {
                discard;
            }
        #else
            discard;
        #endif
    }
    out.object = vec4<f32>(packIntToRGB(f32(uObject.uObjectId)), 1.0);
    out.instance = vec4<f32>(packIntToRGB(in.vPickIds.x), 1.0);
    out.group = vec4<f32>(packIntToRGB(in.vPickIds.y), 1.0);
    out.depth = packDepthToRGBA(fragmentDepth);
`;
