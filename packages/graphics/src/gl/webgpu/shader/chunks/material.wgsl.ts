/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { texture_sample_wgsl } from './quad.wgsl.js';

/** Functions of the fragment chunks below */
export const material_functions_wgsl = /* wgsl */ `
#if defined(dUsePalette)
    ${texture_sample_wgsl('tPalette')}
#endif
#if defined(dXrayShaded)
    fn calcXrayShadedAlpha(alpha: f32, normal: vec3<f32>) -> f32 {
        let d = abs(dot(normal, vec3<f32>(0.0, 0.0, 1.0)));
        #if defined(dXrayShaded_on)
            let a = alpha * (1.0 - pow(d, uGlobal.uXrayEdgeFalloff));
        #else
            let a = alpha * pow(d, uGlobal.uXrayEdgeFalloff);
        #endif
        return clamp(a, 0.001, 0.999);
    }
#endif

/** Dithered fade between levels of detail, coordinates with a bottom-left origin as on WebGL */
fn fadeLod(position: vec4<f32>, modelPosition: vec3<f32>) -> bool {
    let lod = uDraw.uLod;
    if (lod.w != 0.0 || (lod.x == 0.0 && lod.y == 0.0)) {
        return false;
    }
    let plane = uGlobal.uCameraPlane;
    let d = (dot(plane.xyz, modelPosition) + plane.w) / uGlobal.uModelScale;
    let ta = min(smoothstep(lod.x, lod.x + lod.z, d), 1.0 - smoothstep(lod.y - lod.z, lod.y, d));
    let coord = vec2<f32>(position.x, uGlobal.uDrawingBufferSize.y - position.y) + uGlobal.uViewOffset * 0.25;
    let thresholds = mat4x4<f32>(
        1.0 / 17.0, 9.0 / 17.0, 3.0 / 17.0, 11.0 / 17.0,
        13.0 / 17.0, 5.0 / 17.0, 15.0 / 17.0, 7.0 / 17.0,
        4.0 / 17.0, 12.0 / 17.0, 2.0 / 17.0, 10.0 / 17.0,
        16.0 / 17.0, 8.0 / 17.0, 14.0 / 17.0, 6.0 / 17.0
    );
    let at = thresholds[i32(coord.x) % 4][i32(coord.y) % 4];
    return ta < 0.99 && (ta < 0.01 || ta < at);
}
`;

/**
 * \`marker\`, \`material\`, \`emissive\`, \`metalness\`, \`roughness\` and \`bumpiness\` from the uniforms and
 * varyings, with overpaint, substance, x-ray shading and transparency applied; expects \`in\` and, for x-ray
 * shading, \`normal\`
 */
export const assign_material_wgsl = /* wgsl */ `
    #if defined(dNeedsMarker)
        var marker = uObject.uMarker;
        if (uObject.uMarker == -1.0) {
            marker = floor(in.vMarker * 255.0 + 0.5);
        }
    #endif

    #if defined(dUsePalette)
        var material = vec4<f32>(sample_tPalette(vec2<f32>(in.vPaletteV, 0.5)).rgb, uMaterial.uAlpha);
    #elif defined(dColorType_uniform)
        var material = vec4<f32>(uMaterial.uColor, uMaterial.uAlpha);
    #else
        var material = vec4<f32>(in.vColor.rgb, uMaterial.uAlpha);
    #endif

    #if defined(dOverpaint)
        material = vec4<f32>(mix(material.rgb, in.vOverpaint.rgb, in.vOverpaint.a), material.a);
    #endif

    var emissive = uMaterial.uEmissive;
    #if defined(dEmissive)
        emissive += in.vEmissive;
    #endif

    var metalness = uMaterial.uMetalness;
    var roughness = uMaterial.uRoughness;
    var bumpiness = uMaterial.uBumpiness;
    #if defined(dSubstance)
        let sf = clamp(in.vSubstance.a, 0.0, 0.99); // clamp to avoid artifacts
        metalness = mix(metalness, in.vSubstance.r, sf);
        roughness = mix(roughness, in.vSubstance.g, sf);
        bumpiness = mix(bumpiness, in.vSubstance.b, sf);
    #endif

    #if defined(dXrayShaded)
        material.a = calcXrayShadedAlpha(material.a, normal);
    #endif

    #if defined(dTransparency)
        var ta = 1.0 - in.vTransparency;
        if (in.vTransparency < 0.09) {
            ta = 1.0; // hard cutoff looks better
        }
        #if defined(dRenderVariant_pick)
            if (ta * uMaterial.uAlpha < uGlobal.uPickingAlphaThreshold) {
                discard; // ignore so the element below can be picked
            }
        #else
            material.a *= ta;
        #endif
    #endif
`;

/** Discards fragments of the other render mask and transparent back faces; expects \`interior\` */
export const check_transparency_wgsl = /* wgsl */ `
#if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
    #if defined(dTransparentBackfaces_off)
        if (interior && material.a < 1.0) {
            discard;
        }
    #elif defined(dTransparentBackfaces_opaque)
        if (interior) {
            material.a = 1.0;
        }
    #endif

    #if !defined(dXrayShaded)
        let renderMask = uGlobal.uRenderMask;
        if ((renderMask == MaskOpaque && material.a < 1.0) || (renderMask == MaskTransparent && material.a == 1.0)) {
            discard;
        }
    #endif
#elif defined(dRenderVariant_depth)
    #if defined(dTransparentBackfaces_off)
        if (interior && uGlobal.uDepthBack == 0u && uGlobal.uRenderMask != MaskOpaque) {
            discard;
        }
    #endif
#elif defined(dRenderVariant_marking)
    #if defined(dTransparentBackfaces_off)
        if (interior && material.a < 1.0) {
            material.a = 0.0;
        }
    #elif defined(dTransparentBackfaces_opaque)
        if (interior) {
            material.a = 1.0;
        }
    #endif
#endif
`;

/** Interior color and substance of back faces */
export const apply_interior_wgsl = /* wgsl */ `
    if (interior) {
        let interiorColor = uObject.uInteriorColor;
        material = vec4<f32>(mix(material.rgb, interiorColor.rgb, interiorColor.a), material.a);
        let interiorSubstance = uObject.uInteriorSubstance;
        let isf = clamp(interiorSubstance.a, 0.0, 0.99); // clamp to avoid artifacts
        metalness = mix(metalness, interiorSubstance.r, isf);
        roughness = mix(roughness, interiorSubstance.g, isf);
        bumpiness = mix(bumpiness, interiorSubstance.b, isf);
        #if defined(dTransparentBackfaces_opaque)
            material.a = 1.0;
        #endif
    }
`;

/**
 * Lit (or unlit with \`dIgnoreLight\`) \`color\` with bump mapping and exposure; expects \`normal\`,
 * \`modelPosition\` and \`viewPosition\`
 */
export const apply_light_wgsl = /* wgsl */ `
    let bumpFrequency = uMaterial.uBumpFrequency;
    let bumpAmplitude = uMaterial.uBumpAmplitude;
    let hasBump = bumpFrequency > 0.0 && bumpAmplitude > 0.0 && bumpiness > 0.0;
    #if defined(dIgnoreLight)
        if (hasBump) {
            let bump = fbm(modelPosition * bumpFrequency) * bumpAmplitude * bumpiness - 0.5 * bumpAmplitude * bumpiness;
            material = vec4<f32>(material.rgb + bump, material.a);
        }
        var color = vec4<f32>(material.rgb + material.rgb * emissive, material.a);
    #else
        if (hasBump) {
            normal = perturbNormal(modelPosition, viewPosition, normal, bumpFrequency, (bumpAmplitude * bumpiness) / bumpFrequency);
        }
        var color = applyLight(material, normal, viewPosition, metalness, roughness, emissive);
    #endif
    color = vec4<f32>(color.rgb * uGlobal.uExposure, color.a);
`;

/** Highlight, select and dim markers on \`color\`; expects \`marker\` */
export const apply_marker_wgsl = /* wgsl */ `
    #if defined(dColorMarker)
        if (marker > 0.0) {
            let priority = uGlobal.uMarkerPriority;
            if ((priority == 1 && marker != 2.0) || (priority != 1 && marker == 1.0)) {
                color = vec4<f32>(mix(color.rgb, uGlobal.uHighlightColor, uGlobal.uHighlightStrength), color.a);
            } else {
                color = vec4<f32>(mix(color.rgb, uGlobal.uSelectColor, uGlobal.uSelectStrength), color.a);
            }
        } else if (uGlobal.uMarkerAverage > 0.0) {
            color = vec4<f32>(mix(color.rgb, uGlobal.uDimColor, uGlobal.uDimStrength), color.a);
        }
    #endif
`;

/** Fog on \`color\`, premultiplied for transparent backgrounds; expects \`fragmentDepth\` */
export const apply_fog_wgsl = /* wgsl */ `
    let preFogAlpha = color.a;
    let transparentBackground = uGlobal.uTransparentBackground != 0u;
    if (uGlobal.uFog != 0u) {
        let viewZ = depthToViewZ(uGlobal.uIsOrtho, fragmentDepth, uGlobal.uNear, uGlobal.uFar);
        let fogFactor = smoothstep(uGlobal.uFogNear, uGlobal.uFogFar, abs(viewZ));
        let fogAlpha = (1.0 - fogFactor) * color.a;
        if (!transparentBackground) {
            if (color.a < 1.0) {
                // transparent objects are blended with background color
                color.a = fogAlpha;
            } else {
                // mix opaque objects with background color
                color = vec4<f32>(mix(color.rgb, uGlobal.uFogColor, fogFactor), color.a);
            }
        } else {
            #if defined(dRenderVariant_colorDpoit)
                if (color.a < 1.0) {
                    // transparent objects are blended with background color
                    color.a = fogAlpha;
                } else {
                    // opaque objects need to be pre-multiplied alpha
                    color = vec4<f32>(color.rgb * fogAlpha, fogAlpha);
                }
            #else
                // pre-multiplied alpha expected for transparent background
                color = vec4<f32>(color.rgb * fogAlpha, fogAlpha);
            #endif
        }
    } else if (transparentBackground) {
        #if !defined(dRenderVariant_colorDpoit)
            // pre-multiplied alpha expected for transparent background
            color = vec4<f32>(color.rgb * color.a, color.a);
        #endif
    }
`;

/** Highlight, select and dim markers and fog on \`color\`; expects \`marker\` and \`fragmentDepth\` */
export const apply_marker_fog_wgsl = apply_marker_wgsl + apply_fog_wgsl;
