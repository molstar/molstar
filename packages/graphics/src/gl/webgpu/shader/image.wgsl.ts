/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clip_functions_wgsl, clip_pixel_wgsl } from './chunks/clip.wgsl.js';
import { common_wgsl } from './chunks/common.wgsl.js';
import { apply_fog_wgsl, apply_marker_wgsl, material_functions_wgsl } from './chunks/material.wgsl.js';
import {
  dpoit_output_wgsl,
  fragment_output_wgsl,
  marking_output_wgsl,
  tracing_output_wgsl,
  wboit_output_wgsl,
} from './chunks/output.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { texture_sample_wgsl } from './chunks/quad.wgsl.js';
import { read_from_texture_wgsl } from './chunks/read-from-texture.wgsl.js';

/**
 * Images: a textured quad per instance with colors (\`tImageTex\`), groups (\`tGroupTex\`) and values (\`tValueTex\`,
 * for iso-level cut-offs) per texel, optionally cubic interpolation, trimmed by a shape
 */
export const image_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${material_functions_wgsl}
${clip_functions_wgsl}
${packing_wgsl}
${texture_sample_wgsl('tImageTex')}
${texture_sample_wgsl('tGroupTex')}
${texture_sample_wgsl('tValueTex')}

const PaletteScale = 16777214.0; // (1 << 24) - 2

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) vModelPosition: vec3<f32>,
    @location(1) vViewPosition: vec3<f32>,
    @location(2) vUv: vec2<f32>,
    @location(3) @interpolate(flat) vInstance: f32,
    @location(4) vPosition: vec3<f32>,
};

${fragment_output_wgsl()}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    let transform = mat4x4<f32>(in.aTransform0, in.aTransform1, in.aTransform2, in.aTransform3);
    let model = uGlobal.uModel * transform;
    let position4 = vec4<f32>(in.aPosition, 1.0);
    out.vModelPosition = (model * position4).xyz;
    let mvPosition = uGlobal.uView * model * position4;
    out.vViewPosition = mvPosition.xyz;
    out.position = uGlobal.uProjection * mvPosition;
    out.vUv = in.aUv;
    out.vInstance = in.aInstance;
    out.vPosition = in.aPosition;
    return out;
}

#if defined(dInterpolation_catmulrom) || defined(dInterpolation_mitchell) || defined(dInterpolation_bspline)
    #define dInterpolation_cubic
#endif

#if defined(dInterpolation_cubic)
    #if defined(dInterpolation_catmulrom) || defined(dInterpolation_mitchell)
        #if defined(dInterpolation_catmulrom)
            const B = 0.0;
            const C = 0.5;
        #else
            const B = 0.333;
            const C = 0.333;
        #endif

        fn cubicFilter(x: f32) -> f32 {
            let f = abs(x);
            if (f < 1.0) {
                return ((12.0 - 9.0 * B - 6.0 * C) * (f * f * f) + (-18.0 + 12.0 * B + 6.0 * C) * (f * f) + (6.0 - 2.0 * B)) / 6.0;
            } else if (f >= 1.0 && f < 2.0) {
                return ((-B - 6.0 * C) * (f * f * f) + (6.0 * B + 30.0 * C) * (f * f) + (-(12.0 * B) - 48.0 * C) * f + 8.0 * B + 24.0 * C) / 6.0;
            }
            return 0.0;
        }
    #else
        fn cubicFilter(x: f32) -> f32 {
            let f = abs(x);
            if (f >= 0.0 && f <= 1.0) {
                return (2.0 / 3.0) + 0.5 * (f * f * f) - (f * f);
            } else if (f > 1.0 && f <= 2.0) {
                return 1.0 / 6.0 * pow(2.0 - f, 3.0);
            }
            return 1.0;
        }
    #endif

    fn biCubic(texCoordIn: vec2<f32>) -> vec4<f32> {
        let texDim = uObject.uImageTexDim;
        let texelSize = 1.0 / texDim;
        let texCoord = texCoordIn - texelSize / 2.0;
        var nSum = vec4<f32>(0.0);
        var nDenom = 0.0;
        let cell = fract(texCoord * texDim);
        for (var m = -1.0; m <= 2.0; m += 1.0) {
            for (var n = -1.0; n <= 2.0; n += 1.0) {
                let vecData = sample_tImageTex(texCoord + texelSize * vec2<f32>(m, n));
                let c = abs(cubicFilter(m - cell.x) * cubicFilter(-n + cell.y));
                nSum += vecData * c;
                nDenom += c;
            }
        }
        return nSum / nDenom;
    }
#endif

@fragment
fn fragmentMain(in: VertexOutput) -> FragmentOutput {
    var out: FragmentOutput;
    // images have no per-group clipping flags
    let clipFlag = 0;

    let trimType = uObject.uTrimType;
    if (trimType != 0 && getShapeSignedDistance(in.vPosition, trimType, uObject.uTrimCenter, uObject.uTrimRotation, uObject.uTrimScale, uObject.uTrimTransform) > 0.0) {
        discard;
    }

    if (fadeLod(in.position, in.vModelPosition)) {
        discard;
    }
    ${clip_pixel_wgsl('in.vModelPosition')}

    let uv = in.vUv;
    #if defined(dInterpolation_cubic)
        #if defined(dUsePalette)
            var material = sample_tImageTex(uv);
            if (any(material.rgb != vec3<f32>(1.0))) {
                material = biCubic(uv);
            }
        #else
            var material = biCubic(uv);
        #endif
    #else
        var material = sample_tImageTex(uv);
    #endif

    let alpha = uMaterial.uAlpha;
    if (uObject.uIsoLevel >= 0.0) {
        if (sample_tValueTex(uv).r < uObject.uIsoLevel) {
            discard;
        }
        material.a = alpha;
    } else {
        if (material.a == 0.0) {
            discard;
        }
        material.a *= alpha;
    }

    let fragmentDepth = in.position.z;
    let interior = false;

    let packedGroup = sample_tGroupTex(uv).rgb;
    let group = select(unpackRGBToInt(packedGroup), -1.0, all(packedGroup == vec3<f32>(0.0)));
    let instance = in.vInstance;
    let groupCount = f32(uObject.uGroupCount);
    // the marker and layer data per instance or group instance
    let dataIndex = select(i32(instance * groupCount + group), i32(instance), uObject.uInstanceGranularity != 0u);

    // per-group transparency
    #if defined(dTransparency) && !defined(dRenderVariant_depth)
        #if defined(dTransparencyType_texture)
            // an empty pixel (no group) has no data
            var transparency = select(read_tTransparency(dataIndex).r, 0.0, group == -1.0);
        #else
            var transparency = 0.0;
        #endif
        transparency *= uMaterial.uTransparencyStrength;

        var ta = 1.0 - transparency;
        if (transparency < 0.09) {
            ta = 1.0; // hard cutoff looks better
        }

        #if defined(dRenderVariant_pick)
            if (ta * alpha < uGlobal.uPickingAlphaThreshold) {
                discard; // ignore so the element below can be picked
            }
        #elif defined(dRenderVariant_emissive)
            if (ta < 1.0) {
                discard; // emissive not supported with transparency
            }
        #else
            material.a *= ta;
        #endif
    #endif

    let renderMask = uGlobal.uRenderMask;
    if ((renderMask == MaskOpaque && material.a < 1.0) || (renderMask == MaskTransparent && material.a == 1.0)) {
        discard;
    }

    #if defined(dNeedsMarker)
        var marker = uObject.uMarker;
        if (group == -1.0) {
            marker = 0.0;
        } else if (uObject.uMarker == -1.0) {
            marker = read_tMarker(dataIndex).r;
            marker = floor(marker * 255.0 + 0.5);
        }
    #endif

    #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
        var emissive = uMaterial.uEmissive;
        if (group == -1.0) {
            emissive = 0.0;
        } else {
            #if defined(dEmissive)
                #if defined(dEmissiveType_texture)
                    emissive += read_tEmissive(dataIndex).r * uMaterial.uEmissiveStrength;
                #endif
            #endif
        }
    #endif

    #if defined(dRenderVariant_pick)
        if (group == -1.0) {
            discard;
        }
        let pickViewZ = depthToViewZ(uGlobal.uIsOrtho, fragmentDepth, uGlobal.uNear, uGlobal.uFar);
        let pickFogAlpha = (1.0 - smoothstep(uGlobal.uFogNear, uGlobal.uFogFar, abs(pickViewZ))) * alpha;
        // if not opaque enough ignore so the element below can be picked
        if (alpha < uGlobal.uPickingAlphaThreshold || pickFogAlpha < 0.1) {
            discard;
        }
        out.object = vec4<f32>(packIntToRGB(f32(uObject.uObjectId)), 1.0);
        out.instance = vec4<f32>(packIntToRGB(instance), 1.0);
        out.group = vec4<f32>(packIntToRGB(group), 1.0);
        out.depth = packDepthToRGBA(fragmentDepth);
    #elif defined(dRenderVariant_depth)
        if (renderMask == MaskTransparent) {
            out.color = packDepthWithAlphaToRGBA(fragmentDepth, material.a);
        } else {
            out.color = packDepthToRGBA(fragmentDepth);
        }
    #elif defined(dRenderVariant_marking)
        ${marking_output_wgsl}
    #elif defined(dRenderVariant_emissive)
        out.color = vec4<f32>(emissive);
    #else
        #if defined(dUsePalette)
            if (all(material.rgb == vec3<f32>(1.0))) {
                material = vec4<f32>(uObject.uPaletteDefault, material.a);
            } else {
                let v = ((material.r * 256.0 * 256.0 * 255.0 + material.g * 256.0 * 255.0 + material.b * 255.0) - 1.0) / PaletteScale;
                material = vec4<f32>(sample_tPalette(vec2<f32>(v, 0.0)).rgb, material.a);
            }
        #endif

        // mix material with overpaint
        #if defined(dOverpaint)
            var overpaint = vec4<f32>(0.0);
            if (group != -1.0) {
                #if defined(dOverpaintType_texture)
                    overpaint = read_tOverpaint(dataIndex);
                #endif
                overpaint *= uMaterial.uOverpaintStrength;
            }
            material = vec4<f32>(mix(material.rgb, overpaint.rgb, overpaint.a), material.a);
        #endif

        var color = material;
        ${apply_marker_wgsl}
        #if defined(dRenderVariant_tracing)
            ${tracing_output_wgsl('normalize(in.vViewPosition)')}
        #else
            ${apply_fog_wgsl}
            ${wboit_output_wgsl}
            ${dpoit_output_wgsl}
            out.color = color;
        #endif
    #endif
    return out;
}
`;
