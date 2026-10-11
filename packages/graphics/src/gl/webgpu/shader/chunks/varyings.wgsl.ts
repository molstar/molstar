/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clip_flag_wgsl } from './clip.wgsl.js';
import { assign_data_index_wgsl } from './read-from-texture.wgsl.js';

/**
 * \`VertexOutput\` with the varyings shared by the graphics programs at locations 0 to 8 and 15, and the program's
 * own \`extra\` members (from location 9)
 */
export function vertex_output_wgsl(extra: string) {
  return /* wgsl */ `
struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) vModelPosition: vec3<f32>,
    @location(1) vViewPosition: vec3<f32>,
#if defined(dRenderVariant_pick)
    @location(2) @interpolate(flat) vPickIds: vec2<f32>,
#endif
#if defined(dColorType_varying)
    @location(3) vColor: vec4<f32>,
#endif
#if defined(dNeedsMarker)
    #if defined(dMarkerType_instance) || !defined(dVaryingGroup)
        @location(4) vMarker: f32,
    #else
        @location(4) @interpolate(flat) vMarker: f32,
    #endif
#endif
#if defined(dOverpaint)
    @location(5) vOverpaint: vec4<f32>,
#endif
#if defined(dTransparency)
    @location(6) vTransparency: f32,
#endif
#if defined(dEmissive)
    @location(7) vEmissive: f32,
#endif
#if defined(dSubstance)
    @location(8) vSubstance: vec4<f32>,
#endif
#if dClipObjectCount != 0 && defined(dClipping)
    #if defined(dClippingType_instance) || !defined(dVaryingGroup)
        @location(15) vClipping: f32,
    #else
        @location(15) @interpolate(flat) vClipping: f32,
    #endif
#endif
${extra}
};
`;
}

/**
 * Item `vertexId` of the quad from `vertex_index` (6 vertices, two triangles per item) with its `position` and
 * `group` from `tPositionGroup`, and the quad corner `mapping` (-1 or 1 per axis)
 */
export const assign_quad_position_group_wgsl = /* wgsl */ `
    let m = i32(in.vertexIndex) % 6;
    var mapping = vec2<f32>(1.0, 1.0); // vertices 2 and 5
    if (m == 0) {
        mapping = vec2<f32>(-1.0, 1.0);
    } else if (m == 1 || m == 3) {
        mapping = vec2<f32>(-1.0, -1.0);
    } else if (m == 4) {
        mapping = vec2<f32>(1.0, -1.0);
    }
    let vertexId = i32(in.vertexIndex) / 6;
    let positionGroup = readFromTexture(tPositionGroup, vertexId);
    var position = positionGroup.xyz;
    let group = positionGroup.w;
`;

/**
 * Data texture indices and the color, marker and transparency varyings; expects `in.aInstance`, `vertexId`,
 * `group` and `out`
 */
export const assign_varyings_wgsl = /* wgsl */ `
    let instanceIndex = i32(in.aInstance);
    let groupInstanceIndex = i32(in.aInstance * f32(uObject.uGroupCount) + group);
    let vertexInstanceIndex = instanceIndex * uObject.uVertexCount + vertexId;

    #if defined(dRenderVariant_pick)
        out.vPickIds = vec2<f32>(in.aInstance, group);
    #endif

    #if defined(dNeedsMarker)
        ${assign_data_index_wgsl('dMarkerType', 'markerIndex')}
        out.vMarker = readFromTexture(tMarker, markerIndex).r;
    #endif

    #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
        #if defined(dColorType_instance)
            out.vColor = vec4<f32>(readFromTexture(tColor, instanceIndex).rgb, 1.0);
        #elif defined(dColorType_group) || defined(dColorType_groupInstance)
            #if defined(dColorType_group)
                let colorIndex = i32(group);
            #else
                let colorIndex = groupInstanceIndex;
            #endif
            #if defined(dDualColor)
                // two colors per group, except for color mode 2 (cylinders)
                var color2 = vec3<f32>(0.0);
                if (in.aColorMode == 2.0) {
                    out.vColor = vec4<f32>(readFromTexture(tColor, colorIndex).rgb, 1.0);
                } else {
                    out.vColor = vec4<f32>(readFromTexture(tColor, colorIndex * 2).rgb, 1.0);
                    color2 = readFromTexture(tColor, colorIndex * 2 + 1).rgb;
                }
            #else
                out.vColor = vec4<f32>(readFromTexture(tColor, colorIndex).rgb, 1.0);
            #endif
        #elif defined(dColorType_vertex)
            out.vColor = vec4<f32>(readFromTexture(tColor, vertexId).rgb, 1.0);
        #elif defined(dColorType_vertexInstance)
            out.vColor = vec4<f32>(readFromTexture(tColor, vertexInstanceIndex).rgb, 1.0);
        #endif

        #if defined(dOverpaint)
            ${assign_data_index_wgsl('dOverpaintType', 'overpaintIndex')}
            let overpaint = readFromTexture(tOverpaint, overpaintIndex);
            // pre-mix to avoid darkening due to empty overpaint
            #if defined(dColorType_uniform)
                let baseColor = uMaterial.uColor;
            #else
                let baseColor = out.vColor.rgb;
            #endif
            out.vOverpaint = vec4<f32>(mix(baseColor, overpaint.rgb, overpaint.a), overpaint.a) * uMaterial.uOverpaintStrength;
        #endif

        #if defined(dEmissive)
            ${assign_data_index_wgsl('dEmissiveType', 'emissiveIndex')}
            out.vEmissive = readFromTexture(tEmissive, emissiveIndex).r * uMaterial.uEmissiveStrength;
        #endif

        #if defined(dSubstance)
            ${assign_data_index_wgsl('dSubstanceType', 'substanceIndex')}
            let substance = readFromTexture(tSubstance, substanceIndex);
            // pre-mix to avoid artifacts due to empty substance
            let baseSubstance = vec3<f32>(uMaterial.uMetalness, uMaterial.uRoughness, uMaterial.uBumpiness);
            out.vSubstance = vec4<f32>(mix(baseSubstance, substance.rgb, substance.a), substance.a) * uMaterial.uSubstanceStrength;
        #endif
    #endif

    #if defined(dTransparency)
        ${assign_data_index_wgsl('dTransparencyType', 'transparencyIndex')}
        out.vTransparency = readFromTexture(tTransparency, transparencyIndex).r * uMaterial.uTransparencyStrength;
    #endif

    #if dClipObjectCount != 0 && defined(dClipping)
        ${assign_data_index_wgsl('dClippingType', 'clippingIndex')}
        out.vClipping = readFromTexture(tClipping, clippingIndex).r;
    #endif
    ${clip_flag_wgsl('out.vClipping')}
`;

/** \`size\` from the size data; expects the indices of \`assign_varyings_wgsl\` */
export const assign_size_wgsl = /* wgsl */ `
    #if defined(dSizeType_uniform)
        var size = uMaterial.uSize;
    #else
        #if defined(dSizeType_instance)
            let sizeIndex = instanceIndex;
        #elif defined(dSizeType_group)
            let sizeIndex = i32(group);
        #elif defined(dSizeType_groupInstance)
            let sizeIndex = groupInstanceIndex;
        #elif defined(dSizeType_vertex)
            let sizeIndex = vertexId;
        #elif defined(dSizeType_vertexInstance)
            let sizeIndex = vertexInstanceIndex;
        #endif
        // sizes are stored times 100
        var size = unpackRGBToInt(readFromTexture(tSize, sizeIndex).rgb) / 100.0;
    #endif
    size *= uMaterial.uSizeFactor;
`;
