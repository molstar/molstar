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
    #if !defined(dVaryingGroup)
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
#if defined(dUsePalette)
    @location(14) vPaletteV: f32,
#endif
#if defined(dClipObjects) && defined(dClipping)
    #if !defined(dVaryingGroup)
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
    let positionGroup = read_tPositionGroup(vertexId);
    var position = positionGroup.xyz;
    let group = positionGroup.w;
`;

/**
 * Grid textures (smoothed colors etc.) by the prefix of their uniforms and the type define (the color type's grid value
 * is 'grid', see `getColorTypeDefine`; the layers' are their grid types), see `grid_value_wgsl`
 */
const GridTextures = [
  { texture: 'tColorGrid', prefix: 'Color', type: 'dColorType' },
  { texture: 'tOverpaintGrid', prefix: 'Overpaint', type: 'dOverpaintType' },
  { texture: 'tEmissiveGrid', prefix: 'Emissive', type: 'dEmissiveType' },
  { texture: 'tSubstanceGrid', prefix: 'Substance', type: 'dSubstanceType' },
  { texture: 'tTransparencyGrid', prefix: 'Transparency', type: 'dTransparencyType' },
];

/**
 * `texture3dFrom2dLinear_<texture>(pos)`: a grid stored as z-slices in a 2d atlas (rows in upload order, as the CPU
 * and the GPU color smoothing write them), linear within a slice and between slices, as `texture3dFrom2dLinear` on
 * WebGL; for the mesh programs, which can have grid data
 */
export const grid_functions_wgsl = GridTextures.map(
  ({ texture, prefix, type }) => /* wgsl */ `
#if defined(${type}_volume) || defined(${type}_volumeInstance) || defined(${type}_grid)
    fn texture3dFrom2dLinear_${texture}(pos: vec3<f32>) -> vec4<f32> {
        let gridDim = uObject.u${prefix}GridDim;
        let texDim = uObject.u${prefix}GridTexDim;
        let zSlice0 = floor(pos.z * gridDim.z);
        // the tile of a slice in integers: f32 division may be off by an ULP, which put the slices starting an atlas
        // row (x / x just below 1) past the end of the row above
        let width = i32(texDim.x);
        let start0 = i32(zSlice0) * i32(gridDim.x);
        let start1 = start0 + i32(gridDim.x);
        let column0 = f32(start0 % width);
        let row0 = f32(start0 / width);
        let column1 = f32(start1 % width);
        let row1 = f32(start1 / width);
        let offset = pos.xy * gridDim.xy;
        let color0 = textureSampleLevel(${texture}, sLinear, (vec2<f32>(column0, row0 * gridDim.y) + offset) / texDim, 0.0);
        let color1 = textureSampleLevel(${texture}, sLinear, (vec2<f32>(column1, row1 * gridDim.y) + offset) / texDim, 0.0);
        return mix(color0, color1, abs(pos.z * gridDim.z - zSlice0));
    }
#endif
`,
).join('');

/** The value of a grid texture at the model space `position` */
function grid_value_wgsl(texture: string, position: string) {
  const prefix = GridTextures.find((g) => g.texture === texture)!.prefix;
  const transform = `uObject.u${prefix}GridTransform`;
  return `texture3dFrom2dLinear_${texture}((${transform}.w * (${position} - ${transform}.xyz)) / uObject.u${prefix}GridDim)`;
}

/** grid positions: the local position for 'volume', the model position for 'volumeInstance' */
const GridPosition = 'gridPosition';
const GridModelPosition = '(out.vModelPosition / uGlobal.uModelScale)';

/**
 * Data texture indices and the color, marker and transparency varyings; expects `in.aInstance`, `vertexId`,
 * `group` and `out`, and for grid data (mesh programs) `gridPosition` and `out.vModelPosition`
 */
export const assign_varyings_wgsl = /* wgsl */ `
    let instanceIndex = i32(in.aInstance);
    #if defined(dSegmented)
        // the data of merged render items, offset per segment (see \`SegmentSchema\`)
        let groupInstanceIndex = i32(in.aSegment.y + in.aInstance * in.aSegment.x + group);
        let groupIndex = i32(in.aSegment.z + group);
    #else
        let groupInstanceIndex = i32(in.aInstance * f32(uObject.uGroupCount) + group);
        let groupIndex = i32(group);
    #endif
    let vertexInstanceIndex = instanceIndex * uObject.uVertexCount + vertexId;

    #if defined(dRenderVariant_pick)
        out.vPickIds = vec2<f32>(in.aInstance, group);
    #endif

    #if defined(dNeedsMarker)
        ${assign_data_index_wgsl(undefined, 'markerIndex')}
        out.vMarker = read_tMarker(markerIndex).r;
    #endif

    #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
        #if defined(dColorType_texture)
            // the texture's index by the color type (\`ColorTypeIndex\`): instance, group, groupInstance, vertex,
            // vertexInstance
            var colorIndex = instanceIndex;
            switch uObject.uColorIndexType {
                case 1: { colorIndex = groupIndex; }
                case 2: { colorIndex = groupInstanceIndex; }
                case 3: { colorIndex = vertexId; }
                case 4: { colorIndex = vertexInstanceIndex; }
                default: {}
            }
            #if defined(dDualColor)
                // two colors per group, except for color mode 2 (cylinders)
                let dualColor = uObject.uColorIndexType == 1 || uObject.uColorIndexType == 2;
                var color2 = vec3<f32>(0.0);
                if (dualColor && in.aColorMode != 2.0) {
                    out.vColor = vec4<f32>(read_tColor(colorIndex * 2).rgb, 1.0);
                    color2 = read_tColor(colorIndex * 2 + 1).rgb;
                } else {
                    out.vColor = vec4<f32>(read_tColor(colorIndex).rgb, 1.0);
                }
            #else
                out.vColor = vec4<f32>(read_tColor(colorIndex).rgb, 1.0);
            #endif
        #elif defined(dColorType_grid)
            // the local position (volume) or the model position (volumeInstance), see \`ColorTypeIndex\`
            let colorGridPosition = select(${GridPosition}, ${GridModelPosition}, uObject.uColorIndexType == 1);
            out.vColor = vec4<f32>(${grid_value_wgsl('tColorGrid', 'colorGridPosition')}.rgb, 1.0);
        #endif

        #if defined(dUsePalette)
            // the color data hold palette indices, see \`ColorTheme.Palette\`
            out.vPaletteV = ((out.vColor.r * 256.0 * 256.0 * 255.0 + out.vColor.g * 256.0 * 255.0 + out.vColor.b * 255.0) - 1.0) / 16777214.0;
        #endif

        #if defined(dOverpaint)
            #if defined(dOverpaintType_volumeInstance)
                let overpaint = ${grid_value_wgsl('tOverpaintGrid', GridModelPosition)};
            #else
                ${assign_data_index_wgsl('dOverpaintType', 'overpaintIndex')}
                let overpaint = read_tOverpaint(overpaintIndex);
            #endif
            // pre-mix to avoid darkening due to empty overpaint
            #if defined(dColorType_uniform)
                let baseColor = uMaterial.uColor;
            #else
                let baseColor = out.vColor.rgb;
            #endif
            out.vOverpaint = vec4<f32>(mix(baseColor, overpaint.rgb, overpaint.a), overpaint.a) * uMaterial.uOverpaintStrength;
        #endif

        #if defined(dEmissive)
            #if defined(dEmissiveType_volumeInstance)
                out.vEmissive = ${grid_value_wgsl('tEmissiveGrid', GridModelPosition)}.r * uMaterial.uEmissiveStrength;
            #else
                ${assign_data_index_wgsl('dEmissiveType', 'emissiveIndex')}
                out.vEmissive = read_tEmissive(emissiveIndex).r * uMaterial.uEmissiveStrength;
            #endif
        #endif

        #if defined(dSubstance)
            #if defined(dSubstanceType_volumeInstance)
                let substance = ${grid_value_wgsl('tSubstanceGrid', GridModelPosition)};
            #else
                ${assign_data_index_wgsl('dSubstanceType', 'substanceIndex')}
                let substance = read_tSubstance(substanceIndex);
            #endif
            // pre-mix to avoid artifacts due to empty substance
            let baseSubstance = vec3<f32>(uMaterial.uMetalness, uMaterial.uRoughness, uMaterial.uBumpiness);
            out.vSubstance = vec4<f32>(mix(baseSubstance, substance.rgb, substance.a), substance.a) * uMaterial.uSubstanceStrength;
        #endif
    #endif

    #if defined(dTransparency)
        #if defined(dTransparencyType_volumeInstance)
            out.vTransparency = ${grid_value_wgsl('tTransparencyGrid', GridModelPosition)}.r * uMaterial.uTransparencyStrength;
        #else
            ${assign_data_index_wgsl('dTransparencyType', 'transparencyIndex')}
            out.vTransparency = read_tTransparency(transparencyIndex).r * uMaterial.uTransparencyStrength;
        #endif
    #endif

    #if defined(dClipObjects) && defined(dClipping)
        ${assign_data_index_wgsl(undefined, 'clippingIndex')}
        out.vClipping = read_tClipping(clippingIndex).r;
    #endif
    ${clip_flag_wgsl('out.vClipping')}
`;

/** \`size\` from the size data; expects the indices of \`assign_varyings_wgsl\` */
export const assign_size_wgsl = /* wgsl */ `
    #if defined(dSizeType_uniform)
        var size = uMaterial.uSize;
    #else
        // the texture's index by the size type (\`SizeTypeIndex\`): instance, group, groupInstance, vertex,
        // vertexInstance
        var sizeIndex = instanceIndex;
        switch uObject.uSizeIndexType {
            case 1: { sizeIndex = groupIndex; }
            case 2: { sizeIndex = groupInstanceIndex; }
            case 3: { sizeIndex = vertexId; }
            case 4: { sizeIndex = vertexInstanceIndex; }
            default: {}
        }
        // sizes are stored times 100
        var size = unpackRGBToInt(read_tSize(sizeIndex).rgb) / 100.0;
    #endif
    size *= uMaterial.uSizeFactor;
`;
