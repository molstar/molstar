/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clip_flag_wgsl } from './clip.wgsl.js';
import { assign_data_index_wgsl } from './read-from-texture.wgsl.js';
import { CopyMarkerBlock, PagedCopyShift } from '../../../culled-data.js';

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
    @location(2) @interpolate(flat) vPickIds: vec2<u32>,
#endif
// scalars share locations (each counts against the 16 inter-stage variables, whatever its size): the palette
// coordinate in the color's alpha (unread, materials take \`uAlpha\`), the marker with the clipping flags (both flat by
// group), the transparency with the emissive
#if defined(dColorType_varying)
    @location(3) vColor: vec4<f32>,
#endif
#if defined(dNeedsMarker) || (defined(dClipObjects) && defined(dClipping))
    #if !defined(dVaryingGroup)
        @location(4) vMarkerClipping: vec2<f32>,
    #else
        @location(4) @interpolate(flat) vMarkerClipping: vec2<f32>,
    #endif
#endif
#if defined(dOverpaint)
    @location(5) vOverpaint: vec4<f32>,
#endif
#if defined(dTransparency) || defined(dEmissive)
    @location(6) vTransparencyEmissive: vec2<f32>,
#endif
#if defined(dSubstance)
    @location(8) vSubstance: vec4<f32>,
#endif
${extra}
};
`;
}

/**
 * Item `vertexId` of the quad from `vertex_index` (6 vertices, two triangles per item) with its `position` and
 * `group` from `tPositionGroup`, and the quad corner `mapping` (-1 or 1 per axis): quad `p` of a level of detail draws
 * item `p * stride` (the items in their order; merged, from the member's first item, `segmentItem`); quantized centers
 * (`dSites`) from `tSite` by the object's origin and scale (merged, the member's, `segmentSites`), the group the item's
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
    var vertexId = i32(in.vertexIndex) / 6;
    #if defined(dSegmented)
        // the member's stride at the level drawn (1 without levels, or drawn whole) from its first item
        let itemBase = segmentItem(in);
        vertexId = itemBase + (vertexId - itemBase) * i32(segmentStride(in)[clamp(uDraw.uLodLevel, 0, 3)]);
    #else
        // at least 1: a draw before its uniforms are written draws every item
        vertexId *= max(1, i32(uDraw.uLodStride));
    #endif
    #if defined(dSites)
        #if defined(dSegmented)
            let siteOriginScale = segmentSites(in);
            let group = vertexId - itemBase;
        #else
            let siteOriginScale = vec4<f32>(uObject.uSiteOrigin, uObject.uSiteScale);
            let group = vertexId;
        #endif
        let site = read_tSite(vertexId);
        var position = siteOriginScale.xyz + site.xyz * siteOriginScale.w;
        // the index of the site's type in the size table (\`siteType\` sizes; merged, from the member's table)
        #if defined(dSegmented)
            let siteSizeIndex = segmentSizeBase(in) + i32(site.w);
        #else
            let siteSizeIndex = i32(site.w);
        #endif
    #else
        var position = read_tPositionGroup(vertexId).xyz;
        // an integer (\`rgb32f-i32\`)
        let group = readGroup_tPositionGroup(vertexId);
    #endif
`;

/**
 * Levels of detail of the items drawn as quads (spheres, points), after `out.vModelPosition`: `lodFactor` the level's
 * scale (a coarser level's fewer items larger) faded at the level's range (1 without levels), `lodCulled` beyond the
 * range
 */
export const assign_lod_factor_wgsl = /* wgsl */ `
    var lodFactor = 1.0;
    var lodCulled = false;
    let lod = uDraw.uLod;
    if (lod.w != 0.0 && (lod.x != 0.0 || lod.y != 0.0)) {
        #if defined(dSegmented)
            // members of a merged item can have different scales for the same lod level (item-count dependent)
            let lodScale = segmentLod(in)[clamp(uDraw.uLodLevel, 0, 3)];
        #else
            let lodScale = lod.w;
        #endif
        // the distance unscaled, as the instances are culled (any camera scale)
        let d = (dot(uGlobal.uCameraPlane.xyz, out.vModelPosition) + uGlobal.uCameraPlane.w) / uGlobal.uModelScale;
        let level = lodAt(d, lod.xyz, lodScale);
        lodFactor = level.factor;
        lodCulled = level.culled;
    }
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
const GridModelPosition = 'gridModelPosition';

/**
 * Data texture indices and the color, marker and transparency varyings; expects `in.aInstance`, `vertexId`,
 * `group` and `out`, and for grid data (mesh programs) `gridPosition` and `gridModelPosition` (where the instance is
 * placed, without the transient object transform)
 */
export const assign_varyings_wgsl = /* wgsl */ `
    let instanceIndex = i32(in.aInstance);
    #if defined(dCopyLayers)
        // a page store's copy's instance: its copy plus one (none 0) and its template particle's store index, by which
        // its copy's own layers replace its store index's (\`copyLayerOr\`)
        #if defined(dSphereResolve)
            let layerWords = rasterPagedWords(in.instanceIndex);
        #else
            let layerWords = culledPagedWords(in.instanceIndex);
        #endif
        let layerTag = layerWords.w >> ${PagedCopyShift}u;
        let layerIndex = layerWords.z;
    #endif
    #if defined(dSegmented)
        // the data of merged render items, offset per segment (see \`SegmentSchema\`)
        let segment = segmentGroups(in);
        // integers: the instance ids and offsets past 2^24 exact
        let groupInstanceIndex = segment.y + i32(in.aInstance) * segment.x + i32(group);
        let groupIndex = segment.z + i32(group);
    #else
        let groupInstanceIndex = i32(in.aInstance) * i32(uObject.uGroupCount) + i32(group);
        let groupIndex = i32(group);
    #endif
    #if defined(dSegmented)
        // by the segment's offset (a read of its data) where read: the layers of that type; a color of it reads its own
        // (\`uColorIndexType\` 4); sizes of it do not merge (\`canMergeValues\`)
        #if defined(dOverpaintType_vertexInstance) || defined(dTransparencyType_vertexInstance) || defined(dEmissiveType_vertexInstance) || defined(dSubstanceType_vertexInstance)
            let vertexInstanceIndex = segmentVertexInstanceIndex(in, instanceIndex, vertexId);
        #else
            let vertexInstanceIndex = 0;
        #endif
    #else
        let vertexInstanceIndex = instanceIndex * uObject.uVertexCount + vertexId;
    #endif

    #if defined(dRenderVariant_pick)
        #if defined(dCulledPaged)
            // a page store's: its composed pose's slot (a copy's instances' the frame's), its record's identity (store
            // index, copy) gathered before the pick is read (\`paged-pick.ts\`)
            out.vPickIds = vec2<u32>(in.instanceIndex, u32(group));
        #else
            out.vPickIds = vec2<u32>(in.aInstance, u32(group));
        #endif
    #endif

    #if defined(dNeedsMarker)
        ${assign_data_index_wgsl(undefined, 'markerIndex')}
        out.vMarkerClipping.x = read_tMarker(markerIndex).r;
        #if defined(dCopyMarkers)
            // a copy's instance (a page store's template particle, its copy plus one above its pose's segment): its
            // copy's marker (\`sCopyMarkers\`: two words a copy, its marker or its block's offset by store index) or-ed
            // with its store index's (its template's: every copy's, a selection by entity)
            #if defined(dSphereResolve)
                let copyWords = rasterPagedWords(in.instanceIndex);
            #else
                let copyWords = culledPagedWords(in.instanceIndex);
            #endif
            let copyTag = copyWords.w >> ${PagedCopyShift}u;
            if (copyTag != 0u) {
                let copyWord = sCopyMarkers[2u * (copyTag - 1u)];
                var copyMarker = copyWord & 0xffu;
                if ((copyWord & ${CopyMarkerBlock}u) != 0u) {
                    let b = u32(bitcast<i32>(sCopyMarkers[2u * (copyTag - 1u) + 1u]) + i32(copyWords.z));
                    copyMarker = (sCopyMarkers[b / 4u] >> (8u * (b % 4u))) & 0xffu;
                }
                let templateMarker = u32(out.vMarkerClipping.x * 255.0 + 0.5);
                out.vMarkerClipping.x = f32(copyMarker | templateMarker) / 255.0;
            }
        #endif
    #endif

    #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
        #if defined(dColorType_texture)
            // the texture's index by the color type (\`ColorTypeIndex\`): instance, group, groupInstance, vertex,
            // vertexInstance, object (the object's texel; merged, its segment's)
            var colorIndex = instanceIndex;
            switch uObject.uColorIndexType {
                case 1: { colorIndex = groupIndex; }
                case 2: { colorIndex = groupInstanceIndex; }
                case 3: { colorIndex = vertexId; }
                #if defined(dSegmented)
                    case 4: { colorIndex = segmentVertexInstanceIndex(in, instanceIndex, vertexId); }
                #else
                    case 4: { colorIndex = vertexInstanceIndex; }
                #endif
                #if defined(dSegmented)
                    case 5: { colorIndex = i32(in.aSegment.w); }
                #else
                    case 5: { colorIndex = 0; }
                #endif
                default: {}
            }
            #if defined(dDualColor)
                // two colors per group, except for color mode 2 (cylinders)
                let dualColor = uObject.uColorIndexType == 1 || uObject.uColorIndexType == 2;
                var color2 = vec3<f32>(0.0);
                if (dualColor && colorMode != 2.0) {
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
            out.vColor.w = ((out.vColor.r * 256.0 * 256.0 * 255.0 + out.vColor.g * 256.0 * 255.0 + out.vColor.b * 255.0) - 1.0) / 16777214.0;
        #endif

        #if defined(dOverpaint)
            #if defined(dOverpaintType_volumeInstance)
                var overpaint = ${grid_value_wgsl('tOverpaintGrid', GridModelPosition)};
            #else
                ${assign_data_index_wgsl('dOverpaintType', 'overpaintIndex')}
                var overpaint = read_tOverpaint(overpaintIndex);
            #endif
            #if defined(dCopyLayers)
                overpaint = copyLayerOr(0u, 4u, layerTag, layerIndex, overpaint);
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
                out.vTransparencyEmissive.y = ${grid_value_wgsl('tEmissiveGrid', GridModelPosition)}.r * uMaterial.uEmissiveStrength;
            #else
                ${assign_data_index_wgsl('dEmissiveType', 'emissiveIndex')}
                var emissive = read_tEmissive(emissiveIndex).r;
                #if defined(dCopyLayers)
                    emissive = copyLayerOr(3u, 1u, layerTag, layerIndex, vec4<f32>(emissive)).r;
                #endif
                out.vTransparencyEmissive.y = emissive * uMaterial.uEmissiveStrength;
            #endif
        #endif

        #if defined(dSubstance)
            #if defined(dSubstanceType_volumeInstance)
                var substance = ${grid_value_wgsl('tSubstanceGrid', GridModelPosition)};
            #else
                ${assign_data_index_wgsl('dSubstanceType', 'substanceIndex')}
                var substance = read_tSubstance(substanceIndex);
            #endif
            #if defined(dCopyLayers)
                substance = copyLayerOr(4u, 4u, layerTag, layerIndex, substance);
            #endif
            // pre-mix to avoid artifacts due to empty substance
            let baseSubstance = vec3<f32>(uMaterial.uMetalness, uMaterial.uRoughness, uMaterial.uBumpiness);
            out.vSubstance = vec4<f32>(mix(baseSubstance, substance.rgb, substance.a), substance.a) * uMaterial.uSubstanceStrength;
        #endif
    #endif

    #if defined(dTransparency)
        #if defined(dTransparencyType_volumeInstance)
            out.vTransparencyEmissive.x = ${grid_value_wgsl('tTransparencyGrid', GridModelPosition)}.r * uMaterial.uTransparencyStrength;
        #else
            ${assign_data_index_wgsl('dTransparencyType', 'transparencyIndex')}
            var transparency = read_tTransparency(transparencyIndex).r;
            #if defined(dCopyLayers)
                transparency = copyLayerOr(1u, 1u, layerTag, layerIndex, vec4<f32>(transparency)).r;
            #endif
            out.vTransparencyEmissive.x = transparency * uMaterial.uTransparencyStrength;
        #endif
    #endif

    #if defined(dClipObjects) && defined(dClipping)
        ${assign_data_index_wgsl(undefined, 'clippingIndex')}
        out.vMarkerClipping.y = read_tClipping(clippingIndex).r;
        #if defined(dCopyLayers)
            out.vMarkerClipping.y = copyLayerOr(2u, 1u, layerTag, layerIndex, vec4<f32>(out.vMarkerClipping.y)).r;
        #endif
    #endif
    ${clip_flag_wgsl('out.vMarkerClipping.y')}
`;

/** \`size\` from the size data; expects the indices of \`assign_varyings_wgsl\` */
export const assign_size_wgsl = /* wgsl */ `
    #if defined(dSizeType_uniform)
        var size = uMaterial.uSize;
    #else
        var sizeDataIndex = sizeIndex(
            uObject.uSizeIndexType, instanceIndex, groupIndex, groupInstanceIndex, vertexId, vertexInstanceIndex
        );
        #if defined(dSites)
            if (uObject.uSizeIndexType == 5) { sizeDataIndex = siteSizeIndex; }
        #endif
        var size = unpackSize(read_tSize(sizeDataIndex));
    #endif
    size *= uMaterial.uSizeFactor;
`;
