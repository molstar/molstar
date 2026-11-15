/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Items drawn as a quad of two triangles (6 vertices, no elements) with the item data pulled
 * from `tPositionGroup`: sets `mapping` (the quad corner, -1 or 1 per axis), `vertexId` (the item: quad `p` of a level
 * of detail draws item `p * stride`, the items in their order; merged, from the member's first item, `aSegmentItem`),
 * `position` and `group`; quantized centers (`dSites`) from `tSite` by the object's origin and scale (merged, the member's,
 * `aSegmentSites`), the group the item's (the member's), `siteSizeIndex` the index of its type in the size table
 * (`siteType` sizes; merged, from the member's table, `aSegmentItem.z`)
 */
export const assign_quad_position_group = `
#if defined(dInstanceData)
    readInstanceData();
#endif
vec2 mapping = vec2(1.0, 1.0); // vertices 2 and 5
#if __VERSION__ == 100
    int m = imod(VertexID, 6);
#else
    int m = VertexID % 6;
#endif
if (m == 0) {
    mapping = vec2(-1.0, 1.0);
} else if (m == 1 || m == 3) {
    mapping = vec2(-1.0, -1.0);
} else if (m == 4) {
    mapping = vec2(1.0, -1.0);
}

int vertexId = VertexID / 6;
#if defined(dSegmented)
    // the member's stride at the level drawn (1 without levels, or drawn whole) from its first item
    vec4 segmentStrides = aSegmentStride;
    float segmentStride = uLodLevel == 0 ? segmentStrides.x : (uLodLevel == 1 ? segmentStrides.y : (uLodLevel == 2 ? segmentStrides.z : segmentStrides.w));
    vec4 segmentItem = aSegmentItem;
    int itemBase = int(segmentItem.x) * 65536 + int(segmentItem.y);
    vertexId = itemBase + (vertexId - itemBase) * int(segmentStride);
#else
    // at least 1: a draw before its uniforms are written draws every item
    vertexId *= int(max(1.0, uLodStride));
#endif

#if defined(dSites)
    #if defined(dSegmented)
        vec4 siteOriginScale = aSegmentSites;
        GroupId group = GroupId(vertexId - itemBase);
    #else
        vec4 siteOriginScale = vec4(uSiteOrigin, uSiteScale);
        GroupId group = GroupId(vertexId);
    #endif
    vec4 site = read_tSite(vertexId);
    vec3 position = siteOriginScale.xyz + site.xyz * siteOriginScale.w;
    #if defined(dSegmented)
        int siteSizeIndex = int(aSegmentItem.z) + int(site.w);
    #else
        int siteSizeIndex = int(site.w);
    #endif
#else
    vec3 position = read_tPositionGroup(vertexId).rgb;
    // an integer (\`rgb32f-i32\`; WebGL1 a float, exact to 2^24)
    GroupId group = readGroup_tPositionGroup(vertexId);
#endif
`;
