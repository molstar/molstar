/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const assign_size = `
#if defined(dSizeType_uniform)
    float size = uSize;
#elif defined(dSizeType_texture)
    // the texture's index by the size type (SizeTypeIndex): instance, group, groupInstance, vertex, vertexInstance,
    // siteType (quantized centers' types, siteSizeIndex of the quad chunk)
    int sizeIndex = int(aInstance);
    if (uSizeIndexType == 1) {
        sizeIndex = groupIndex(group);
    } else if (uSizeIndexType == 2) {
        sizeIndex = groupInstanceIndex(aInstance, group);
    } else if (uSizeIndexType == 3) {
        sizeIndex = vertexId;
    } else if (uSizeIndexType == 4) {
        sizeIndex = int(aInstance) * uVertexCount + vertexId;
    }
    #if defined(dSites)
        if (uSizeIndexType == 5) sizeIndex = siteSizeIndex;
    #endif
    float size = unpackRGBToInt(read_tSize(sizeIndex).rgb);
    size /= 100.0; // NOTE factor also set in TypeScript
#endif

size *= uSizeFactor;
`;
