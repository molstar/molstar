export const assign_size = `
#if defined(dSizeType_uniform)
    float size = uSize;
#elif defined(dSizeType_texture)
    // the texture's index by the size type (SizeTypeIndex): instance, group, groupInstance, vertex, vertexInstance
    float sizeIndex = aInstance;
    if (uSizeIndexType == 1) {
        sizeIndex = groupIndex(group);
    } else if (uSizeIndexType == 2) {
        sizeIndex = groupInstanceIndex(aInstance, group);
    } else if (uSizeIndexType == 3) {
        sizeIndex = float(vertexId);
    } else if (uSizeIndexType == 4) {
        sizeIndex = float(int(aInstance) * uVertexCount + vertexId);
    }
    float size = unpackRGBToInt(readFromTexture(tSize, sizeIndex, uSizeTexDim).rgb);
    size /= 100.0; // NOTE factor also set in TypeScript
#endif

size *= uSizeFactor;
`;
