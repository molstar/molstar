export const assign_size = `
#if defined(dSizeType_uniform)
    float size = uSize;
#elif defined(dSizeType_attribute)
    float size = aSize;
#elif defined(dSizeType_instance)
    float size = unpackRGBToInt(readFromTexture(tSize, aInstance, uSizeTexDim).rgb);
#elif defined(dSizeType_group)
    float size = unpackRGBToInt(readFromTexture(tSize, group, uSizeTexDim).rgb);
#elif defined(dSizeType_groupInstance)
    float size = unpackRGBToInt(readFromTexture(tSize, aInstance * float(uGroupCount) + group, uSizeTexDim).rgb);
#elif defined(dSizeType_vertex)
    float size = unpackRGBToInt(readFromTexture(tSize, vertexId, uSizeTexDim).rgb);
#elif defined(dSizeType_vertexInstance)
    float size = unpackRGBToInt(readFromTexture(tSize, int(aInstance) * uVertexCount + vertexId, uSizeTexDim).rgb);
#endif

#if defined(dSizeType_instance) || defined(dSizeType_group) || defined(dSizeType_groupInstance) || defined(dSizeType_vertex) || defined(dSizeType_vertexInstance)
    size /= 100.0; // NOTE factor also set in TypeScript
#endif

size *= uSizeFactor;
`;