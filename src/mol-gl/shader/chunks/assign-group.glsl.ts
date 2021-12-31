export const assign_group = `
#ifdef dGeometryType_textureMesh
    float group = unpackRGBToInt(readFromTexture(tGroup, VertexID, uGeoTexDim).rgb);
#else
    float group = aGroup;
#endif
`;