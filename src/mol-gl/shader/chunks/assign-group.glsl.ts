export default `
#ifdef dGeoTexture
    float group = decodeFloatRGB(readFromTexture(tGroup, VertexID, uGeoTexDim).rgb);
#else
    float group = aGroup;
#endif
`;