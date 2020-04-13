export default `
#ifdef dGeoTexture
    // aGroup is used as a vertex index here and the group id is retirieved from tPositionGroup
    float group = readFromTexture(tPositionGroup, aGroup, uGeoTexDim).w;
#else
    float group = aGroup;
#endif
`;