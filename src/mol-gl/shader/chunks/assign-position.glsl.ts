export default `
mat4 modelView = uView * uModel * aTransform;
#ifdef dGeoTexture
    vec3 position = readFromTexture(tPositionGroup, aGroup, uGeoTexDim).xyz;
#else
    vec3 position = aPosition;
#endif
vec4 mvPosition = modelView * vec4(position, 1.0);
vViewPosition = mvPosition.xyz;
gl_Position = uProjection * mvPosition;
`