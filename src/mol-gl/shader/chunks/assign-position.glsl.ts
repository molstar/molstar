export default `
mat4 model = uModel * aTransform;
mat4 modelView = uView * model;
#ifdef dGeoTexture
    vec3 position = readFromTexture(tPositionGroup, aGroup, uGeoTexDim).xyz;
#else
    vec3 position = aPosition;
#endif
vec4 position4 = vec4(position, 1.0);
vModelPosition = (model * position4).xyz; // for clipping in frag shader
vec4 mvPosition = modelView * position4;
vViewPosition = mvPosition.xyz;
gl_Position = uProjection * mvPosition;
`;