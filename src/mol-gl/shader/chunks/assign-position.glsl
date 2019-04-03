mat4 modelView = uView * uModel * aTransform;
#ifdef dPositionTexture
    vec3 position = readFromTexture(tPosition, aIndex, uPositionTexDim).xyz;
#else
    vec3 position = aPosition;
#endif
vec4 mvPosition = modelView * vec4(position, 1.0);
vViewPosition = mvPosition.xyz;
gl_Position = uProjection * mvPosition;