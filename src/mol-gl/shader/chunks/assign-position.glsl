mat4 modelView = uView * uModel * aTransform;
vec4 mvPosition = modelView * vec4(aPosition, 1.0);
vViewPosition = mvPosition.xyz;
gl_Position = uProjection * mvPosition;