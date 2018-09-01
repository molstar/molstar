uniform mat4 uProjection, uModel, uView;

uniform int uObjectId;
uniform int uInstanceCount;
uniform int uGroupCount;

uniform vec2 uMarkerTexDim;
uniform sampler2D tMarker;
varying float vMarker;

varying vec3 vViewPosition;

#pragma glslify: readFromTexture = require(../utils/read-from-texture.glsl)