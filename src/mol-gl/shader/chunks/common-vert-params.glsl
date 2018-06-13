uniform mat4 uProjection, uModel, uView;

uniform int uObjectId;
uniform int uInstanceCount;
uniform int uElementCount;

uniform vec2 uMarkerTexSize;
uniform sampler2D tMarker;
varying float vMarker;
#pragma glslify: readFromTexture = require(../utils/read-from-texture.glsl)