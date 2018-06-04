uniform mat4 uProjection, uModel, uView;

uniform int uObjectId;
uniform int uInstanceCount;
uniform int uElementCount;

uniform vec2 uFlagTexSize;
uniform sampler2D tFlag;
varying float vFlag;
#pragma glslify: readFromTexture = require(../utils/read-from-texture.glsl)