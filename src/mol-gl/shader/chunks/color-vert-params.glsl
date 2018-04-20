#if defined(UNIFORM_COLOR)
    uniform vec3 color;
#elif defined(ATTRIBUTE_COLOR)
    varying vec3 vColor;
    attribute vec3 color;
#elif defined(INSTANCE_COLOR) || defined(ELEMENT_COLOR) || defined(ELEMENT_INSTANCE_COLOR)
    varying vec3 vColor;
    uniform vec2 colorTexSize;
    uniform sampler2D colorTex;
#endif

#pragma glslify: read_vec3 = require(../utils/read-vec3.glsl)