#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_attribute)
    varying vec3 vColor;
    attribute vec3 aColor;
#elif defined(dColorType_instance) || defined(dColorType_element) || defined(dColorType_elementInstance)
    varying vec3 vColor;
    uniform vec2 uColorTexSize;
    uniform sampler2D tColor;
#endif

#pragma glslify: read_vec3 = require(../utils/read-from-texture.glsl)