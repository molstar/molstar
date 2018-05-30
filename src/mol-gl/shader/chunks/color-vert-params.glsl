#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_attribute)
    varying vec4 vColor;
    attribute vec3 aColor;
#elif defined(dColorType_instance) || defined(dColorType_element) || defined(dColorType_elementInstance)
    varying vec4 vColor;
    uniform vec2 uColorTexSize;
    uniform sampler2D tColor;
#elif defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_elementPicking)
    varying vec4 vColor;
    #pragma glslify: encodeIdRGBA = require(../utils/encode-id-rgba.glsl)
#endif

#pragma glslify: read_vec3 = require(../utils/read-from-texture.glsl)