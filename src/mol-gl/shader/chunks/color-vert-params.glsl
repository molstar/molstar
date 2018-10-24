#if defined(dColorType_uniform)
    uniform vec3 uColor;
#elif defined(dColorType_attribute)
    varying vec4 vColor;
    attribute vec3 aColor;
#elif defined(dColorType_instance) || defined(dColorType_group) || defined(dColorType_groupInstance)
    varying vec4 vColor;
    uniform vec2 uColorTexDim;
    uniform sampler2D tColor;
#elif defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
    varying vec4 vColor;
    #pragma glslify: encodeIdRGB = require(../utils/encode-id-rgb.glsl)
#endif