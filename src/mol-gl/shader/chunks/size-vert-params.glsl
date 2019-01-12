#if defined(dSizeType_uniform)
    uniform float uSize;
#elif defined(dSizeType_attribute)
    attribute float aSize;
#elif defined(dSizeType_instance) || defined(dSizeType_group) || defined(dSizeType_groupInstance)
    uniform vec2 uSizeTexDim;
    uniform sampler2D tSize;
    #pragma glslify: decodeFloatLog = require(../utils/decode-float-log.glsl)
#endif