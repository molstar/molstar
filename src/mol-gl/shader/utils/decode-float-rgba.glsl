float decodeFloatRGBA(const in vec4 rgba) {
    return dot(rgba, vec4(1.0, 1/255.0, 1/65025.0, 1/16581375.0));
}

#pragma glslify: export(decodeFloatRGBA)