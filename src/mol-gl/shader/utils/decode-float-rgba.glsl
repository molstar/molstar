float decodeFloatRGBA(const in vec4 rgba) {
    return dot(rgba, vec4(1.0, 1.0/255.0, 1.0/65025.0, 1.0/16581375.0));
}

#pragma glslify: export(decodeFloatRGBA)