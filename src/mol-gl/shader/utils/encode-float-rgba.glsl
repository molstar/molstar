vec4 encodeFloatRGBA(in float value) {
    value = clamp(value, 0., 16777216.);
    vec3 c = vec3(0.);
    c.b = mod(value, 256.);
    value = floor(value/256.);
    c.g = mod(value, 256.);
    value = floor(value/256.);
    c.r = mod(value, 256.);
    return vec4(c/255., 1.);
}

#pragma glslify: export(encodeFloatRGBA)