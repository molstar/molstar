/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export default `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform vec2 uQuadShift;
uniform vec3 uDimensions;
uniform vec3 uMin;
uniform vec3 uDelta;

uniform sampler2D tCenters;
uniform sampler2D tInfo;
uniform sampler2D tCoeff;
uniform sampler2D tAlpha;

uniform int uNCenters;
uniform int uNCoeff;
uniform int uNAlpha;

uniform int uLittleEndian;

float shiftRight (float v, float amt) {
  v = floor(v) + 0.5;
  return floor(v / exp2(amt));
}
float shiftLeft (float v, float amt) {
    return floor(v * exp2(amt) + 0.5);
}
float maskLast (float v, float bits) {
    return mod(v, shiftLeft(1.0, bits));
}
float extractBits (float num, float from, float to) {
    from = floor(from + 0.5); to = floor(to + 0.5);
    return maskLast(shiftRight(num, from), to - from);
}
// Adapted from https://github.com/equinor/glsl-float-to-rgba
// MIT License, Copyright (c) 2020 Equinor
vec4 floatToRgba(float texelFloat) {
    if (texelFloat == 0.0) return vec4(0, 0, 0, 0);
    float sign = texelFloat > 0.0 ? 0.0 : 1.0;
    texelFloat = abs(texelFloat);
    float exponent = floor(log2(texelFloat));
    float biased_exponent = exponent + 127.0;
    float fraction = ((texelFloat / exp2(exponent)) - 1.0) * 8388608.0;
    float t = biased_exponent / 2.0;
    float last_bit_of_biased_exponent = fract(t) * 2.0;
    float remaining_bits_of_biased_exponent = floor(t);
    float byte4 = extractBits(fraction, 0.0, 8.0) / 255.0;
    float byte3 = extractBits(fraction, 8.0, 16.0) / 255.0;
    float byte2 = (last_bit_of_biased_exponent * 128.0 + extractBits(fraction, 16.0, 23.0)) / 255.0;
    float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0;
    return (
      uLittleEndian > 0
      ? vec4(byte4, byte3, byte2, byte1)
      : vec4(byte1, byte2, byte3, byte4)
    );
}

float L1(vec3 p, float a0, float a1, float a2) {
    return a0 * p.z + a1 * p.x + a2 * p.y;
}

float L2(vec3 p, float a0, float a1, float a2, float a3, float a4) {
    float x = p.x, y = p.y, z = p.z;
    float xx = x * x, yy = y * y, zz = z * z;
    return (
        a0 * (-0.5 * xx - 0.5 * yy + zz) +
        a1 * (1.7320508075688772 * x * z) +
        a2 * (1.7320508075688772 * y * z) +
        a3 * (0.8660254037844386 * xx - 0.8660254037844386 * yy) +
        a4 * (1.7320508075688772 * x * y)
    );
}

float L3(vec3 p, float a0, float a1, float a2, float a3, float a4, float a5, float a6) {
    float x = p.x, y = p.y, z = p.z;
    float xx = x * x, yy = y * y, zz = z * z;
    float xxx = xx * x, yyy = yy * y, zzz = zz * z;
    return (
        a0 * (-1.5 * xx * z - 1.5 * yy * z + zzz) +
        a1 * (-0.6123724356957945 * xxx - 0.6123724356957945 * x * yy + 2.449489742783178 * x * zz) +
        a2 * (-0.6123724356957945 * xx * y - 0.6123724356957945 * yyy + 2.449489742783178 * y * zz) +
        a3 * (1.9364916731037085 * xx * z - 1.9364916731037085 * yy * z) +
        a4 * (3.872983346207417 * x * y * z) +
        a5 * (0.7905694150420949 * xxx - 2.3717082451262845 * x * yy) +
        a6 * (2.3717082451262845 * xx * y - 0.7905694150420949 * yyy)
    );
}

void main(void) {
    // TODO: use "integer division" to avoid rounding errors for non power of 2 dimensions
    vec3 idx = vec3(gl_FragCoord.x - 0.5, mod(gl_FragCoord.y - 0.5, uDimensions.z), floor((gl_FragCoord.y - 0.5) / uDimensions.y));
    float offset = idx.x + idx.y * uDimensions.x + idx.z * uDimensions.x * uDimensions.y;

    // TODO: is there render target ordering that does not require the remapping of the coords?
    float k = mod(offset, uDimensions.z), kk = floor(offset / uDimensions.z);
    float j = mod(kk, uDimensions.y);
    float i = floor(kk / uDimensions.y);
    
    vec3 xyz = uMin + uDelta * vec3(i, j, k);

    float fCenter = 1.0 / float(uNCenters - 1);
    float fCoeff = 1.0 / float(uNCoeff - 1);
    float fAlpha = 1.0 / float(uNAlpha - 1);

    float v = 0.0;

    for (int i = 0; i < uNCenters; i++) {
        vec2 cCoord = vec2(float(i) * fCenter, 0.5);

        vec4 center = texture2D(tCenters, cCoord);
        vec3 X = xyz - center.xyz;
        float R2 = dot(X, X);

        // center.w is squared cutoff radius
        if (R2 > center.w) {
            continue;
        }

        vec4 info = texture2D(tInfo, cCoord);

        int L = int(info.x);
        float alphaOffset = info.y;
        int coeffStart = int(info.z);
        int coeffEnd = int(info.w);


        float gauss = 0.0;
        for (int j = coeffStart; j < coeffEnd; j++) {
            vec2 c = texture2D(tCoeff, vec2(float(j) * fCoeff, 0.5)).xy;
            gauss += c.x * exp(-c.y * R2);
        }

        float spherical = 0.0;
        if (L == 0) {
            spherical = texture2D(tAlpha, vec2(alphaOffset * fAlpha, 0.5)).x;
        } else if (L == 1) {
            spherical = L1(
                X,
                texture2D(tAlpha, vec2(alphaOffset * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 1.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 2.0) * fAlpha, 0.5)).x
            );
        } else if (L == 2) {
            spherical = L2(
                X,
                texture2D(tAlpha, vec2(alphaOffset * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 1.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 2.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 3.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 4.0) * fAlpha, 0.5)).x
            );
        } else if (L == 3) {
            spherical = L3(
                X,
                texture2D(tAlpha, vec2(alphaOffset * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 1.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 2.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 3.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 4.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 5.0) * fAlpha, 0.5)).x,
                texture2D(tAlpha, vec2((alphaOffset + 6.0) * fAlpha, 0.5)).x
            );
        } 
        // TODO: do we need L > 3?

        v += gauss * spherical;
    }

    // TODO: render to single component float32 texture in WebGL2
    gl_FragColor = floatToRgba(v);
}
`;