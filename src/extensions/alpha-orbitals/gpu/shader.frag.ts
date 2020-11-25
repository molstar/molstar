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

uniform float uWidth;

#ifndef uNCenters
    uniform int uNCenters;
#endif

uniform int uNCoeff;
uniform int uNAlpha;

uniform bool uDensity;
uniform float uOccupancy;
uniform sampler2D tCumulativeSum;

uniform bool uLittleEndian;

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
        uLittleEndian
            ? vec4(byte4, byte3, byte2, byte1)
            : vec4(byte1, byte2, byte3, byte4)
    );
}

///////////////////////////////////////////////////////

ivec4 floatsToBytes(vec4 inputFloats) {
  ivec4 bytes = ivec4(inputFloats * 255.0);
  return (
    uLittleEndian
    ? bytes.abgr
    : bytes
  );
}

// Break the four bytes down into an array of 32 bits.
void bytesToBits(const in ivec4 bytes, out bool bits[32]) {
  for (int channelIndex = 0; channelIndex < 4; ++channelIndex) {
    float acc = float(bytes[channelIndex]);
    for (int indexInByte = 7; indexInByte >= 0; --indexInByte) {
      float powerOfTwo = exp2(float(indexInByte));
      bool bit = acc >= powerOfTwo;
      bits[channelIndex * 8 + (7 - indexInByte)] = bit;
      acc = mod(acc, powerOfTwo);
    }
  }
}

// Compute the exponent of the 32-bit float.
float getExponent(bool bits[32]) {
  const int startIndex = 1;
  const int bitStringLength = 8;
  const int endBeforeIndex = startIndex + bitStringLength;
  float acc = 0.0;
  int pow2 = bitStringLength - 1;
  for (int bitIndex = startIndex; bitIndex < endBeforeIndex; ++bitIndex) {
    acc += float(bits[bitIndex]) * exp2(float(pow2--));
  }
  return acc;
}

// Compute the mantissa of the 32-bit float.
float getMantissa(bool bits[32], bool subnormal) {
  const int startIndex = 9;
  const int bitStringLength = 23;
  const int endBeforeIndex = startIndex + bitStringLength;
  // Leading/implicit/hidden bit convention:
  // If the number is not subnormal (with exponent 0), we add a leading 1 digit.
  float acc = float(!subnormal) * exp2(float(bitStringLength));
  int pow2 = bitStringLength - 1;
  for (int bitIndex = startIndex; bitIndex < endBeforeIndex; ++bitIndex) {
    acc += float(bits[bitIndex]) * exp2(float(pow2--));
  }
  return acc;
}

// Parse the float from its 32 bits.
float bitsToFloat(bool bits[32]) {
  float signBit = float(bits[0]) * -2.0 + 1.0;
  float exponent = getExponent(bits);
  bool subnormal = abs(exponent - 0.0) < 0.01;
  float mantissa = getMantissa(bits, subnormal);
  float exponentBias = 127.0;
  return signBit * mantissa * exp2(exponent - exponentBias - 23.0);
}

// Decode a 32-bit float from the RGBA color channels of a texel.
float rgbaToFloat(vec4 texelRGBA) {
  ivec4 rgbaBytes = floatsToBytes(texelRGBA);
  bool bits[32];
  bytesToBits(rgbaBytes, bits);
  return bitsToFloat(bits);
}

///////////////////////////////////////////////////////

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

float L4(vec3 p, float a0, float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8) {
    float x = p.x, y = p.y, z = p.z;
    float xx = x * x, yy = y * y, zz = z * z;
    float xxx = xx * x, yyy = yy * y, zzz = zz * z;
    float xxxx = xxx * x, yyyy = yyy * y, zzzz = zzz * z;
    return (
        a0 * (0.375 * xxxx + 0.75 * xx * yy + 0.375 * yyyy - 3.0 * xx * zz - 3.0 * yy * zz + zzzz) +
        a1 * (-2.3717082451262845 * xxx * z - 2.3717082451262845 * x * yy * z + 3.1622776601683795 * x * zzz) +
        a2 * (-2.3717082451262845 * xx * y * z - 2.3717082451262845 * yyy * z + 3.1622776601683795 * y * zzz) +
        a3 * (-0.5590169943749475 * xxxx + 0.5590169943749475 * yyyy + 3.3541019662496847 * xx * zz - 3.3541019662496847 * yy * zz) +
        a4 * (-1.118033988749895 * xxx * y - 1.118033988749895 * x * yyy + 6.708203932499369 * x * y * zz) +
        a5 * (2.091650066335189 * xxx * z + -6.274950199005566 * x * yy * z) +
        a6 * (6.274950199005566 * xx * y * z + -2.091650066335189 * yyy * z) +
        a7 * (0.739509972887452 * xxxx - 4.437059837324712 * xx * yy + 0.739509972887452 * yyyy) +
        a8 * (2.958039891549808 * xxx * y + -2.958039891549808 * x * yyy)
    );
}

float alpha(float offset, float f) {
    #ifdef uMaxCoeffs
        // in webgl1, the value is in the alpha channel!
        return texture2D(tAlpha, vec2(offset * f, 0.5)).a;
    #endif

    #ifndef uMaxCoeffs
        return texture2D(tAlpha, vec2(offset * f, 0.5)).x;
    #endif
}

float Y(int L, vec3 X, float aO, float fA) {
    if (L == 0) {
        return alpha(aO, fA);
    } else if (L == 1) {
        return L1(X,
            alpha(aO, fA), alpha(aO + 1.0, fA), alpha(aO + 2.0, fA)
        );
    } else if (L == 2) {
        return L2(X,
            alpha(aO, fA), alpha(aO + 1.0, fA), alpha(aO + 2.0, fA), alpha(aO + 3.0, fA), alpha(aO + 4.0, fA)
        );
    } else if (L == 3) {
        return L3(X,
            alpha(aO, fA), alpha(aO + 1.0, fA), alpha(aO + 2.0, fA), alpha(aO + 3.0, fA), alpha(aO + 4.0, fA),
            alpha(aO + 5.0, fA), alpha(aO + 6.0, fA)
        );
    } else if (L == 4) {
        return L4(X,
            alpha(aO, fA), alpha(aO + 1.0, fA), alpha(aO + 2.0, fA), alpha(aO + 3.0, fA), alpha(aO + 4.0, fA),
            alpha(aO + 5.0, fA), alpha(aO + 6.0, fA), alpha(aO + 7.0, fA), alpha(aO + 8.0, fA)
        );
    }
    // TODO: do we need L > 4?
    return 0.0;
}

#ifndef uMaxCoeffs
    float R(float R2, int start, int end, float fCoeff) {
        float gauss = 0.0;
        for (int i = start; i < end; i++) {
            vec2 c = texture2D(tCoeff, vec2(float(i) * fCoeff, 0.5)).xy;
            gauss += c.x * exp(-c.y * R2);
        }
        return gauss;
    }
#endif

#ifdef uMaxCoeffs
    float R(float R2, int start, int end, float fCoeff) {
        float gauss = 0.0;
        int o = start;
        for (int i = 0; i < uMaxCoeffs; i++) {
            if (o >= end) break;

            vec2 c = texture2D(tCoeff, vec2(float(o) * fCoeff, 0.5)).xy;
            gauss += c.x * exp(-c.y * R2);
            o++;
        }
        return gauss;
    }
#endif

float intDiv(float a, float b) { return float(int(a) / int(b)); }
float intMod(float a, float b) { return a - b * float(int(a) / int(b)); }

void main(void) {
    float offset = floor(gl_FragCoord.x) + floor(gl_FragCoord.y) * uWidth;

    // axis order fast to slow Z, Y, X
    // TODO: support arbitrary axis orders?
    float k = intMod(offset, uDimensions.z), kk = intDiv(offset, uDimensions.z);
    float j = intMod(kk, uDimensions.y);
    float i = intDiv(kk, uDimensions.y);

    vec3 xyz = uMin + uDelta * vec3(i, j, k);

    float fCenter = 1.0 / float(uNCenters - 1);
    float fCoeff = 1.0 / float(uNCoeff - 1);
    float fA = 1.0 / float(uNAlpha - 1);

    // gl_FragColor = floatToRgba(offset);
    // return;

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
        float aO = info.y;
        int coeffStart = int(info.z);
        int coeffEnd = int(info.w);

        v += R(R2, coeffStart, coeffEnd, fCoeff) * Y(L, X, aO, fA);
    }

    
    if (uDensity) {
        float current = rgbaToFloat(texture2D(tCumulativeSum, gl_FragCoord.xy / vec2(uWidth, uWidth)));
        gl_FragColor = floatToRgba(current + uOccupancy * v * v);
    } else {
         gl_FragColor = floatToRgba(v);
    }
}
`;