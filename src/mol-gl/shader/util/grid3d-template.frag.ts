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
uniform bool uLittleEndian;
uniform float uWidth;

#ifdef CUMULATIVE
uniform sampler2D tCumulativeSum;
#endif

{UNIFORMS}

{UTILS}

//////////////////////////////////////////////////////////

// floatToRgba adapted from https://github.com/equinor/glsl-float-to-rgba
// MIT License, Copyright (c) 2020 Equinor

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

// rgbaToFloat adapted from https://github.com/ihmeuw/glsl-rgba-to-float
// BSD 3-Clause License
//
// Copyright (c) 2019, Institute for Health Metrics and Evaluation All rights reserved.
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//  - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//  - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//  - Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
// OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.

#ifdef CUMULATIVE

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

float rgbaToFloat(vec4 texelRGBA) {
  ivec4 rgbaBytes = floatsToBytes(texelRGBA);
  bool bits[32];
  bytesToBits(rgbaBytes, bits);
  return bitsToFloat(bits);
}

#endif

///////////////////////////////////////////////////////

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

    {MAIN}

    #ifdef CUMULATIVE
        float current = rgbaToFloat(texture2D(tCumulativeSum, gl_FragCoord.xy / vec2(uWidth, uWidth)));
    #endif
    gl_FragColor = floatToRgba({RETURN});
}
`;