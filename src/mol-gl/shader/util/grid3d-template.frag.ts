/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export const grid3dTemplate_frag = `
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

#include float_to_rgba
#ifdef CUMULATIVE
    #include rgba_to_float
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

    {MAIN}

    #ifdef CUMULATIVE
        float current = rgbaToFloat(texture2D(tCumulativeSum, gl_FragCoord.xy / vec2(uWidth, uWidth)), uLittleEndian);
    #endif
    gl_FragColor = floatToRgba({RETURN}, uLittleEndian);
}
`;