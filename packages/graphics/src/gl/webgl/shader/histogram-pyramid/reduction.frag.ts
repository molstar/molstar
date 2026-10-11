import { common } from '../chunks/common.glsl.js';

export const reduction_frag = `
precision highp float;
precision highp int;

${common}

void main(void) {
    float k = 0.5 * uSize;
    vec2 position = floor((gl_FragCoord.xy / uTexSize) / uSize) * uSize;

    #if __VERSION__ == 100
        float a, b, c, d;

        if (uFirst) {
            a = texture2D(tInputLevel, position).r * 255.0;
            b = texture2D(tInputLevel, position + vec2(k, 0.0)).r * 255.0;
            c = texture2D(tInputLevel, position + vec2(0.0, k)).r * 255.0;
            d = texture2D(tInputLevel, position + vec2(k, k)).r * 255.0;
        } else {
            a = unpackRGBToInt(texture2D(tPreviousLevel, position).rgb);
            b = unpackRGBToInt(texture2D(tPreviousLevel, position + vec2(k, 0.0)).rgb);
            c = unpackRGBToInt(texture2D(tPreviousLevel, position + vec2(0.0, k)).rgb);
            d = unpackRGBToInt(texture2D(tPreviousLevel, position + vec2(k, k)).rgb);
        }
        gl_FragColor = vec4(packIntToRGB(a + b + c + d), 1.0);
    #else
        int a, b, c, d;

        if (uFirst) {
            a = int(texture2D(tInputLevel, position).r * 255.0);
            b = int(texture2D(tInputLevel, position + vec2(k, 0.0)).r * 255.0);
            c = int(texture2D(tInputLevel, position + vec2(0.0, k)).r * 255.0);
            d = int(texture2D(tInputLevel, position + vec2(k, k)).r * 255.0);
        } else {
            a = texture2D(tPreviousLevel, position).r;
            b = texture2D(tPreviousLevel, position + vec2(k, 0.0)).r;
            c = texture2D(tPreviousLevel, position + vec2(0.0, k)).r;
            d = texture2D(tPreviousLevel, position + vec2(k, k)).r;
        }
        gl_FragColor = ivec4(a + b + c + d);
    #endif
}
`;
