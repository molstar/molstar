/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
export default `
precision highp float;
precision highp int;

#include common
#include read_from_texture
#include common_frag_params

uniform vec2 uImageTexDim;
uniform sampler2D tImageTex;
uniform sampler2D tGroupTex;

uniform vec2 uMarkerTexDim;
uniform sampler2D tMarker;

varying vec2 vUv;
varying float vInstance;

#if defined(dInterpolation_catmulrom) || defined(dInterpolation_mitchell) || defined(dInterpolation_bspline)
    #define dInterpolation_cubic
#endif

#if defined(dInterpolation_cubic)
    #if defined(dInterpolation_catmulrom) || defined(dInterpolation_mitchell)
        #if defined(dInterpolation_catmulrom)
            const float B = 0.0;
            const float C = 0.5;
        #elif defined(dInterpolation_mitchell)
            const float B = 0.333;
            const float C = 0.333;
        #endif

        float cubicFilter(float x){
            float f = x;
            if (f < 0.0) {
                f = -f;
            }
            if (f < 1.0) {
                return ((12.0 - 9.0 * B - 6.0 * C) * (f * f * f) +
                    (-18.0 + 12.0 * B + 6.0 * C) * (f * f) +
                    (6.0 - 2.0 * B)) / 6.0;
            }else if (f >= 1.0 && f < 2.0){
                return ((-B - 6.0 * C) * ( f * f * f)
                    + (6.0 * B + 30.0 * C) * (f * f) +
                    (-(12.0 * B) - 48.0 * C) * f +
                    8.0 * B + 24.0 * C) / 6.0;
            }else{
                return 0.0;
            }
        }
    #elif defined(dInterpolation_bspline)
        float cubicFilter(float x) {
            float f = x;
            if (f < 0.0) {
                f = -f;
            }
            if (f >= 0.0 && f <= 1.0){
                return (2.0 / 3.0) + (0.5) * (f * f * f) - (f * f);
            } else if (f > 1.0 && f <= 2.0) {
                return 1.0 / 6.0 * pow((2.0 - f), 3.0);
            }
            return 1.0;
        }
    #endif

    vec4 biCubic(sampler2D tex, vec2 texCoord) {
        vec2 texelSize = 1.0 / uImageTexDim;
        texCoord -= texelSize / 2.0;
        vec4 nSum = vec4(0.0);
        float nDenom = 0.0;
        vec2 cell = fract(texCoord * uImageTexDim);
        for (float m = -1.0; m <= 2.0; ++m) {
            for (float n = -1.0; n <= 2.0; ++n) {
                vec4 vecData = texture2D(tex, texCoord + texelSize * vec2(m, n));
                float c = cubicFilter(m - cell.x) * cubicFilter(-n + cell.y);
                nSum += vecData * c;
                nDenom += c;
            }
        }
        return nSum / nDenom;
    }
#endif

void main() {
    #if defined(dInterpolation_cubic)
        vec4 imageData = biCubic(tImageTex, vUv);
    #else
        vec4 imageData = texture2D(tImageTex, vUv);
    #endif
    imageData.a = clamp(imageData.a, 0.0, 1.0);
    if (imageData.a > 0.9) imageData.a = 1.0;

    #if defined(dRenderVariant_pick)
        if (imageData.a < 0.3)
            discard;

        #if defined(dRenderVariant_pickObject)
            gl_FragColor = vec4(encodeFloatRGB(float(uObjectId)), 1.0);
        #elif defined(dRenderVariant_pickInstance)
            gl_FragColor = vec4(encodeFloatRGB(vInstance), 1.0);
        #elif defined(dRenderVariant_pickGroup)
            float group = texture2D(tGroupTex, vUv).r;
            gl_FragColor = vec4(encodeFloatRGB(group), 1.0);
        #endif
    #elif defined(dRenderVariant_depth)
        if (imageData.a < 0.05)
            discard;

        #ifdef enabledFragDepth
            gl_FragColor = packDepthToRGBA(gl_FragDepthEXT);
        #else
            gl_FragColor = packDepthToRGBA(gl_FragCoord.z);
        #endif
    #elif defined(dRenderVariant_color)
        if (imageData.a < 0.05)
            discard;

        gl_FragColor = imageData;
        gl_FragColor.a *= uAlpha;

        float group = texture2D(tGroupTex, vUv).r;
        float vMarker = readFromTexture(tMarker, vInstance * float(uGroupCount) + group, uMarkerTexDim).a;
        #include apply_marker_color
        #include apply_fog
    #endif
}
`;