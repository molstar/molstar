/**
 * Copyright (c) 2020-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
export const image_frag = `
precision highp float;
precision highp int;

#include common
#include read_from_texture
#include common_frag_params
#include common_clip

// Density value to estimate object thickness
uniform float uDensity;

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
                float c = abs(cubicFilter(m - cell.x) * cubicFilter(-n + cell.y));
                nSum += vecData * c;
                nDenom += c;
            }
        }
        return nSum / nDenom;
    }
#endif

void main() {
    #include fade_lod
    #include clip_pixel

    #if defined(dInterpolation_cubic)
        vec4 imageData = biCubic(tImageTex, vUv);
    #else
        vec4 imageData = texture2D(tImageTex, vUv);
    #endif
    imageData.a = clamp(imageData.a, 0.0, 1.0);
    if (imageData.a > 0.9) imageData.a = 1.0;

    imageData.a *= uAlpha;
    if (imageData.a < 0.05)
        discard;

    float fragmentDepth = gl_FragCoord.z;

    if ((uRenderMask == MaskOpaque && imageData.a < 1.0) ||
        (uRenderMask == MaskTransparent && imageData.a == 1.0)
    ) {
        discard;
    }

    #if defined(dRenderVariant_pick)
        if (imageData.a < 0.3)
            discard;
        #ifdef requiredDrawBuffers
            gl_FragColor = vec4(packIntToRGB(float(uObjectId)), 1.0);
            gl_FragData[1] = vec4(packIntToRGB(vInstance), 1.0);
            gl_FragData[2] = vec4(texture2D(tGroupTex, vUv).rgb, 1.0);
            gl_FragData[3] = packDepthToRGBA(fragmentDepth);
        #else
            gl_FragColor = vColor;
            if (uPickType == 1) {
                gl_FragColor = vec4(packIntToRGB(float(uObjectId)), 1.0);
            } else if (uPickType == 2) {
                gl_FragColor = vec4(packIntToRGB(vInstance), 1.0);
            } else {
                gl_FragColor = vec4(texture2D(tGroupTex, vUv).rgb, 1.0);
            }
        #endif
    #elif defined(dRenderVariant_depth)
        if (imageData.a < 0.05)
            discard;
        if (uRenderMask == MaskOpaque) {
            gl_FragColor = packDepthToRGBA(fragmentDepth);
        } else if (uRenderMask == MaskTransparent) {
            gl_FragColor = packDepthWithAlphaToRGBA(fragmentDepth, imageData.a);
        }
    #elif defined(dRenderVariant_marking)
        float marker = uMarker;
        if (uMarker == -1.0) {
            float group = unpackRGBToInt(texture2D(tGroupTex, vUv).rgb);
            marker = readFromTexture(tMarker, vInstance * float(uGroupCount) + group, uMarkerTexDim).a;
            marker = floor(marker * 255.0 + 0.5); // rounding required to work on some cards on win
        }
        if (uMarkingType == 1) {
            if (marker > 0.0 || imageData.a < 0.05)
                discard;
            gl_FragColor = packDepthToRGBA(fragmentDepth);
        } else {
            if (marker == 0.0 || imageData.a < 0.05)
                discard;
            float depthTest = 1.0;
            if (uMarkingDepthTest) {
                depthTest = (fragmentDepth >= getDepthPacked(gl_FragCoord.xy / uDrawingBufferSize)) ? 1.0 : 0.0;
            }
            bool isHighlight = intMod(marker, 2.0) > 0.1;
            gl_FragColor = vec4(0.0, depthTest, isHighlight ? 1.0 : 0.0, 1.0);
        }
    #elif defined(dRenderVariant_emissive)
        gl_FragColor = vec4(0.0);
    #elif defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
        gl_FragColor = imageData;

        float marker = uMarker;
        if (uMarker == -1.0) {
            float group = unpackRGBToInt(texture2D(tGroupTex, vUv).rgb);
            marker = readFromTexture(tMarker, vInstance * float(uGroupCount) + group, uMarkerTexDim).a;
            marker = floor(marker * 255.0 + 0.5); // rounding required to work on some cards on win
        }

        #include apply_marker_color

        #if defined(dRenderVariant_color)
            #include apply_fog
            #include wboit_write
            #include dpoit_write
        #elif defined(dRenderVariant_tracing)
            gl_FragData[1] = vec4(normalize(vViewPosition), 0.0);
            gl_FragData[2] = vec4(imageData.rgb, uDensity);
        #endif
    #endif
}
`;
