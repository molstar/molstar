/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';

export const dof_frag = `
precision highp float;
precision highp int;

${common}

// Function to convert depth value from depth buffer to view space Z
float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

// Retrieve depth from opaque depth texture
float getDepthOpaque(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepthOpaque, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepthOpaque, coords));
    #endif
}

// Retrieve depth from transparent depth texture
float getDepthTransparent(const in vec2 coords) {
    return unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords)).x;
}

bool isBackground(const in float depth) {
    return depth == 0.0;
}

float getDepth(const in vec2 coords) {
    return max(getDepthOpaque(coords), getDepthTransparent(coords));
}

float getCOC(vec2 uv) {
    float depth = getDepth(uv);
    float viewDist = getViewZ(depth);
    vec3 aposition = screenSpaceToViewSpace(vec3(uv.xy, depth), uInvProjection);
    float focusDist = length(aposition - uCenter);
    float coc = 0.0;  // Circle of Confusion
    if (uMode == 0) { // planar Depth of field
        coc = (abs(viewDist) - uInFocus) / uPPM;  //focus distance, focus range
    } else if(uMode == 1) { // spherical Depth of field
        coc = focusDist / uPPM ;
    }
    coc = clamp(coc, -1.0, 1.0);
    return coc;
}

// the color weighted by the circle of confusion at uv (rgb), and the weight (a)
vec4 getWeightedColor(const in vec2 uv) {
    float coc = smoothstep(0.0, 1.0, abs(getCOC(uv)));
    return vec4(texture2D(tColor, uv).rgb * coc, coc);
}

// the offset of sample i of a row or column, in texels, centered (texel centers for odd sizes, not their edges, where
// rounding picks a neighbour)
float sampleOffset(const in int i) {
    return float(i) - float(dBlurSize - 1) / 2.0;
}

// Box blur weighted by the circle of confusion: weighted color sums (rgb) and weights (a). The weights are per sample,
// so the blur separates exactly: dDofPass 'horizontal' writes the sums of a row, 'vertical' adds those of the column
vec4 getBlurSums(const in vec2 coords) {
    vec2 texelSize = vec2(1.0 / uTexSize.x, 1.0 / uTexSize.y);
    vec4 sums = vec4(0.0);
    #if defined(dDofPass_horizontal)
        for (int x = 0; x < int(dBlurSize); x++) {
            sums += getWeightedColor(coords + vec2(sampleOffset(x), 0.0) * texelSize * uBlurSpread);
        }
    #elif defined(dDofPass_vertical)
        for (int y = 0; y < int(dBlurSize); y++) {
            sums += texture2D(tBlur, coords + vec2(0.0, sampleOffset(y)) * texelSize * uBlurSpread);
        }
    #else
        for (int x = 0; x < int(dBlurSize); x++) {
            for (int y = 0; y < int(dBlurSize); y++) {
                sums += getWeightedColor(coords + vec2(sampleOffset(x), sampleOffset(y)) * texelSize * uBlurSpread);
            }
        }
    #endif
    return sums;
}

// simplification from https://catlikecoding.com/unity/tutorials/advanced-rendering/depth-of-field/
void main() {
    vec2 uv = gl_FragCoord.xy / uTexSize;
    #if defined(dDofPass_horizontal)
        gl_FragColor = getBlurSums(uv);
    #else
        vec4 color = texture2D(tColor, uv);
        vec4 sums = getBlurSums(uv);
        // no sample out of focus: the pixel itself
        vec3 blurColor = sums.a > 0.0 ? sums.rgb / sums.a : color.rgb;
        float coc = getCOC(uv); // Circle of Confusion
        color.rgb = mix(color.rgb, blurColor, smoothstep(0.0, 1.0, abs(coc))); // Smooth blending based on CoC
        gl_FragColor = color;
    #endif
}
`;
