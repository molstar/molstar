/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';

/**
 * Direct volumes rendered at a lower resolution (`tVolume`, premultiplied) upsampled to the target: the four nearest
 * texels weighted bilinearly and by how close the opaque depth their rays stopped at is to the pixel's, so that the
 * volume does not bleed across the edges of opaque objects
 */
export const volumeComposite_frag = `
precision highp float;
precision highp int;

${common}

float viewZAt(const in vec2 pixel) {
    float depth = texture2D(tDepth, (pixel + 0.5) / uTexSize).r;
    return depthToViewZ(uIsOrtho, depth, uNear, uFar);
}

vec4 volumeTexel(vec2 texel, const in float bilinear, const in float viewZ, inout float weightSum) {
    texel = clamp(texel, vec2(0.0), uVolumeSize - 1.0);
    // the pixel whose opaque depth the texel's ray read (\`gl_FragCoord.xy / uDrawingBufferSize\` at the lower size)
    float texelViewZ = viewZAt(floor((texel + 0.5) / uVolumeSize * uTexSize));
    float weight = bilinear / (0.001 + abs(texelViewZ - viewZ));
    weightSum += weight;
    return texture2D(tVolume, (texel + 0.5) / uVolumeSize) * weight;
}

void main(void) {
    float viewZ = viewZAt(floor(gl_FragCoord.xy));
    vec2 p = gl_FragCoord.xy / uTexSize * uVolumeSize - 0.5;
    vec2 i = floor(p);
    vec2 f = p - i;
    float weightSum = 0.0;
    vec4 sum = volumeTexel(i, (1.0 - f.x) * (1.0 - f.y), viewZ, weightSum);
    sum += volumeTexel(i + vec2(1.0, 0.0), f.x * (1.0 - f.y), viewZ, weightSum);
    sum += volumeTexel(i + vec2(0.0, 1.0), (1.0 - f.x) * f.y, viewZ, weightSum);
    sum += volumeTexel(i + vec2(1.0, 1.0), f.x * f.y, viewZ, weightSum);
    gl_FragColor = weightSum > 0.0 ? sum / weightSum : vec4(0.0);
}
`;
