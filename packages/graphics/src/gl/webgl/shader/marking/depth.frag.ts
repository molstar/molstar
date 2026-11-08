/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from '../chunks/common.glsl.js';

/**
 * The marking depth from the scene's: its opaque depth, or its front transparent layer where that is nearer and less
 * transparent than \`uTransparentThreshold\`; packed with the opacity as the marking depth pass writes it, cleared to 0.0
 * (the far plane) where there is nothing. Read \`uShift\` pixels away where the scene was drawn with another camera
 * jitter than the sample: the farthest of the four texels around that point, so that a marked surface is not hidden by
 * itself drawn half a pixel away (occluders' edges recede by up to a pixel, the jitter's own uncertainty)
 */
export const markingDepth_frag = `
precision highp float;

${common}

vec2 getDepth(const in vec2 coords) {
    float depth = texture2D(tDepth, coords).r;
    float alpha = 1.0;
    if (uHasTransparentDepth) {
        vec2 layer = unpackRGBAToDepthWithAlpha(texture2D(tTransparentDepth, coords));
        // depth is reversed: nearer is larger
        if (layer.x > depth && 1.0 - layer.y < uTransparentThreshold) {
            depth = layer.x;
            alpha = layer.y;
        }
    }
    return vec2(depth, alpha);
}

void main() {
    vec2 d;
    if (uShift == vec2(0.0)) {
        d = getDepth(gl_FragCoord.xy * uTexSizeInv);
    } else {
        // the centers of the four texels around the point
        vec2 p = floor(gl_FragCoord.xy + uShift - 0.5) + 0.5;
        d = getDepth(p * uTexSizeInv);
        vec2 d1 = getDepth((p + vec2(1.0, 0.0)) * uTexSizeInv);
        if (d1.x < d.x) d = d1;
        vec2 d2 = getDepth((p + vec2(0.0, 1.0)) * uTexSizeInv);
        if (d2.x < d.x) d = d2;
        vec2 d3 = getDepth((p + vec2(1.0, 1.0)) * uTexSizeInv);
        if (d3.x < d.x) d = d3;
    }
    if (d.x <= 0.0)
        discard;
    gl_FragColor = packDepthWithAlphaToRGBA(d.x, d.y);
}
`;
