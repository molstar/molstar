/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';

/**
 * The ghost of what the clip objects remove over the image: a flat tint at `uAlpha` and its edges
 * at `uEdgeAlpha` where the ghost's depth (`tGhostDepth`, packed, at a lower resolution) is nearer than the scene's:
 * its opaque depth, or its front transparent layer where that is less transparent than `uTransparentThreshold` (as
 * the marking depth reads it). Edges: a ghost texel next to an empty one or one farther apart in view space than
 * `uEdgeDistance` of its distance. Premultiplied
 */
export const clipGhost_frag = `
precision highp float;

${common}

float getGhostDepth(const in vec2 coords) {
    return unpackRGBAToDepth(texture2D(tGhostDepth, coords));
}

float getSceneDepth(const in vec2 coords) {
    #if defined(dPackedDepth)
        float depth = unpackRGBAToDepth(texture2D(tDepth, coords));
    #else
        float depth = texture2D(tDepth, coords).r;
    #endif
    if (uHasTransparentDepth) {
        vec2 layer = unpackRGBAToDepthWithAlpha(texture2D(tTransparentDepth, coords));
        // depth is reversed: nearer is larger
        if (layer.x > depth && 1.0 - layer.y < uTransparentThreshold) depth = layer.x;
    }
    return depth;
}

bool isEdge(const in vec2 coords, const in float viewZ) {
    float n = getGhostDepth(coords);
    return n <= 0.0 || abs(depthToViewZ(uIsOrtho, n, uNear, uFar) - viewZ) > uEdgeDistance * abs(viewZ);
}

void main() {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;
    float depth = getGhostDepth(coords);
    // nothing removed here (cleared to 0.0, the far plane), or behind the scene (depth is reversed)
    if (depth <= 0.0 || depth < getSceneDepth(coords)) discard;

    float viewZ = depthToViewZ(uIsOrtho, depth, uNear, uFar);
    vec2 step = uGhostTexSizeInv;
    bool edge = isEdge(coords + vec2(step.x, 0.0), viewZ) || isEdge(coords - vec2(step.x, 0.0), viewZ)
        || isEdge(coords + vec2(0.0, step.y), viewZ) || isEdge(coords - vec2(0.0, step.y), viewZ);
    float alpha = edge ? uEdgeAlpha : uAlpha;
    gl_FragColor = vec4(uColor * alpha, alpha);
}
`;
