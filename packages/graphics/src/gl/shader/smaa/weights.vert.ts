/**
 * Slightly adapted from https://github.com/mrdoob/three.js
 * MIT License Copyright (c) 2010-2020 three.js authors
 *
 * WebGL port of Subpixel Morphological Antialiasing (SMAA) v2.8
 * Preset: SMAA 1x Medium (with color edge detection)
 * https://github.com/iryoku/smaa/releases/tag/v2.8
 */

export const weights_vert = `
precision highp float;

attribute vec2 aPosition;
uniform vec2 uQuadScale;

uniform vec2 uTexSizeInv;
uniform vec4 uViewport;

varying vec2 vUv;
varying vec4 vOffset[3];
varying vec2 vPixCoord;

void SMAABlendingWeightCalculationVS(vec2 texCoord) {
    vPixCoord = texCoord / uTexSizeInv;

    // We will use these offsets for the searches later on (see @PSEUDO_GATHER4):
    vOffset[0] = texCoord.xyxy + uTexSizeInv.xyxy * vec4(-0.25, 0.125, 1.25, 0.125); // WebGL port note: Changed sign in Y and W components
    vOffset[1] = texCoord.xyxy + uTexSizeInv.xyxy * vec4(-0.125, 0.25, -0.125, -1.25); // WebGL port note: Changed sign in Y and W components

    // And these for the searches, they indicate the ends of the loops:
    vOffset[2] = vec4(vOffset[0].xz, vOffset[1].yw) + vec4(-2.0, 2.0, -2.0, 2.0) * uTexSizeInv.xxyy * float(dMaxSearchSteps);
}

void main() {
    vec2 scale = uViewport.zw * uTexSizeInv;
    vec2 shift = uViewport.xy * uTexSizeInv;
    vUv = (aPosition + 1.0) * 0.5 * scale + shift;
    SMAABlendingWeightCalculationVS(vUv);
    vec2 position = aPosition * uQuadScale - vec2(1.0, 1.0) + uQuadScale;
    gl_Position = vec4(position, 0.0, 1.0);
}
`;