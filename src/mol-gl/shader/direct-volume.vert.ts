/**
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export const directVolume_vert = `
precision highp float;

attribute vec3 aPosition;
attribute mat4 aTransform;
attribute float aInstance;

uniform mat4 uModelView;
uniform mat4 uProjection;
uniform vec4 uInvariantBoundingSphere;

varying vec3 vOrigPos;
varying float vInstance;
varying vec4 vBoundingSphere;
varying mat4 vTransform;

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform mat4 uTransform;

uniform mat4 uUnitToCartn;

void main() {
    vec4 unitCoord = vec4(aPosition + vec3(0.5), 1.0);
    vec4 mvPosition = uModelView * aTransform * uUnitToCartn * unitCoord;

    vOrigPos = (aTransform * uUnitToCartn * unitCoord).xyz;
    vInstance = aInstance;
    vBoundingSphere = vec4(
        (aTransform * vec4(uInvariantBoundingSphere.xyz, 1.0)).xyz,
        uInvariantBoundingSphere.w
    );
    vTransform = aTransform;

    gl_Position = uProjection * mvPosition;

    // move z position to near clip plane (but not too close to get precision issues)
    gl_Position.z = gl_Position.w - 0.01;
}
`;