/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export default `
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

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform mat4 uTransform;

uniform mat4 uUnitToCartn;
uniform mat4 uCartnToUnit;

void main() {
    vec3 unitCoord = aPosition + vec3(0.5);
    vec4 mvPosition = uModelView * uUnitToCartn * vec4(unitCoord, 1.0);

    vOrigPos = (uUnitToCartn * vec4(unitCoord, 1.0)).xyz;
    vInstance = aInstance;
    vBoundingSphere = uInvariantBoundingSphere;

    gl_Position = uProjection * mvPosition;

    // move z position to near clip plane
    gl_Position.z = gl_Position.w - 0.0001;
}
`;