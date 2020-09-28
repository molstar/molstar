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

varying vec3 unitCoord;
varying vec3 origPos;
varying float instance;

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform mat4 uTransform;

uniform mat4 uUnitToCartn;
uniform mat4 uCartnToUnit;

uniform mat4 uModelView;
uniform mat4 uProjection;

void main() {
    unitCoord = aPosition + vec3(0.5);
    vec4 mvPosition = uModelView * uUnitToCartn * vec4(unitCoord, 1.0);
    origPos = (uUnitToCartn * vec4(unitCoord, 1.0)).xyz;
    instance = aInstance;
    gl_Position = uProjection * mvPosition;

    // clamp z position to clip space
    if(gl_Position.z > gl_Position.w) {
        gl_Position.z = gl_Position.w - 0.0001;
    } else if(gl_Position.z < -gl_Position.w) {
        gl_Position.z = -gl_Position.w + 0.0001;
    }
}
`;