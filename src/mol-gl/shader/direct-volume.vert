/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

precision highp float;

attribute vec3 aPosition;

varying vec3 unitCoord;
varying vec3 origPos;

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform mat4 uTransform;

uniform mat4 uInvView;
uniform mat4 uModelView;
uniform mat4 uInvModelView;
uniform mat4 uProjection, uView, uModel;

varying vec4 vNearPos;
varying vec4 vFarPos;
varying vec3 vPosition;

void main() {
    unitCoord = aPosition + vec3(0.5);
    vec4 mvPosition = uView * uModel * vec4(unitCoord * uBboxSize + uBboxMin, 1.0);
    // vec4 mvPosition = vec4(unitCoord * uBboxSize + uBboxMin, 1.0);
    // origPos = mvPosition.xyz;
    origPos = unitCoord * uBboxSize + uBboxMin;
    gl_Position = uProjection * mvPosition;
}

// void main() {
//     // Project local vertex coordinate to camera position. Then do a step
//     // backward (in cam coords) to the near clipping plane, and project back. Do
//     // the same for the far clipping plane. This gives us all the information we
//     // need to calculate the ray and truncate it to the viewing cone.
//     vec3 position = aPosition * uBboxSize + uBboxMin;
//     vec4 position4 = vec4(position, 1.0);
//     vec4 posInCam = uView * position4;
//     // Intersection of ray and near clipping plane (z = -1 in clip coords)
//     posInCam.z = -posInCam.w;
//     vNearPos = uInvView * posInCam;
//     // Intersection of ray and far clipping plane (z = +1 in clip coords)
//     posInCam.z = posInCam.w;
//     vFarPos = uInvView * posInCam;
//     // Set varyings and output pos
//     vPosition = position;
//     gl_Position = uProjection * uModelView * position4;
// }