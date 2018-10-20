/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

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

uniform mat4 uModelView;
uniform mat4 uProjection;

void main() {
    unitCoord = aPosition + vec3(0.5);
    vec4 mvPosition = uModelView * vec4(unitCoord * uBboxSize + uBboxMin, 1.0);
    origPos = unitCoord * uBboxSize + uBboxMin;
    instance = aInstance;
    gl_Position = uProjection * mvPosition;
}