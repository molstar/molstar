/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

precision highp float;

attribute vec3 aPosition;
attribute float aRadius;

varying vec3 position;
varying float radius;

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform vec3 uBboxMax;
uniform vec3 uGridDim;
uniform float uCurrentSlice;

void main() {
    radius = aRadius;
    float scale = max(uBboxSize.z, max(uBboxSize.x, uBboxSize.y));
    gl_PointSize = (radius / scale) * max(uGridDim.x, uGridDim.y) * 6.0;
    position = (aPosition - uBboxMin) / uBboxSize;
    gl_Position = vec4(position * 2.0 - 1.0, 1.0);
}