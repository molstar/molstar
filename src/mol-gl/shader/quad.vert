/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;

attribute vec2 aPosition;

// the output UV coordinate for each fragment
varying vec2 vCoordinate;

uniform vec2 uScale;

const vec2 scale = vec2(0.5, 0.5);

void main(void) {
    vec2 s = scale * uScale;
    vec2 position = aPosition * uScale - vec2(1.0, 1.0) + uScale;
    vCoordinate = position * s + s; // scale vertex attribute to [0,1] range
    gl_Position = vec4(position, 0.0, 1.0);
}