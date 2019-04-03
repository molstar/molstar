/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision mediump float;

attribute vec2 aPosition;

// the output UV coordinate for each fragment
varying vec2 vCoordinate;

const vec2 scale = vec2(0.5, 0.5);

void main(void) {
    vCoordinate  = aPosition * scale + scale; // scale vertex attribute to [0,1] range
    gl_Position = vec4(aPosition, 0.0, 1.0);
}