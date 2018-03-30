/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision mediump float;
uniform mat4 projection, model, view;
attribute vec3 position;

// instanced
// attribute mat4 transform;
uniform mat4 transform;

varying vec3 vPosition;

void main(){
    gl_Position = projection * view * model * transform * vec4(position, 1.0);
}