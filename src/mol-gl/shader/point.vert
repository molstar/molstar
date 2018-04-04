/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision mediump float;

uniform mat4 projection, model, view;

attribute vec3 position; //, color;
attribute vec4 transformColumn0, transformColumn1, transformColumn2, transformColumn3;
// attribute int instanceId;

// instanced
// attribute mat4 transform;
// uniform mat4 transform;

// varying vec3 vColor;

void main(){
    mat4 transform = mat4(transformColumn0, transformColumn1, transformColumn2, transformColumn3);
    // vColor = color;
    gl_PointSize = 20.0;
    gl_Position = projection * view * model * transform * vec4(position, 1.0);
}