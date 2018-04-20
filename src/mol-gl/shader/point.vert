/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;

uniform mat4 projection, model, view;

uniform int objectId;
uniform int instanceCount;
uniform int elementCount;

#pragma glslify: import('./chunks/color-vert-params.glsl')

attribute vec3 position;
attribute vec4 transformColumn0, transformColumn1, transformColumn2, transformColumn3;
attribute float instanceId;
attribute float elementId;

void main(){
    #pragma glslify: import('./chunks/color-assign-varying.glsl')

    mat4 transform = mat4(transformColumn0, transformColumn1, transformColumn2, transformColumn3);
    mat4 modelView = view * model * transform;

    gl_PointSize = 1.0;
    gl_Position = projection * modelView * vec4(position, 1.0);
}