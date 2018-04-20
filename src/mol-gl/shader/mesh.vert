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

attribute vec3 normal;

varying vec3 vNormal;
varying vec3 vViewPosition;

#pragma glslify: inverse = require(./utils/inverse.glsl)
#pragma glslify: transpose = require(./utils/transpose.glsl)

void main(){
    #pragma glslify: import('./chunks/color-assign-varying.glsl')

    mat4 transform = mat4(transformColumn0, transformColumn1, transformColumn2, transformColumn3);
    mat4 modelView = view * model * transform;

    vec4 mvPosition = modelView * vec4(position, 1.0);
    vViewPosition = mvPosition.xyz;
    gl_Position = projection * mvPosition;

    mat3 normalMatrix = transpose(inverse(mat3(modelView)));
    vNormal = normalize(normalMatrix * normal);
}