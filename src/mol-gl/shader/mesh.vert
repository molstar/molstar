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
attribute mat4 transform;
attribute float instanceId;
attribute float elementId;

#ifndef FLAT_SHADED
    attribute vec3 normal;
    varying vec3 vNormal;
#endif

varying vec3 vViewPosition;

#pragma glslify: inverse = require(./utils/inverse.glsl)
#pragma glslify: transpose = require(./utils/transpose.glsl)

void main(){
    #pragma glslify: import('./chunks/color-assign-varying.glsl')

    mat4 modelView = view * model * transform;
    vec4 mvPosition = modelView * vec4(position, 1.0);
    vViewPosition = mvPosition.xyz;
    gl_Position = projection * mvPosition;

    #ifndef FLAT_SHADED
        mat3 normalMatrix = transpose(inverse(mat3(modelView)));
        vec3 transformedNormal = normalize(normalMatrix * normal);
        #ifdef FLIP_SIDED
            transformedNormal = -transformedNormal;
        #endif
        vNormal = transformedNormal;
    #endif
}