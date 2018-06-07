/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-vert-params.glsl')
#pragma glslify: import('./chunks/color-vert-params.glsl')

attribute vec3 aPosition;
attribute mat4 aTransform;
attribute float aInstanceId;
attribute float aElementId;

#ifndef dFlatShaded
    attribute vec3 aNormal;
    varying vec3 vNormal;
#endif

varying vec3 vViewPosition;

#pragma glslify: inverse = require(./utils/inverse.glsl)
#pragma glslify: transpose = require(./utils/transpose.glsl)

void main(){
    #pragma glslify: import('./chunks/assign-color-varying.glsl')
    #pragma glslify: import('./chunks/assign-marker-varying.glsl')

    mat4 modelView = uView * uModel * aTransform;
    vec4 mvPosition = modelView * vec4(aPosition, 1.0);
    vViewPosition = mvPosition.xyz;
    gl_Position = uProjection * mvPosition;

    #ifndef dFlatShaded
        mat3 normalMatrix = transpose(inverse(mat3(modelView)));
        vec3 transformedNormal = normalize(normalMatrix * normalize(aNormal));
        #if defined(dFlipSided) && !defined(dDoubleSided) // TODO checking dDoubleSided should not be required, ASR
            transformedNormal = -transformedNormal;
        #endif
        vNormal = transformedNormal;
    #endif
}