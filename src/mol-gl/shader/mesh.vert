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
attribute float aInstance;
attribute float aGroup;

#ifndef dFlatShaded
    attribute vec3 aNormal;
    varying vec3 vNormal;
#endif

void main(){
    #pragma glslify: import('./chunks/assign-color-varying.glsl')
    #pragma glslify: import('./chunks/assign-marker-varying.glsl')
    #pragma glslify: import('./chunks/assign-position.glsl')

    #ifndef dFlatShaded
        mat3 normalMatrix = transpose(inverse(mat3(modelView)));
        vec3 transformedNormal = normalize(normalMatrix * normalize(aNormal));
        #if defined(dFlipSided) && !defined(dDoubleSided) // TODO checking dDoubleSided should not be required, ASR
            transformedNormal = -transformedNormal;
        #endif
        vNormal = transformedNormal;
    #endif
}