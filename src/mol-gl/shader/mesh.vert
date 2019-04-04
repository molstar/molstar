/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-vert-params.glsl')
#pragma glslify: import('./chunks/color-vert-params.glsl')

#ifdef dGeoTexture
    uniform vec2 uGeoTexDim;
    uniform sampler2D tPositionGroup;
#else
    attribute vec3 aPosition;
#endif
attribute mat4 aTransform;
attribute float aInstance;
attribute float aGroup;

#ifndef dFlatShaded
    #ifdef dGeoTexture
        uniform sampler2D tNormal;
    #else
        attribute vec3 aNormal;
    #endif
    varying vec3 vNormal;
#endif

void main(){
    #pragma glslify: import('./chunks/assign-group.glsl')
    #pragma glslify: import('./chunks/assign-color-varying.glsl')
    #pragma glslify: import('./chunks/assign-marker-varying.glsl')
    #pragma glslify: import('./chunks/assign-position.glsl')

    #ifndef dFlatShaded
        #ifdef dGeoTexture
            vec3 normal = readFromTexture(tNormal, aGroup, uGeoTexDim).xyz;
        #else
            vec3 normal = aNormal;
        #endif
        mat3 normalMatrix = transpose(inverse(mat3(modelView)));
        vec3 transformedNormal = normalize(normalMatrix * normalize(normal));
        #if defined(dFlipSided) && !defined(dDoubleSided) // TODO checking dDoubleSided should not be required, ASR
            transformedNormal = -transformedNormal;
        #endif
        vNormal = transformedNormal;
    #endif
}