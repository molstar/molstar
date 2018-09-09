/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-vert-params.glsl')
#pragma glslify: import('./chunks/color-vert-params.glsl')

uniform float uPixelRatio;
uniform float uViewportHeight;

#if defined(dSizeType_uniform)
    uniform float uSize;
#elif defined(dSizeType_attribute)
    attribute float aSize;
#elif defined(dSizeType_instance) || defined(dSizeType_group) || defined(dSizeType_groupInstance)
    varying vec4 vSize;
    uniform vec2 uSizeTexDim;
    uniform sampler2D tSize;
#endif

attribute vec3 aPosition;
attribute mat4 aTransform;
attribute float aInstance;
attribute float aGroup;

void main(){
    #pragma glslify: import('./chunks/assign-color-varying.glsl')
    #pragma glslify: import('./chunks/assign-position.glsl')

    #if defined(dSizeType_uniform)
        float size = uSize;
    #elif defined(dSizeType_attribute)
        float size = aSize;
    #elif defined(dSizeType_instance)
        float size = readFromTexture(tSize, aInstance, uSizeTexDim).r;
    #elif defined(dSizeType_group)
        float size = readFromTexture(tSize, aGroup, uSizeTexDim).r;
    #elif defined(dSizeType_groupInstance)
        float size = readFromTexture(tSize, aInstance * float(uGroupCount) + aGroup, uSizeTexDim).r;
    #endif

    #ifdef dPointSizeAttenuation
        gl_PointSize = size * uPixelRatio * ((uViewportHeight / 2.0) / -mvPosition.z) * 5.0;
    #else
        gl_PointSize = size * uPixelRatio;
    #endif

    gl_Position = uProjection * mvPosition;
}