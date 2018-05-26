/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;

uniform mat4 uProjection, uModel, uView;

uniform int uObjectId;
uniform int uInstanceCount;
uniform int uElementCount;

uniform float uPixelRatio;
uniform float uViewportHeight;

#pragma glslify: import('./chunks/color-vert-params.glsl')

#if defined(dSizeType_uniform)
    uniform float uSize;
#elif defined(dSizeType_attribute)
    attribute float aSize;
#endif

attribute vec3 aPosition;
attribute mat4 aTransform;
attribute float aInstanceId;
attribute float aElementId;

void main(){
    #pragma glslify: import('./chunks/color-assign-varying.glsl')

    mat4 modelView = uView * uModel * aTransform;
    vec4 mvPosition = modelView * vec4(aPosition, 1.0);

    #if defined(dSizeType_uniform)
        float size = uSize;
    #elif defined(dSizeType_attribute)
        float size = aSize;
    #endif

    #ifdef dPointSizeAttenuation
        gl_PointSize = size * uPixelRatio * ((uViewportHeight / 2.0) / -mvPosition.z) * 5.0;
    #else
        gl_PointSize = size * uPixelRatio;
    #endif

    gl_Position = uProjection * mvPosition;
}