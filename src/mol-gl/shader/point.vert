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
#endif

attribute vec3 aPosition;
attribute mat4 aTransform;
attribute float aInstanceId;
attribute float aElementId;

void main(){
    #pragma glslify: import('./chunks/assign-color-varying.glsl')

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