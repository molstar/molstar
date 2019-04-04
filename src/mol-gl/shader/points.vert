/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-vert-params.glsl')
#pragma glslify: import('./chunks/color-vert-params.glsl')
#pragma glslify: import('./chunks/size-vert-params.glsl')

uniform float uPixelRatio;
uniform float uViewportHeight;

attribute vec3 aPosition;
attribute mat4 aTransform;
attribute float aInstance;
attribute float aGroup;

void main(){
    #pragma glslify: import('./chunks/assign-group.glsl')
    #pragma glslify: import('./chunks/assign-color-varying.glsl')
    #pragma glslify: import('./chunks/assign-marker-varying.glsl')
    #pragma glslify: import('./chunks/assign-position.glsl')
    #pragma glslify: import('./chunks/assign-size.glsl')

    #ifdef dPointSizeAttenuation
        gl_PointSize = size * uPixelRatio * ((uViewportHeight / 2.0) / -mvPosition.z) * 5.0;
    #else
        gl_PointSize = size * uPixelRatio;
    #endif

    gl_Position = uProjection * mvPosition;
}