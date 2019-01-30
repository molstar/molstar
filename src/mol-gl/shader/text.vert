/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-vert-params.glsl')
#pragma glslify: import('./chunks/color-vert-params.glsl')
#pragma glslify: import('./chunks/size-vert-params.glsl')

uniform mat4 uModelView;

attribute vec3 aPosition;
attribute vec2 aMapping;
attribute float aDepth;
attribute vec2 aTexCoord;
attribute mat4 aTransform;
attribute float aInstance;
attribute float aGroup;

uniform float uOffsetX;
uniform float uOffsetY;
uniform float uOffsetZ;

// uniform bool ortho;
uniform float uPixelRatio;
uniform float uViewportHeight;

varying vec2 vTexCoord;

#pragma glslify: matrixScale = require(./utils/matrix-scale.glsl)

void main(void){
    #pragma glslify: import('./chunks/assign-color-varying.glsl')
    #pragma glslify: import('./chunks/assign-marker-varying.glsl')
    #pragma glslify: import('./chunks/assign-size.glsl')

    vTexCoord = aTexCoord;

    float scale = matrixScale(uModelView);

    float offsetX = uOffsetX * scale;
    float offsetY = uOffsetY * scale;
    float offsetZ = (uOffsetZ + aDepth * 0.95) * scale;
    if (vTexCoord.x == 10.0) {
        offsetZ -= 0.01;
    }

    vec4 mvPosition = uModelView * aTransform * vec4(aPosition, 1.0);

    // #ifdef FIXED_SIZE
    //     if (ortho) {
    //         scale /= pixelRatio * ((uViewportHeight / 2.0) / -uCameraPosition.z) * 0.1;
    //     } else {
    //         scale /= pixelRatio * ((uViewportHeight / 2.0) / -mvPosition.z) * 0.1;
    //     }
    // #endif

    vec4 mvCorner = vec4(mvPosition.xyz, 1.0);
    mvCorner.xy += aMapping * size * scale;
    mvCorner.x += offsetX;
    mvCorner.y += offsetY;
    // if(ortho){
    //     mvCorner.xyz += normalize(-uCameraPosition) * offsetZ;
    // } else {
    //     mvCorner.xyz += normalize(-mvCorner.xyz) * offsetZ;
    // }
    mvCorner.xyz += normalize(-mvCorner.xyz) * offsetZ;

    gl_Position = uProjection * mvCorner;

    vViewPosition = -mvCorner.xyz;
}