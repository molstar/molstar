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

uniform float pixelRatio;
uniform float viewportHeight;

#pragma glslify: import('./chunks/color-vert-params.glsl')

#if defined(UNIFORM_SIZE)
    uniform float size;
#elif defined(ATTRIBUTE_SIZE)
    attribute float size;
#endif

attribute vec3 position;
attribute vec4 transformColumn0, transformColumn1, transformColumn2, transformColumn3;
attribute float instanceId;
attribute float elementId;

void main(){
    #pragma glslify: import('./chunks/color-assign-varying.glsl')

    mat4 transform = mat4(transformColumn0, transformColumn1, transformColumn2, transformColumn3);
    mat4 modelView = view * model * transform;
    vec4 mvPosition = modelView * vec4(position, 1.0);

    #ifdef POINT_SIZE_ATTENUATION
        gl_PointSize = size * pixelRatio * ((viewportHeight / 2.0) / -mvPosition.z) * 5.0;
    #else
        gl_PointSize = size * pixelRatio;
    #endif

    gl_Position = projection * mvPosition;
}