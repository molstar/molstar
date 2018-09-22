/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * heavily based on code by WestLangley from https://github.com/WestLangley/three.js/blob/af28b2fb706ac109771ecad0a7447fad90ab3210/examples/js/lines/LineMaterial.js
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

attribute vec2 aMapping;
attribute vec3 aStart;
attribute vec3 aEnd;

void trimSegment(const in vec4 start, inout vec4 end) {
    // trim end segment so it terminates between the camera plane and the near plane
    // conservative estimate of the near plane
    float a = uProjection[2][2];  // 3nd entry in 3th column
    float b = uProjection[3][2];  // 3nd entry in 4th column
    float nearEstimate = -0.5 * b / a;
    float alpha = (nearEstimate - start.z) / (end.z - start.z);
    end.xyz = mix(start.xyz, end.xyz, alpha);
}

void main(){
    #pragma glslify: import('./chunks/assign-color-varying.glsl')

    // TODO move to chunk (also in point.vert)
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

    float linewidth = 3.0; // size;

    mat4 modelView = uView * uModel * aTransform;

    // camera space
    vec4 start = modelView * vec4(aStart, 1.0);
    vec4 end = modelView * vec4(aEnd, 1.0);

    // special case for perspective projection, and segments that terminate either in, or behind, the camera plane
    // clearly the gpu firmware has a way of addressing this issue when projecting into ndc space
    // but we need to perform ndc-space calculations in the shader, so we must address this issue directly
    // perhaps there is a more elegant solution -- WestLangley
    bool perspective = (uProjection[2][3] == -1.0); // 4th entry in the 3rd column
    if (perspective) {
        if (start.z < 0.0 && end.z >= 0.0) {
            trimSegment(start, end);
        } else if (end.z < 0.0 && start.z >= 0.0) {
            trimSegment(end, start);
        }
    }

    // clip space
    vec4 clipStart = uProjection * start;
    vec4 clipEnd = uProjection * end;

    // ndc space
    vec2 ndcStart = clipStart.xy / clipStart.w;
    vec2 ndcEnd = clipEnd.xy / clipEnd.w;

    // direction
    vec2 dir = ndcEnd - ndcStart;

    // account for clip-space aspect ratio
    dir.x *= uPixelRatio;
    dir = normalize(dir);

    // perpendicular to dir
    vec2 offset = vec2(dir.y, - dir.x);

    // undo aspect ratio adjustment
    dir.x /= uPixelRatio;
    offset.x /= uPixelRatio;

    // sign flip
    if (aMapping.x < 0.0) offset *= -1.0;

    // adjust for linewidth
    offset *= linewidth;

    // adjust for clip-space to screen-space conversion
    offset /= uViewportHeight;

    // select end
    vec4 clip = (aMapping.y < 0.5) ? clipStart : clipEnd;

    // back to clip space
    offset *= clip.w;
    clip.xy += offset;
    gl_Position = clip;

    // gl_Position = uProjection * (modelView * vec4(aEnd.x * 5.0 - 5.0, aMapping.y * 5.0, 2.0, 1.0));

    // TODO
    // vViewPosition = (projectionMatrixInverse * clip).xyz;
}