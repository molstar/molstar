/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

#define ATTRIBUTE_COLOR
// #define INSTANCE_COLOR

precision highp float;

uniform mat4 projection, model, view;

uniform int objectId;
uniform int instanceCount;

#if defined( ATTRIBUTE_COLOR )
    attribute vec3 color;
#elif defined( INSTANCE_COLOR ) || defined( ELEMENT_COLOR )
    uniform vec2 colorTexSize;
    uniform sampler2D colorTex;
#endif

attribute vec3 position;
attribute vec3 normal;
attribute vec4 transformColumn0, transformColumn1, transformColumn2, transformColumn3;
attribute float instanceId;
// attribute int elementId;

varying vec3 vColor;
varying vec3 vNormal;
varying vec3 vViewPosition;

#pragma glslify: inverse = require(./inverse.glsl)
#pragma glslify: read_vec3 = require(./read-vec3.glsl)
#pragma glslify: transpose = require(./transpose.glsl)

void main(){
    #if defined( ATTRIBUTE_COLOR )
        vColor = color;
    #elif defined( INSTANCE_COLOR )
        vColor = read_vec3(colorTex, instanceId, colorTexSize);
    // #elif defined( ELEMENT_COLOR )
    //     vColor = read_vec3(colorTex, instanceId * instanceCount + elementId, colorTexSize);
    #else
        vColor = vec3(0.0, 1.0, 0.0);
    #endif

    mat4 transform = mat4(transformColumn0, transformColumn1, transformColumn2, transformColumn3);
    mat4 modelView = view * model * transform;

    vec4 mvPosition = modelView * vec4(position, 1.0);
    vViewPosition = mvPosition.xyz;
    gl_Position = projection * mvPosition;

    mat3 normalMatrix = transpose(inverse(mat3(modelView)));
    vNormal = normalize(normalMatrix * normal);
}