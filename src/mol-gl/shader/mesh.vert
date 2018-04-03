/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

#define INSTANCE_COLOR

precision mediump float;

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
attribute vec4 transformColumn0, transformColumn1, transformColumn2, transformColumn3;
attribute float instanceId;
// attribute int elementId;

varying vec3 vColor;

#pragma glslify: read_vec3 = require(./read-vec3.glsl)

void main(){
    mat4 transform = mat4(transformColumn0, transformColumn1, transformColumn2, transformColumn3);
    #if defined( ATTRIBUTE_COLOR )
        vColor = color;
    #elif defined( INSTANCE_COLOR )
        vColor = read_vec3(colorTex, instanceId, colorTexSize);
    // #elif defined( ELEMENT_COLOR )
    //     vColor = read_vec3(colorTex, instanceId * instanceCount + elementId, colorTexSize);
    #else
        vColor = vec3(0.0, 1.0, 0.0);
    #endif

    gl_Position = projection * view * model * transform * vec4(position, 1.0);
}