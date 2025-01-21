/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
export const image_vert = `
precision highp float;
precision highp int;

#include common
#include common_vert_params

attribute vec3 aPosition;
attribute vec2 aUv;
attribute mat4 aTransform;
attribute float aInstance;

varying vec2 vUv;
varying float vInstance;

void main() {
    int vertexId = VertexID;

    #include assign_position

    vUv = aUv;
    vInstance = aInstance;
}
`;