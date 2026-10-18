/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
import { common } from './chunks/common.glsl.js';
import { common_vert_params } from './chunks/common-vert-params.glsl.js';
import { assign_position } from './chunks/assign-position.glsl.js';

export const image_vert = `
precision highp float;
precision highp int;

${common}
${common_vert_params}

varying vec2 vUv;
varying float vInstance;
varying vec3 vPosition;

void main() {
    int vertexId = VertexID;

    ${assign_position}

    vUv = aUv;
    vInstance = aInstance;
    vPosition = aPosition;
}
`;
