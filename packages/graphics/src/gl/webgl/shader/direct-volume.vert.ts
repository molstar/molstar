/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import { common } from './chunks/common.glsl.js';

export const directVolume_vert = `
precision highp float;

${common}

varying vec3 vModelPosition;
varying float vInstance;
varying vec4 vBoundingSphere;
varying mat4 vTransform;
// the normal matrix per instance: per fragment it costs as if it were computed per lit sample (WebGL2: up to 25 %)
varying mat3 vNormalMatrix;

void main() {
    vec4 unitCoord = vec4(aPosition + vec3(0.5), 1.0);
    vec4 mvPosition = uModelView * aTransform * uUnitToCartn * unitCoord;

    vModelPosition = (uModel * aTransform * uUnitToCartn * unitCoord).xyz;
    vInstance = aInstance;
    vBoundingSphere = vec4(
        (uModel * aTransform * vec4(uInvariantBoundingSphere.xyz, 1.0)).xyz,
        uModelScale * uInvariantBoundingSphere.w
    );
    vTransform = aTransform;
    vNormalMatrix = adjoint(uModelView * aTransform);

    gl_Position = uProjection * mvPosition;

    // move z position to the far clip plane (but not too close to get precision issues; clip z is reversed)
    #ifdef clipZeroToOne
        gl_Position.z = 0.005;
    #else
        gl_Position.z = -gl_Position.w + 0.01;
    #endif
}
`;
