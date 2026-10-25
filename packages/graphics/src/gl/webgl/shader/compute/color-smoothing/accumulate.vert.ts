/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from '../../chunks/common.glsl.js';
import { read_from_texture } from '../../chunks/read-from-texture.glsl.js';

export const accumulate_vert = `
precision highp float;

${common}
${read_from_texture}

#define SampleID int(aSample)

varying vec3 vPosition;
varying vec4 vColor;

void main() {
    vec3 position = readFromTexture(tPosition, SampleID, uGeoTexDim).xyz;
    vec3 packedGroup = readFromTexture(tGroup, SampleID, uGeoTexDim).rgb;
    // a vertex of a cell no corner of which has a group: no color (outside the target)
    if (packedGroup == vec3(0.0)) {
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
        gl_PointSize = 0.0;
        return;
    }
    float group = unpackRGBToInt(packedGroup);

    position = (aTransform * vec4(position, 1.0)).xyz;
    gl_PointSize = 7.0;
    vPosition = (position - uBboxMin) / uResolution;
    gl_Position = vec4(((position - uBboxMin) / uBboxSize) * 2.0 - 1.0, 1.0);
    #ifdef clipZeroToOne
        // z inside the clip depth range, as with [-1, 1]
        gl_Position.z = gl_Position.z * 0.5 + 0.5;
    #endif

    #if defined(dColorType_group)
        vColor = readFromTexture(tColor, group, uColorTexDim);
    #elif defined(dColorType_groupInstance)
        vColor = readFromTexture(tColor, aInstance * float(uGroupCount) + group, uColorTexDim);
    #endif
}
`;
