/**
 * Copyright (c) 2021-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    float group = unpackRGBToInt(readFromTexture(tGroup, SampleID, uGeoTexDim).rgb);

    position = (aTransform * vec4(position, 1.0)).xyz;
    gl_PointSize = 7.0;
    vPosition = (position - uBboxMin) / uResolution;
    gl_Position = vec4(((position - uBboxMin) / uBboxSize) * 2.0 - 1.0, 1.0);

    #if defined(dColorType_group)
        vColor = readFromTexture(tColor, group, uColorTexDim);
    #elif defined(dColorType_groupInstance)
        vColor = readFromTexture(tColor, aInstance * float(uGroupCount) + group, uColorTexDim);
    #endif
}
`;
