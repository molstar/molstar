/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export const gaussianDensity_vert = `
precision highp float;

attribute vec3 aPosition;
attribute float aRadius;

varying vec3 vPosition;
varying float vRadiusSqInv;

#if defined(dCalcType_groupId)
    attribute float aGroup;
    varying float vGroup;
#endif

uniform vec3 uBboxSize;
uniform vec3 uBboxMin;
uniform float uResolution;

void main() {
    vRadiusSqInv = 1.0 / (aRadius * aRadius);
    #if defined(dCalcType_groupId)
        vGroup = aGroup;
    #endif
    gl_PointSize = ceil(((aRadius * 3.0) / uResolution) + uResolution);
    vPosition = (aPosition - uBboxMin) / uResolution;
    gl_Position = vec4(((aPosition - uBboxMin) / uBboxSize) * 2.0 - 1.0, 1.0);
}
`;