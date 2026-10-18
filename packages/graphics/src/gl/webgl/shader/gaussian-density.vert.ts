/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export const gaussianDensity_vert = `
precision highp float;

varying vec3 vPosition;
varying float vRadiusSqInv;

#if defined(dCalcType_groupId)
    varying float vGroup;
#endif

void main() {
    vRadiusSqInv = 1.0 / (aRadius * aRadius);
    #if defined(dCalcType_groupId)
        vGroup = aGroup;
    #endif
    // a quad around the atom covering its splat, in pixels of the slice (uGridDim.xy)
    float size = ceil(((aRadius * 3.0) / uResolution) + uResolution);
    vPosition = (aPosition - uBboxMin) / uResolution;
    gl_Position = vec4(((aPosition - uBboxMin) / uBboxSize) * 2.0 - 1.0, 1.0);
    #ifdef clipZeroToOne
        // z inside the clip depth range, as with [-1, 1]
        gl_Position.z = gl_Position.z * 0.5 + 0.5;
    #endif
    gl_Position.xy += aMapping * size / uGridDim.xy;
}
`;
