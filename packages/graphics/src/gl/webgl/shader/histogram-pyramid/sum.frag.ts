/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const sum_frag = `
precision highp float;
precision highp int;

void main(void) {
    #if __VERSION__ == 100
        gl_FragColor = texture2D(tTexture, vec2(0.5));
    #else
        gl_FragColor = ivec4(texture2D(tTexture, vec2(0.5)).r);
    #endif
}
`;
