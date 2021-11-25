/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const normalize_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tColor;
uniform sampler2D tCount;
uniform vec2 uTexSize;

void main(void) {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    vec4 color = texture2D(tColor, coords);
    float count = texture2D(tCount, coords).r;

    gl_FragColor = color / count;
}
`;