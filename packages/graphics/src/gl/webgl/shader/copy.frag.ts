import { common } from './chunks/common.glsl.js';

export const copy_frag = `
precision highp float;

${common}

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    #ifdef dPackDepth
        gl_FragColor = packDepthToRGBA(texture2D(tColor, coords).r);
    #else
        gl_FragColor = texture2D(tColor, coords);
    #endif
}
`;
