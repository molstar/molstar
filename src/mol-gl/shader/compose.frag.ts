

import { loadCheckpoint } from '../../mol-util/debug';
loadCheckpoint(`mol-gl/shader/compose.frag.ts::start`);
export const compose_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tColor;
uniform vec2 uTexSize;
uniform float uWeight;

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    gl_FragColor = texture2D(tColor, coords) * uWeight;
}
`;
loadCheckpoint(`mol-gl/shader/compose.frag.ts::end`);
