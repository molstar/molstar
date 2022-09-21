

import { loadCheckpoint } from '../../mol-util/debug';
loadCheckpoint(`mol-gl/shader/background.vert.ts::start`);
export const background_vert = `
precision mediump float;

attribute vec2 aPosition;

varying vec4 vPosition;

void main() {
    vPosition = vec4(aPosition, 1.0, 1.0);
    gl_Position = vec4(aPosition, 1.0, 1.0);
}
`;

loadCheckpoint(`mol-gl/shader/background.vert.ts::end`);
