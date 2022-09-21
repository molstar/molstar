

import { loadCheckpoint } from '../../../mol-util/debug';
loadCheckpoint(`mol-gl/shader/chunks/apply-interior-color.glsl.ts::start`);
export const apply_interior_color = `
if (interior) {
    if (uInteriorColorFlag) {
        gl_FragColor.rgb = uInteriorColor;
    } else {
        gl_FragColor.rgb *= 1.0 - uInteriorDarkening;
    }

    #ifdef dTransparentBackfaces_opaque
        gl_FragColor.a = 1.0;
    #endif
}
`;
loadCheckpoint(`mol-gl/shader/chunks/apply-interior-color.glsl.ts::end`);
