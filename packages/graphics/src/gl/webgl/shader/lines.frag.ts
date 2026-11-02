/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';
import { common_frag_params } from './chunks/common-frag-params.glsl.js';
import { color_frag_params } from './chunks/color-frag-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { fade_lod } from './chunks/fade-lod.glsl.js';
import { clip_pixel } from './chunks/clip-pixel.glsl.js';
import { assign_material_color } from './chunks/assign-material-color.glsl.js';
import { check_transparency } from './chunks/check-transparency.glsl.js';
import { check_picking_alpha } from './chunks/check-picking-alpha.glsl.js';
import { apply_marking } from './chunks/apply-marking.glsl.js';
import { apply_marker_color } from './chunks/apply-marker-color.glsl.js';
import { apply_fog } from './chunks/apply-fog.glsl.js';
import { wboit_write } from './chunks/wboit-write.glsl.js';
import { dpoit_write } from './chunks/dpoit-write.glsl.js';

export const lines_frag = `
precision highp float;
precision highp int;

${common}
${common_frag_params}
${color_frag_params}
${common_clip}

#if defined(dClipVariant_pixel) && defined(dClipObjects)
    // 1 where the clip objects may cut the line (\`clipSphereClass\` of its bounding sphere): its clip tests skipped if 0
    #if __VERSION__ == 100
        varying float vClipSome;
    #else
        flat in float vClipSome;
    #endif
#endif

void main(){
    ${fade_lod}
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        if (vClipSome > 0.5) {
            ${clip_pixel}
        }
    #endif

    float fragmentDepth = gl_FragCoord.z;
    ${assign_material_color}
    ${check_transparency}

    #if defined(dRenderVariant_pick)
        ${check_picking_alpha}
        #ifdef requiredDrawBuffers
            gl_FragColor = vObject;
            gl_FragData[1] = vInstance;
            gl_FragData[2] = vGroup;
            gl_FragData[3] = packDepthToRGBA(fragmentDepth);
        #else
            gl_FragColor = vColor;
        #endif
    #elif defined(dRenderVariant_depth)
        gl_FragColor = material;
    #elif defined(dRenderVariant_marking)
        ${apply_marking}
        gl_FragColor = material;
    #elif defined(dRenderVariant_emissive)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
        gl_FragColor = material;
        ${apply_marker_color}

        #if defined(dRenderVariant_color)
            ${apply_fog}
            ${wboit_write}
            ${dpoit_write}
        #elif defined(dRenderVariant_tracing)
            gl_FragData[1] = packNormalToRGBA(normalize(vViewPosition), emissive);
            gl_FragData[2] = vec4(material.rgb, uDensity);
        #endif
    #endif
}
`;
