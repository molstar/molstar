/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { common } from './chunks/common.glsl.js';
import { common_frag_params } from './chunks/common-frag-params.glsl.js';
import { color_frag_params } from './chunks/color-frag-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { fade_lod } from './chunks/fade-lod.glsl.js';
import { assign_material_color } from './chunks/assign-material-color.glsl.js';
import { check_transparency } from './chunks/check-transparency.glsl.js';
import { check_picking_alpha } from './chunks/check-picking-alpha.glsl.js';
import { apply_marking } from './chunks/apply-marking.glsl.js';
import { apply_marker_color } from './chunks/apply-marker-color.glsl.js';
import { apply_fog } from './chunks/apply-fog.glsl.js';
import { wboit_write } from './chunks/wboit-write.glsl.js';
import { dpoit_write } from './chunks/dpoit-write.glsl.js';

export const gaussianSplats_frag = `
precision highp float;
precision highp int;

${common}
#if defined(dColorType_direct)
    // SOG's base color, decoded in the vertex stage (\`vColor\`)
    #define dColorType_varying
#endif
${common_frag_params}
${color_frag_params}
${common_clip}

/** the fragment's offset from the center in standard deviations along the 2D Gaussian's axes */
varying vec2 vQuad;
varying float vSplatOpacity;
varying vec3 vSplatSh;

void main(){
    ${fade_lod}
    // pixel clipping in the vertex stage: a splat tested at its center

    float fragmentDepth = gl_FragCoord.z;
    ${assign_material_color}
    #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
        // the view-dependent color over the base color, not below 0 (as 3D Gaussian splatting clamps)
        material.rgb = max(material.rgb + vSplatSh, 0.0);
    #endif

    // the Gaussian's falloff times the splat's opacity, in every variant (picked where it is seen)
    float coverage = exp(-0.5 * dot(vQuad, vQuad)) * vSplatOpacity;
    if (coverage < uSplatAlphaCutoff) discard;
    #if !defined(dRenderVariant_pick)
        // the pick variant has no material (it writes ids)
        material.a *= coverage;
    #endif
    #if defined(dRenderVariant_pick) || defined(dRenderVariant_marking)
        // picked and outlined where at least as opaque as the picking threshold: not by the faint tails of splats in
        // front of what is seen
        if (coverage < uPickingAlphaThreshold) discard;
    #endif

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
