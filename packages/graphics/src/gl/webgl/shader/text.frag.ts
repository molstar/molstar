/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
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

export const text_frag = `
precision highp float;
precision highp int;

${common}
${common_frag_params}
${color_frag_params}
${common_clip}

varying vec2 vTexCoord;

void main(){
    // anti-aliased SDF text (.v6/designs/labels.md §5)
    // sampled and differentiated before anything can discard: derivatives need every fragment of the 2x2 block
    float rawSdf = texture2D(tFont, vTexCoord).r;
    // how much the distance field changes across one pixel, at any label size
    float sdfPerPixel = max(length(vec2(dFdx(rawSdf), dFdy(rawSdf))), 1e-4);

    ${fade_lod}
    ${clip_pixel}

    float fragmentDepth = gl_FragCoord.z;

    // determine if this is a background or glyph fragment
    bool isBackground = vTexCoord.x > 1.0;

    // discard background for pick/marking/emissive, but keep it in depth
    #if !defined(dRenderVariant_color) && !defined(dRenderVariant_tracing) && !defined(dRenderVariant_depth)
        if (isBackground) discard;
    #endif

    // glyph + border coverage: in color and tracing, anti-aliased over one pixel (the fractional edge is drawn in the
    // transparent pass, check_transparency); in the other variants, the edge at the field's middle as before
    float edgeAlpha = 1.0;
    float fillAlpha = 1.0;
    if (!isBackground) {
        float sdf = rawSdf + min(uBorderWidth, 0.49); // clamp to avoid exceeding max SDF range
        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            // edges a quarter pixel out: keeps one-pixel stems of small text solid, invisible on large text
            edgeAlpha = clamp((sdf - 0.5) / sdfPerPixel + 0.75, 0.0, 1.0);
            fillAlpha = clamp((rawSdf - 0.5) / sdfPerPixel + 0.75, 0.0, 1.0);
            if (edgeAlpha < 1.0 / 255.0) discard;
        #else
            if (sdf < 0.5) discard;
        #endif
    }

    #ifdef enabledFragDepth
        gl_FragDepthEXT = fragmentDepth;
    #endif

    // background writes its own opacity into depth alpha; glyphs stay opaque
    #define dHasMaterialOpacity
    float materialOpacity = isBackground ? uBackgroundOpacity : 1.0;
    ${assign_material_color}

    if (isBackground) {
        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            material = vec4(uBackgroundColor, uBackgroundOpacity * material.a);
        #endif
    } else {
        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            material.a *= edgeAlpha;
            if (uBorderWidth > 0.0) material.xyz = mix(uBorderColor, material.xyz, fillAlpha);
            if (uBorderWidth == 0.0 || fillAlpha >= 0.5) {
                // push text fragments forward in depth so they render in front of border (depth is reversed)
                #ifdef enabledFragDepth
                    gl_FragDepthEXT = fragmentDepth + 0.0001;
                #endif
            }
        #endif
    }

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
            gl_FragData[1] = packNormalToRGBA(-normalize(vViewPosition), emissive);
            gl_FragData[2] = vec4(material.rgb, uDensity);
        #endif
    #endif
}
`;
