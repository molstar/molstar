/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { common } from './chunks/common.glsl.js';
import { common_frag_params } from './chunks/common-frag-params.glsl.js';
import { color_frag_params } from './chunks/color-frag-params.glsl.js';
import { light_frag_params } from './chunks/light-frag-params.glsl.js';
import { normal_frag_params } from './chunks/normal-frag-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { fade_lod } from './chunks/fade-lod.glsl.js';
import { clip_pixel } from './chunks/clip-pixel.glsl.js';
import { assign_material_color } from './chunks/assign-material-color.glsl.js';
import { check_transparency } from './chunks/check-transparency.glsl.js';
import { check_picking_alpha } from './chunks/check-picking-alpha.glsl.js';
import { apply_marking } from './chunks/apply-marking.glsl.js';
import { apply_interior_color } from './chunks/apply-interior-color.glsl.js';
import { apply_light_color } from './chunks/apply-light-color.glsl.js';
import { apply_marker_color } from './chunks/apply-marker-color.glsl.js';
import { apply_fog } from './chunks/apply-fog.glsl.js';
import { wboit_write } from './chunks/wboit-write.glsl.js';
import { dpoit_write } from './chunks/dpoit-write.glsl.js';

export const mesh_frag = `
precision highp float;
precision highp int;

#define bumpEnabled

${common}
${common_frag_params}
${color_frag_params}
${light_frag_params}
${normal_frag_params}
${common_clip}

void main() {
    #ifdef dSolidInterior
        float fragmentDepth = gl_FragCoord.z;
        bool capPass = uSolidInteriorPass != 0;
        vec3 viewPosition = vViewPosition;
        vec3 modelPosition = vModelPosition;
        #ifdef enabledFragDepth
            if (capPass) {
                float nearZ = -uNear * 1.0001;
                vec3 nearPosition = vec3(mix(vViewPosition.xy * (nearZ / vViewPosition.z), vViewPosition.xy, uIsOrtho), nearZ);
                float s = 0.0;
                #if defined(dClipObjects)
                    if (uSolidInteriorClip >= 0) {
                        // exempt parts are not cut by this clip object, so they neither count nor fill
                        if (clipIgnored(uSolidInteriorClip)) discard;
                        s = clipCapExit((uInvView * vec4(nearPosition, 1.0)).xyz / uModelScale, vModelPosition / uModelScale, uDepthBack);
                        if (s < 0.0 || s > 1.0) discard;
                    }
                #endif
                if (uSolidInteriorPass == 2) {
                    gl_FragColor = vec4(0.0);
                    gl_FragDepthEXT = fragmentDepth;
                    return;
                }
                viewPosition = mix(nearPosition, vViewPosition, s);
                modelPosition = (uInvView * vec4(viewPosition, 1.0)).xyz;
                fragmentDepth = mix(calcDepth(viewPosition), gl_FragCoord.z, 0.0001);
                // beyond the far plane (depth is reversed)
                if (fragmentDepth < 0.0) discard;
            }
            gl_FragDepthEXT = fragmentDepth;
        #endif
        #if defined(dClipVariant_pixel) && defined(dClipObjects)
            if (clipTest(modelPosition / uModelScale)) discard;
        #endif
        vec3 vViewPosition = viewPosition;
        vec3 vModelPosition = modelPosition;
        // after the mark pass returned (it needs every face) and at the cap position
        ${fade_lod}
    #else
        ${fade_lod}
        ${clip_pixel}
        float fragmentDepth = gl_FragCoord.z;
    #endif

    interior = !gl_FrontFacing;
    #ifdef dSolidInterior
        // the fill may draw front faces (reflected instances), a cap is always interior
        #if defined(dRenderVariant_depth)
            // a cap is the visible surface, so the depth passes must not discard it as interior
            if (capPass) interior = false;
        #else
            if (capPass) interior = true;
        #endif
    #endif

    #ifdef dNeedsNormal
        #if defined(dFlatShaded)
            vec3 fdx = dFdx(vViewPosition);
            vec3 fdy = dFdy(vViewPosition);
            vec3 normal = -normalize(cross(fdx,fdy));
        #else
            vec3 normal = -normalize(vNormal);
            if (uDoubleSided) normal *= float(gl_FrontFacing) * 2.0 - 1.0;
        #endif

        #if defined(dFlipSided)
            normal *= -1.0;
        #endif

        #ifdef dSolidInterior
            if (capPass) {
                normal = vec3(0.0, 0.0, -1.0);
                #if defined(dClipObjects)
                    if (uSolidInteriorClip >= 0) normal = normalize(clipCapNormal(vModelPosition / uModelScale) * mat3(uInvView));
                #endif
            }
        #endif
    #endif

    ${assign_material_color}
    ${check_transparency}

    #if defined(dRenderVariant_pick)
        ${check_picking_alpha}
        #ifdef requiredDrawBuffers
            gl_FragColor = vObject;
            gl_FragData[1] = vInstance;
            #ifdef dSolidInterior
                // 16777214 is PickingId.Null, the cap picks the whole instance
                gl_FragData[2] = capPass ? vec4(packIntToRGB(16777214.0), 1.0) : vGroup;
            #else
                gl_FragData[2] = vGroup;
            #endif
            gl_FragData[3] = packDepthToRGBA(fragmentDepth);
        #else
            gl_FragColor = vColor;
            #ifdef dSolidInterior
                if (capPass && uPickType == 3) gl_FragColor = vec4(packIntToRGB(16777214.0), 1.0);
            #endif
        #endif
    #elif defined(dRenderVariant_depth)
        gl_FragColor = material;
    #elif defined(dRenderVariant_marking)
        ${apply_marking}
        gl_FragColor = material;
    #elif defined(dRenderVariant_emissive)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
        ${apply_interior_color}
        ${apply_light_color}
        ${apply_marker_color}

        #if defined(dRenderVariant_color)
            ${apply_fog}
            ${wboit_write}
            ${dpoit_write}
        #elif defined(dRenderVariant_tracing)
            gl_FragData[1] = packNormalToRGBA(normal, emissive);
            gl_FragData[2] = vec4(material.rgb, uDensity);
        #endif
    #endif
}
`;
