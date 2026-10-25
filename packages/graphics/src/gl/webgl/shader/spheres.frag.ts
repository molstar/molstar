/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { common } from './chunks/common.glsl.js';
import { common_frag_params } from './chunks/common-frag-params.glsl.js';
import { color_frag_params } from './chunks/color-frag-params.glsl.js';
import { light_frag_params } from './chunks/light-frag-params.glsl.js';
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

export const spheres_frag = `
precision highp float;
precision highp int;

#define bumpEnabled

${common}
${common_frag_params}
${color_frag_params}
${light_frag_params}
${common_clip}

varying float vRadius;
varying vec3 vPoint;
varying vec3 vPointViewPosition;

#ifdef dSolidInterior
    const bool solidInterior = true;
#else
    const bool solidInterior = false;
#endif

bool isCap = false;

bool SphereImpostor(out vec3 modelPos, out vec3 cameraPos, out vec3 cameraNormal, out bool interior, out float fragmentDepth){
    vec3 cameraSpherePos = -vPointViewPosition;

    vec3 rayOrigin = mix(vec3(0.0, 0.0, 0.0), vPoint, uIsOrtho);
    vec3 rayDirection = mix(normalize(vPoint), vec3(0.0, 0.0, 1.0), uIsOrtho);
    vec3 cameraSphereDir = mix(cameraSpherePos, rayOrigin - cameraSpherePos, uIsOrtho);

    float B = dot(rayDirection, cameraSphereDir);
    float det = B * B + vRadius * vRadius - dot(cameraSphereDir, cameraSphereDir);

    if (det < 0.0) return false;

    float sqrtDet = sqrt(det);
    float posT = mix(B + sqrtDet, B - sqrtDet, uIsOrtho);
    float negT = mix(B - sqrtDet, B + sqrtDet, uIsOrtho);

    #if defined(dRenderVariant_depth)
        if (uDepthBack) {
            // back-depth pass wants the farthest point that is not clipped
            cameraPos = rayDirection * posT + rayOrigin;
            modelPos = (uInvView * vec4(cameraPos, 1.0)).xyz;
            cameraNormal = -normalize(cameraPos - cameraSpherePos);
            interior = true;
            #if defined(dSolidInterior) && !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
                vec3 frontModelPos = (uInvView * vec4(rayDirection * negT + rayOrigin, 1.0)).xyz;
                float s = clipExit(modelPos / uModelScale, frontModelPos / uModelScale, 0.0);
                if (s < 0.0) return false;
                if (s > 0.0) {
                    cameraPos = rayDirection * mix(posT, negT, s) + rayOrigin;
                    modelPos = (uInvView * vec4(cameraPos, 1.0)).xyz;
                    cameraNormal = -normalize(clipNormal(modelPos / uModelScale) * mat3(uInvView));
                    isCap = true;
                }
            #endif
            fragmentDepth = calcDepth(cameraPos);
            // depth is reversed: in front of the near plane above 1
            return fragmentDepth < 1.0;
        }
    #endif

    cameraPos = rayDirection * negT + rayOrigin;
    modelPos = (uInvView * vec4(cameraPos, 1.0)).xyz;
    fragmentDepth = calcDepth(cameraPos);
    float frontDepth = fragmentDepth;
    vec3 frontModelPos = modelPos;

    bool objectClipped = false;

    #if !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
        if (clipTest(modelPos / uModelScale)) {
            objectClipped = true;
            fragmentDepth = 2.0;
        }
    #endif

    if (fragmentDepth < 1.0) {
        cameraNormal = normalize(cameraPos - cameraSpherePos);
        interior = false;
        return true;
    } else if (uDoubleSided || solidInterior) {
        cameraPos = rayDirection * posT + rayOrigin;
        modelPos = (uInvView * vec4(cameraPos, 1.0)).xyz;
        fragmentDepth = calcDepth(cameraPos);
        cameraNormal = -normalize(cameraPos - cameraSpherePos);
        interior = true;
        if (fragmentDepth < 1.0) {
            #ifdef dSolidInterior
                float nearT = - (uNear + rayOrigin.z) / rayDirection.z;
                float sNear = frontDepth < 1.0 ? 0.0 : clamp((nearT - negT) / (posT - negT), 0.0, 1.0);
                #if !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
                    float sCap = clipExit(frontModelPos / uModelScale, modelPos / uModelScale, sNear);
                    if (sCap < 0.0) return false;
                #else
                    float sCap = sNear;
                #endif
                if (sCap == sNear && frontDepth >= 1.0) {
                    fragmentDepth = 1.0 - (0.0000001 / vRadius);
                    cameraNormal = -mix(normalize(vPoint), vec3(0.0, 0.0, -1.0), uIsOrtho);
                    cameraPos = rayDirection * nearT + rayOrigin;
                    modelPos = (uInvView * vec4(cameraPos, 1.0)).xyz;
                }
                #if !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
                    else {
                        cameraPos = rayDirection * mix(negT, posT, sCap) + rayOrigin;
                        modelPos = (uInvView * vec4(cameraPos, 1.0)).xyz;
                        fragmentDepth = calcDepth(cameraPos) - (0.0000001 / vRadius);
                        cameraNormal = -normalize(clipNormal(modelPos / uModelScale) * mat3(uInvView));
                        isCap = true;
                    }
                #endif
            #endif
            return true;
        }
    }

    return false;
}

void main(void){
    vec3 cameraNormal;
    float fragmentDepth;

    #ifdef dApproximate
        vec3 pointDir = -vPointViewPosition - vPoint;
        if (dot(pointDir, pointDir) > vRadius * vRadius) discard;
        vec3 vViewPosition = -vPointViewPosition;
        fragmentDepth = gl_FragCoord.z;
        #if !defined(dIgnoreLight) || defined(dXrayShaded) || defined(dRenderVariant_tracing)
            pointDir.z -= cos(length(pointDir)) * vRadius * 0.5;
            cameraNormal = -normalize(pointDir);
        #endif
        interior = false;
    #else
        vec3 modelPos;
        vec3 cameraPos;
        bool hit = SphereImpostor(modelPos, cameraPos, cameraNormal, interior, fragmentDepth);
        if (!hit) discard;

        if (fragmentDepth < 0.0) discard;
        if (fragmentDepth > 1.0) discard;

        gl_FragDepthEXT = fragmentDepth;

        vec3 vModelPosition = modelPos;
        vec3 vViewPosition = cameraPos;
    #endif

    ${fade_lod}
    #if !defined(dClipPrimitive) && defined(dClipVariant_pixel) && defined(dClipObjects)
        if (!isCap) {
            ${clip_pixel}
        }
    #endif

    #ifdef dNeedsNormal
        vec3 normal = -cameraNormal;
    #endif

    ${assign_material_color}

    #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing) || defined(dRenderVariant_marking)
        #if defined(dRenderVariant_marking)
        if (material.a < 1.0 && uAlphaThickness > 0.0) {
        #else
        if (uRenderMask == MaskTransparent && uAlphaThickness > 0.0) {
        #endif
            material.a *= min(1.0, vRadius / (uAlphaThickness * uModelScale));
        }
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
