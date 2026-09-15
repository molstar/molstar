/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

export const mesh_frag = `
precision highp float;
precision highp int;

#define bumpEnabled

#include common
#include common_frag_params
#include color_frag_params
#include light_frag_params
#include normal_frag_params
#include common_clip

uniform vec4 uInteriorColor;
uniform vec4 uInteriorSubstance;

#ifdef dSolidInterior
    uniform int uSolidInteriorPass;
    uniform vec4 uSolidInteriorPlane;
    varying vec4 vCapPosition;
#endif

void main() {
    #include fade_lod

    #ifdef dSolidInterior
        if (uSolidInteriorPass == 2) {
            if (dot(uSolidInteriorPlane.xyz, vViewPosition) + uSolidInteriorPlane.w > 0.0) discard;
            gl_FragColor = vec4(0.0);
            return;
        }
        bool capPass = uSolidInteriorPass != 0;
        vec3 viewPosition = vViewPosition;
        vec3 modelPosition = vModelPosition;
        if (capPass) {
            viewPosition = vCapPosition.xyz / vCapPosition.w;
            modelPosition = (uInvView * vec4(viewPosition, 1.0)).xyz;
        }
        #if defined(dClipVariant_pixel) && dClipObjectCount != 0
            if (clipTest(modelPosition / uModelScale)) discard;
        #endif
        vec3 vViewPosition = viewPosition;
        vec3 vModelPosition = modelPosition;
    #else
        #include clip_pixel
    #endif

    interior = !gl_FrontFacing;

    float fragmentDepth = gl_FragCoord.z;

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
                normal = -uSolidInteriorPlane.xyz;
            }
        #endif
    #endif

    #include assign_material_color
    #include check_transparency

    #if defined(dRenderVariant_pick)
        #include check_picking_alpha
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
        gl_FragColor = material;
    #elif defined(dRenderVariant_emissive)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
        #include apply_interior_color
        #include apply_light_color
        #include apply_marker_color

        #if defined(dRenderVariant_color)
            #include apply_fog
            #ifdef dSolidInterior
                if (uSolidInteriorPass == 3) {
                    if (fragmentDepth >= getDepth(gl_FragCoord.xy / uDrawingBufferSize)) discard;
                    return;
                }
            #endif
            #include wboit_write
            #include dpoit_write
        #elif defined(dRenderVariant_tracing)
            gl_FragData[1] = vec4(normal, emissive);
            gl_FragData[2] = vec4(material.rgb, uDensity);
        #endif
    #endif
}
`;
