/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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

void main() {
    #include fade_lod
    #include clip_pixel

    // Workaround for buggy gl_FrontFacing (e.g. on some integrated Intel GPUs)
    vec3 fdx = dFdx(vViewPosition);
    vec3 fdy = dFdy(vViewPosition);
    vec3 faceNormal = normalize(cross(fdx,fdy));
    bool frontFacing = dot(vNormal, faceNormal) > 0.0;

    #if defined(dFlipSided)
        interior = frontFacing;
    #else
        interior = !frontFacing;
    #endif

    float fragmentDepth = gl_FragCoord.z;

    #ifdef dNeedsNormal
        #if defined(dFlatShaded)
            vec3 normal = -faceNormal;
        #else
            vec3 normal = -normalize(vNormal);
            if (uDoubleSided) normal *= float(frontFacing) * 2.0 - 1.0;
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
        #include apply_light_color
        #include apply_interior_color
        #include apply_marker_color

        #if defined(dRenderVariant_color)
            #include apply_fog
            #include wboit_write
            #include dpoit_write
        #elif defined(dRenderVariant_tracing)
            gl_FragData[1] = vec4(normal, emissive);
            gl_FragData[2] = vec4(material.rgb, uDensity);
        #endif
    #endif
}
`;
