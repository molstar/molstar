/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    #include assign_material_color

    #if defined(dRenderVariant_pick)
        #include check_picking_alpha
        gl_FragColor = material;
    #elif defined(dRenderVariant_depth)
        gl_FragColor = material;
    #elif defined(dRenderVariant_marking)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color)
        #ifdef dIgnoreLight
            gl_FragColor = material;
        #else
            #if defined(dFlatShaded)
                vec3 normal = -faceNormal;
            #else
                vec3 normal = -normalize(vNormal);
                #ifdef dDoubleSided
                    normal = normal * (float(frontFacing) * 2.0 - 1.0);
                #endif
            #endif
            #include apply_light_color
        #endif

        #include apply_interior_color
        #include apply_marker_color
        #include apply_fog
        #include wboit_write
    #endif
}
`;