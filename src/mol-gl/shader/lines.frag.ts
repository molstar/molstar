/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const lines_frag = `
precision highp float;
precision highp int;

#include common
#include common_frag_params
#include color_frag_params
#include common_clip

void main(){
    #include fade_lod
    #include clip_pixel

    float fragmentDepth = gl_FragCoord.z;
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
    #elif defined(dRenderVariant_color)
        gl_FragColor = material;

        #include apply_marker_color
        #include apply_fog
        #include wboit_write
        #include dpoit_write
    #endif
}
`;
