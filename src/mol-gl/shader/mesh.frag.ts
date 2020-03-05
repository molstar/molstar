/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export default `
precision highp float;
precision highp int;

#include common
#include common_frag_params
#include color_frag_params
#include light_frag_params
#include normal_frag_params

void main() {
    // Workaround for buggy gl_FrontFacing (e.g. on some integrated Intel GPUs)
    #if defined(enabledStandardDerivatives)
        vec3 fdx = dFdx(vViewPosition);
        vec3 fdy = dFdy(vViewPosition);
        vec3 faceNormal = normalize(cross(fdx,fdy));
        bool frontFacing = dot(vNormal, faceNormal) > 0.0;
    #else
        bool frontFacing = dot(vNormal, vViewPosition) < 0.0;
    #endif

    interior = !frontFacing; // TODO take dFlipSided into account

    #include assign_material_color

    #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
        #include check_picking_alpha
        gl_FragColor = material;
    #elif defined(dColorType_depth)
        gl_FragColor = material;
    #else
        #ifdef dIgnoreLight
            gl_FragColor = material;
        #else
            #include assign_normal
            #include apply_light_color
        #endif

        #include apply_interior_color
        #include apply_marker_color
        #include apply_fog
    #endif
}
`