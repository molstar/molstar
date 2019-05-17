/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    #include assign_material_color

    #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
        if (uAlpha < uPickingAlphaThreshold)
            discard; // ignore so the element below can be picked
        gl_FragColor = material;
    #elif defined(dColorType_depth)
        gl_FragColor = material;
    #else
        #include assign_normal
        #include apply_light_color
        #include apply_marker_color
        #include apply_fog
    #endif
}
`