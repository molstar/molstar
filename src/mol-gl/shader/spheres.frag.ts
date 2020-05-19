/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
#include common_clip

uniform mat4 uProjection;

uniform float uClipNear;
uniform float uIsOrtho;

varying float vRadius;
varying float vRadiusSq;
varying vec3 vPoint;
varying vec3 vPointViewPosition;

vec3 cameraPos;
vec3 cameraNormal;

// Calculate depth based on the given camera position.
float calcDepth(const in vec3 cameraPos){
    vec2 clipZW = cameraPos.z * uProjection[2].zw + uProjection[3].zw;
    return 0.5 + 0.5 * clipZW.x / clipZW.y;
}

float calcClip(const in vec3 cameraPos) {
    return dot(vec4(cameraPos, 1.0), vec4(0.0, 0.0, 1.0, uClipNear - 0.5));
}

bool Impostor(out vec3 cameraPos, out vec3 cameraNormal){
    vec3 cameraSpherePos = -vPointViewPosition;
    cameraSpherePos.z += vRadius;

    vec3 rayOrigin = mix(vec3(0.0, 0.0, 0.0), vPoint, uIsOrtho);
    vec3 rayDirection = mix(normalize(vPoint), vec3(0.0, 0.0, 1.0), uIsOrtho);
    vec3 cameraSphereDir = mix(cameraSpherePos, rayOrigin - cameraSpherePos, uIsOrtho);

    float B = dot(rayDirection, cameraSphereDir);
    float det = B * B + vRadiusSq - dot(cameraSphereDir, cameraSphereDir);

    if (det < 0.0){
        discard;
        return false;
    }

    float sqrtDet = sqrt(det);
    float posT = mix(B + sqrtDet, B + sqrtDet, uIsOrtho);
    float negT = mix(B - sqrtDet, sqrtDet - B, uIsOrtho);

    cameraPos = rayDirection * negT + rayOrigin;


    if (calcDepth(cameraPos) <= 0.0) {
        cameraPos = rayDirection * posT + rayOrigin;
        interior = true;
    } else {
        interior = false;
    }

    cameraNormal = normalize(cameraPos - cameraSpherePos);
    cameraNormal *= float(!interior) * 2.0 - 1.0;

    return !interior;
}

void main(void){
    #include clip_pixel

    bool flag = Impostor(cameraPos, cameraNormal);
    #ifndef dDoubleSided
        if (interior)
            discard;
    #endif

    // FIXME not compatible with custom clipping plane
    // Set the depth based on the new cameraPos.
    gl_FragDepthEXT = calcDepth(cameraPos);
    if (!flag && gl_FragDepthEXT >= 0.0) {
        gl_FragDepthEXT = 0.0 + (0.0000001 / vRadius);
    }

    // bugfix (mac only?)
    if (gl_FragDepthEXT < 0.0)
        discard;
    if (gl_FragDepthEXT > 1.0)
        discard;

    #include assign_material_color

    #if defined(dRenderVariant_pick)
        #include check_picking_alpha
        gl_FragColor = material;
    #elif defined(dRenderVariant_depth)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color)
        #ifdef dIgnoreLight
            gl_FragColor = material;
        #else
            vec3 normal = -cameraNormal;
            vec3 vViewPosition = cameraPos;
            #include apply_light_color
        #endif

        #include apply_interior_color
        #include apply_marker_color
        #include apply_fog
    #endif
}
`;