/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const cylinders_frag = `
precision highp float;
precision highp int;

#define bumpEnabled

uniform mat4 uView;

varying mat4 vTransform;
varying vec3 vStart;
varying vec3 vEnd;
varying float vSize;
varying float vCap;

uniform vec3 uCameraDir;
uniform vec3 uCameraPosition;
uniform mat4 uInvView;

#include common
#include common_frag_params
#include color_frag_params
#include light_frag_params
#include common_clip

// adapted from https://www.shadertoy.com/view/4lcSRn
// The MIT License, Copyright 2016 Inigo Quilez
bool CylinderImpostor(
    in vec3 rayOrigin, in vec3 rayDir,
    in vec3 start, in vec3 end, in float radius,
    out vec4 intersection, out bool interior
){
    vec3 ba = end - start;
    vec3 oc = rayOrigin - start;

    float baba = dot(ba, ba);
    float bard = dot(ba, rayDir);
    float baoc = dot(ba, oc);

    float k2 = baba - bard*bard;
    float k1 = baba * dot(oc, rayDir) - baoc * bard;
    float k0 = baba * dot(oc, oc) - baoc * baoc - radius * radius * baba;

    float h = k1 * k1 - k2 * k0;
    if (h < 0.0) return false;

    bool topCap = (vCap > 0.9 && vCap < 1.1) || vCap >= 2.9;
    bool bottomCap = (vCap > 1.9 && vCap < 2.1) || vCap >= 2.9;

    // body outside
    h = sqrt(h);
    float t = (-k1 - h) / k2;
    float y = baoc + t * bard;
    if (y > 0.0 && y < baba) {
        interior = false;
        intersection = vec4(t, (oc + t * rayDir - ba * y / baba) / radius);
        return true;
    }

    if (topCap && y < 0.0) {
        // top cap
        t = -baoc / bard;
        if (abs(k1 + k2 * t) < h) {
            interior = false;
            intersection = vec4(t, ba * sign(y) / baba);
            return true;
        }
    } else if(bottomCap && y >= 0.0) {
        // bottom cap
        t = (baba - baoc) / bard;
        if (abs(k1 + k2 * t) < h) {
            interior = false;
            intersection = vec4(t, ba * sign(y) / baba);
            return true;
        }
    }

    #ifdef dDoubleSided
        // body inside
        h = -h;
        t = (-k1 - h) / k2;
        y = baoc + t * bard;
        if (y > 0.0 && y < baba) {
            interior = true;
            intersection = vec4(t, (oc + t * rayDir - ba * y / baba) / radius);
            return true;
        }

        // TODO: handle inside caps???
    #endif

    return false;
}

void main() {
    #include clip_pixel

    vec3 rayDir = mix(normalize(vModelPosition - uCameraPosition), uCameraDir, uIsOrtho);

    vec4 intersection;
    bool interior;
    bool hit = CylinderImpostor(vModelPosition, rayDir, vStart, vEnd, vSize, intersection, interior);
    if (!hit) discard;

    vec3 vViewPosition = vModelPosition + intersection.x * rayDir;
    vViewPosition = (uView * vec4(vViewPosition, 1.0)).xyz;
    gl_FragDepthEXT = calcDepth(vViewPosition);

    vec3 vModelPosition = (uInvView * vec4(vViewPosition, 1.0)).xyz;

    if (gl_FragDepthEXT < 0.0) discard;
    if (gl_FragDepthEXT > 1.0) discard;

    float fragmentDepth = gl_FragDepthEXT;
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
            mat3 normalMatrix = transpose3(inverse3(mat3(uView)));
            vec3 normal = normalize(normalMatrix * -normalize(intersection.yzw));
            #include apply_light_color
        #endif

        #include apply_interior_color
        #include apply_marker_color
        #include apply_fog
        #include wboit_write
    #endif
}
`;