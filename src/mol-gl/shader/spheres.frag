/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-frag-params.glsl')
#pragma glslify: import('./chunks/color-frag-params.glsl')
#pragma glslify: import('./chunks/light-frag-params.glsl')

uniform mat4 uProjection;
// uniform vec3 uInteriorColor;
// uniform float uInteriorDarkening;
vec3 uInteriorColor = vec3(0.5, 0.5, 0.5);
float uInteriorDarkening = 0.0;

uniform float uClipNear;
uniform float uIsOrtho;

varying float vRadius;
varying float vRadiusSq;
varying vec3 vPoint;
varying vec3 vPointViewPosition;

bool flag2 = false;
bool interior = false;
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

    if(det < 0.0){
        discard;
        return false;
    }

    float sqrtDet = sqrt(det);
    float posT = mix(B + sqrtDet, B + sqrtDet, uIsOrtho);
    float negT = mix(B - sqrtDet, sqrtDet - B, uIsOrtho);

    cameraPos = rayDirection * negT + rayOrigin;

    #ifdef NEAR_CLIP
        if(calcDepth(cameraPos) <= 0.0){
            cameraPos = rayDirection * posT + rayOrigin;
            interior = true;
        }else if(calcClip(cameraPos) > 0.0){
            cameraPos = rayDirection * posT + rayOrigin;
            interior = true;
            flag2 = true;
        }
    #else
        if(calcDepth(cameraPos) <= 0.0){
            cameraPos = rayDirection * posT + rayOrigin;
            interior = true;
        }
    #endif

    cameraNormal = normalize(cameraPos - cameraSpherePos);
    cameraNormal *= float(!interior) * 2.0 - 1.0;

    return !interior;
}

void main(void){
    bool flag = Impostor(cameraPos, cameraNormal);

    #ifdef NEAR_CLIP
        if(calcClip(cameraPos) > 0.0)
            discard;
    #endif

    // FIXME not compatible with custom clipping plane
    // Set the depth based on the new cameraPos.
    gl_FragDepthEXT = calcDepth(cameraPos);
    if(!flag){
        // clamp to near clipping plane and add a tiny value to
        // make spheres with a greater radius occlude smaller ones
        #ifdef NEAR_CLIP
            if( flag2 ){
                gl_FragDepthEXT = max(0.0, calcDepth(vec3(-(uClipNear - 0.5))) + (0.0000001 / vRadius));
            }else if(gl_FragDepthEXT >= 0.0){
                gl_FragDepthEXT = 0.0 + (0.0000001 / vRadius);
            }
        #else
            if(gl_FragDepthEXT >= 0.0){
                gl_FragDepthEXT = 0.0 + (0.0000001 / vRadius);
            }
        #endif
    }

    // bugfix (mac only?)
    if (gl_FragDepthEXT < 0.0)
        discard;
    if (gl_FragDepthEXT > 1.0)
        discard;

    // material color
    #pragma glslify: import('./chunks/assign-material-color.glsl')

    #if defined(dColorType_objectPicking) || defined(dColorType_instancePicking) || defined(dColorType_groupPicking)
        if (uAlpha < uPickingAlphaThreshold)
            discard; // ignore so the element below can be picked
        gl_FragColor = material;
    #else
        vec3 normal = cameraNormal;
        vec3 vViewPosition = -cameraPos;
        #pragma glslify: import('./chunks/apply-light-color.glsl')

        if(interior){
            #ifdef USE_INTERIOR_COLOR
                gl_FragColor.rgb = uInteriorColor;
            #endif
            gl_FragColor.rgb *= 1.0 - uInteriorDarkening;
        }

        #pragma glslify: import('./chunks/apply-marker-color.glsl')
        #pragma glslify: import('./chunks/apply-fog.glsl')
    #endif
}