/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;
precision highp int;

#pragma glslify: import('./chunks/common-frag-params.glsl')
#pragma glslify: import('./chunks/color-frag-params.glsl')

// uniform vec3 uLightPosition;
uniform vec3 uLightColor;
uniform vec3 uLightAmbient;
uniform mat4 uView;

uniform mat4 uProjection;
// uniform vec3 interiorColor;
// uniform float interiorDarkening;
vec3 interiorColor = vec3(1.0, 0.5, 0.5);
float interiorDarkening = 0.0;

uniform float clipNear;
// uniform float ortho;
float ortho = 0.0;

varying float vRadius;
varying float vRadiusSq;
varying vec3 vPoint;
varying vec3 vPointViewPosition;

#pragma glslify: attenuation = require(./utils/attenuation.glsl)
#pragma glslify: calculateSpecular = require(./utils/phong-specular.glsl)
#pragma glslify: calculateDiffuse = require(./utils/oren-nayar-diffuse.glsl)

const float specularScale = 0.15;
const float shininess = 200.0;
const float roughness = 100.0;
const float albedo = 0.95;

bool flag2 = false;
bool interior = false;
vec3 cameraPos;
vec3 cameraNormal;

// Calculate depth based on the given camera position.
float calcDepth(in vec3 cameraPos){
    vec2 clipZW = cameraPos.z * uProjection[2].zw + uProjection[3].zw;
    return 0.5 + 0.5 * clipZW.x / clipZW.y;
}

float calcClip(in vec3 cameraPos) {
    return dot(vec4(cameraPos, 1.0), vec4(0.0, 0.0, 1.0, clipNear - 0.5));
}

bool Impostor(out vec3 cameraPos, out vec3 cameraNormal){
    vec3 cameraSpherePos = -vPointViewPosition;
    cameraSpherePos.z += vRadius;

    vec3 rayOrigin = mix(vec3(0.0, 0.0, 0.0), vPoint, ortho);
    vec3 rayDirection = mix(normalize(vPoint), vec3(0.0, 0.0, 1.0), ortho);
    vec3 cameraSphereDir = mix(cameraSpherePos, rayOrigin - cameraSpherePos, ortho);

    float B = dot(rayDirection, cameraSphereDir);
    float det = B * B + vRadiusSq - dot(cameraSphereDir, cameraSphereDir);

    if(det < 0.0){
        discard;
        return false;
    }

    float sqrtDet = sqrt(det);
    float posT = mix(B + sqrtDet, B + sqrtDet, ortho);
    float negT = mix(B - sqrtDet, sqrtDet - B, ortho);

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

void main2(void){
    gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
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
                gl_FragDepthEXT = max(0.0, calcDepth(vec3(-(clipNear - 0.5))) + (0.0000001 / vRadius));
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

        vec3 vNormal = cameraNormal;
        vec3 vViewPosition = -cameraPos;

        // determine surface to light direction
        // vec4 viewLightPosition = view * vec4(lightPosition, 1.0);
        // vec3 lightVector = viewLightPosition.xyz - vViewPosition;
        vec3 lightVector = vViewPosition;

        vec3 L = normalize(lightVector); // light direction
        vec3 V = normalize(vViewPosition); // eye direction

        vec3 N = normalize(vNormal);
        #ifdef dDoubleSided
            N = N * (float(gl_FrontFacing) * 2.0 - 1.0);
        #endif

        // compute our diffuse & specular terms
        float specular = calculateSpecular(L, V, N, shininess) * specularScale;
        vec3 diffuse = uLightColor * calculateDiffuse(L, V, N, roughness, albedo);
        vec3 ambient = uLightAmbient;

        // add the lighting
        vec3 finalColor = material.rgb * (diffuse + ambient) + specular;

        // gl_FragColor.rgb = N;
        // gl_FragColor.a = 1.0;
        // gl_FragColor.rgb = vec3(1.0, 0.0, 0.0);
        gl_FragColor.rgb = finalColor;
        gl_FragColor.a = material.a;

        if(interior){
            #ifdef USE_INTERIOR_COLOR
                gl_FragColor.rgb = interiorColor;
            #endif
            gl_FragColor.rgb *= 1.0 - interiorDarkening;
        }

        #pragma glslify: import('./chunks/apply-marker-color.glsl')
        #pragma glslify: import('./chunks/apply-fog.glsl')
    #endif
}