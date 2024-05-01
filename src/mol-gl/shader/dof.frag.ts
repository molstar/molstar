/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
*/

export const dof_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

#include common

uniform sampler2D tColor;
uniform sampler2D tDepth;

uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform int blurSize;
uniform float blurSpread;
uniform float inFocus;
uniform float PPM;// = 20.0;

uniform float uNear;// = 0.1; // uNear plane
uniform float uFar;// = 100.0; // uFar plane

uniform mat4 uInvProjection; // Inverse projection
uniform mat4 uProjection; // projection
uniform mat4 uModelView; // Model view matrix

uniform int uMode;  // 0-planar,  1-spherical
uniform vec3 uCenter;  // Center of focus sphere in view space

// Function to convert depth value from depth buffer to view space Z
float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

// Retrieve depth from depth texture
float getDepthOpaque(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepth, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #endif
}

float liuNearizeDepth(vec2 uv)
{
    float z = texture2D(tDepth, uv).x;
    return (2.0 * uNear) / (uFar + uNear - z * (uFar - uNear));   
}

vec3 getBlurredImage1(vec2 coords){
    vec4 blurColor      = vec4(0);
    vec2 texelSize =vec2(1.0/uTexSize.x, 1.0/uTexSize.y);
    for (int x = 0; x < int(blurSize); x++) 
    {
        for (int y = 0; y < int(blurSize); y++) 
        {
            vec2 offset = vec2(float(x) - float(blurSize) / 2.0, float(y) - float(blurSize) / 2.0);
            // float offset    = i - blurD / 2.0;
            vec2 uvPixel    = coords.xy + offset * texelSize * blurSpread;
            blurColor   += vec4(texture2D(tColor, uvPixel).xyz, 1.0);
        }
    }
    blurColor  = blurColor / (float(blurSize)*float(blurSize));

    return blurColor.rgb;
}

// Simple box blur for blurring the image
vec3 getBlurredImage(vec2 coords, int blurSize, float blurSpread) {
    vec4 blurColor = vec4(0.0);
    vec2 texelSize = vec2(1.0 / uTexSize.x, 1.0 / uTexSize.y);
    for (int x = -blurSize / 2; x <= blurSize / 2; ++x) {
        for (int y = -blurSize / 2; y <= blurSize / 2; ++y) {
            vec2 offset = vec2(float(x), float(y)) * texelSize * blurSpread;
            blurColor += texture2D(tColor, coords + offset);
        }
    }
    blurColor /= float((blurSize + 1) * (blurSize + 1));
    return blurColor.rgb;
}

// simplification from https://catlikecoding.com/unity/tutorials/advanced-rendering/depth-of-field/
void main()
{
    vec2 uv = gl_FragCoord.xy / uTexSize;
    vec4 color = texture2D(tColor, uv);
    float depth = getDepthOpaque(uv);
    float viewDist = getViewZ(depth);

    vec4 center = uProjection * vec4(uCenter, 1.0);
    center.xyz = (center.xyz / center.w) * 0.5 + 0.5;
    float cdepth = getDepthOpaque(center.xy);
    float cview = getViewZ(cdepth);
    
    vec3 aposition = screenSpaceToViewSpace(vec3(uv.xy, depth), uInvProjection);
    vec4 focusPoint = uModelView * vec4(uCenter, 1.0);
    float focusDist = length(aposition - focusPoint.xyz); // Adjust the center point if necessary
    // vec3 blurColor = getBlurredImage(uv, blurSize, blurSpread);
    vec3 blurColor = getBlurredImage1(uv);
    float coc = 0.0;  // Circle of Confusion

    // planar Depth of field
    if (uMode == 0)
    {
        coc = (abs(viewDist) - inFocus) / PPM;  //focus distance, focus range
    }
    // spherical Depth of field
    else if(uMode == 1)
    {
        coc = focusDist / PPM ;
    }
    coc = clamp(coc, -1.0, 1.0);
    // for debugging the coc
    // color.rgb = (coc < 0.0) ? (1.0 - abs(coc)) * vec3(1.0,0.0,0.0) : vec3(0.0, 1.0 - coc, 0.0) ;//mix(color.rgb, blurColor.rgb, abs(coc));
    color.rgb = mix(color.rgb, blurColor, smoothstep(0.0, 1.0, abs(coc)));  // Smooth blending based on CoC
    gl_FragColor  = color;
}   
`;
