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


float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}


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
// simplification from https://catlikecoding.com/unity/tutorials/advanced-rendering/depth-of-field/
void main()
{
    float u = gl_FragCoord.x/uTexSize.x;
    float v = gl_FragCoord.y/uTexSize.y;
    vec2 uv = vec2(u, v);
    vec2 offset = vec2(0.5, 0.5);
    vec4 color = texture2D(tColor, uv);

    float z = getDepthOpaque(uv);
    float viewDist = abs(getViewZ(z));
    // if we want to use a focuspoint.
    //vec3 aposition = screenSpaceToViewSpace(vec3(uv.xy, z), uInvProjection);
    //float focusdepth = getDepthOpaque(offset);
    //vec3 focusPoint = screenSpaceToViewSpace(vec3(offset.xy, focusdepth), uInvProjection);
    //float focusDist = distance(aposition, focusPoint);

    //Bluring
    vec3 blurColor = getBlurredImage1(uv.xy);
    float coc = (viewDist - inFocus) / PPM;  //focus distance, focus range

    coc = clamp(coc, -1.0, 1.0);
    // for debugging the coc
    // color.rgb = (coc < 0.0) ? abs(coc) * vec3(1.0,0.0,0.0) : vec3(coc, coc, coc) ;//mix(color.rgb, blurColor.rgb, abs(coc));
    color.rgb = mix(color.rgb, blurColor.rgb, abs(coc));
    gl_FragColor  = color;
}   
`;
