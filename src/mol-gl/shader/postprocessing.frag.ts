/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

export const postprocessing_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tSsaoDepth;
uniform sampler2D tSsaoDepthTransparent;
uniform sampler2D tColor;
uniform sampler2D tTransparentColor;
uniform sampler2D tDepthOpaque;
uniform sampler2D tDepthTransparent;
uniform sampler2D tShadows;
uniform sampler2D tOutlines;
uniform vec2 uTexSize;

uniform float uNear;
uniform float uFar;
uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;
uniform vec3 uOutlineColor;
uniform vec3 uOcclusionColor;
uniform bool uTransparentBackground;
uniform vec2 uOcclusionOffset;

#include common

float getViewZ(const in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

float getDepthOpaque(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepthOpaque, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepthOpaque, coords));
    #endif
}

float getDepthTransparent(const in vec2 coords) {
    return unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords)).x;
}

bool isBackground(const in float depth) {
    return depth > 0.999; // handle depth packing precision issues
}

float getOutline(const in vec2 coords, const in float opaqueDepth, out float closestTexel, out float isTransparent) {
    float backgroundViewZ = 2.0 * uFar;
    vec2 invTexSize = 1.0 / uTexSize;

    float outline = 1.0;
    closestTexel = 1.0;
    isTransparent = 0.0;
    for (int y = -dOutlineScale; y <= dOutlineScale; y++) {
        for (int x = -dOutlineScale; x <= dOutlineScale; x++) {
            if (x * x + y * y > dOutlineScale * dOutlineScale) {
                continue;
            }

            vec2 sampleCoords = coords + vec2(float(x), float(y)) * invTexSize;

            vec4 sampleOutlineCombined = texture2D(tOutlines, sampleCoords);
            float sampleOutline = sampleOutlineCombined.r;
            float sampleOutlineDepth = unpackRGToUnitInterval(sampleOutlineCombined.gb);

            if (sampleOutline == 0.0 && sampleOutlineDepth < closestTexel) {
                outline = 0.0;
                closestTexel = sampleOutlineDepth;
                isTransparent = sampleOutlineCombined.a;
            }
        }
    }
    return closestTexel < opaqueDepth ? outline : 1.0;
}

float getSsao(vec2 coords) {
    float rawSsao = unpackRGToUnitInterval(texture2D(tSsaoDepth, coords).xy);
    if (rawSsao > 0.999) {
        return 1.0;
    } else if (rawSsao > 0.001) {
        return rawSsao;
    }
    // treat values close to 0.0 as errors and return no occlusion
    return 1.0;
}

float getSsaoTransparent(vec2 coords) {
    float rawSsao = unpackRGToUnitInterval(texture2D(tSsaoDepthTransparent, coords).xy);
    if (rawSsao > 0.999) {
        return 1.0;
    } else if (rawSsao > 0.001) {
        return rawSsao;
    }
    // treat values close to 0.0 as errors and return no occlusion
    return 1.0;
}

void main(void) {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    vec4 color = texture2D(tColor, coords);

    float viewDist;
    float fogFactor;
    float opaqueDepth = getDepthOpaque(coords);

    vec4 transparentColor = texture2D(tTransparentColor, coords);
    float alpha = transparentColor.a;

    #if defined(dOcclusionEnable) && (defined(dOcclusionIncludeOpacity) || defined(dOcclusionIncludeTransparency))
        #ifdef dOcclusionSingleDepth
            float occlusionFactor = getSsaoTransparent(coords + uOcclusionOffset);
            float transparentDepth = getDepthTransparent(coords);
            #ifdef dOcclusionIncludeOpacity
                if (!isBackground(opaqueDepth)) {
                    viewDist = abs(getViewZ(opaqueDepth));
                    fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
                    /*float opaqueOcclusionFactor = occlusionFactor;
                    if (alpha != 0.0) {
                        //when transparency is on top of opaque, use alpha as occlusion factor
                        opaqueOcclusionFactor = alpha;
                    }*/
                    if (!uTransparentBackground) {
                        color.rgb = mix(mix(uOcclusionColor, uFogColor, fogFactor), color.rgb, occlusionFactor);
                    } else {
                        color.rgb = mix(uOcclusionColor * (1.0 - fogFactor), color.rgb, occlusionFactor);
                    }
                }
            #endif
            #ifdef dOcclusionIncludeTransparency                
                if (!isBackground(transparentDepth)) {
                    viewDist = abs(getViewZ(transparentDepth));
                    fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
                    transparentColor.rgb = mix(uOcclusionColor * (1.0 - fogFactor), transparentColor.rgb, occlusionFactor);
                }
            #endif
        #else
            #ifdef dOcclusionIncludeOpacity
                if (!isBackground(opaqueDepth)) {
                    viewDist = abs(getViewZ(opaqueDepth));
                    fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
                    float occlusionFactor = getSsao(coords + uOcclusionOffset);
                    
                    if (!uTransparentBackground) {
                        color.rgb = mix(mix(uOcclusionColor, uFogColor, fogFactor), color.rgb, occlusionFactor);
                    } else {
                        color.rgb = mix(uOcclusionColor * (1.0 - fogFactor), color.rgb, occlusionFactor);
                    }
                }
            #endif
            #ifdef dOcclusionIncludeTransparency
                float transparentDepth = getDepthTransparent(coords);
                if (!isBackground(transparentDepth)) {
                    viewDist = abs(getViewZ(transparentDepth));
                    fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
                    float occlusionFactor = getSsaoTransparent(coords + uOcclusionOffset);
                    transparentColor.rgb = mix(uOcclusionColor * (1.0 - fogFactor), transparentColor.rgb, occlusionFactor);
                }
            #endif
        #endif        
    #endif

    #ifdef dShadowEnable
        if (!isBackground(opaqueDepth)) {
            viewDist = abs(getViewZ(opaqueDepth));
            fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
            vec4 shadow = texture2D(tShadows, coords);
            if (!uTransparentBackground) {
                color.rgb = mix(mix(vec3(0), uFogColor, fogFactor), color.rgb, shadow.a);
            } else {
                color.rgb = mix(vec3(0) * (1.0 - fogFactor), color.rgb, shadow.a);
            }
        }
    #endif

    // outline needs to be handled after occlusion and shadow to keep them clean
    #ifdef dOutlineEnable
        float closestTexel;
        float isTransparentOutline;
        float outline = getOutline(coords, opaqueDepth, closestTexel, isTransparentOutline);
        if (outline == 0.0) {
            viewDist = abs(getViewZ(closestTexel));
            fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
            if (isTransparentOutline == 0.0){
                if (!uTransparentBackground) {
                    color.rgb = mix(uOutlineColor, uFogColor, fogFactor);
                } else {
                    color.a = 1.0 - fogFactor;
                    color.rgb = mix(uOutlineColor, color.rgb, fogFactor);
                }
            }            
        }
    #endif

    if (alpha != 0.0) {
        //blending
        color = transparentColor + color * (1.0 -  alpha);
    }
    
    //apply transparent outlines after blending transparent color
    #ifdef dOutlineEnable
        if (isTransparentOutline == 1.0){
            color.a = 1.0 - fogFactor;
            color.rgb = mix(uOutlineColor, color.rgb, fogFactor);
        }
    #endif

    gl_FragColor = color;
}
`;