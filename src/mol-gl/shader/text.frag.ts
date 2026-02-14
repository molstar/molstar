/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

export const text_frag = `
precision highp float;
precision highp int;

#include common
#include common_frag_params
#include color_frag_params
#include common_clip

uniform sampler2D tFont;

uniform vec3 uBorderColor;
uniform float uBorderWidth;
uniform vec3 uBackgroundColor;
uniform float uBackgroundOpacity;

varying vec2 vTexCoord;

void main(){
    #include fade_lod
    #include clip_pixel

    // discard label when anchor approaches the near clip plane
    float anchorDepth = vViewPosition.z;
    float nearFadeEnd = uNear + (uFar - uNear) * 0.03;
    if (anchorDepth < nearFadeEnd) discard;

    float fragmentDepth = gl_FragCoord.z;

    // handle depth variant before assign_material_color,
    // since the chunk classifies opaque/transparent based on uAlpha alone
    // and doesn't know about uBackgroundOpacity
    #if defined(dRenderVariant_depth)
    {
        if (fragmentDepth > getDepth(gl_FragCoord.xy / uDrawingBufferSize)) {
            discard;
        }
        if (vTexCoord.x > 1.0) {
            discard; // background is cosmetic, skip depth to avoid outline artifacts
        }
        float rawSdf = texture2D(tFont, vTexCoord).a;
        float sdf = rawSdf + min(uBorderWidth, 0.49); // clamp to avoid exceeding max SDF range
        if (sdf < 0.5) discard;

        vec4 material;
        if (uRenderMask == MaskOpaque) {
            if (uAlpha < 1.0) {
                discard;
            }
            material = packDepthToRGBA(fragmentDepth);
        } else if (uRenderMask == MaskTransparent) {
            if (uAlpha == 1.0) {
                discard;
            }
            material = packDepthWithAlphaToRGBA(fragmentDepth, uAlpha);
        }
        gl_FragColor = material;
        return;
    }
    #endif

    #ifdef enabledFragDepth
        gl_FragDepthEXT = fragmentDepth;
    #endif

    #include assign_material_color

    if (vTexCoord.x > 1.0) {
        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            material = vec4(uBackgroundColor, uBackgroundOpacity * material.a);
        #endif
    } else {
        // retrieve signed distance
        float rawSdf = texture2D(tFont, vTexCoord).a;
        float borderWidth = min(uBorderWidth, 0.49);
        float sdf = rawSdf + borderWidth;

        if (sdf < 0.5) discard;

        #if defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
            if (uBorderWidth > 0.0 && rawSdf < 0.5) {
                material.xyz = uBorderColor;
            } else {
                // push text fragments forward in depth so they render in front of border
                #ifdef enabledFragDepth
                    gl_FragDepthEXT = fragmentDepth - 0.0001;
                #endif
            }
        #endif
    }

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
    #elif defined(dRenderVariant_emissive)
        gl_FragColor = material;
    #elif defined(dRenderVariant_color) || defined(dRenderVariant_tracing)
        gl_FragColor = material;
        #include apply_marker_color

        #if defined(dRenderVariant_color)
            #include apply_fog
            #include wboit_write
            #include dpoit_write
        #elif defined(dRenderVariant_tracing)
            gl_FragData[1] = vec4(-normalize(vViewPosition), emissive);
            gl_FragData[2] = vec4(material.rgb, uDensity);
        #endif
    #endif
}
`;
