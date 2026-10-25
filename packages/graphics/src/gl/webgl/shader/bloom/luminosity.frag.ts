import { common } from '../chunks/common.glsl.js';

export const luminosity_frag = `
precision highp float;
precision highp int;

${common}

float getDepthOpaque(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepthOpaque, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepthOpaque, coords));
    #endif
}

float getDepth(const in vec2 coords) {
    #ifdef dBlendTransparency
        return max(getDepthOpaque(coords), unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords)).x);
    #else
        return getDepthOpaque(coords);
    #endif
}

bool isBackground(const in float depth) {
    // the far plane (depth is reversed), packed or raw
    return depth == 0.0;
}

void main(void) {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;
    vec4 opaqueTexel = texture2D(tColorOpaque, coords);
    #ifdef dBlendTransparency
        vec4 transparentTexel = texture2D(tColorTransparent, coords);
    #else
        vec4 transparentTexel = vec4(0.0);
    #endif
    // fog illumination's un-fogged opaque seed by the opaque depth, matching standard's per-layer fog
    if (!uOpaqueFogged) {
        opaqueTexel.rgb *= 1.0 - smoothstep(uFogNear, uFogFar, abs(depthToViewZ(uIsOrtho, getDepthOpaque(coords), uNear, uFar)));
    }
    // PMA OVER composite, matches postprocessing's final transparency blend.
    vec4 texel = transparentTexel + opaqueTexel * (1.0 - transparentTexel.a);
    float emissive = texture2D(tEmissive, coords).a;
    float depth = getDepth(coords);

    if (isBackground(depth)) {
        gl_FragColor = vec4(0.0, 0.0, 0.0, 0.0);
        return;
    }

    vec4 outputColor = vec4(uDefaultColor.rgb, uDefaultOpacity);

    #if defined(dMode_luminosity)
        vec3 luma = vec3(0.299, 0.587, 0.114);
        float v = dot(texel.xyz, luma);
        float alpha = smoothstep(uLuminosityThreshold, uLuminosityThreshold + uSmoothWidth, v);

        gl_FragColor = mix(outputColor, texel, alpha);
    #elif defined(dMode_emissive)
        // the prepass already holds the fogged emissive color, free of lighting/exposure
        gl_FragColor = texture2D(tEmissive, coords);
    #endif
}
`;
