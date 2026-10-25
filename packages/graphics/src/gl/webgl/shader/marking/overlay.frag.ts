import { common } from '../chunks/common.glsl.js';

export const overlay_frag = `
precision highp float;

${common}

// occlusion at \`coords\`, depth-aware upsampled from a lower SSAO resolution with the scene's opaque depth
float getOcclusion(const in vec2 coords) {
    #ifdef dOcclusionUpsample
        #ifdef dPackedDepth
            float depth = unpackRGBAToDepth(texture2D(tDepth, coords));
        #else
            float depth = texture2D(tDepth, coords).r;
        #endif
        float viewZ = depthToViewZ(uIsOrtho, depth, uNear, uFar);
        return getSsaoUpsampled(tSsaoDepth, uOcclusionTexSize, coords + uOcclusionOffset, viewZ, uNear, uFar);
    #else
        return unpackSsao(texture2D(tSsaoDepth, coords + uOcclusionOffset));
    #endif
}

// mixes toward the SSAO occlusion color, like the postprocessing pass
vec3 applyOcclusion(const in vec3 color, const in float occlusion, const in float fogAlpha) {
    vec3 occlusionColor = uTransparentBackground ? uOcclusionColor * fogAlpha : mix(uFogColor, uOcclusionColor, fogAlpha);
    return mix(occlusionColor, color, occlusion);
}

// the path traced occlusion and shadows, i.e. the base image relative to the direct shading
float getTracedShading(const in float baseLuminance, const in float shadedLuminance, const in float fogAlpha) {
    float a = max(fogAlpha, 0.001);
    // undo the fog of the base image
    float l = uTransparentBackground ? baseLuminance / a : (baseLuminance - luminance(uFogColor) * (1.0 - fogAlpha)) / a;
    return clamp(l / max(shadedLuminance, 0.001), 0.0, 1.0);
}

// edge of the marked regions at \`coords\` (alpha 0 if none): the mask is antialiased, so the coverage can be taken from
// single texels, \`uEdgeScale\` texels apart; rgb: hidden, outer (the mask), highlight
vec4 getEdge(const in vec2 coords) {
    vec4 offset = vec4(uEdgeScale, 0.0, 0.0, uEdgeScale) * vec4(uTexSizeInv, uTexSizeInv);

    vec4 m1 = texture2D(tMaskTexture, coords + offset.xy);
    vec4 m2 = texture2D(tMaskTexture, coords - offset.xy);
    vec4 m3 = texture2D(tMaskTexture, coords + offset.yw);
    vec4 m4 = texture2D(tMaskTexture, coords - offset.yw);

    float edge = max(abs(m1.r - m2.r), abs(m3.r - m4.r));
    if (edge <= 0.0) return vec4(0.0);

    // the sample covering the most of the marked region, mask values are 0.0 where marked
    vec4 mx = m1.r < m2.r ? m1 : m2;
    vec4 my = m3.r < m4.r ? m3 : m4;
    vec4 m = mx.r < my.r ? mx : my;

    // the samples average marked and unmarked texels and unmarked ones are 1.0 in every channel,
    // so the values of the marked texels can be recovered from how much they are covered
    vec3 marked = clamp((m.gba - m.r) / max(1.0 - m.r, 0.001), 0.0, 1.0);

    float visibility = marked.x > 0.5 ? 1.0 : 0.0;
    float mask = texture2D(tMaskTexture, coords).r;
    float marker = marked.y > 0.5 ? 1.0 : 0.0;
    float fogAlpha = marked.z;
    return vec4(visibility, mask, marker, edge * fogAlpha);
}

void main() {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;

    // coverage of the visible unmarked geometry, see dim.frag
    float dimAlpha = 0.0;
    float dimFogAlpha = 1.0;
    if (uDimStrength > 0.0) {
        vec4 d = texture2D(tDimTexture, coords);
        dimAlpha = d.r * uDimStrength;
        dimFogAlpha = d.r / max(d.g, 0.001);
    }

    // solid fill covering the interior of marked regions, sampled directly from the mask
    vec3 fillRgb = vec3(0.0);
    float fillAlpha = 0.0;
    float fillFogAlpha = 1.0;
    if (uHighlightFillStrength > 0.0 || uSelectFillStrength > 0.0) {
        vec4 m = texture2D(tMaskTexture, coords);
        float coverage = 1.0 - m.r;
        if (coverage > 0.0) {
            vec3 marked = clamp((m.gba - m.r) / max(coverage, 0.001), 0.0, 1.0);
            bool isHighlight = marked.y > 0.5;
            fillRgb = isHighlight ? uHighlightFillColor : uSelectFillColor;
            float fillStrength = isHighlight ? uHighlightFillStrength : uSelectFillStrength;
            fillFogAlpha = marked.z;
            // marked.x: 1.0 = hidden, 0.0 = visible; without depth test every marked texel reads as hidden
            float visibility = uDepthTest && marked.x > 0.5 ? 0.0 : 1.0;
            fillAlpha = coverage * fillStrength * visibility * fillFogAlpha * unpackMarkingOpacity(marked.x);
        }
    }

    // keep the shading of the image, so that blending toward the color does not wash it out
    vec3 dimRgb = uDimColor;
    if (dimAlpha > 0.0 || fillAlpha > 0.0) {
        #if defined(dMarkingShading_ssao)
            float occlusion = getOcclusion(coords);
            dimRgb = applyOcclusion(dimRgb, occlusion, dimFogAlpha);
            fillRgb = applyOcclusion(fillRgb, occlusion, fillFogAlpha);
        #elif defined(dMarkingShading_traced)
            float baseLuminance = luminance(texture2D(tBase, coords).rgb);
            float shadedLuminance = luminance(texture2D(tShaded, coords).rgb);
            dimRgb *= getTracedShading(baseLuminance, shadedLuminance, dimFogAlpha);
            fillRgb *= getTracedShading(baseLuminance, shadedLuminance, fillFogAlpha);
        #endif
    }

    vec4 e = getEdge(coords);
    vec3 edgeRgb = vec3(0.0);
    float edgeAlpha = 0.0;
    if (e.a > 0.0) {
        vec3 edgeValue = e.rgb;
        bool isHighlight = edgeValue.z > 0.5;
        edgeRgb = isHighlight ? uHighlightEdgeColor : uSelectEdgeColor;
        edgeRgb = edgeValue.y > 0.5 ? edgeRgb : edgeRgb * uInnerEdgeFactor;
        edgeAlpha = (edgeValue.x > 0.5 ? uGhostEdgeStrength : 1.0) * e.a;
        edgeAlpha *= isHighlight ? uHighlightEdgeStrength : uSelectEdgeStrength;
    }

    // premultiplied, so that partially covered pixels are not weighted by alpha twice; dim, fill, edge from bottom to top
    vec3 rgb = fillRgb * fillAlpha + dimRgb * dimAlpha * (1.0 - fillAlpha);
    float alpha = fillAlpha + dimAlpha * (1.0 - fillAlpha);
    gl_FragColor.rgb = edgeRgb * edgeAlpha + rgb * (1.0 - edgeAlpha);
    gl_FragColor.a = edgeAlpha + alpha * (1.0 - edgeAlpha);
}
`;
