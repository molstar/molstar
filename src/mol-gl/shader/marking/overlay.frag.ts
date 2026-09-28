export const overlay_frag = `
precision highp float;
precision highp sampler2D;

uniform vec2 uTexSizeInv;
uniform sampler2D tEdgeTexture;
uniform sampler2D tMaskTexture;
uniform vec3 uHighlightEdgeColor;
uniform vec3 uSelectEdgeColor;
uniform float uHighlightEdgeStrength;
uniform float uSelectEdgeStrength;
uniform float uGhostEdgeStrength;
uniform float uInnerEdgeFactor;
uniform vec3 uHighlightFillColor;
uniform vec3 uSelectFillColor;
uniform float uHighlightFillStrength;
uniform float uSelectFillStrength;

void main() {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;

    // solid fill covering the interior of marked regions, sampled directly from the mask
    vec3 fillRgb = vec3(0.0);
    float fillAlpha = 0.0;
    if (uHighlightFillStrength > 0.0 || uSelectFillStrength > 0.0) {
        vec4 m = texture2D(tMaskTexture, coords);
        float coverage = 1.0 - m.r;
        if (coverage > 0.0) {
            vec3 marked = clamp((m.gba - m.r) / max(coverage, 0.001), 0.0, 1.0);
            bool isHighlight = marked.y > 0.5;
            fillRgb = isHighlight ? uHighlightFillColor : uSelectFillColor;
            float fillStrength = isHighlight ? uHighlightFillStrength : uSelectFillStrength;
            // marked.x: 1.0 = occluded (hidden), 0.0 = visible
            fillAlpha = coverage * fillStrength * (marked.x > 0.5 ? uGhostEdgeStrength : 1.0) * marked.z;
        }
    }

    vec4 e = texture2D(tEdgeTexture, coords);
    vec3 edgeRgb = vec3(0.0);
    float edgeAlpha = 0.0;
    if (e.a > 0.0) {
        // the edge is stored premultiplied, see edge.frag
        vec3 edgeValue = e.rgb / e.a;
        bool isHighlight = edgeValue.z > 0.5;
        edgeRgb = isHighlight ? uHighlightEdgeColor : uSelectEdgeColor;
        edgeRgb = edgeValue.y > 0.5 ? edgeRgb : edgeRgb * uInnerEdgeFactor;
        edgeAlpha = (edgeValue.x > 0.5 ? uGhostEdgeStrength : 1.0) * e.a;
        edgeAlpha *= isHighlight ? uHighlightEdgeStrength : uSelectEdgeStrength;
    }

    // premultiplied, so that partially covered edge and fill pixels are not weighted by alpha twice
    gl_FragColor.rgb = edgeRgb * edgeAlpha + fillRgb * fillAlpha * (1.0 - edgeAlpha);
    gl_FragColor.a = edgeAlpha + fillAlpha * (1.0 - edgeAlpha);
}
`;