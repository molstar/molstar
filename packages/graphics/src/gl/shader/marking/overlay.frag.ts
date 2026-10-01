export const overlay_frag = `
precision highp float;
precision highp sampler2D;

uniform vec2 uTexSizeInv;
uniform sampler2D tEdgeTexture;
uniform vec3 uHighlightEdgeColor;
uniform vec3 uSelectEdgeColor;
uniform float uHighlightEdgeStrength;
uniform float uSelectEdgeStrength;
uniform float uGhostEdgeStrength;
uniform float uInnerEdgeFactor;

void main() {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;
    vec4 edgeValue = texture2D(tEdgeTexture, coords);
    if (edgeValue.a > 0.0) {
        vec3 edgeColor = edgeValue.b == 1.0 ? uHighlightEdgeColor : uSelectEdgeColor;
        gl_FragColor.rgb = edgeValue.g > 0.0 ? edgeColor : edgeColor * uInnerEdgeFactor;
        gl_FragColor.a = (edgeValue.r == 1.0 ? uGhostEdgeStrength : 1.0) * edgeValue.a;
        float edgeStrength = edgeValue.b == 1.0 ? uHighlightEdgeStrength : uSelectEdgeStrength;
        gl_FragColor.a *= edgeStrength;
    } else {
        gl_FragColor = vec4(0.0);
    }
}
`;