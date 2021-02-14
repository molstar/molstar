export const depthMerge_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tDepthPrimitives;
uniform sampler2D tDepthVolumes;
uniform vec2 uTexSize;

#include common

float getDepth(const in vec2 coords, sampler2D tDepth) {
    #ifdef dPackedDepth
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #else
        return texture2D(tDepth, coords).r;
    #endif
}

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;
    float depth = min(getDepth(coords, tDepthPrimitives), getDepth(coords, tDepthVolumes));
    gl_FragColor = packDepthToRGBA(depth);
}
`;