export const ssaoDepthMerge_frag = `
precision highp float;
precision highp sampler2D;

uniform sampler2D tDepth;
uniform sampler2D tDepthTransparent;
uniform vec2 uTexSize;

#include common

float getDepth(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepth, coords).r;        
    #else
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #endif
}

bool isBackground(const in float depth) {
    return depth > 0.999; // handle precision issues with packed depth
}

vec2 getDepthTransparentWithAlpha(const in vec2 coords) {
    return unpackRGBAToDepthWithAlpha(texture2D(tDepthTransparent, coords));
}

void main() {
    vec2 coords = gl_FragCoord.xy / uTexSize;

    float opaqueDepth = getDepth(coords);
    vec2 transparentDepthWithAlpha = getDepthTransparentWithAlpha(coords);

    bool isBackgroundOpaque = isBackground(opaqueDepth);
    bool isBackgroundTransparent = isBackground(transparentDepthWithAlpha.x);

    float depth = min(opaqueDepth, transparentDepthWithAlpha.x);

    float alpha = transparentDepthWithAlpha.y;
    if (!isBackground(opaqueDepth)) alpha = 1.0;

    gl_FragColor = packDepthWithAlphaToRGBA(depth, alpha);
}
`;