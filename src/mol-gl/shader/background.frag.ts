export const background_frag = `
precision mediump float;
precision mediump samplerCube;
precision mediump sampler2D;

#if defined(dVariant_skybox)
    uniform samplerCube tSkybox;
    uniform mat4 uViewDirectionProjectionInverse;
#elif defined(dVariant_image)
    uniform sampler2D tImage;
    uniform vec2 uImageScale;
    uniform vec2 uImageOffset;
#elif defined(dVariant_horizontalGradient) || defined(dVariant_radialGradient)
    uniform vec3 uGradientColorA;
    uniform vec3 uGradientColorB;
    uniform float uGradientRatio;
#endif

uniform vec2 uTexSize;
uniform float uOpacity;
varying vec4 vPosition;

// TODO: add as general pp option to remove banding?
// Iestyn's RGB dither from http://alex.vlachos.com/graphics/Alex_Vlachos_Advanced_VR_Rendering_GDC2015.pdf
vec3 ScreenSpaceDither(vec2 vScreenPos) {
    vec3 vDither = vec3(dot(vec2(171.0, 231.0), vScreenPos.xy));
    vDither.rgb = fract(vDither.rgb / vec3(103.0, 71.0, 97.0));
    return vDither.rgb / 255.0;
}

void main() {
    #if defined(dVariant_skybox)
        vec4 t = uViewDirectionProjectionInverse * vPosition;
        gl_FragColor = textureCube(tSkybox, normalize(t.xyz / t.w));
        gl_FragColor.a = uOpacity;
    #elif defined(dVariant_image)
        vec2 coords = (gl_FragCoord.xy / uImageScale) + uImageOffset;
        gl_FragColor = texture2D(tImage, vec2(coords.x, 1.0 - coords.y));
        gl_FragColor.a = uOpacity;
    #elif defined(dVariant_horizontalGradient)
        float d = (gl_FragCoord.y / uTexSize.y) + 1.0 - (uGradientRatio * 2.0);
        gl_FragColor = vec4(mix(uGradientColorB, uGradientColorA, clamp(d, 0.0, 1.0)), uOpacity);
        gl_FragColor.rgb += ScreenSpaceDither(gl_FragCoord.xy);
    #elif defined(dVariant_radialGradient)
        float d = distance(vec2(0.5), gl_FragCoord.xy / uTexSize) + uGradientRatio - 0.5;
        gl_FragColor = vec4(mix(uGradientColorB, uGradientColorA, 1.0 - clamp(d, 0.0, 1.0)), uOpacity);
        gl_FragColor.rgb += ScreenSpaceDither(gl_FragCoord.xy);
    #endif
}
`;
