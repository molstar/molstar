export const blur_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tInput;
uniform vec2 uTexSizeInv;

uniform vec2 uDirection;
uniform float uGaussianCoefficients[dKernelRadius];

void main(void) {
    vec2 coords = gl_FragCoord.xy * uTexSizeInv;
    float weightSum = uGaussianCoefficients[0];
    vec4 diffuseSum = texture2D(tInput, coords) * weightSum;

    for(int i = 1; i < dKernelRadius; ++i) {
        float x = float(i);
        float w = uGaussianCoefficients[i];
        vec2 offset = uDirection * uTexSizeInv * x;
        vec4 sample1 = texture2D(tInput, coords + offset);
        vec4 sample2 = texture2D(tInput, coords - offset);
        diffuseSum += (sample1 + sample2) * w;
        weightSum += 2.0 * w;
    }

    gl_FragColor = diffuseSum / weightSum;
}
`;