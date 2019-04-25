precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tColor;
uniform sampler2D tDepth;
uniform vec2 uTexSize;

uniform float uOcclusionBias;
uniform float uOcclusionRadius;

uniform float uOutlineScale;
uniform float uOutlineThreshold;

const float noiseAmount = 0.0002;

float noise(vec2 coords) {
	float a = 12.9898;
	float b = 78.233;
	float c = 43758.5453;
	float dt = dot(coords, vec2(a,b));
	float sn = mod(dt, 3.14159);

	return fract(sin(sn) * c);
}

float calcSSAO(in vec2 coords, in float depth) {
	float occlusionFactor = 0.0;

	for (int i = -dOcclusionKernelSize; i <= dOcclusionKernelSize; i++) {
		for (int j = -dOcclusionKernelSize; j <= dOcclusionKernelSize; j++) {
			vec2 coordsDelta = coords + uOcclusionRadius / float(dOcclusionKernelSize) * vec2(float(i) / uTexSize.x, float(j) / uTexSize.y);
            coordsDelta += noiseAmount * (noise(coordsDelta) - 0.5) / uTexSize;
            coordsDelta = clamp(coordsDelta, 0.5 / uTexSize, 1.0 - 1.0 / uTexSize);
			if (texture2D(tDepth, coordsDelta).r < depth) occlusionFactor += 1.0;
		}
	}

	return occlusionFactor / float((2 * dOcclusionKernelSize + 1) * (2 * dOcclusionKernelSize + 1));
}

float calcEdgeDepth(in vec2 coords) {
    vec2 invTexSize = 1.0 / uTexSize;
    float halfScaleFloor = floor(uOutlineScale * 0.5);
    float halfScaleCeil = ceil(uOutlineScale * 0.5);

    vec2 bottomLeftUV = coords - invTexSize * halfScaleFloor;
    vec2 topRightUV = coords + invTexSize * halfScaleCeil;
    vec2 bottomRightUV = coords + vec2(invTexSize.x * halfScaleCeil, -invTexSize.y * halfScaleFloor);
    vec2 topLeftUV = coords + vec2(-invTexSize.x * halfScaleFloor, invTexSize.y * halfScaleCeil);

    float depth0 = texture2D(tDepth, bottomLeftUV).r;
    float depth1 = texture2D(tDepth, topRightUV).r;
    float depth2 = texture2D(tDepth, bottomRightUV).r;
    float depth3 = texture2D(tDepth, topLeftUV).r;

    float depthFiniteDifference0 = depth1 - depth0;
    float depthFiniteDifference1 = depth3 - depth2;

    return sqrt(pow(depthFiniteDifference0, 2.0) + pow(depthFiniteDifference1, 2.0)) * 100.0;
}

void main(void) {
	vec2 coords = gl_FragCoord.xy / uTexSize;
	vec4 color = texture2D(tColor, coords);

	#ifdef dOcclusionEnable
		float depth = texture2D(tDepth, coords).r;
		if (depth != 1.0) {
			float occlusionFactor = calcSSAO(coords, depth);
			color = mix(color, vec4(0.0, 0.0, 0.0, 1.0), uOcclusionBias * occlusionFactor);
		}
	#endif

	#ifdef dOutlineEnable
    	color.rgb *= (step(calcEdgeDepth(coords), uOutlineThreshold));
	#endif

	gl_FragColor = color;
}