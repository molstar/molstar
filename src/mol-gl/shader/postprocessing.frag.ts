export default `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tSsaoDepth;
uniform sampler2D tColor;
uniform sampler2D tDepth;
uniform sampler2D tOutlines;
uniform vec2 uTexSize;

uniform float uNear;
uniform float uFar;
uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

uniform float uOcclusionBias;
uniform float uOcclusionRadius;

uniform float uOutlineScale;
uniform float uOutlineThreshold;

uniform float uMaxPossibleViewZDiff;

const vec4 occlusionColor = vec4(0.0, 0.0, 0.0, 1.0);

#include common

float perspectiveDepthToViewZ(const in float invClipZ, const in float near, const in float far) {
	return (near * far) / ((far - near) * invClipZ - far);
}

float orthographicDepthToViewZ(const in float linearClipZ, const in float near, const in float far) {
	return linearClipZ * (near - far) - near;
}

float getViewZ(const in float depth) {
	#if dOrthographic == 1
		return orthographicDepthToViewZ(depth, uNear, uFar);
	#else
		return perspectiveDepthToViewZ(depth, uNear, uFar);
	#endif
}

float getDepth(const in vec2 coords) {
	return unpackRGBAToDepth(texture2D(tDepth, coords));
}

bool isBackground(const in float depth) {
    return depth >= 0.99;
}

float getOutline(const in vec2 coords, out float closestTexel) {
	float backgroundViewZ = uFar + 3.0 * uMaxPossibleViewZDiff;
	vec2 invTexSize = 1.0 / uTexSize;

	float selfDepth = getDepth(coords);
	float selfViewZ = isBackground(selfDepth) ? backgroundViewZ : getViewZ(getDepth(coords));

	float outline = 1.0;
	closestTexel = 1.0;
	for (float y = -uOutlineScale; y <= uOutlineScale; y++) {
		for (float x = -uOutlineScale; x <= uOutlineScale; x++) {
			if (x * x + y * y > uOutlineScale * uOutlineScale) {
				continue;
			}

			vec2 sampleCoords = coords + vec2(x, y) * invTexSize;

			vec4 sampleOutlineCombined = texture2D(tOutlines, sampleCoords);
			float sampleOutline = sampleOutlineCombined.r;
			float sampleOutlineDepth = unpackRGToUnitInterval(sampleOutlineCombined.gb);

			if (sampleOutline == 0.0 && sampleOutlineDepth < closestTexel && abs(selfViewZ - sampleOutlineDepth) > uMaxPossibleViewZDiff) {
				outline = 0.0;
				closestTexel = sampleOutlineDepth;
			}
		}
	}
	return outline;
}

float getSsao(vec2 coords) {
	float rawSsao = unpackRGToUnitInterval(texture(tSsaoDepth, coords).xy);
	if (rawSsao > 0.999) {
		return 1.0;
	} else if (rawSsao > 0.001) {
		return rawSsao;
	}
	return 0.0;
}

void main(void) {
	vec2 coords = gl_FragCoord.xy / uTexSize;
	vec4 color = texture(tColor, coords);

	#ifdef dOutlineEnable
		float closestTexel;
		float outline = getOutline(coords, closestTexel);
    	
		if (outline == 0.0) {
			color.rgb *= outline;
			float viewDist = abs(getViewZ(closestTexel));
			float fogFactor = smoothstep(uFogNear, uFogFar, viewDist);
			if (color.a != 1.0) {
				color.a = 1.0 - fogFactor;
			}
			color.rgb = mix(color.rgb, vec3(1.0), fogFactor);
		}
	#endif

	// occlusion needs to be handled after outline to darken them properly
	#ifdef dOcclusionEnable
		float depth = getDepth(coords);
		if (!isBackground(depth)) {
			float occlusionFactor = getSsao(coords);
			color = mix(occlusionColor, color, occlusionFactor);
		}
	#endif

	gl_FragColor = color;
}
`;