/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * adapted from three.js (https://github.com/mrdoob/three.js/)
 * which under the MIT License, Copyright © 2010-2019 three.js authors
 */

export default `
uniform float uLightIntensity;
uniform float uAmbientIntensity;
uniform float uReflectivity;
uniform float uMetalness;
uniform float uRoughness;

struct PhysicalMaterial {
	vec3 diffuseColor;
	float specularRoughness;
	vec3 specularColor;
};

struct IncidentLight {
	vec3 color;
	vec3 direction;
};

struct ReflectedLight {
	vec3 directDiffuse;
	vec3 directSpecular;
	vec3 indirectDiffuse;
};

struct GeometricContext {
	vec3 position;
	vec3 normal;
	vec3 viewDir;
};

vec3 F_Schlick(const in vec3 specularColor, const in float dotLH) {
	// Original approximation by Christophe Schlick '94
	// float fresnel = pow( 1.0 - dotLH, 5.0 );
	// Optimized variant (presented by Epic at SIGGRAPH '13)
	// https://cdn2.unrealengine.com/Resources/files/2013SiggraphPresentationsNotes-26915738.pdf
	float fresnel = exp2((-5.55473 * dotLH - 6.98316) * dotLH);
	return (1.0 - specularColor) * fresnel + specularColor;
}

// Moving Frostbite to Physically Based Rendering 3.0 - page 12, listing 2
// https://seblagarde.files.wordpress.com/2015/07/course_notes_moving_frostbite_to_pbr_v32.pdf
float G_GGX_SmithCorrelated(const in float alpha, const in float dotNL, const in float dotNV) {
	float a2 = pow2(alpha);
	// dotNL and dotNV are explicitly swapped. This is not a mistake.
	float gv = dotNL * sqrt(a2 + (1.0 - a2) * pow2(dotNV));
	float gl = dotNV * sqrt(a2 + (1.0 - a2) * pow2(dotNL));
	return 0.5 / max(gv + gl, EPSILON);
}

// Microfacet Models for Refraction through Rough Surfaces - equation (33)
// http://graphicrants.blogspot.com/2013/08/specular-brdf-reference.html
// alpha is "roughness squared" in Disney’s reparameterization
float D_GGX(const in float alpha, const in float dotNH) {
	float a2 = pow2(alpha);
	float denom = pow2(dotNH) * (a2 - 1.0) + 1.0; // avoid alpha = 0 with dotNH = 1
	return RECIPROCAL_PI * a2 / pow2(denom);
}

vec3 BRDF_Diffuse_Lambert(const in vec3 diffuseColor) {
	return RECIPROCAL_PI * diffuseColor;
}

// GGX Distribution, Schlick Fresnel, GGX-Smith Visibility
vec3 BRDF_Specular_GGX(const in IncidentLight incidentLight, const in GeometricContext geometry, const in vec3 specularColor, const in float roughness) {
	float alpha = pow2(roughness); // UE4's roughness
	vec3 halfDir = normalize(incidentLight.direction + geometry.viewDir);

	float dotNL = saturate(dot(geometry.normal, incidentLight.direction));
	float dotNV = saturate(dot(geometry.normal, geometry.viewDir));
	float dotNH = saturate(dot(geometry.normal, halfDir));
	float dotLH = saturate(dot(incidentLight.direction, halfDir));

	vec3 F = F_Schlick(specularColor, dotLH);
	float G = G_GGX_SmithCorrelated(alpha, dotNL, dotNV);
	float D = D_GGX(alpha, dotNH);
	return F * (G * D);
}

// ref: https://www.unrealengine.com/blog/physically-based-shading-on-mobile - environmentBRDF for GGX on mobile
vec3 BRDF_Specular_GGX_Environment(const in GeometricContext geometry, const in vec3 specularColor, const in float roughness) {
	float dotNV = saturate(dot(geometry.normal, geometry.viewDir));
	const vec4 c0 = vec4(-1, -0.0275, -0.572, 0.022);
	const vec4 c1 = vec4(1, 0.0425, 1.04, -0.04);
	vec4 r = roughness * c0 + c1;
	float a004 = min(r.x * r.x, exp2(-9.28 * dotNV)) * r.x + r.y;
	vec2 AB = vec2(-1.04, 1.04) * a004 + r.zw;
	return specularColor * AB.x + AB.y;
}

void RE_Direct_Physical(const in IncidentLight directLight, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight) {
	float dotNL = saturate(dot(geometry.normal, directLight.direction));
    vec3 irradiance = dotNL * directLight.color;
	irradiance *= PI; // punctual light

	reflectedLight.directSpecular += irradiance * BRDF_Specular_GGX(directLight, geometry, material.specularColor, material.specularRoughness);
	reflectedLight.directDiffuse += irradiance * BRDF_Diffuse_Lambert(material.diffuseColor);
}

void RE_IndirectDiffuse_Physical(const in vec3 irradiance, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight) {
	reflectedLight.indirectDiffuse += irradiance * BRDF_Diffuse_Lambert(material.diffuseColor);
}
`;