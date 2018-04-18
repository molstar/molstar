/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

precision highp float;

struct Light {
    vec3 position;
    vec3 color;
    vec3 ambient;
    float falloff;
    float radius;
};

uniform Light light;
uniform mat4 view;

varying vec3 vNormal, vViewPosition, vColor;

float phongSpecular(vec3 lightDirection, vec3 viewDirection, vec3 surfaceNormal, float shininess) {
    //Calculate Phong power
    vec3 R = -reflect(lightDirection, surfaceNormal);
    return pow(max(0.0, dot(viewDirection, R)), shininess);
}

#define PI 3.14159265

float orenNayarDiffuse(vec3 lightDirection, vec3 viewDirection, vec3 surfaceNormal, float roughness, float albedo) {
    float LdotV = dot(lightDirection, viewDirection);
    float NdotL = dot(lightDirection, surfaceNormal);
    float NdotV = dot(surfaceNormal, viewDirection);

    float s = LdotV - NdotL * NdotV;
    float t = mix(1.0, max(NdotL, NdotV), step(0.0, s));

    float sigma2 = roughness * roughness;
    float A = 1.0 + sigma2 * (albedo / (sigma2 + 0.13) + 0.5 / (sigma2 + 0.33));
    float B = 0.45 * sigma2 / (sigma2 + 0.09);

    return albedo * max(0.0, NdotL) * (A + B * s / t) / PI;
}

#pragma glslify: attenuation = require(./attenuation.glsl)

const float specularScale = 0.65;
const float shininess = 30.0;
const float roughness = 5.0;
const float albedo = 0.95;

void main() {
    // determine surface to light direction
    // vec4 lightPosition = view * vec4(light.position, 1.0);
    vec4 lightPosition = vec4(vec3(0.0, 0.0, -10000.0), 1.0);
    vec3 lightVector = lightPosition.xyz - vViewPosition;

    // calculate attenuation
    // float lightDistance = length(lightVector);
    float falloff = 1.0; // attenuation(light.radius, light.falloff, lightDistance);

    vec3 L = normalize(lightVector); // light direction
    vec3 V = normalize(vViewPosition); // eye direction
    vec3 N = normalize(-vNormal); // surface normal

    // compute our diffuse & specular terms
    float specular = phongSpecular(L, V, N, shininess) * specularScale * falloff;
    vec3 diffuse = light.color * orenNayarDiffuse(L, V, N, roughness, albedo) * falloff;
    vec3 ambient = light.ambient;

    // add the lighting
    vec3 color = vColor * (diffuse + ambient) + specular;

    // gl_FragColor.rgb = N;
    // gl_FragColor.rgb = vec3(1.0, 0.0, 0.0);
    gl_FragColor.rgb = color;
    gl_FragColor.a = 1.0;
}