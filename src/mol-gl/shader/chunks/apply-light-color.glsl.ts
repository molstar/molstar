/**
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * adapted from three.js (https://github.com/mrdoob/three.js/)
 * which under the MIT License, Copyright © 2010-2021 three.js authors
 */

export const apply_light_color = `
#ifdef bumpEnabled
    if (uBumpFrequency > 0.0 && uBumpAmplitude > 0.0) {
        vec3 bumpNormal = perturbNormal(-vViewPosition, normal, fbm(vModelPosition * uBumpFrequency), (uBumpAmplitude * bumpiness) / uBumpFrequency);
        #ifdef enabledFragDepth
            if (!any(isNaN(bumpNormal))) normal = bumpNormal;
        #else
            normal = bumpNormal;
        #endif
    }
#endif

vec4 color = material;

ReflectedLight reflectedLight = ReflectedLight(vec3(0.0), vec3(0.0), vec3(0.0), vec3(0.0));

PhysicalMaterial physicalMaterial;
physicalMaterial.diffuseColor = color.rgb * (1.0 - metalness);
#ifdef enabledFragDepth
    physicalMaterial.roughness = min(max(roughness, 0.0525), 1.0);
#else
    vec3 dxy = max(abs(dFdx(normal)), abs(dFdy(normal)));
    float geometryRoughness = max(max(dxy.x, dxy.y), dxy.z);
    physicalMaterial.roughness = min(max(roughness, 0.0525) + geometryRoughness, 1.0);
#endif
physicalMaterial.specularColor = mix(vec3(0.04), color.rgb, metalness);
physicalMaterial.specularF90 = 1.0;

GeometricContext geometry;
geometry.position = -vViewPosition;
geometry.normal = normal;
geometry.viewDir = normalize(vViewPosition);

IncidentLight directLight;
#pragma unroll_loop_start
for (int i = 0; i < dLightCount; ++i) {
    directLight.direction = uLightDirection[i];
    directLight.color = uLightColor[i] * PI; // * PI for punctual light
    RE_Direct_Physical(directLight, geometry, physicalMaterial, reflectedLight);
}
#pragma unroll_loop_end

vec3 irradiance = uAmbientColor * PI; // * PI for punctual light
RE_IndirectDiffuse_Physical(irradiance, geometry, physicalMaterial, reflectedLight);

// indirect specular only metals
vec3 radiance = uAmbientColor * metalness;
vec3 iblIrradiance = uAmbientColor * metalness;
vec3 clearcoatRadiance = vec3(0.0);
RE_IndirectSpecular_Physical(radiance, iblIrradiance, clearcoatRadiance, geometry, physicalMaterial, reflectedLight);

vec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + reflectedLight.directSpecular + reflectedLight.indirectSpecular;

gl_FragColor = vec4(outgoingLight, color.a);

#ifdef dXrayShaded
    gl_FragColor.a *= 1.0 - pow(abs(dot(normal, vec3(0.0, 0.0, 1.0))), uXrayEdgeFalloff);
#endif
`;