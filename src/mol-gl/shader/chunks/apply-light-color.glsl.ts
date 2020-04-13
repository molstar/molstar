/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * adapted from three.js (https://github.com/mrdoob/three.js/)
 * which under the MIT License, Copyright © 2010-2019 three.js authors
 */

export default `
// inputs
// - vec4 material
// - vec3 vViewPosition
// - vec3 normal
// - float uMetalness
// - float uRoughness
// - float uReflectivity
// - float uLightIntensity
// - float uAmbientIntensity

// outputs
// - sets gl_FragColor

vec4 color = material;

ReflectedLight reflectedLight = ReflectedLight(vec3(0.0), vec3(0.0), vec3(0.0));

PhysicalMaterial physicalMaterial;
physicalMaterial.diffuseColor = color.rgb * (1.0 - uMetalness);
physicalMaterial.specularRoughness = clamp(uRoughness, 0.04, 1.0);
physicalMaterial.specularColor = mix(vec3(0.16 * pow2(uReflectivity)), color.rgb, uMetalness);

GeometricContext geometry;
geometry.position = -vViewPosition;
geometry.normal = normal;
geometry.viewDir = normalize(vViewPosition);

IncidentLight directLight;
directLight.direction = vec3(0.0, 0.0, -1.0);
directLight.color = vec3(uLightIntensity);

RE_Direct_Physical(directLight, geometry, physicalMaterial, reflectedLight);

vec3 irradiance = vec3(uAmbientIntensity) * PI;
RE_IndirectDiffuse_Physical(irradiance, geometry, physicalMaterial, reflectedLight);

vec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + reflectedLight.directSpecular;

gl_FragColor = vec4(outgoingLight, color.a);
`;