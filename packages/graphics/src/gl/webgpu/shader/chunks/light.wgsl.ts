/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Physically based lighting of the color variants (direct lights, ambient, multi-scattering
 * approximation for indirect specular), cel shading with \`dCelShaded\`; loops over the lights
 * bounded by \`oLightCount\`
 */
export const light_wgsl = /* wgsl */ `
override oLightCount: i32 = 1;

struct PhysicalMaterial {
    diffuseColor: vec3<f32>,
    roughness: f32,
    specularColor: vec3<f32>,
    specularF90: f32,
};

fn BRDF_Lambert(diffuseColor: vec3<f32>) -> vec3<f32> {
    return RECIPROCAL_PI * diffuseColor;
}

fn F_Schlick(f0: vec3<f32>, f90: f32, dotVH: f32) -> vec3<f32> {
    // optimized variant (presented by Epic at SIGGRAPH '13)
    let fresnel = exp2((-5.55473 * dotVH - 6.98316) * dotVH);
    return f0 * (1.0 - fresnel) + (f90 * fresnel);
}

// Moving Frostbite to Physically Based Rendering 3.0 - page 12, listing 2
fn V_GGX_SmithCorrelated(alpha: f32, dotNL: f32, dotNV: f32) -> f32 {
    let a2 = pow2(alpha);
    let gV = dotNL * sqrt(a2 + (1.0 - a2) * pow2(dotNV));
    let gL = dotNV * sqrt(a2 + (1.0 - a2) * pow2(dotNL));
    return 0.5 / max(gV + gL, EPSILON);
}

// Microfacet Models for Refraction through Rough Surfaces - equation (33)
fn D_GGX(alpha: f32, dotNH: f32) -> f32 {
    let a2 = pow2(alpha);
    let denom = pow2(dotNH) * (a2 - 1.0) + 1.0;
    return RECIPROCAL_PI * a2 / pow2(denom);
}

// GGX distribution, Schlick Fresnel, GGX-Smith correlated visibility
fn BRDF_GGX(lightDir: vec3<f32>, viewDir: vec3<f32>, normal: vec3<f32>, f0: vec3<f32>, f90: f32, roughness: f32) -> vec3<f32> {
    let alpha = pow2(roughness);
    let halfDir = normalize(lightDir + viewDir);
    let dotNL = saturate(dot(normal, lightDir));
    let dotNV = saturate(dot(normal, viewDir));
    let dotNH = saturate(dot(normal, halfDir));
    let dotVH = saturate(dot(viewDir, halfDir));
    let F = F_Schlick(f0, f90, dotVH);
    let V = V_GGX_SmithCorrelated(alpha, dotNL, dotNV);
    let D = D_GGX(alpha, dotNH);
    return F * (V * D);
}

// analytical approximation of the DFG LUT, see 'environmentBRDF' of "Physically Based Shading on Mobile"
fn DFGApprox(normal: vec3<f32>, viewDir: vec3<f32>, roughness: f32) -> vec2<f32> {
    let dotNV = saturate(dot(normal, viewDir));
    let c0 = vec4<f32>(-1.0, -0.0275, -0.572, 0.022);
    let c1 = vec4<f32>(1.0, 0.0425, 1.04, -0.04);
    let r = roughness * c0 + c1;
    let a004 = min(r.x * r.x, exp2(-9.28 * dotNV)) * r.x + r.y;
    return vec2<f32>(-1.04, 1.04) * a004 + r.zw;
}

fn perturbNormal(modelPosition: vec3<f32>, viewPosition: vec3<f32>, normal: vec3<f32>, frequency: f32, scale: f32) -> vec3<f32> {
    let t1 = normalize(cross(normal, select(vec3<f32>(0.0, 1.0, 0.0), vec3<f32>(1.0, 0.0, 0.0), abs(normal.x) < 0.9)));
    let t2 = cross(normal, t1);
    let d = select(2.0 * abs(viewPosition.z), 2.0, uGlobal.uIsOrtho == 1.0);
    let e = max(d / (uGlobal.uProjection[1][1] * uGlobal.uDrawingBufferSize.y), 0.01 / frequency);
    let p = modelPosition * frequency;
    let invView = mat3From4(uGlobal.uInvView);
    let h0 = fbm(p);
    let h1 = fbm(p + (invView * t1) * (e * frequency));
    let h2 = fbm(p + (invView * t2) * (e * frequency));
    return normalize(normal - scale * ((h1 - h0) * t1 + (h2 - h0) * t2) / e);
}

/** Lit color of \`material\`, before the exposure */
fn applyLight(material: vec4<f32>, normal: vec3<f32>, viewPosition: vec3<f32>, metalnessIn: f32, roughnessIn: f32, emissive: f32) -> vec4<f32> {
    var metalness = metalnessIn;
    var roughness = roughnessIn;
    #if defined(dCelShaded)
        // clamp to avoid artifacts
        metalness = clamp(metalness, 0.0, 0.99);
        roughness = clamp(roughness, 0.05, 1.0);
    #endif

    let viewDir = normalize(viewPosition);
    var pm: PhysicalMaterial;
    pm.diffuseColor = material.rgb * (1.0 - metalness);
    pm.roughness = min(max(roughness, 0.0525), 1.0);
    pm.specularColor = mix(vec3<f32>(0.04), material.rgb, metalness);
    pm.specularF90 = 1.0;

    var outgoingLight = vec3<f32>(0.0);

    #if defined(dCelShaded)
        for (var i = 0; i < oLightCount; i++) {
            let direction = uGlobal.uLightDirection[i];
            let color = uGlobal.uLightColor[i] * PI; // * PI for punctual light

            let celDiffuse = RECIPROCAL_PI * max(dot(normal, direction), 0.0) * (1.0 - metalness);
            let celSpecular = luminance(saturate(dot(normal, direction)) * BRDF_GGX(direction, viewDir, normal, pm.specularColor, pm.specularF90, roughness));
            let celIntensity = ceil((celDiffuse + celSpecular) * uGlobal.uCelSteps) / uGlobal.uCelSteps;
            outgoingLight += material.rgb * color * celIntensity;
        }
        outgoingLight += pm.diffuseColor * uGlobal.uAmbientColor;
    #else
        var directDiffuse = vec3<f32>(0.0);
        var directSpecular = vec3<f32>(0.0);
        for (var i = 0; i < oLightCount; i++) {
            let direction = uGlobal.uLightDirection[i];
            let color = uGlobal.uLightColor[i] * PI; // * PI for punctual light
            let irradiance = saturate(dot(normal, direction)) * color;
            directSpecular += irradiance * BRDF_GGX(direction, viewDir, normal, pm.specularColor, pm.specularF90, pm.roughness);
            directDiffuse += irradiance * BRDF_Lambert(pm.diffuseColor);
        }

        let irradiance = uGlobal.uAmbientColor * PI; // * PI for punctual light
        var indirectDiffuse = irradiance * BRDF_Lambert(pm.diffuseColor);

        // indirect specular only for metals, multi-scattering (Fdez-Aguera) to preserve energy
        let radiance = uGlobal.uAmbientColor * metalness;
        let iblIrradiance = uGlobal.uAmbientColor * metalness;
        let fab = DFGApprox(normal, viewDir, pm.roughness);
        let singleScattering = pm.specularColor * fab.x + pm.specularF90 * fab.y;
        let Ess = fab.x + fab.y;
        let Ems = 1.0 - Ess;
        let Favg = pm.specularColor + (1.0 - pm.specularColor) * 0.047619; // 1/21
        let Fms = singleScattering * Favg / (1.0 - Ems * Favg);
        let multiScattering = Fms * Ems;
        let cosineWeightedIrradiance = iblIrradiance * RECIPROCAL_PI;
        let diffuse = pm.diffuseColor * (1.0 - (singleScattering + multiScattering));
        let indirectSpecular = radiance * singleScattering + multiScattering * cosineWeightedIrradiance;
        indirectDiffuse += diffuse * cosineWeightedIrradiance;

        outgoingLight = directDiffuse + indirectDiffuse + directSpecular + indirectSpecular;
    #endif
    // prevents black artifacts on specular highlight with transparent background
    outgoingLight = clamp(outgoingLight, vec3<f32>(0.01), vec3<f32>(0.99));
    outgoingLight += material.rgb * emissive;
    return vec4<f32>(outgoingLight, material.a);
}
`;
