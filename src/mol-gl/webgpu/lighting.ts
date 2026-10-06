/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * Physical lighting adapted from Mol*'s MIT-licensed three.js-derived shader chunks.
 */
import { Camera } from '../../mol-canvas3d/camera';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Color } from '../../mol-util/color';
import { getLight, getTransformedLightDirection, RendererProps } from '../renderer';

// Keep the lighting block below WebGPU's portable 64 KiB uniform-buffer limit
// while leaving the eight storage bindings available to volume themes. Each
// light occupies two vec4 values, so 1024 lights use about 32 KiB.
export const MaxWebGPULights = 1024;
export const WebGPULightingByteSize = (12 + MaxWebGPULights * 8) * 4;

export function createLightingData(camera: Camera, props: RendererProps, width: number, height: number) {
    if (props.light.length > MaxWebGPULights) throw new Error(`WebGPU supports at most ${MaxWebGPULights} directional lights per frame.`);
    const data = new Float32Array(WebGPULightingByteSize / 4);
    data.set(Color.toRgbNormalized(props.ambientColor).map(c => c * props.ambientIntensity));
    data.set([props.light.length, props.exposure, props.celSteps, props.xrayEdgeFalloff], 4);
    data.set([width, height, props.pickingAlphaThreshold, 0], 8);
    const light = getLight(props.light);
    const direction = Mat4.isZero(camera.headRotation) ? light.direction : getTransformedLightDirection(light, Mat4.invert(Mat4(), camera.headRotation));
    for (let i = 0; i < light.count; i++) {
        data.set(direction.slice(i * 3, i * 3 + 3), 12 + i * 8);
        data.set(light.color.slice(i * 3, i * 3 + 3), 16 + i * 8);
    }
    return data;
}

export const lightingShader = /* wgsl */ `
struct DirectionalLight { direction: vec4f, color: vec4f };
struct Lighting { ambient: vec4f, config: vec4f, dimensions: vec4f, lights: array<DirectionalLight, ${MaxWebGPULights}> };
@group(0) @binding(1) var<uniform> lighting: Lighting;
const PI = 3.141592653589793;
fn normalized(v: vec3f) -> vec3f { return v / max(length(v), 0.000001); }
fn schlick(f0: vec3f, dotVH: f32) -> vec3f {
    let fresnel = exp2((-5.55473 * dotVH - 6.98316) * dotVH);
    return f0 * (1.0 - fresnel) + fresnel;
}
fn ggx(light: vec3f, view: vec3f, normal: vec3f, f0: vec3f, roughness: f32) -> vec3f {
    let alpha = roughness * roughness; let a2 = alpha * alpha;
    let half = normalized(light + view);
    let nl = clamp(dot(normal, light), 0.0, 1.0); let nv = clamp(dot(normal, view), 0.0, 1.0);
    let nh = clamp(dot(normal, half), 0.0, 1.0); let vh = clamp(dot(view, half), 0.0, 1.0);
    let gv = nl * sqrt(a2 + (1.0 - a2) * nv * nv); let gl = nv * sqrt(a2 + (1.0 - a2) * nl * nl);
    let visibility = 0.5 / max(gv + gl, 0.000001);
    let denominator = nh * nh * (a2 - 1.0) + 1.0;
    let distribution = a2 / (PI * denominator * denominator);
    return schlick(f0, vh) * (visibility * distribution);
}
fn dfg(normal: vec3f, view: vec3f, roughness: f32) -> vec2f {
    let nv = clamp(dot(normal, view), 0.0, 1.0);
    let r = roughness * vec4f(-1.0, -0.0275, -0.572, 0.022) + vec4f(1.0, 0.0425, 1.04, -0.04);
    let a004 = min(r.x * r.x, exp2(-9.28 * nv)) * r.x + r.y;
    return vec2f(-1.04, 1.04) * a004 + r.zw;
}
// Normals and view directions follow the original Mol* inward-facing convention.
fn shadeMaterial(base: vec3f, normal: vec3f, viewPosition: vec3f, material: vec4f, emission: f32, ignoreLight: bool, geometryRoughness: f32) -> vec3f {
    if (ignoreLight) { return base * (1.0 + emission) * lighting.config.y; }
    let view = normalized(viewPosition);
    let cel = material.w > 0.0;
    let metalness = clamp(material.x, 0.0, select(1.0, 0.99, cel));
    let roughness = clamp(material.y, select(0.0525, 0.05, cel), 1.0);
    let physicalRoughness = min(roughness + geometryRoughness, 1.0);
    let diffuse = base * (1.0 - metalness); let specular = mix(vec3f(0.04), base, metalness);
    var outgoing = diffuse * lighting.ambient.rgb;
    for (var i = 0u; i < u32(lighting.config.x); i++) {
        let light = lighting.lights[i]; let nl = clamp(dot(normal, light.direction.xyz), 0.0, 1.0);
        let brdf = ggx(light.direction.xyz, view, normal, specular, select(physicalRoughness, roughness, cel));
        if (cel) {
            let diffuseIntensity = nl * (1.0 - metalness) / PI;
            let specularIntensity = dot(nl * brdf, vec3f(0.2125, 0.7154, 0.0721));
            let steps = max(lighting.config.z, 1.0);
            let intensity = ceil((diffuseIntensity + specularIntensity) * steps) / steps;
            outgoing += base * light.color.rgb * PI * intensity;
        } else { outgoing += nl * light.color.rgb * (diffuse + PI * brdf); }
    }
    if (!cel) {
        let fab = dfg(normal, view, physicalRoughness);
        let singleScatter = specular * fab.x + fab.y;
        let ems = 1.0 - fab.x - fab.y;
        let average = specular + (1.0 - specular) * 0.047619;
        let multiScatter = singleScatter * average / max(1.0 - ems * average, vec3f(0.000001)) * ems;
        let irradiance = lighting.ambient.rgb * metalness / PI;
        outgoing += lighting.ambient.rgb * metalness * singleScatter + multiScatter * irradiance;
        outgoing += diffuse * (1.0 - singleScatter - multiScatter) * irradiance;
    }
    return (clamp(outgoing, vec3f(0.01), vec3f(0.99)) + base * emission) * lighting.config.y;
}
`;
