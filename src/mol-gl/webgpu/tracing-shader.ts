/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { lightingShader } from './lighting';

/** Native port of shader/illumination/trace.frag.ts, using WebGPU depth and top-origin pixels. */
export const tracingShader = lightingShader + /* wgsl */ `
struct Settings {
    projection: mat4x4f, inverseProjection: mat4x4f,
    viewport: vec4f, dimensions: vec4f, tracing: vec4f, thickness: vec4f, options: vec4f,
};
@group(0) @binding(0) var<uniform> settings: Settings;
@group(0) @binding(2) var shaded: texture_2d<f32>;
@group(0) @binding(3) var normals: texture_2d<f32>;
@group(0) @binding(4) var albedo: texture_2d<f32>;
@group(0) @binding(5) var depth: texture_depth_2d;
@group(0) @binding(6) var backDepth: texture_depth_2d;
@group(0) @binding(7) var previous: texture_2d<f32>;
@group(0) @binding(8) var output: texture_storage_2d<rgba32float, write>;
fn pcg(state: ptr<function, u32>) -> u32 {
    *state = *state * 747796405u + 2891336453u;
    let word = ((*state >> ((*state >> 28u) + 4u)) ^ *state) * 277803737u;
    return (word >> 22u) ^ word;
}
fn random(state: ptr<function, u32>) -> f32 { return f32(pcg(state)) / 4294967296.0; }
fn randomUnit(state: ptr<function, u32>) -> vec3f {
    let z = random(state) * 2.0 - 1.0; let angle = random(state) * 2.0 * PI;
    let radius = sqrt(max(0.0, 1.0 - z * z));
    return vec3f(radius * cos(angle), radius * sin(angle), z);
}
fn coord(p: vec2f) -> vec2i { return vec2i(clamp(p, settings.viewport.xy + 0.5, settings.viewport.xy + settings.viewport.zw - 0.5)); }
fn viewPosition(p: vec2f, d: f32) -> vec3f {
    let uv = (p - settings.viewport.xy) / settings.viewport.zw;
    let v = settings.inverseProjection * vec4f(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0, d * 2.0 - 1.0, 1.0);
    return v.xyz / v.w;
}
fn screenPosition(v: vec3f) -> vec2f {
    let p = settings.projection * vec4f(v, 1.0); let uv = p.xy / p.w;
    return settings.viewport.xy + vec2f(uv.x * 0.5 + 0.5, 0.5 - uv.y * 0.5) * settings.viewport.zw;
}
fn rayOffset(v: vec3f) -> f32 {
    let clip = settings.projection * vec4f(v, 1.0);
    return max(distance(v, viewPosition(screenPosition(v) + vec2f(1.0, 0.0), (clip.z / clip.w) * 0.5 + 0.5)) * 0.1, 0.001 * settings.dimensions.z);
}
struct March { position: vec3f, pixel: vec2i, missed: bool };
fn march(direction: vec3f, thickness: f32, start: vec3f) -> March {
    var position = start; let begin = rayOffset(start); var step = direction * begin;
    let growth = pow(settings.thickness.x / begin, 1.0 / settings.tracing.y);
    var pixel = coord(screenPosition(start));
    for (var i = 1u; i < u32(settings.tracing.y); i++) {
        position += step; let stepZ = abs(step.z); step *= growth;
        pixel = coord(screenPosition(position));
        let d = textureLoad(depth, pixel, 0);
        let z = viewPosition(vec2f(pixel) + 0.5, d).z; let difference = z - position.z;
        if (difference >= 0.0) {
            var t = thickness;
            if (t == 0.0) {
                t = settings.thickness.w;
                if (settings.options.x > 0.0) {
                    let backZ = viewPosition(vec2f(pixel) + 0.5, textureLoad(backDepth, pixel, 0)).z;
                    t = max(settings.thickness.y, (z - backZ) * settings.thickness.z * textureLoad(albedo, pixel, 0).a);
                }
            }
            if (difference < max(t, stepZ)) {
                if (settings.tracing.z > 0.0) {
                    step *= 0.5; position -= step;
                    for (var j = 0u; j < u32(settings.tracing.z); j++) {
                        pixel = coord(screenPosition(position));
                        let sampleZ = viewPosition(vec2f(pixel) + 0.5, textureLoad(depth, pixel, 0)).z;
                        step *= 0.5; position += select(step, -step, sampleZ - position.z >= 0.0);
                    }
                    pixel = coord(screenPosition(position));
                }
                return March(position, pixel, false);
            }
        }
    }
    return March(position, pixel, true);
}
struct Hit { position: vec3f, normal: vec3f, color: vec3f, emission: vec3f, missed: bool };
fn traceRay(position: vec3f, direction: vec3f) -> Hit {
    let hit = march(direction, 0.0, position); let n = textureLoad(normals, hit.pixel, 0);
    let color = textureLoad(albedo, hit.pixel, 0).rgb;
    let emission = select(color * n.a * 2.0, vec3f(0.0), textureLoad(depth, hit.pixel, 0) == 1.0);
    return Hit(hit.position, -n.xyz, color, emission, hit.missed);
}
fn colorForRay(pixel: vec2i, state: ptr<function, u32>) -> vec3f {
    let n = textureLoad(normals, pixel, 0); let material = textureLoad(albedo, pixel, 0).rgb;
    let view = viewPosition(vec2f(pixel) + 0.5, textureLoad(depth, pixel, 0));
    var hit = Hit(view, -n.xyz, textureLoad(shaded, pixel, 0).rgb, material * n.a, all(n.xyz == vec3f(0.0)));
    var previousHit = hit; var throughput = vec3f(1.0); var result = vec3f(0.0);
    var position = view; var direction = normalized(view);
    for (var bounce = 0u; bounce <= u32(settings.tracing.w); bounce++) {
        if (bounce > 0u) { previousHit = hit; hit = traceRay(position, direction); }
        else if (settings.options.y > 0.0) {
            var direct = lighting.ambient.rgb; var full = direct;
            for (var i = 0u; i < u32(lighting.config.x); i++) {
                let light = lighting.lights[i]; let ndotl = clamp(dot(hit.normal, -light.direction.xyz), 0.0, 1.0);
                let irradiance = ndotl * light.color.rgb; full += irradiance;
                if (ndotl > 0.0) {
                    let start = view + hit.normal * (0.0001 * settings.dimensions.z) - light.direction.xyz * random(state) * settings.dimensions.z;
                    let shadow = march(-light.direction.xyz + randomUnit(state) * settings.options.z, settings.options.w, start);
                    if (shadow.missed) { direct += irradiance; }
                }
            }
            hit.color *= direct / max(full, vec3f(0.0001));
        }
        if (hit.missed) {
            var escape = previousHit.color;
            if (bounce > 1u) {
                var irradiance = lighting.ambient.rgb;
                for (var i = 0u; i < u32(lighting.config.x); i++) {
                    let light = lighting.lights[i]; irradiance += clamp(dot(previousHit.normal, -light.direction.xyz), 0.0, 1.0) * light.color.rgb;
                }
                escape = min(previousHit.color * irradiance, vec3f(0.99)) * lighting.config.y;
            }
            result += escape * throughput; break;
        }
        result += hit.emission * throughput;
        position = hit.position + hit.normal * (0.0001 * settings.dimensions.z);
        direction = normalized(hit.normal + randomUnit(state));
        if (bounce == 0u) { continue; }
        throughput *= hit.color;
        let probability = max(throughput.r, max(throughput.g, throughput.b));
        if (random(state) > probability || probability <= 0.0) { break; }
        throughput /= probability;
    }
    return result;
}
@compute @workgroup_size(8, 8) fn trace(@builtin(global_invocation_id) id: vec3u) {
    if (any(id.xy >= vec2u(settings.dimensions.xy))) { return; }
    let pixel = vec2i(id.xy); let p = vec2f(pixel) + 0.5;
    if (any(p < settings.viewport.xy) || any(p >= settings.viewport.xy + settings.viewport.zw) || textureLoad(depth, pixel, 0) == 1.0) {
        textureStore(output, pixel, vec4f(0.0)); return;
    }
    var state = (id.x * 1973u + id.y * 9277u + u32(settings.dimensions.w) * 26699u) | 1u;
    var color = vec3f(0.0);
    for (var ray = 0u; ray < u32(settings.tracing.x); ray++) { color += colorForRay(pixel, &state) / settings.tracing.x; }
    let last = textureLoad(previous, pixel, 0);
    let weight = select(1.0 / (1.0 + 1.0 / max(last.a, 0.000001)), 1.0, settings.dimensions.w < 1.0 || last.a == 0.0);
    textureStore(output, pixel, vec4f(mix(last.rgb, color, weight), weight));
}
`;
