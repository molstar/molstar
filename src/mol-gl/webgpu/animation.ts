/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
export const animationShader = /* wgsl */ `
@group(1) @binding(12) var<storage, read> wiggleWeights: array<u32>;
fn noiseHash(h: f32) -> f32 { return fract(sin(h) * 43758.5453123); }
fn noise(p: vec3f) -> f32 {
    let cell = floor(p); let a = fract(p); let f = a * a * (3.0 - 2.0 * a);
    let n = cell.x + cell.y * 157.0 + cell.z * 113.0;
    return mix(mix(mix(noiseHash(n), noiseHash(n + 1.0), f.x), mix(noiseHash(n + 157.0), noiseHash(n + 158.0), f.x), f.y),
        mix(mix(noiseHash(n + 113.0), noiseHash(n + 114.0), f.x), mix(noiseHash(n + 270.0), noiseHash(n + 271.0), f.x), f.y), f.z);
}
fn fbm(p: vec3f) -> f32 { return 0.5 * noise(p) + 0.25 * noise(p * 2.01) + 0.125 * noise(p * 2.01 * 2.02); }
fn wiggle(position: vec3f, group: f32, instance: u32) -> vec3f {
    if (camera.animation.y == 0.0) { return position; }
    var amplitude = object.wiggle.y;
    if (object.wiggleOverlay.x > 0.0) {
        let index = select(instance * u32(object.clipConfig.w) + u32(group), instance, object.wiggleOverlay.y > 0.0);
        if (index / 4u < arrayLength(&wiggleWeights)) { amplitude += f32((wiggleWeights[index / 4u] >> (index % 4u * 8u)) & 255u) / 255.0 * object.wiggleOverlay.z; }
    }
    if (amplitude <= 0.0 || object.wiggle.x <= 0.0 || object.wiggle.z <= 0.0) { return position; }
    let t = camera.animation.x * object.wiggle.x;
    var seed = position;
    if (object.wiggle.w > 0.0) { seed = fract(sin(group * vec3f(127.1, 269.5, 419.2)) * vec3f(43758.5453, 21639.7182, 32517.3926)) * 1000.0; }
    seed *= object.wiggle.z;
    let delta = vec3f(fbm(vec3f(seed.x, seed.y + t, seed.z)), fbm(vec3f(seed.x + 37.0, seed.y, seed.z + t)), fbm(vec3f(seed.x + t, seed.y + 73.0, seed.z)));
    return position + (delta / 0.4375 - 1.0) * amplitude;
}
fn tumble(transform: mat4x4f, instance: u32) -> mat4x4f {
    if (camera.animation.y == 0.0 || object.tumble.x <= 0.0 || object.tumble.y <= 0.0 || object.tumble.z <= 0.0) { return transform; }
    let amplitude = object.tumble.y / max(object.sphere.w, 1.0);
    let t = camera.animation.x * object.tumble.x;
    let seed = (f32(instance) * 127.1 + f32(object.info.x) * 311.7) * object.tumble.z;
    let angles = (vec3f(fbm(vec3f(seed, t, 0.0)), fbm(vec3f(seed, 0.0, t)), fbm(vec3f(0.0, seed, t))) / 0.4375 - 1.0) * amplitude;
    let c = cos(angles); let s = sin(angles);
    let rotation = mat3x3f(vec3f(c.y*c.z, c.x*s.z+s.x*s.y*c.z, s.x*s.z-c.x*s.y*c.z),
        vec3f(-c.y*s.z, c.x*c.z-s.x*s.y*s.z, s.x*c.z+c.x*s.y*s.z), vec3f(s.y, -s.x*c.y, c.x*c.y));
    let shifted = seed + 31.7;
    let offset = (vec3f(fbm(vec3f(shifted, t, 0.0)), fbm(vec3f(shifted, 0.0, t)), fbm(vec3f(0.0, shifted, t))) / 0.4375 - 1.0) * amplitude;
    let center = mat3x3f(transform[0].xyz, transform[1].xyz, transform[2].xyz) * object.sphere.xyz;
    return mat4x4f(vec4f(rotation * transform[0].xyz, 0.0), vec4f(rotation * transform[1].xyz, 0.0),
        vec4f(rotation * transform[2].xyz, 0.0), vec4f(transform[3].xyz + center - rotation * center + offset, 1.0));
}
`;
