/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Wiggle (per position or group) and tumble (per instance) animations; expects `fbm` */
export const animation_wgsl = /* wgsl */ `
fn applyWiggle(pos: vec3<f32>, groupId: f32, wiggle: f32) -> vec3<f32> {
    if (uGlobal.uEnableAnimation == 0u) {
        return pos;
    }
    let amplitude = uMaterial.uWiggleAmplitude + wiggle;
    if (amplitude > 0.0 && uMaterial.uWiggleSpeed > 0.0 && uMaterial.uWiggleFrequency > 0.0) {
        let t = uGlobal.uTime * uMaterial.uWiggleSpeed;
        var s: vec3<f32>;
        if (uMaterial.uWiggleMode == 0) {
            // position mode: spatial position correlates nearby atoms
            s = pos;
        } else {
            // group mode: per-group independent noise
            s = vec3<f32>(
                fract(sin(groupId * 127.1) * 43758.5453) * 1000.0,
                fract(sin(groupId * 269.5) * 21639.7182) * 1000.0,
                fract(sin(groupId * 419.2) * 32517.3926) * 1000.0
            );
        }
        s *= uMaterial.uWiggleFrequency;
        return pos + vec3<f32>(
            fbm(vec3<f32>(s.x, s.y + t, s.z)) / 0.4375 - 1.0,
            fbm(vec3<f32>(s.x + 37.0, s.y, s.z + t)) / 0.4375 - 1.0,
            fbm(vec3<f32>(s.x + t, s.y + 73.0, s.z)) / 0.4375 - 1.0
        ) * amplitude;
    }
    return pos;
}

fn applyTumble(transform: mat4x4<f32>, instanceIndex: f32, objectId: f32, sphere: vec4<f32>) -> mat4x4<f32> {
    if (uGlobal.uEnableAnimation == 0u) {
        return transform;
    }
    if (uMaterial.uTumbleAmplitude > 0.0 && uMaterial.uTumbleSpeed > 0.0 && uMaterial.uTumbleFrequency > 0.0) {
        // scale amplitude inversely with bounding-sphere radius (Stokes-Einstein: D ~ 1/r)
        let amplitude = uMaterial.uTumbleAmplitude / max(sphere.w, 1.0);
        let t = uGlobal.uTime * uMaterial.uTumbleSpeed;
        let seed = (instanceIndex * 127.1 + objectId * 311.7) * uMaterial.uTumbleFrequency;

        // per-instance rotation angles from layered noise (Brownian-like)
        let angleX = (fbm(vec3<f32>(seed, t, 0.0)) / 0.4375 - 1.0) * amplitude;
        let angleY = (fbm(vec3<f32>(seed, 0.0, t)) / 0.4375 - 1.0) * amplitude;
        let angleZ = (fbm(vec3<f32>(0.0, seed, t)) / 0.4375 - 1.0) * amplitude;

        let cx = cos(angleX); let sx = sin(angleX);
        let cy = cos(angleY); let sy = sin(angleY);
        let cz = cos(angleZ); let sz = sin(angleZ);

        // combined rotation matrix (Rz * Ry * Rx)
        let rot = mat3x3<f32>(
            cy * cz, cx * sz + sx * sy * cz, sx * sz - cx * sy * cz,
            -cy * sz, cx * cz - sx * sy * sz, sx * cz + cx * sy * sz,
            sy, -sx * cy, cx * cy
        );

        // per-instance translation offset from layered noise (Brownian-like)
        let offset = vec3<f32>(
            (fbm(vec3<f32>(seed + 31.7, t, 0.0)) / 0.4375 - 1.0),
            (fbm(vec3<f32>(seed + 31.7, 0.0, t)) / 0.4375 - 1.0),
            (fbm(vec3<f32>(0.0, seed + 31.7, t)) / 0.4375 - 1.0)
        ) * amplitude;

        // bounding-sphere center transformed by the linear part only (no translation)
        let localCenter = mat3From4(transform) * sphere.xyz;

        // rotate basis vectors, pivot around the transformed center
        return mat4x4<f32>(
            vec4<f32>(rot * transform[0].xyz, transform[0].w),
            vec4<f32>(rot * transform[1].xyz, transform[1].w),
            vec4<f32>(rot * transform[2].xyz, transform[2].w),
            vec4<f32>(transform[3].xyz + localCenter - rot * localCenter + offset, transform[3].w)
        );
    }
    return transform;
}
`;

/**
 * \`wiggle\`, the amplitude of the wiggle data (\`tWiggle\`, per instance or group instance, offset per segment of merged
 * items) for \`applyWiggle\`; expects \`in.aInstance\` and \`group\`
 */
export const assign_wiggle_wgsl = /* wgsl */ `
    #if defined(dWiggle)
        #if defined(dSegmented)
            let wiggleGroupInstance = i32(in.aSegment.y + in.aInstance * in.aSegment.x + group);
        #else
            let wiggleGroupInstance = i32(in.aInstance * f32(uObject.uGroupCount) + group);
        #endif
        let wiggleIndex = select(wiggleGroupInstance, i32(in.aInstance), uObject.uInstanceGranularity != 0u);
        let wiggle = read_tWiggle(wiggleIndex).r * uMaterial.uWiggleStrength;
    #else
        let wiggle = 0.0;
    #endif
`;
