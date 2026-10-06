/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { lightingShader } from './lighting';
import { cameraShader } from './shader';
import { clipShader } from './clip';
import { weightedTransparencyShader } from './transparency';
import { depthPeelingShader } from './depth-peeling';

export const volumeShader = cameraShader + lightingShader + clipShader + weightedTransparencyShader + depthPeelingShader + /* wgsl */ `
struct Volume {
    inverseCamera: mat4x4f,
    dimensions: vec4f,
    strides: vec4f,
    options: vec4f,
    color: vec4f,
    domain: vec4f,
    info: vec4u,
    styles: vec4f,
    granularities: vec4u,
    clipConfig: vec4u,
    emission: vec4f,
    material: vec4f,
};
struct Instance { unitToWorld: mat4x4f, worldToUnit: mat4x4f, info: vec4u };
@group(1) @binding(0) var<uniform> volume: Volume;
@group(1) @binding(1) var<storage, read> instances: array<Instance>;
@group(1) @binding(2) var grid: texture_3d<f32>;
@group(1) @binding(3) var transfer: texture_2d<f32>;
@group(1) @binding(4) var palette: texture_2d<f32>;
@group(1) @binding(5) var linearSampler: sampler;
@group(1) @binding(6) var opaqueDepth: texture_depth_2d;
@group(1) @binding(7) var<storage, read> markers: array<u32>;
@group(1) @binding(8) var<storage, read> transparencies: array<u32>;
@group(1) @binding(9) var<storage, read> overpaints: array<u32>;
@group(1) @binding(10) var<storage, read> emissives: array<u32>;
@group(1) @binding(11) var<storage, read> colors: array<u32>;

@group(1) @binding(12) var<storage, read> clipObjects: array<ClipObject>;
@group(1) @binding(13) var<storage, read> clipMasks: array<u32>;
fn volumeClipped(point: vec3f, group: u32, instance: u32) -> bool {
    var mask = 0u;
    if (volume.clipConfig.w > 0u) {
        let index = select(instance * volume.info.y + group, instance, volume.clipConfig.z > 0u);
        if (index / 4u < arrayLength(&clipMasks)) { mask = (clipMasks[index / 4u] >> (index % 4u * 8u)) & 255u; }
    }
    for (var i = 0u; i < volume.clipConfig.x; i++) {
        if ((mask & (i + 1u)) != 0u) { continue; }
        let clip = clipObjects[i];
        if (clipInside(point, clip) != (clip.info.y > 0.0)) { return true; }
    }
    return false;
}
struct Varying {
    @builtin(position) position: vec4f,
    @location(0) @interpolate(flat) instance: u32,
};
@vertex fn vs(@builtin(vertex_index) vertex: u32, @builtin(instance_index) instance: u32) -> Varying {
    let positions = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    var out: Varying;
    out.position = vec4f(positions[vertex], 0.0, 1.0);
    out.instance = instance;
    return out;
}
fn voxel(cell: vec3i) -> vec4f {
    return textureLoad(grid, clamp(cell, vec3i(0), vec3i(volume.dimensions.xyz) - vec3i(1)), 0);
}
fn sampleGrid(unit: vec3f) -> vec4f {
    let position = unit * volume.dimensions.xyz;
    let cell = vec3i(floor(position));
    let f = fract(position);
    let a = mix(voxel(cell), voxel(cell + vec3i(1, 0, 0)), f.x);
    let b = mix(voxel(cell + vec3i(0, 1, 0)), voxel(cell + vec3i(1, 1, 0)), f.x);
    let c = mix(voxel(cell + vec3i(0, 0, 1)), voxel(cell + vec3i(1, 0, 1)), f.x);
    let d = mix(voxel(cell + vec3i(0, 1, 1)), voxel(cell + vec3i(1, 1, 1)), f.x);
    return mix(mix(a, b, f.y), mix(c, d, f.y), f.z);
}
fn markerAt(index: u32) -> u32 {
    if (index / 4u >= arrayLength(&markers)) { return 0u; }
    return (markers[index / 4u] >> ((index % 4u) * 8u)) & 255u;
}
fn transparencyAt(index: u32) -> f32 {
    if (index / 4u >= arrayLength(&transparencies)) { return 0.0; }
    return f32((transparencies[index / 4u] >> ((index % 4u) * 8u)) & 255u) / 255.0;
}
fn emissiveAt(index: u32) -> f32 {
    if (index / 4u >= arrayLength(&emissives)) { return 0.0; }
    return f32((emissives[index / 4u] >> ((index % 4u) * 8u)) & 255u) / 255.0;
}
fn colorByte(index: u32) -> f32 {
    if (index / 4u >= arrayLength(&colors)) { return 1.0; }
    return f32((colors[index / 4u] >> ((index % 4u) * 8u)) & 255u) / 255.0;
}
fn locationIndex(granularity: u32, group: u32, instance: u32) -> u32 {
    if (granularity == 0u) { return instance; }
    return instance * volume.info.y + group;
}
// Position themes use the position iterator's canonical z-fast order, independently
// of the source tensor's group order. Clamp neighbors and use fract at exact voxels.
fn vertexTheme(unit: vec3f, offset: u32, overpaint: bool) -> vec4f {
    let position = clamp(unit * volume.dimensions.xyz, vec3f(0.0), volume.dimensions.xyz - vec3f(1.0));
    let cell = vec3u(floor(position)); let fraction = fract(position);
    let dimensions = vec3u(volume.dimensions.xyz);
    var result = vec4f(0.0);
    for (var z = 0u; z < 2u; z++) {
        for (var y = 0u; y < 2u; y++) {
            for (var x = 0u; x < 2u; x++) {
                let neighbor = min(cell + vec3u(x, y, z), dimensions - vec3u(1u));
                let index = offset + neighbor.z + neighbor.y * dimensions.z + neighbor.x * dimensions.y * dimensions.z;
                var color = vec4f(colorByte(index * 3u), colorByte(index * 3u + 1u), colorByte(index * 3u + 2u), 1.0);
                if (overpaint) {
                    color = vec4f(0.0);
                    if (index < arrayLength(&overpaints)) { color = unpack4x8unorm(overpaints[index]); }
                }
                let weight = select(1.0 - fraction, fraction, vec3u(x, y, z) > vec3u(0u));
                result += color * weight.x * weight.y * weight.z;
            }
        }
    }
    return result;
}
struct Fragment {
    @location(0) color: vec4f,
    @location(1) pick: vec4u,
    @location(2) emissive: vec4f,
    @builtin(frag_depth) depth: f32,
};
fn volumeFragment(in: Varying, colorPass: bool) -> Fragment {
    let instance = instances[in.instance];
    if (volume.clipConfig.y > 0u && volumeClipped((instance.unitToWorld * vec4f(0.5, 0.5, 0.5, 1.0)).xyz, 0u, instance.info.x)) { discard; }
    let pixel = in.position.xy;
    let local = (pixel - camera.viewportRect.xy) / camera.viewportRect.zw;
    let ndc = vec2f(local.x * 2.0 - 1.0, 1.0 - local.y * 2.0);
    let nearClip = volume.inverseCamera * vec4f(ndc, -1.0, 1.0);
    let farClip = volume.inverseCamera * vec4f(ndc, 1.0, 1.0);
    let nearWorld = nearClip.xyz / nearClip.w;
    let farWorld = farClip.xyz / farClip.w;
    let origin = (instance.worldToUnit * vec4f(nearWorld, 1.0)).xyz;
    let end = (instance.worldToUnit * vec4f(farWorld, 1.0)).xyz;
    let direction = normalize(end - origin);
    let safeDirection = select(direction, vec3f(0.0000001), abs(direction) < vec3f(0.0000001));
    let distancesA = -origin / safeDirection;
    let distancesB = (vec3f(1.0) - origin) / safeDirection;
    let nearDistances = min(distancesA, distancesB);
    let farDistances = max(distancesA, distancesB);
    let start = max(0.0, max(nearDistances.x, max(nearDistances.y, nearDistances.z)));
    let stop = min(farDistances.x, min(farDistances.y, farDistances.z));
    if (stop <= start) { discard; }
    let sceneDepth = textureLoad(opaqueDepth, vec2i(pixel), 0);
    let step = volume.options.y / max(length(direction * volume.dimensions.xyz), 0.000001);
    var distance = start + step * 0.5;
    var accumulated = vec4f(0.0);
    var emitted = vec3f(0.0);
    var firstDepth = 1.0;
    var firstGroup = 0u;
    let modelView = camera.view * instance.unitToWorld;
    let m = mat3x3f(modelView[0].xyz, modelView[1].xyz, modelView[2].xyz);
    let normalMatrix = mat3x3f(cross(m[1], m[2]), cross(m[2], m[0]), cross(m[0], m[1]));
    for (var i = 0u; i < u32(volume.options.x); i++) {
        if (distance > stop || accumulated.a > 0.995) { break; }
        let unit = origin + direction * distance;
        let world = instance.unitToWorld * vec4f(unit, 1.0);
        let view = camera.view * world;
        let clip = camera.projection * view;
        let depth = (clip.z / clip.w + 1.0) * 0.5;
        if (depth > sceneDepth || depth > 1.0) { break; }
        let cell = sampleGrid(unit);
        let transferAlpha = textureSampleLevel(transfer, linearSampler, vec2f(cell.a, 0.5), 0.0).r;
        let coord = clamp(floor(unit * volume.dimensions.xyz + vec3f(0.5)), vec3f(0.0), volume.dimensions.xyz - vec3f(1.0));
        let group = u32(dot(coord, volume.strides.xyz));
        let logicalInstance = instance.info.x;
        if (volume.clipConfig.y == 0u && volumeClipped(world.xyz, group, logicalInstance)) { distance += step; continue; }
        var color = volume.color.rgb;
        if (volume.info.z == 1u) {
            let paletteValue = clamp((cell.a - volume.domain.x) / max(volume.domain.y - volume.domain.x, 0.000001), 0.0, 1.0);
            if (volume.domain.z > 0.0) { color = textureSampleLevel(palette, linearSampler, vec2f(paletteValue, 0.5), 0.0).rgb; }
            else {
                let size = i32(textureDimensions(palette).x);
                color = textureLoad(palette, vec2i(clamp(i32(floor(paletteValue * f32(size))), 0, size - 1), 0), 0).rgb;
            }
        } else if (volume.info.z >= 5u) {
            let offset = select(0u, logicalInstance * volume.info.w, volume.info.z == 6u);
            color = vertexTheme(unit, offset, false).rgb;
        } else if (volume.info.z > 1u) {
            var index = group;
            if (volume.info.z == 2u) { index = logicalInstance; }
            if (volume.info.z == 4u) { index += logicalInstance * volume.info.y; }
            color = vec3f(colorByte(index * 3u), colorByte(index * 3u + 1u), colorByte(index * 3u + 2u));
        }
        if (volume.granularities.y == 2u) {
            let overlay = vertexTheme(unit, logicalInstance * volume.info.w, true);
            color = mix(color, overlay.rgb, overlay.a * volume.styles.y);
        } else {
            let overlayIndex = locationIndex(volume.granularities.y, group, logicalInstance);
            var overlay = vec4f(0.0);
            if (overlayIndex < arrayLength(&overpaints)) { overlay = unpack4x8unorm(overpaints[overlayIndex]); }
            color = mix(color, overlay.rgb, overlay.a * volume.styles.y);
        }
        let transparencyIndex = locationIndex(volume.granularities.z, group, logicalInstance);
        var alpha = transferAlpha * volume.options.y * volume.options.z * (1.0 - transparencyAt(transparencyIndex) * volume.styles.z);
        if (alpha > 0.001) {
            if (firstDepth == 1.0) { firstDepth = depth; firstGroup = group; }
            var emissiveColor = color * (volume.emission.x + emissiveAt(locationIndex(volume.granularities.w, group, logicalInstance)) * volume.styles.w);
            let gradient = normalMatrix * (cell.xyz * 2.0 - vec3f(1.0));
            let normal = -normalized(gradient);
            let emission = volume.emission.x + emissiveAt(locationIndex(volume.granularities.w, group, logicalInstance)) * volume.styles.w;
            color = shadeMaterial(color, normal, view.xyz, volume.material, emission, volume.options.w > 0.0, 0.0);
            var marker = u32(max(volume.styles.x, 0.0));
            if (volume.styles.x < 0.0) { marker = markerAt(locationIndex(volume.granularities.x, group, logicalInstance)); }
            if (marker > 0u) {
                let marking = select(camera.selection, camera.highlight, (marker & 1u) == 1u);
                color = mix(color, marking.rgb, marking.a);
            }
            if (camera.fogLight.y > camera.fogLight.x) { let fog = smoothstep(camera.fogLight.x, camera.fogLight.y, -view.z); if (colorPass && camera.animation.w > 0.0) { alpha *= 1.0 - fog; } else { color = mix(color, camera.background.rgb, fog); emissiveColor *= 1.0 - fog; } }
            emitted += emissiveColor * alpha * (1.0 - accumulated.a);
            accumulated += vec4f(color * alpha, alpha) * (1.0 - accumulated.a);
        }
        distance += step;
    }
    if (accumulated.a <= 0.001) { discard; }
    var out: Fragment;
    out.color = accumulated;
    out.emissive = vec4f(emitted, accumulated.a);
    out.pick = vec4u(volume.info.x + 1u, instance.info.x, firstGroup, bitcast<u32>(firstDepth));
    out.depth = firstDepth;
    return out;
}
@fragment fn fs(in: Varying) -> Fragment { return volumeFragment(in, true); }
struct WeightedVolumeFragment { @location(0) color: vec4f, @location(1) weight: vec4f, @location(2) emission: vec4f, @builtin(frag_depth) depth: f32 };
@fragment fn weighted(in: Varying) -> WeightedVolumeFragment {
    let result = volumeFragment(in, true); let accum = weightedFragment(result.color, result.emissive, result.depth);
    return WeightedVolumeFragment(accum.color, accum.weight, accum.emission, result.depth);
}
struct PickFragment { @location(0) ids: vec4u, @builtin(frag_depth) depth: f32 };
@fragment fn colorDepthId(in: Varying) -> PickFragment {
    let result = volumeFragment(in, true); return PickFragment(result.pick, result.depth);
}
@fragment fn peelDepth(in: Varying) -> @builtin(frag_depth) f32 {
    let result = volumeFragment(in, true); return peelSearch(result.depth, in.position.xy);
}
@fragment fn peelFront(in: Varying) -> PeelFragment {
    let result = volumeFragment(in, true); return peelLayer(result.color, result.emissive, result.depth, in.position.xy, true);
}
@fragment fn peelBack(in: Varying) -> PeelFragment {
    let result = volumeFragment(in, true); return peelLayer(result.color, result.emissive, result.depth, in.position.xy, false);
}
@fragment fn transparentColor(in: Varying) -> @location(0) vec4f { return volumeFragment(in, true).color; }
struct OutlineFragment { @location(0) depthAlpha: vec4f, @builtin(frag_depth) depth: f32 };
@fragment fn outlineDepth(in: Varying) -> OutlineFragment {
    let result = volumeFragment(in, false);
    var out: OutlineFragment; out.depthAlpha = vec4f(result.depth, result.color.a, 0.0, 0.0); out.depth = result.depth; return out;
}
@fragment fn pick(in: Varying) -> PickFragment {
    let result = volumeFragment(in, false);
    if (result.color.a < lighting.dimensions.z) { discard; }
    var out: PickFragment; out.ids = result.pick; out.depth = result.depth;
    return out;
}

`;
