/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { lightingShader } from './lighting';
import { imageShader } from './image';
import { clipShader } from './clip';
import { animationShader } from './animation';
import { weightedTransparencyShader } from './transparency';
import { depthPeelingShader } from './depth-peeling';

export const cameraShader = /* wgsl */ `
struct Camera {
    projection: mat4x4f,
    view: mat4x4f,
    background: vec4f,
    viewport: vec4f,
    highlight: vec4f,
    selection: vec4f,
    fogLight: vec4f,
    viewportRect: vec4f,
    animation: vec4f,
    lodRaster: vec4f,
};
@group(0) @binding(0) var<uniform> camera: Camera;
`;

export const geometryShader = cameraShader + lightingShader + clipShader + weightedTransparencyShader + depthPeelingShader + /* wgsl */ `
struct Object {
    info: vec4u,
    options: vec4f,
    textOffset: vec4f,
    border: vec4f,
    textBackground: vec4f,
    image: vec4f,
    paletteDefault: vec4f,
    trimCenter: vec4f,
    trimRotation: vec4f,
    trimScale: vec4f,
    trimTransform: mat4x4f,
    clipConfig: vec4f,
    sphere: vec4f,
    wiggle: vec4f,
    tumble: vec4f,
    wiggleOverlay: vec4f,
    material: vec4f,
    shading: vec4f,
    interiorColor: vec4f,
    interiorSubstance: vec4f,
    surface: vec4f,
    substanceGrid: vec4f,
    substanceTransform: vec4f,
    overpaintGrid: vec4f,
    overpaintTransform: vec4f,
    transparencyGrid: vec4f,
    transparencyTransform: vec4f,
    emissiveGrid: vec4f,
    emissiveTransform: vec4f,
    colorGrid: vec4f,
    colorTransform: vec4f,
    colorPalette: vec4f,
    lod: vec4f,
};
struct Vertex { position: vec4f, normal: vec4f, end: vec4f, mapping: vec4f, group: vec4f };
struct Theme { color: vec4f, style: vec4f, substance: vec4f };
@group(1) @binding(0) var<uniform> object: Object;
@group(1) @binding(1) var<storage, read> vertices: array<Vertex>;
@group(1) @binding(2) var<storage, read> transforms: array<mat4x4f>;
@group(1) @binding(3) var<storage, read> themes: array<Theme>;
@group(1) @binding(4) var atlas: texture_2d<f32>;
@group(1) @binding(5) var atlasSampler: sampler;
@group(1) @binding(13) var substanceGrid: texture_2d<f32>;
@group(1) @binding(14) var overpaintGrid: texture_2d<f32>;
@group(1) @binding(15) var transparencyGrid: texture_2d<f32>;
@group(1) @binding(16) var emissiveGrid: texture_2d<f32>;
@group(1) @binding(17) var colorGrid: texture_2d<f32>;
@group(1) @binding(18) var<storage, read> overpaintThemes: array<vec4f>;
fn sampleThemeGrid(tex: texture_2d<f32>, point: vec3f, dim: vec3f, transform: vec4f) -> vec4f {
    let p = (point - transform.xyz) * transform.w;
    let texDim = vec2f(textureDimensions(tex));
    let z = floor(p.z);
    let origin0 = vec2f((z * dim.x) % texDim.x, floor(z * dim.x / texDim.x) * dim.y);
    let origin1 = vec2f(((z + 1.0) * dim.x) % texDim.x, floor((z + 1.0) * dim.x / texDim.x) * dim.y);
    return mix(textureSampleLevel(tex, atlasSampler, (origin0 + p.xy) / texDim, 0.0),
        textureSampleLevel(tex, atlasSampler, (origin1 + p.xy) / texDim, 0.0), fract(p.z));
}
` + imageShader + animationShader + /* wgsl */ `
@group(1) @binding(10) var<storage, read> clipObjects: array<ClipObject>;
@group(1) @binding(11) var<storage, read> clipMasks: array<u32>;
fn clippingMask(group: u32, instance: u32) -> u32 {
    if (object.clipConfig.z == 0.0) { return 0u; }
    let index = select(instance * u32(object.clipConfig.w) + group, instance, object.clipConfig.y > 0.0);
    if (index / 4u >= arrayLength(&clipMasks)) { return 0u; }
    return (clipMasks[index / 4u] >> (index % 4u * 8u)) & 255u;
}
fn clipped(point: vec3f, mask: u32) -> bool {
    for (var i = 0u; i < object.info.w; i++) {
        // Preserve Mol* clip-group identifiers, which start at 1.
        if ((mask & (i + 1u)) != 0u) { continue; }
        let clip = clipObjects[i];
        if (clipInside(point, clip) != (clip.info.y > 0.0)) { return true; }
    }
    return false;
}
struct Varying {
    @builtin(position) position: vec4f,
    @location(0) color: vec4f,
    @location(1) normal: vec3f,
    @location(2) uv: vec2f,
    @location(3) viewDepth: f32,
    @location(4) @interpolate(flat) ids: vec4u,
    @location(5) @interpolate(flat) marker: u32,
    @location(6) emissive: f32,
    @location(7) localPosition: vec3f,
    @location(8) worldPosition: vec3f,
    @location(9) viewPosition: vec4f,
    @location(11) substance: vec4f,
    @location(12) overpaint: vec4f,
    @location(13) @interpolate(flat) primitiveStart: vec4f,
    @location(14) @interpolate(flat) primitiveEnd: vec4f,
    @location(15) @interpolate(flat) primitiveEye: vec4f,
};
fn inverseAffinePoint(m: mat4x4f, p: vec3f) -> vec3f {
    let a = m[0].xyz; let b = m[1].xyz; let c = m[2].xyz;
    let r = p - m[3].xyz;
    return vec3f(dot(cross(b, c), r), dot(cross(c, a), r), dot(cross(a, b), r)) / dot(a, cross(b, c));
}
fn project(position: vec4f) -> vec4f {
    var clip = camera.projection * position;
    // Mol* cameras use OpenGL's -1..1 clip depth; WebGPU requires 0..1.
    clip.z = (clip.z + clip.w) * 0.5;
    return clip;
}
fn distanceCoverage(lod: vec4f, distance: f32) -> f32 {
    if (lod.z <= 0.0) { return select(0.0, 1.0, distance >= lod.x && distance <= lod.y); }
    return min(smoothstep(lod.x, lod.x + lod.z, distance), 1.0 - smoothstep(lod.y - lod.z, lod.y, distance));
}
fn distanceFade(distance: f32) -> f32 { return distanceCoverage(object.lod, distance); }
@vertex fn vs(@builtin(vertex_index) vid: u32, @builtin(instance_index) iid: u32) -> Varying {
    let vertex = vertices[vid];
    let theme = themes[iid * object.info.y + vid];
    let kind = object.info.z;
    var transform = transforms[iid];
    var position = vertex.position;
    var lineEnd = vertex.end.xyz;
    if (kind != 4u) {
        transform = tumble(transform, u32(theme.style.w));
        position = vec4f(wiggle(position.xyz, vertex.group.x, u32(theme.style.w)), 1.0);
        if (kind == 1u) { lineEnd = wiggle(lineEnd, vertex.group.x, u32(theme.style.w)); }
    }
    let modelView = camera.view * transform;
    var size = theme.style.x;
    var lodHidden = false;
    if (kind == 5u && vertex.normal.w > 0.0) {
        let distance = -(modelView * position).z / camera.viewport.w;
        var factor = vertex.normal.w;
        if (camera.viewport.w == 1.0) {
            lodHidden = distance < vertex.normal.x || distance > vertex.normal.y;
            factor *= distanceCoverage(vertex.normal, distance);
        }
        size *= factor;
    }
    let nearPatch = (kind == 5u || kind == 6u) && vertex.mapping.w < 0.0;
    let modelScale = max(length(modelView[0].xyz), max(length(modelView[1].xyz), length(modelView[2].xyz)));
    var primitiveStart = vec4f(0.0); var primitiveEnd = vec4f(0.0);
    var capMask = 0.0;
    var localStart = position.xyz; var localEnd = localStart;
    if (kind == 5u) {
        primitiveStart = vec4f((modelView * position).xyz, size * vertex.mapping.x * modelScale);
        primitiveEnd.w = select(0.0, -1.0, nearPatch);
    } else if (kind == 6u) {
        var start = position; var end = vec4f(vertex.end.xyz, 1.0);
        let first = u32(select(vertex.mapping.w, vertex.normal.w, nearPatch)) - 1u;
        capMask = vertices[first].normal.w;
        if (!nearPatch) {
            start = vec4f(wiggle(vertices[first].position.xyz, vertex.group.x, u32(theme.style.w)), 1.0);
            end = vec4f(wiggle(vertices[first + 17u].position.xyz, vertex.group.x, u32(theme.style.w)), 1.0);
        } else { end = vec4f(wiggle(end.xyz, vertex.group.x, u32(theme.style.w)), 1.0); }
        primitiveStart = vec4f((modelView * start).xyz, theme.style.x * select(vertex.mapping.x, vertex.mapping.z, nearPatch) * modelScale);
        primitiveEnd = vec4f((modelView * end).xyz, select(select(0.0, -2.0, vertex.normal.w < 0.0), -1.0, nearPatch));
        localStart = start.xyz; localEnd = end.xyz;
    }
    if ((kind == 5u || kind == 6u) && !nearPatch) { position = vec4f(position.xyz + vertex.end.xyz * size * vertex.mapping.x, 1.0); }
    var view = modelView * position;
    if (nearPatch) {
        let nearZ = -(camera.projection[3][2] + camera.projection[3][3]) / (camera.projection[2][2] + camera.projection[2][3]);
        let start = primitiveStart.xyz; let end = select(start, primitiveEnd.xyz, kind == 6u);
        let minimum = min(start.xy, end.xy) - primitiveStart.w;
        let maximum = max(start.xy, end.xy) + primitiveStart.w;
        let corner = select(vertex.end.xy, vertex.mapping.xy, kind == 6u);
        view = vec4f(mix(minimum, maximum, corner * 0.5 + 0.5), nearZ, 1.0);
        position = vec4f(inverseAffinePoint(modelView, view.xyz), 1.0);
    }
    var clip = project(view);
    if (nearPatch) { clip.z = max(0.0000001 / max(primitiveStart.w, 0.00001), 0.00000001) * clip.w; }
    var uv = vertex.group.zw;
    if (kind == 1u) {
        var start = modelView * position;
        var end = modelView * vec4f(lineEnd, 1.0);
        // Trim segments crossing the eye plane before dividing by w.
        let near = -0.5 * camera.projection[3][2] / camera.projection[2][2];
        if (camera.projection[2][3] == -1.0) {
            if (start.z >= 0.0 && end.z < 0.0) { start = mix(end, start, (near - end.z) / (start.z - end.z)); }
            if (end.z >= 0.0 && start.z < 0.0) { end = mix(start, end, (near - start.z) / (end.z - start.z)); }
        }
        let a = project(start); let b = project(end);
        let delta = (b.xy / b.w - a.xy / a.w) * camera.viewport.xy;
        let direction = delta / max(length(delta), 0.00001);
        let perpendicular = vec2f(direction.y, -direction.x);
        var size = theme.style.x * camera.viewport.z;
        if (object.options.y > 0.0) { size *= camera.viewport.y / max(-start.z, 0.01) * 2.5; }
        clip = select(a, b, vertex.mapping.y >= 0.5);
        clip = vec4f(clip.xy + perpendicular * vertex.mapping.x * max(1.0, size) / camera.viewport.xy * clip.w, clip.zw);
    } else if (kind == 2u) {
        var size = theme.style.x * camera.viewport.z;
        if (object.options.y > 0.0) { size *= camera.viewport.y / max(-view.z, 0.01); }
        clip = vec4f(clip.xy + vertex.mapping.xy * max(1.0, size) / camera.viewport.xy * clip.w, clip.zw);
        uv = vertex.mapping.xy;
    } else if (kind == 3u) {
        var center = view;
        center.z += (object.textOffset.z + vertex.end.w * 0.95) * camera.viewport.w;
        clip = project(center);
        let offset = (vertex.mapping.xy * theme.style.x + object.textOffset.xy) * camera.viewport.w;
        clip = vec4f(clip.xy + vec2f(camera.projection[0][0], camera.projection[1][1]) * offset, clip.zw);
    }
    let m = mat3x3f(modelView[0].xyz, modelView[1].xyz, modelView[2].xyz);
    let normalMatrix = mat3x3f(cross(m[1], m[2]), cross(m[2], m[0]), cross(m[0], m[1]));
    var out: Varying;
    if (object.clipConfig.x > 0.0 && clipped((transform * vec4f(object.sphere.xyz, 1.0)).xyz, clippingMask(u32(vertex.group.x), u32(theme.style.w)))) {
        clip.z = 2.0 * clip.w;
    }
    if (lodHidden || (kind == 5u && size <= 0.0)) { clip.z = 2.0 * clip.w; }
    out.position = clip;
    out.color = theme.color;
    if (nearPatch && kind == 6u) {
        let first = u32(vertex.normal.w) - 1u;
        let axis = primitiveEnd.xyz - primitiveStart.xyz;
        let fraction = clamp(dot(view.xyz - primitiveStart.xyz, axis) / dot(axis, axis), 0.0, 1.0);
        out.color = mix(themes[iid * object.info.y + first].color, themes[iid * object.info.y + first + 17u].color, fraction);
    }
    out.normal = normalMatrix * vertex.normal.xyz * select(1.0, -1.0, object.options.w > 0.0);
    out.uv = uv;
    out.viewDepth = -view.z;
    out.viewPosition = vec4f(view.xyz, -select(view.z, primitiveStart.z, kind == 5u) / camera.viewport.w);
    let primitiveRadius = size * select(vertex.mapping.x, vertex.mapping.z, nearPatch && kind == 6u);
    out.ids = vec4u(object.info.x + 1u, u32(theme.style.w), u32(vertex.group.x), bitcast<u32>(primitiveRadius));
    out.marker = u32(theme.style.y);
    out.emissive = theme.style.z;
    out.primitiveStart = vec4f(localStart, primitiveRadius);
    if (kind == 6u) { out.primitiveStart.x = capMask; }
    out.primitiveEnd = vec4f(localEnd - localStart, primitiveEnd.w);
    out.primitiveEye = vec4f(0.0);
    if (kind == 5u) { out.primitiveEnd.x = f32(iid); out.primitiveEnd.y = f32(vid); }
    if (kind == 5u || kind == 6u) {
        let eye = inverseAffinePoint(modelView, vec3f(0.0));
        let orthoDirection = inverseAffinePoint(modelView, vec3f(0.0, 0.0, -1.0)) - eye;
        out.primitiveEye = select(vec4f(eye - localStart, -1.0), vec4f(orthoDirection, length(orthoDirection)), camera.projection[2][3] == 0.0);
    }
    out.substance = theme.substance;
    out.localPosition = vertex.position.xyz;
    if (kind == 5u || kind == 6u) { out.localPosition = position.xyz - localStart; }
    var modelPosition = position;
    if (kind == 1u) { modelPosition = mix(position, vec4f(lineEnd, 1.0), vertex.mapping.y); }
    out.worldPosition = (transform * modelPosition).xyz;
    if (kind != 5u && object.colorGrid.w >= 0.0) {
        let point = select(modelPosition.xyz, out.worldPosition, object.colorGrid.w > 1.0);
        out.color = vec4f(sampleThemeGrid(colorGrid, point, object.colorGrid.xyz, object.colorTransform).rgb, out.color.a);
    }
    if (kind != 5u && object.substanceGrid.w >= 0.0) {
        let s = sampleThemeGrid(substanceGrid, out.worldPosition, object.substanceGrid.xyz, object.substanceTransform);
        out.substance = vec4f(mix(object.material.xyz, s.rgb, s.a), s.a) * object.substanceGrid.w;
    }
    out.overpaint = vec4f(0.0);
    if (kind != 5u && object.overpaintGrid.w >= 0.0) {
        let o = sampleThemeGrid(overpaintGrid, out.worldPosition, object.overpaintGrid.xyz, object.overpaintTransform);
        out.overpaint = vec4f(mix(out.color.rgb, o.rgb, o.a), o.a) * object.overpaintGrid.w;
    } else if (kind != 4u && object.colorPalette.w > 0.0) {
        let o = overpaintThemes[iid * object.info.y + vid];
        out.overpaint = vec4f(mix(out.color.rgb, o.rgb, o.a), o.a) * object.colorPalette.z;
    }
    if (kind != 5u && object.transparencyGrid.w >= 0.0) {
        let t = sampleThemeGrid(transparencyGrid, out.worldPosition, object.transparencyGrid.xyz, object.transparencyTransform).a;
        out.color.a *= 1.0 - t * object.transparencyGrid.w;
    }
    if (kind != 5u && object.emissiveGrid.w >= 0.0) {
        out.emissive += sampleThemeGrid(emissiveGrid, out.worldPosition, object.emissiveGrid.xyz, object.emissiveTransform).a * object.emissiveGrid.w;
    }
    return out;
}
struct Fragment {
    @location(0) color: vec4f,
    @location(1) pick: vec4u,
    @location(2) emissive: vec4f,
    @builtin(frag_depth) depth: f32,
};
struct SurfaceFragment { color: vec4f, pick: vec4u, emissive: vec4f, normal: vec4f, albedo: vec4f, depth: f32 };
fn geometryFragment(input: Varying, front: bool, applyThickness: bool, tracing: bool, backDepth: bool) -> SurfaceFragment {
    var in = input;
    let kind = object.info.z;
    let nearPatch = in.primitiveEnd.w == -1.0;
    let artificialCap = in.primitiveEnd.w == -2.0;
    var interior = !front || nearPatch || artificialCap;
    var normal = -normalized(in.normal) * select(-1.0, 1.0, front);
    let flatNormal = normalized(cross(dpdx(in.viewPosition.xyz), dpdy(in.viewPosition.xyz)));
    if (object.shading.z > 0.0) { normal = flatNormal * select(1.0, -1.0, object.options.w > 0.0); }
    let dxy = max(abs(dpdx(normal)), abs(dpdy(normal)));
    let geometryRoughness = select(max(dxy.x, max(dxy.y, dxy.z)), 0.0, kind != 0u);
    if (kind == 5u && !nearPatch) {
        let ortho = in.primitiveEye.w > 0.0;
        let direction = select(normalized(in.localPosition - in.primitiveEye.xyz), normalized(in.primitiveEye.xyz), ortho);
        let viewDirection = select(normalized(in.viewPosition.xyz), vec3f(0.0, 0.0, -1.0), ortho);
        let ratio = select(length(in.viewPosition.xyz) / length(in.localPosition - in.primitiveEye.xyz), 1.0 / in.primitiveEye.w, ortho);
        let b = dot(in.localPosition, direction);
        let determinant = b * b + in.primitiveStart.w * in.primitiveStart.w - dot(in.localPosition, in.localPosition);
        if (determinant < 0.0 || in.primitiveStart.w <= 0.0) { discard; }
        let root = sqrt(max(0.0, determinant));
        var t = -b - root;
        var point = in.viewPosition.xyz + viewDirection * t * ratio;
        var world = inverseAffinePoint(camera.view, point);
        var projected = project(vec4f(point, 1.0));
        let clippedFront = object.clipConfig.x == 0.0 && clipped(world, clippingMask(in.ids.z, in.ids.y));
        var useBack = backDepth || clippedFront || projected.z < 0.0;
        if (useBack && !backDepth && !clippedFront && object.surface.w > 0.0) { discard; }
        if (!front && !useBack) { discard; }
        if (front && useBack && !backDepth) { discard; }
        if (useBack) {
            t = -b + root; point = in.viewPosition.xyz + viewDirection * t * ratio;
            world = inverseAffinePoint(camera.view, point); projected = project(vec4f(point, 1.0));
        }
        let depth = projected.z / projected.w;
        if (depth < 0.0 || depth > 1.0 || projected.w <= 0.0) { discard; }
        let localNormal = normalized(in.localPosition + direction * t);
        let transform = tumble(transforms[u32(in.primitiveEnd.x)], in.ids.y);
        let mv = camera.view * transform;
        let m = mat3x3f(mv[0].xyz, mv[1].xyz, mv[2].xyz);
        let normalMatrix = mat3x3f(cross(m[1], m[2]), cross(m[2], m[0]), cross(m[0], m[1]));
        normal = normalized(normalMatrix * localNormal) * select(-1.0, 1.0, useBack);
        interior = useBack;
        in.localPosition = in.primitiveStart.xyz + in.localPosition + direction * t;
        in.worldPosition = world; in.viewPosition = vec4f(point, in.viewPosition.w);
        in.viewDepth = -point.z; in.position.z = depth;
    }
    if (object.lod.w == 0.0 && (object.lod.x != 0.0 || object.lod.y != 0.0)) {
        let fade = distanceFade(in.viewPosition.w);
        if (applyThickness || tracing) {
            // The same ordered coverage mask as the molecular GLSL renderer.
            let thresholds = array<f32, 16>(1, 9, 3, 11, 13, 5, 15, 7, 4, 12, 2, 10, 16, 8, 14, 6);
            let pixel = vec2i(floor(vec2f(in.position.x, camera.lodRaster.x - in.position.y) + camera.lodRaster.yz));
            let column = u32((pixel.x % 4 + 4) % 4); let row = u32((pixel.y % 4 + 4) % 4);
            let threshold = thresholds[column * 4u + row] / 17.0;
            if (fade < 0.99 && (fade < 0.01 || fade < threshold)) { discard; }
        } else if (fade < lighting.dimensions.z) { discard; }
    }
    if ((kind == 5u && nearPatch) || kind == 6u) {
        if (artificialCap && !backDepth) { normal = select(normalized(in.viewPosition.xyz), vec3f(0.0, 0.0, -1.0), camera.projection[2][3] == 0.0); }
        if (nearPatch) {
            if (object.surface.w == 0.0 || backDepth) { discard; }
            var distance = length(in.localPosition);
            if (kind == 6u) {
                let axis = in.primitiveEnd.xyz; let fraction = dot(in.localPosition, axis) / dot(axis, axis);
                if (fraction < 0.0 || fraction > 1.0) { discard; }
                distance = length(in.localPosition - axis * fraction);
            }
            if (distance > in.primitiveStart.w) { discard; }
            normal = select(normalized(in.viewPosition.xyz), vec3f(0.0, 0.0, -1.0), camera.projection[2][3] == 0.0);
        } else if ((!front || artificialCap) && !backDepth) {
            // Impostors emit one intersection per ray. Do not blend a second
            // interior layer through a visible front surface.
            let ortho = in.primitiveEye.w > 0.0;
            let direction = select(normalized(in.localPosition - in.primitiveEye.xyz), normalized(in.primitiveEye.xyz), ortho);
            let viewDirection = select(normalized(in.viewPosition.xyz), vec3f(0.0, 0.0, -1.0), ortho);
            let ratio = select(length(in.viewPosition.xyz) / length(in.localPosition - in.primitiveEye.xyz), 1.0 / in.primitiveEye.w, ortho);
            let offset = in.localPosition; let radius = in.primitiveStart.w;
            var t = 0.0; var hit = false;
            if (kind == 5u) {
                let b = dot(offset, direction); let det = b * b + radius * radius - dot(offset, offset);
                if (det >= 0.0) { t = -b - sqrt(det); hit = true; }
            } else {
                let axis = in.primitiveEnd.xyz;
                let length2 = dot(axis, axis); let ad = dot(axis, direction); let ao = dot(axis, offset);
                let a = length2 - ad * ad; let b = length2 * dot(offset, direction) - ao * ad;
                let c = length2 * dot(offset, offset) - ao * ao - radius * radius * length2;
                let det = b * b - a * c;
                if (det >= 0.0 && abs(a) > 0.000001) {
                    let intersection = (-b - sqrt(det)) / a; let y = ao + intersection * ad;
                    if (y >= 0.0 && y <= length2) { t = intersection; hit = true; }
                }
                if (abs(ad) > 0.000001) {
                    for (var cap = 0u; cap < 2u; cap++) {
                        if (object.surface.w == 0.0 && (u32(in.primitiveStart.x) & (1u << cap)) == 0u) { continue; }
                        let intersection = (f32(cap) * length2 - ao) / ad;
                        let radial = offset + direction * intersection - axis * f32(cap);
                        if (dot(radial, radial) <= radius * radius && (!hit || intersection < t)) { t = intersection; hit = true; }
                    }
                }
            }
            if (hit) {
                let point = in.viewPosition.xyz + viewDirection * t * ratio; let projected = project(vec4f(point, 1.0));
                let world = inverseAffinePoint(camera.view, point);
                let frontClipped = object.clipConfig.x == 0.0 && clipped(world, clippingMask(in.ids.z, in.ids.y));
                if (artificialCap) {
                    if (frontClipped || t < -0.0001 || projected.z <= 0.0) { discard; }
                } else if (!frontClipped && (projected.z > 0.0 || object.surface.w > 0.0)) { discard; }
            }
        }
    }
    var color = in.color;
    if (kind == 5u) {
        if (object.colorGrid.w >= 0.0) {
            let localPoint = select(in.localPosition, in.primitiveStart.xyz + in.localPosition, nearPatch);
            let point = select(localPoint, in.worldPosition, object.colorGrid.w > 1.0);
            color = vec4f(sampleThemeGrid(colorGrid, point, object.colorGrid.xyz, object.colorTransform).rgb, color.a);
        }
        if (object.substanceGrid.w >= 0.0) {
            let sample = sampleThemeGrid(substanceGrid, in.worldPosition, object.substanceGrid.xyz, object.substanceTransform);
            in.substance = vec4f(mix(object.material.xyz, sample.rgb, sample.a) * object.substanceGrid.w, sample.a * object.substanceGrid.w);
        }
        if (object.overpaintGrid.w >= 0.0) {
            let sample = sampleThemeGrid(overpaintGrid, in.worldPosition, object.overpaintGrid.xyz, object.overpaintTransform);
            in.overpaint = vec4f(mix(color.rgb, sample.rgb, sample.a), sample.a) * object.overpaintGrid.w;
        } else if (object.colorPalette.w > 0.0) {
            let overlay = overpaintThemes[u32(in.primitiveEnd.x) * object.info.y + u32(in.primitiveEnd.y)];
            in.overpaint = vec4f(mix(color.rgb, overlay.rgb, overlay.a), overlay.a) * object.colorPalette.z;
        }
        if (object.transparencyGrid.w >= 0.0) { color.a *= 1.0 - sampleThemeGrid(transparencyGrid, in.worldPosition, object.transparencyGrid.xyz, object.transparencyTransform).a * object.transparencyGrid.w; }
        if (object.emissiveGrid.w >= 0.0) { in.emissive += sampleThemeGrid(emissiveGrid, in.worldPosition, object.emissiveGrid.xyz, object.emissiveTransform).a * object.emissiveGrid.w; }
    }
    if (kind != 4u && object.colorPalette.x > 0.0) {
        let encoded = (dot(color.rgb, vec3f(16711680.0, 65280.0, 255.0)) - 1.0) / 16777214.0;
        color = vec4f(paletteColor(encoded), color.a);
    }
    var ids = in.ids.xyz;
    var marker = in.marker;
    var emission = in.emissive;
    if (kind == 2u && object.options.z > 0.0) {
        let radius = length(in.uv);
        if (radius > 1.0) { discard; }
        if (object.options.z == 2.0) { color.a *= 1.0 - smoothstep(0.0, 1.0, radius); }
    }
    if (kind == 3u) {
        if (in.uv.x > 1.0) { color = vec4f(object.textBackground.rgb, object.textBackground.a * color.a); }
        else {
            let sdf = textureSampleLevel(atlas, atlasSampler, in.uv, 0.0).r;
            if (sdf + min(object.border.a, 0.49) < 0.5) { discard; }
            if (sdf < 0.5) { color = vec4f(object.border.rgb, color.a); }
        }
    } else if (kind == 4u) {
        if (imageOutside(in.localPosition)) { discard; }
        color = imageSample(in.uv);
        if (object.image.x >= 0.0) {
            if (!(textureLoad(imageValues, imageCoord(in.uv), 0).r >= object.image.x)) { discard; }
            color.a = in.color.a;
        } else { color.a *= in.color.a; }
        if (object.image.z > 0.0) {
            if (all(color.rgb == vec3f(1.0))) { color = vec4f(object.paletteDefault.rgb, color.a); }
            else {
                let encoded = (dot(color.rgb, vec3f(16711680.0, 65280.0, 255.0)) - 1.0) / 16777214.0;
                color = vec4f(paletteColor(encoded), color.a);
            }
        }
        let packed = vec3u(round(textureLoad(imageGroups, imageCoord(in.uv), 0).rgb * 255.0));
        let group = packed.x * 65536u + packed.y * 256u + packed.z;
        ids.z = select(0u, group - 1u, group > 0u);
        marker = 0u;
        if (group > 0u && ids.z < u32(object.image.w)) {
            let style = imageStyles[ids.y * u32(object.image.w) + ids.z];
            color = vec4f(mix(color.rgb, style.color.rgb, style.color.a), color.a * (1.0 - style.style.y));
            marker = u32(style.style.x); emission = style.style.z;
        } else { ids = vec3u(0u); }
    }
    color = vec4f(mix(color.rgb, in.overpaint.rgb, in.overpaint.a), color.a);
    if (object.clipConfig.x == 0.0 && clipped(in.worldPosition, clippingMask(ids.z, ids.y))) { discard; }
    if (kind == 0u || kind == 5u || kind == 6u) {
        if (object.shading.w > 0.0) {
            let facing = pow(clamp(abs(normal.z), 0.0, 1.0), lighting.config.w);
            let opacity = select(1.0 - facing, facing, object.shading.w > 1.0);
            color.a = clamp(color.a * opacity, 0.001, 0.999);
        }
        // Model scaling cancels between radius and thickness, as in the sphere impostor shader.
        // Picking deliberately uses the theme alpha, matching the existing selection behavior.
        if (applyThickness && kind == 5u && color.a < 1.0 && object.surface.y > 0.0) {
            color.a *= min(1.0, bitcast<f32>(in.ids.w) / object.surface.y);
        }
        if (interior && (applyThickness || backDepth)) {
            if (object.surface.x == 0.0 && color.a < 1.0 && !backDepth) { discard; }
            if (object.surface.x == 2.0) { color.a = 1.0; }
        }
    }
    if (color.a <= 0.001) { discard; }
    var albedo = color.rgb;
    var emissiveColor = color.rgb * max(emission, 0.0);
    if (kind == 0u || kind == 5u || kind == 6u) {
        var base = color.rgb; var material = object.material;
        material = vec4f(mix(material.xyz, in.substance.rgb, clamp(in.substance.a, 0.0, 0.99)), material.w);
        if (interior) {
            base = mix(base, object.interiorColor.rgb, object.interiorColor.a);
            material = vec4f(mix(material.xyz, object.interiorSubstance.rgb, clamp(object.interiorSubstance.a, 0.0, 0.99)), material.w);
        }
        let frequency = object.shading.x; let amplitude = object.shading.y * material.z;
        if (frequency > 0.0 && amplitude > 0.0) {
            if (object.options.x > 0.0) { base += (fbm(in.worldPosition * frequency) - 0.5) * amplitude; }
            else {
                let t1 = normalized(cross(normal, select(vec3f(0.0, 1.0, 0.0), vec3f(1.0, 0.0, 0.0), abs(normal.x) < 0.9)));
                let t2 = cross(normal, t1);
                let distance = select(2.0 * abs(in.viewPosition.z), 2.0, camera.projection[2][3] == 0.0);
                let e = max(distance / (camera.projection[1][1] * lighting.dimensions.y), 0.01 / frequency);
                let inverseView = transpose(mat3x3f(camera.view[0].xyz, camera.view[1].xyz, camera.view[2].xyz));
                let p = in.worldPosition * frequency; let h0 = fbm(p);
                let h1 = fbm(p + inverseView * t1 * e * frequency); let h2 = fbm(p + inverseView * t2 * e * frequency);
                normal = normalized(normal - (amplitude / frequency) * ((h1 - h0) * t1 + (h2 - h0) * t2) / e);
            }
        }
        albedo = base;
        color = vec4f(shadeMaterial(base, normal, in.viewPosition.xyz, material, emission, object.options.x > 0.0, geometryRoughness), color.a);
    }
    if (marker > 0u) {
        let markerColor = select(camera.selection, camera.highlight, (marker & 1u) == 1u);
        color = vec4f(mix(color.rgb, markerColor.rgb, markerColor.a), color.a);
    }
    let coverageAlpha = color.a;
    if (!tracing && camera.fogLight.y > camera.fogLight.x) {
        let fog = smoothstep(camera.fogLight.x, camera.fogLight.y, in.viewDepth);
        if (applyThickness && camera.animation.w > 0.0) { color.a *= 1.0 - fog; }
        else { color = vec4f(mix(color.rgb, camera.background.rgb, fog), color.a); emissiveColor *= 1.0 - fog; }
    }
    var out: SurfaceFragment;
    out.normal = vec4f(select(normalized(in.viewPosition.xyz), normal, kind == 0u || kind == 5u || kind == 6u), max(emission, 0.0));
    out.albedo = vec4f(albedo, object.surface.z);
    out.color = vec4f(color.rgb * color.a, color.a);
    out.emissive = vec4f(emissiveColor * color.a, coverageAlpha);
    out.pick = vec4u(ids, bitcast<u32>(in.position.z));
    out.depth = in.position.z;
    if (kind == 3u && in.uv.x > 1.0) { out.pick = vec4u(0u); }
    return out;
}
@fragment fn fs(in: Varying, @builtin(front_facing) front: bool) -> Fragment { let result = geometryFragment(in, front, true, false, false); return Fragment(result.color, result.pick, result.emissive, result.depth); }
@fragment fn fsOpaque(in: Varying, @builtin(front_facing) front: bool) -> Fragment {
    let result = geometryFragment(in, front, true, false, false); if (result.emissive.a < 1.0) { discard; } return Fragment(result.color, result.pick, result.emissive, result.depth);
}
@fragment fn fsTransparent(in: Varying, @builtin(front_facing) front: bool) -> Fragment {
    let result = geometryFragment(in, front, true, false, false); if (result.emissive.a >= 1.0) { discard; } return Fragment(result.color, result.pick, result.emissive, result.depth);
}
struct ColorFragment { @location(0) color: vec4f, @builtin(frag_depth) depth: f32 };
struct PickFragment { @location(0) pick: vec4u, @builtin(frag_depth) depth: f32 };
struct WeightedGeometryFragment { @location(0) color: vec4f, @location(1) weight: vec4f, @location(2) emission: vec4f, @builtin(frag_depth) depth: f32 };
@fragment fn transparentColor(in: Varying, @builtin(front_facing) front: bool) -> ColorFragment { let result = geometryFragment(in, front, true, false, false); return ColorFragment(result.color, result.depth); }
@fragment fn transparentSurfaceColor(in: Varying, @builtin(front_facing) front: bool) -> ColorFragment {
    let result = geometryFragment(in, front, true, false, false); if (result.emissive.a >= 1.0) { discard; } return ColorFragment(result.color, result.depth);
}
@fragment fn weighted(in: Varying, @builtin(front_facing) front: bool) -> WeightedGeometryFragment {
    let result = geometryFragment(in, front, true, false, false);
    if ((object.info.z == 0u || object.info.z == 5u || object.info.z == 6u) && result.emissive.a >= 1.0) { discard; }
    let accum = weightedFragment(result.color, result.emissive, result.depth); return WeightedGeometryFragment(accum.color, accum.weight, accum.emission, result.depth);
}
@fragment fn colorDepthId(in: Varying, @builtin(front_facing) front: bool) -> PickFragment {
    let result = geometryFragment(in, front, true, false, false);
    if (result.pick.x == 0u || ((object.info.z == 0u || object.info.z == 5u || object.info.z == 6u) && result.emissive.a >= 1.0)) { discard; }
    return PickFragment(result.pick, result.depth);
}
fn peelSurface(in: Varying, front: bool) -> SurfaceFragment {
    let result = geometryFragment(in, front, true, false, false);
    if ((object.info.z == 0u || object.info.z == 5u || object.info.z == 6u) && result.emissive.a >= 1.0) { discard; }
    return result;
}
@fragment fn peelDepth(in: Varying, @builtin(front_facing) front: bool) -> @builtin(frag_depth) f32 {
    let result = peelSurface(in, front); return peelSearch(result.depth, in.position.xy);
}
@fragment fn peelFront(in: Varying, @builtin(front_facing) front: bool) -> PeelFragment {
    let result = peelSurface(in, front); return peelLayer(result.color, result.emissive, result.depth, in.position.xy, true);
}
@fragment fn peelBack(in: Varying, @builtin(front_facing) front: bool) -> PeelFragment {
    let result = peelSurface(in, front); return peelLayer(result.color, result.emissive, result.depth, in.position.xy, false);
}
@fragment fn outlineDepth(in: Varying, @builtin(front_facing) front: bool) -> ColorFragment {
    let result = geometryFragment(in, front, false, false, false);
    return ColorFragment(vec4f(result.depth, result.color.a, 0.0, 0.0), result.depth);
}
@fragment fn outlineSurfaceDepth(in: Varying, @builtin(front_facing) front: bool) -> ColorFragment {
    let result = geometryFragment(in, front, false, false, false);
    if (result.color.a >= 1.0) { discard; }
    return ColorFragment(vec4f(result.depth, result.color.a, 0.0, 0.0), result.depth);
}
@fragment fn pick(in: Varying, @builtin(front_facing) front: bool) -> PickFragment {
    let result = geometryFragment(in, front, false, false, false);
    if (result.color.a < lighting.dimensions.z || result.pick.x == 0u) { discard; }
    return PickFragment(result.pick, result.depth);
}

struct TracingFragment { @location(0) shaded: vec4f, @location(1) normal: vec4f, @location(2) albedo: vec4f, @builtin(frag_depth) depth: f32 };
@fragment fn tracing(in: Varying, @builtin(front_facing) front: bool) -> TracingFragment {
    let result = geometryFragment(in, front, true, true, false);
    if (result.emissive.a < 1.0) { discard; }
    return TracingFragment(result.color, result.normal, result.albedo, result.depth);
}
@fragment fn tracingBackDepth(in: Varying, @builtin(front_facing) front: bool) -> @builtin(frag_depth) f32 {
    let result = geometryFragment(in, front, false, true, true);
    if (result.emissive.a < 1.0) { discard; }
    return result.depth;
}

@group(2) @binding(0) var unmarkedDepth: texture_depth_2d;
fn fragmentMarker(in: Varying, result: SurfaceFragment) -> u32 {
    if (object.info.z != 4u) { return in.marker; }
    if (result.pick.x == 0u || result.pick.z >= u32(object.image.w)) { return 0u; }
    return u32(imageStyles[result.pick.y * u32(object.image.w) + result.pick.z].style.x);
}
@fragment fn markingDepth(in: Varying, @builtin(front_facing) front: bool) -> @builtin(frag_depth) f32 {
    let result = geometryFragment(in, front, false, false, false);
    if (fragmentMarker(in, result) > 0u) { discard; }
    return result.depth;
}
@fragment fn markingMask(in: Varying, @builtin(front_facing) front: bool) -> ColorFragment {
    let result = geometryFragment(in, front, false, false, false); let marker = fragmentMarker(in, result);
    if (marker == 0u) { discard; }
    var hidden = 1.0;
    if (camera.animation.z > 0.0) { hidden = select(0.0, 1.0, result.depth >= textureLoad(unmarkedDepth, vec2i(in.position.xy), 0)); }
    var fog = 0.0;
    if (camera.fogLight.y > camera.fogLight.x) { fog = smoothstep(camera.fogLight.x, camera.fogLight.y, in.viewDepth); }
    if (fog >= 1.0) { discard; }
    return ColorFragment(vec4f(0.0, hidden, select(0.0, 1.0, (marker & 1u) == 1u), 1.0 - fog), result.depth);
}

`;

/** Raster geometry keeps hardware depth interpolation; spheres write ray depth. */
export const rasterGeometryShader = geometryShader
    .replace('@builtin(frag_depth) depth: f32,', '')
    .replace(/, @builtin\(frag_depth\) depth: f32/g, '')
    .replace(/\b(Fragment|ColorFragment|PickFragment|WeightedGeometryFragment|TracingFragment)\(([^;\n]+), result.depth\)/g, '$1($2)');
