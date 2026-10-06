/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
export const postprocessingShader = /* wgsl */ `
struct Settings {
    inverseProjection: mat4x4f,
    dimensions: vec4f,
    viewport: vec4f,
    fxaa: vec4f,
    outline: vec4f,
    outlineOptions: vec4f,
    projection: mat4x4f,
    occlusion: vec4f,
    occlusionBlur: vec4f,
    occlusionColor: vec4f, // RGB tint, opaque AO enable flag in W.
    occlusionLevels: vec4f,
    sharpening: vec4f,
    dof: vec4f,
    dofCenter: vec4f,
    shadow: vec4f,
    ambient: vec4f,
    fog: vec4f,
    aoDimensions: vec4f,
};
@group(0) @binding(0) var<uniform> settings: Settings;
@group(0) @binding(1) var color: texture_2d<f32>;
@group(0) @binding(2) var linearSampler: sampler;
@group(0) @binding(3) var depth: texture_depth_2d;
@group(0) @binding(4) var picking: texture_2d<u32>;
@group(0) @binding(5) var occlusionTexture: texture_2d<f32>;
@group(0) @binding(6) var<storage, read> samples: array<vec4f>;
@group(0) @binding(7) var<storage, read> levels: array<vec4f>;
@group(0) @binding(8) var outlineTexture: texture_2d<f32>;
@group(0) @binding(9) var bloomTexture: texture_2d<f32>;
struct Light { direction: vec4f, color: vec4f };
@group(0) @binding(10) var<storage, read> lights: array<Light>;
// Geometry coverage is retained in the emissive attachment alpha, independently of the background.
@group(0) @binding(11) var coverageTexture: texture_2d<f32>;
@group(0) @binding(12) var transparentDepth: texture_2d<f32>;
@group(0) @binding(13) var transparentColor: texture_2d<f32>;
@group(0) @binding(14) var aoDepthPyramid: texture_2d<f32>;
@vertex fn vs(@builtin(vertex_index) i: u32) -> @builtin(position) vec4f {
    let positions = array<vec2f, 3>(vec2f(-1.0, -1.0), vec2f(3.0, -1.0), vec2f(-1.0, 3.0));
    return vec4f(positions[i], 0.0, 1.0);
}
fn insideViewport(p: vec2f) -> bool { return all(p >= settings.viewport.xy) && all(p < settings.viewport.xy + settings.viewport.zw); }
fn coord(p: vec2i) -> vec2i { return clamp(p, vec2i(settings.viewport.xy), vec2i(settings.viewport.xy + settings.viewport.zw) - 1); }
fn sampleColor(uv: vec2f) -> vec4f {
    let low = (settings.viewport.xy + 0.5) / settings.dimensions.xy;
    let high = (settings.viewport.xy + settings.viewport.zw - 0.5) / settings.dimensions.xy;
    return textureSampleLevel(color, linearSampler, clamp(uv, low, high), 0.0);
}
fn luma(rgb: vec3f) -> f32 { return sqrt(max(0.0, dot(rgb, vec3f(0.299, 0.587, 0.114)))); }
fn sampleLuma(uv: vec2f) -> f32 { return luma(sampleColor(uv).rgb); }
fn neighborLuma(uv: vec2f, x: f32, y: f32) -> f32 { return sampleLuma(uv + vec2f(x, y) / settings.dimensions.xy); }
fn quality(q: u32) -> f32 {
    if (q < 5u) { return 1.0; } if (q == 5u) { return 1.5; }
    if (q < 10u) { return 2.0; } if (q < 11u) { return 4.0; } return 8.0;
}
// Native port of Mol*'s FXAA shader (MIT, adapted from Rendu / Simon Rodriguez).
@fragment fn fxaa(@builtin(position) position: vec4f) -> @location(0) vec4f {
    if (!insideViewport(position.xy)) { return textureLoad(color, vec2i(position.xy), 0); }
    let uv = position.xy / settings.dimensions.xy; let inverseSize = 1.0 / settings.dimensions.xy;
    let center = sampleColor(uv); let lc = luma(center.rgb);
    let down = neighborLuma(uv, 0.0, -1.0); let up = neighborLuma(uv, 0.0, 1.0);
    let left = neighborLuma(uv, -1.0, 0.0); let right = neighborLuma(uv, 1.0, 0.0);
    let low = min(lc, min(min(down, up), min(left, right))); let high = max(lc, max(max(down, up), max(left, right)));
    let range = high - low;
    if (range < max(settings.fxaa.x, high * settings.fxaa.y) || range == 0.0) { return center; }
    let dl = neighborLuma(uv, -1.0, -1.0); let ur = neighborLuma(uv, 1.0, 1.0);
    let ul = neighborLuma(uv, -1.0, 1.0); let dr = neighborLuma(uv, 1.0, -1.0);
    let du = down + up; let lr = left + right;
    let leftCorners = dl + ul; let downCorners = dl + dr; let rightCorners = dr + ur; let upCorners = ur + ul;
    let horizontal = abs(-2.0 * left + leftCorners) + abs(-2.0 * lc + du) * 2.0 + abs(-2.0 * right + rightCorners);
    let vertical = abs(-2.0 * up + upCorners) + abs(-2.0 * lc + lr) * 2.0 + abs(-2.0 * down + downCorners);
    let isHorizontal = horizontal >= vertical;
    var step = select(inverseSize.x, inverseSize.y, isHorizontal);
    let l1 = select(left, down, isHorizontal); let l2 = select(right, up, isHorizontal);
    let g1 = l1 - lc; let g2 = l2 - lc; let steepest1 = abs(g1) >= abs(g2);
    let gradient = 0.25 * max(abs(g1), abs(g2));
    let average = 0.5 * (select(l2, l1, steepest1) + lc);
    if (steepest1) { step = -step; }
    let across = select(vec2f(step, 0.0), vec2f(0.0, step), isHorizontal);
    let current = uv + across * 0.5;
    let offset = select(vec2f(0.0, inverseSize.y), vec2f(inverseSize.x, 0.0), isHorizontal);
    var uv1 = current - offset; var uv2 = current + offset;
    var end1 = sampleLuma(uv1) - average; var end2 = sampleLuma(uv2) - average;
    var reached1 = abs(end1) >= gradient; var reached2 = abs(end2) >= gradient;
    if (!reached1) { uv1 -= offset; } if (!reached2) { uv2 += offset; }
    for (var i = 2u; i < u32(settings.fxaa.z); i++) {
        if (reached1 && reached2) { break; }
        if (!reached1) { end1 = sampleLuma(uv1) - average; }
        if (!reached2) { end2 = sampleLuma(uv2) - average; }
        reached1 = abs(end1) >= gradient; reached2 = abs(end2) >= gradient;
        if (!reached1) { uv1 -= offset * quality(i); } if (!reached2) { uv2 += offset * quality(i); }
    }
    let distance1 = select(uv.y - uv1.y, uv.x - uv1.x, isHorizontal);
    let distance2 = select(uv2.y - uv.y, uv2.x - uv.x, isHorizontal);
    let closer1 = distance1 < distance2; let centerSmaller = lc < average;
    let correct = select((end2 < 0.0) != centerSmaller, (end1 < 0.0) != centerSmaller, closer1);
    let pixelOffset = 0.5 - min(distance1, distance2) / max(distance1 + distance2, 0.000001);
    var finalOffset = select(0.0, pixelOffset, correct);
    let neighborhoodAverage = (2.0 * (du + lr) + leftCorners + rightCorners) / 12.0;
    let sub = clamp(abs(neighborhoodAverage - lc) / range, 0.0, 1.0);
    let subSmooth = (-2.0 * sub + 3.0) * sub * sub;
    finalOffset = max(finalOffset, subSmooth * subSmooth * settings.fxaa.w);
    return sampleColor(uv + across * finalOffset);
}
fn viewPosition(p: vec2f, d: f32) -> vec3f {
    let uv = (p - settings.viewport.xy) / settings.viewport.zw;
    let view = settings.inverseProjection * vec4f(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0, d * 2.0 - 1.0, 1.0);
    return view.xyz / view.w;
}
fn viewZ(p: vec2i, d: f32) -> f32 { return viewPosition(vec2f(p) + 0.5, d).z; }
fn outlineDepthAt(p: vec2i, transparent: bool) -> vec2f {
    if (transparent) {
        if (settings.outlineOptions.y == 0.0) { return vec2f(1.0, 0.0); }
        return textureLoad(transparentDepth, coord(p), 0).rg;
    }
    return vec2f(textureLoad(depth, coord(p), 0), 1.0);
}
fn outlineZ(p: vec2i, d: f32) -> f32 {
    if (d == 1.0) { return 2.0 * settings.dimensions.w; }
    return viewZ(p, d);
}
fn outlinePixelSize(p: vec2i, d: f32) -> f32 {
    return distance(viewPosition(vec2f(p) + 0.5, d), viewPosition(vec2f(p) + vec2f(1.5, 0.5), d)) * settings.outline.w;
}
fn layerEdge(p: vec2i, transparent: bool) -> vec2f {
    let selfDepth = outlineDepthAt(p, transparent).x; let selfZ = outlineZ(p, selfDepth);
    let size = outlinePixelSize(p, selfDepth);
    var best = vec2f(1.0, 0.0);
    for (var y = -1; y <= 1; y++) { for (var x = -1; x <= 1; x++) {
        let q = p + vec2i(x, y); let sample = outlineDepthAt(q, transparent);
        if (selfDepth > sample.x && abs(selfZ - outlineZ(q, sample.x)) > size && sample.x <= best.x) { best = sample; }
    } }
    if (best.x < 1.0 && selfDepth < 1.0) {
        let zLeft = outlineZ(p, outlineDepthAt(p - vec2i(1, 0), transparent).x); let zRight = outlineZ(p, outlineDepthAt(p + vec2i(1, 0), transparent).x);
        let zUp = outlineZ(p, outlineDepthAt(p - vec2i(0, 1), transparent).x); let zDown = outlineZ(p, outlineDepthAt(p + vec2i(0, 1), transparent).x);
        if (max(abs(zLeft + zRight - 2.0 * selfZ), abs(zUp + zDown - 2.0 * selfZ)) < size * 0.75) { best = vec2f(1.0, 0.0); }
    }
    return best;
}
fn edge(p: vec2i) -> vec4f {
    let opaque = layerEdge(p, false); var transparent = layerEdge(p, true);
    let selfOpaque = outlineDepthAt(p, false).x;
    if (transparent.x < 1.0 && abs(outlineZ(p, selfOpaque) - viewZ(p, transparent.x)) < outlinePixelSize(p, selfOpaque)) { transparent = vec2f(1.0, 0.0); }
    let flags = select(0u, 1u, opaque.x < 1.0) | select(0u, 2u, transparent.x < 1.0);
    return vec4f(transparent.y, abs(viewZ(p, opaque.x)), abs(viewZ(p, transparent.x)), f32(flags));
}
@fragment fn outlineEdges(@builtin(position) position: vec4f) -> @location(0) vec4f {
    if (!insideViewport(position.xy)) { return vec4f(0.0); }
    return edge(vec2i(position.xy));
}
fn outlineFog(distance: f32) -> f32 {
    if (settings.fog.w > settings.ambient.w) { return smoothstep(settings.ambient.w, settings.fog.w, distance); }
    return 0.0;
}
@fragment fn outline(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let p = vec2i(position.xy); let original = textureLoad(color, p, 0);
    if (!insideViewport(position.xy)) { return original; }
    let radius = i32(settings.outlineOptions.x) - 1;
    var opaqueDistance = 1e30; var transparentDistance = 1e30; var sourceAlpha = 0.0;
    for (var y = -radius; y <= radius; y++) { for (var x = -radius; x <= radius; x++) {
        let edge = textureLoad(outlineTexture, coord(p + vec2i(x, y)), 0); let flags = u32(edge.a);
        if ((flags & 1u) != 0u) { opaqueDistance = min(opaqueDistance, edge.g); }
        if ((flags & 2u) != 0u && edge.b < transparentDistance) { transparentDistance = edge.b; sourceAlpha = edge.r; }
    } }
    let hasOpaque = opaqueDistance < 1e30; let hasTransparent = transparentDistance < 1e30;
    if (!hasOpaque && !hasTransparent) { return original; }
    var translucent = textureLoad(transparentColor, p, 0);
    var opaque = clamp((original - translucent) / max(1.0 - translucent.a, 0.00001), vec4f(0.0), vec4f(1.0));
    if (hasOpaque) {
        let alpha = 1.0 - outlineFog(opaqueDistance);
        opaque = select(vec4f(mix(settings.fog.rgb, settings.outline.rgb, alpha), 1.0), vec4f(settings.outline.rgb * alpha, alpha), settings.outlineOptions.z > 0.0);
    }
    if (hasTransparent) {
        if (hasOpaque && opaqueDistance < transparentDistance) { return opaque; }
        let alpha = max(translucent.a, min(1.0, sourceAlpha * 2.0) * (1.0 - outlineFog(transparentDistance)));
        translucent = vec4f(settings.outline.rgb * alpha, alpha);
    }
    return translucent + opaque * (1.0 - translucent.a);
}
// View normals reconstructed from the closest depth slopes, as in the original SSAO pass.
fn opaqueDepthAt(p: vec2i) -> f32 { return textureLoad(depth, coord(p), 0); }
fn opaquePosition(p: vec2i) -> vec3f { return viewPosition(vec2f(p) + 0.5, opaqueDepthAt(p)); }
fn pyramidDepth(p: vec2i, level: i32) -> vec3f {
    let size = vec2f(textureDimensions(aoDepthPyramid, level));
    let low = vec2i(floor(settings.viewport.xy * size / settings.dimensions.xy));
    let high = vec2i(ceil((settings.viewport.xy + settings.viewport.zw) * size / settings.dimensions.xy)) - 1;
    return textureLoad(aoDepthPyramid, clamp(vec2i((vec2f(p) + 0.5) * size / settings.dimensions.xy), low, high), level).rgb;
}
fn aoDepthAt(p: vec2i, transparent: bool) -> f32 {
    if (transparent && settings.dimensions.z == 0.0) { return 1.0; }
    let depths = pyramidDepth(p, 0);
    return select(depths.r, depths.g, transparent);
}
fn aoPosition(p: vec2i, transparent: bool) -> vec3f { return viewPosition(vec2f(p) + 0.5, aoDepthAt(p, transparent)); }
fn viewNormal(p: vec2i, transparent: bool) -> vec3f {
    let d = aoDepthAt(p, transparent); let center = aoPosition(p, transparent);
    let step = max(vec2i(1), vec2i(round(settings.dimensions.xy / settings.aoDimensions.xy)));
    let hl = abs(2.0 * aoDepthAt(p - vec2i(step.x, 0), transparent) - aoDepthAt(p - vec2i(2 * step.x, 0), transparent) - d);
    let hr = abs(2.0 * aoDepthAt(p + vec2i(step.x, 0), transparent) - aoDepthAt(p + vec2i(2 * step.x, 0), transparent) - d);
    let vl = abs(2.0 * aoDepthAt(p - vec2i(0, step.y), transparent) - aoDepthAt(p - vec2i(0, 2 * step.y), transparent) - d);
    let vr = abs(2.0 * aoDepthAt(p + vec2i(0, step.y), transparent) - aoDepthAt(p + vec2i(0, 2 * step.y), transparent) - d);
    let horizontal = select(aoPosition(p + vec2i(step.x, 0), transparent) - center, center - aoPosition(p - vec2i(step.x, 0), transparent), hl < hr);
    let vertical = select(aoPosition(p + vec2i(0, step.y), transparent) - center, center - aoPosition(p - vec2i(0, step.y), transparent), vl < vr);
    let n = -cross(horizontal, vertical);
    return n / max(length(n), 0.000001);
}
fn randomValue(p: vec2f) -> f32 {
    let value = dot(p, vec2f(12.9898, 78.233));
    return abs(fract(sin(value - floor(value / 3.14159265) * 3.14159265) * 43758.5453));
}
fn smootherstep(x: f32) -> f32 { let v = clamp(x, 0.0, 1.0); return v*v*v*(v*(v*6.0-15.0)+10.0); }
fn layerOcclusion(position: vec4f, transparent: bool) -> f32 {
    let p = vec2i(position.xy); let d = aoDepthAt(p, transparent);
    if (!insideViewport(position.xy) || d >= 1.0 || d <= 0.0) { return 1.0; }
    let center = aoPosition(p, transparent); let normal = viewNormal(p, transparent);
    let uv = position.xy / settings.dimensions.xy;
    let random = normalize(vec3f(vec2f(randomValue(uv), randomValue(uv + vec2f(3.14159265, 2.71828))) * 2.0 - 1.0, 0.0));
    let tangentRaw = random - normal * dot(random, normal);
    let tangent = tangentRaw / max(length(tangentRaw), 0.000001);
    let basis = mat3x3f(tangent, cross(normal, tangent), normal);
    let pixelSize = distance(center, viewPosition(position.xy + vec2f(1.0, 0.0), d));
    var occlusion = 0.0;
    for (var level = 0u; level < u32(settings.occlusionLevels.x); level++) {
        let radius = levels[level].x; let bias = levels[level].y;
        if (settings.occlusionLevels.y > 0.0 && (pixelSize * settings.occlusionLevels.z > radius || pixelSize * settings.occlusionLevels.w < radius)) { continue; }
        var sum = 0.0;
        for (var i = 0u; i < u32(settings.occlusion.z); i++) {
            let samplePosition = center + (basis * samples[i].xyz) * radius;
            let projected = settings.projection * vec4f(samplePosition, 1.0);
            let sampleUv = vec2f(projected.x / projected.w * 0.5 + 0.5, 0.5 - projected.y / projected.w * 0.5);
            let q = coord(vec2i(settings.viewport.xy + sampleUv * settings.viewport.zw));
            var mip = 0;
            if (settings.occlusionLevels.y > 0.0) {
                let offset = distance(vec2f(q) / settings.dimensions.xy, uv);
                mip = select(select(0, 1, offset > 0.05), 2, offset > 0.1);
                mip = min(mip, i32(textureNumLevels(aoDepthPyramid)) - 1);
            }
            let depths = pyramidDepth(q, mip);
            let sampleDepth = depths.r;
            var contribution = 0.0;
            if (sampleDepth < 1.0) {
                let z = viewPosition(vec2f(q) + 0.5, sampleDepth).z;
                contribution = select(0.0, smootherstep(radius / max(abs(center.z - z), 0.000001)), z >= samplePosition.z + 0.025) * bias;
            }
            if (settings.dimensions.z > 0.0) {
                let depthAlpha = depths.gb;
                if (depthAlpha.x < 1.0) {
                    let z = viewPosition(vec2f(q) + 0.5, depthAlpha.x).z;
                    let value = select(0.0, smootherstep(radius / max(abs(center.z - z), 0.000001)), z >= samplePosition.z + 0.025) * bias * depthAlpha.y;
                    contribution = max(contribution, value);
                }
            }
            sum += contribution;
        }
        occlusion = max(occlusion, sum / settings.occlusion.z);
    }
    return clamp(1.0 - settings.occlusion.y * occlusion, 0.01, 1.0);
}
fn aoPixel(p: vec2i) -> vec2i {
    let low = vec2i(floor(settings.viewport.xy * settings.aoDimensions.xy / settings.dimensions.xy));
    let high = vec2i(ceil((settings.viewport.xy + settings.viewport.zw) * settings.aoDimensions.xy / settings.dimensions.xy)) - 1;
    return clamp(vec2i((vec2f(p) + 0.5) * settings.aoDimensions.xy / settings.dimensions.xy), low, high);
}
fn sampleOcclusion(p: vec2i) -> vec2f {
    let low = (floor(settings.viewport.xy * settings.aoDimensions.xy / settings.dimensions.xy) + 0.5) / settings.aoDimensions.xy;
    let high = (ceil((settings.viewport.xy + settings.viewport.zw) * settings.aoDimensions.xy / settings.dimensions.xy) - 0.5) / settings.aoDimensions.xy;
    return textureSampleLevel(occlusionTexture, linearSampler, clamp((vec2f(p) + 0.5) / settings.dimensions.xy + settings.aoDimensions.zw, low, high), 0.0).rg;
}
fn blurOcclusion(position: vec4f, transparent: bool) -> f32 {
    let p = vec2i(position.xy); let d = aoDepthAt(p, transparent);
    if (!insideViewport(position.xy) || d >= 1.0 || d <= 0.0) { return 1.0; }
    let selfZ = viewZ(p, d); let pixelSize = distance(viewPosition(position.xy, d), viewPosition(position.xy + vec2f(settings.dimensions.x / settings.aoDimensions.x, 0.0), d));
    let halfKernel = i32(settings.occlusionBlur.x) / 2; let sigma = settings.occlusionBlur.x / 3.0;
    var sum = 0.0; var weights = 0.0;
    for (var i = -halfKernel; i <= halfKernel; i++) {
        if (abs(i) > 1 && f32(abs(i)) * pixelSize > 0.8) { continue; }
        let q = vec2i(position.xy + settings.occlusionBlur.zw * f32(i) * settings.dimensions.xy / settings.aoDimensions.xy);
        if (!insideViewport(vec2f(q) + 0.5)) { continue; }
        let sampleDepth = aoDepthAt(q, transparent);
        if (sampleDepth >= 1.0 || sampleDepth <= 0.0 || abs(selfZ - viewZ(q, sampleDepth)) >= settings.occlusionBlur.y) { continue; }
        let weight = exp(-f32(i * i) / (2.0 * sigma * sigma));
        sum += select(textureLoad(occlusionTexture, aoPixel(q), 0).r, textureLoad(occlusionTexture, aoPixel(q), 0).g, transparent) * weight; weights += weight;
    }
    return select(1.0, sum / max(weights, 0.000001), weights > 0.0);
}
@fragment fn ssao(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let source = vec4f(position.xy * settings.dimensions.xy / settings.aoDimensions.xy, position.zw);
    var opaque = 1.0;
    if (settings.occlusionColor.w > 0.0) { opaque = layerOcclusion(source, false); }
    return vec4f(opaque, layerOcclusion(source, true), 0.0, 1.0);
}
@fragment fn ssaoBlur(@builtin(position) position: vec4f) -> @location(0) vec4f { let source = vec4f(position.xy * settings.dimensions.xy / settings.aoDimensions.xy, position.zw); return vec4f(blurOcclusion(source, false), blurOcclusion(source, true), 0.0, 1.0); }
fn recoverOpaque(original: vec4f, translucent: vec4f) -> vec4f { return clamp((original - translucent) / max(1.0 - translucent.a, 0.00001), vec4f(0.0), vec4f(1.0)); }
fn sourceCoverage(alpha: f32, fog: f32) -> f32 {
    // Environmental fog is already in the geometry alpha; do not apply it twice to AO tint.
    return select(alpha, alpha / max(1.0 - fog, 0.00001), settings.outlineOptions.w > 0.0);
}
fn occludedTransparent(p: vec2i) -> vec4f {
    let original = textureLoad(transparentColor, p, 0); let d = textureLoad(transparentDepth, coord(p), 0).r;
    if (d >= 1.0 || settings.dimensions.z == 0.0) { return original; }
    let fog = outlineFog(abs(viewZ(p, d)));
    let visibility = sampleOcclusion(p).g;
    return vec4f(mix(settings.occlusionColor.rgb * sourceCoverage(original.a, fog) * (1.0 - fog), original.rgb, visibility), original.a);
}
@fragment fn ssaoTransparentColor(@builtin(position) position: vec4f) -> @location(0) vec4f { return occludedTransparent(vec2i(position.xy)); }
@fragment fn ssaoCompose(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let p = vec2i(position.xy); let original = textureLoad(color, p, 0);
    let translucent = textureLoad(transparentColor, p, 0); var opaque = recoverOpaque(original, translucent);
    let d = opaqueDepthAt(p);
    if (d < 1.0 && settings.occlusionColor.w > 0.0) {
        let fog = outlineFog(abs(viewZ(p, d))); let visibility = sampleOcclusion(p).r;
        let tint = select(mix(settings.occlusionColor.rgb, settings.fog.rgb, fog), settings.occlusionColor.rgb * (1.0 - fog), settings.outlineOptions.z > 0.0);
        opaque = vec4f(mix(tint * sourceCoverage(opaque.a, fog), opaque.rgb, visibility), opaque.a);
    }
    let transparent = occludedTransparent(p);
    return vec4f((transparent + opaque * (1.0 - transparent.a)).rgb, original.a);
}

// Screen-space directional shadows, matching the existing depth ray march.
fn shadowVisibility(origin: vec3f, direction: vec3f) -> f32 {
    let step = -direction * settings.shadow.x / settings.shadow.z;
    var ray = origin;
    for (var i = 0u; i < u32(settings.shadow.z); i++) {
        ray += step;
        let clip = settings.projection * vec4f(ray, 1.0);
        if (clip.w <= 0.0) { return 1.0; }
        let uv = vec2f(clip.x / clip.w * 0.5 + 0.5, 0.5 - clip.y / clip.w * 0.5);
        let point = settings.viewport.xy + uv * settings.viewport.zw;
        if (!insideViewport(point)) { return 1.0; }
        let p = coord(vec2i(point)); let d = opaqueDepthAt(p);
        if (d < 1.0 && ray.z - viewZ(p, d) < settings.shadow.y) {
            let fade = max(12.0 * abs(uv - 0.5) - 5.0, vec2f(0.0));
            return 1.0 - clamp(1.0 - dot(fade, fade), 0.0, 1.0);
        }
    }
    return 1.0;
}
@fragment fn shadowCompose(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let p = vec2i(position.xy); let original = textureLoad(color, p, 0);
    let d = opaqueDepthAt(p);
    if (!insideViewport(position.xy) || d >= 1.0 || settings.shadow.x <= 0.0) { return original; }
    let origin = viewPosition(position.xy, d);
    var total = length(settings.ambient.rgb); var visible = total;
    for (var i = 0u; i < u32(settings.shadow.w); i++) {
        let weight = length(lights[i].color.rgb);
        total += weight; visible += weight * shadowVisibility(origin, lights[i].direction.xyz);
    }
    let factor = select(1.0, visible / max(total, 0.000001), total > 0.0);
    var fog = 0.0;
    if (settings.fog.w > settings.ambient.w) { fog = smoothstep(settings.ambient.w, settings.fog.w, abs(origin.z)); }
    // Fog has already been applied to the material; leave its contribution unshadowed.
    let translucent = textureLoad(transparentColor, p, 0); let opaque = recoverOpaque(original, translucent);
    let shaded = vec4f(mix(settings.fog.rgb * fog * opaque.a, opaque.rgb, factor), opaque.a);
    return vec4f((translucent + shaded * (1.0 - translucent.a)).rgb, original.a);
}

// FidelityFX RCAS, matching the existing contrast-adaptive sharpening pass.
@fragment fn sharpen(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let p = vec2i(position.xy); let center = textureLoad(color, p, 0);
    if (!insideViewport(position.xy) || center.a == 0.0) { return center; }
    let b = textureLoad(color, coord(p + vec2i(0, -1)), 0).rgb;
    let d = textureLoad(color, coord(p + vec2i(-1, 0)), 0).rgb;
    let f = textureLoad(color, coord(p + vec2i(1, 0)), 0).rgb;
    let h = textureLoad(color, coord(p + vec2i(0, 1)), 0).rgb;
    let e = center.rgb;
    let low = min(b, min(f, h)); let high = max(b, max(f, h));
    let hitMin = low / max(4.0 * high, vec3f(0.000001));
    let hitMax = (1.0 - high) / min(4.0 * low - 4.0, vec3f(-0.000001));
    let lobeRGB = max(-hitMin, hitMax);
    var lobe = max(-0.1875, min(max(lobeRGB.r, max(lobeRGB.g, lobeRGB.b)), 0.0)) * settings.sharpening.x;
    if (settings.sharpening.y > 0.0) {
        let coefficients = vec3f(0.5, 1.0, 0.5);
        let bl = dot(b, coefficients); let dl = dot(d, coefficients); let el = dot(e, coefficients);
        let fl = dot(f, coefficients); let hl = dot(h, coefficients);
        let range = max(max(bl, dl), max(el, max(fl, hl))) - min(min(bl, dl), min(el, min(fl, hl)));
        let noise = clamp(abs(0.25 * (bl + dl + fl + hl) - el) / max(range, 0.000001), 0.0, 1.0);
        lobe *= 1.0 - 0.5 * noise;
    }
    let sharpened = (lobe * (b + d + h + f) + e) / (4.0 * lobe + 1.0);
    return vec4f(clamp(sharpened, vec3f(0.0), vec3f(center.a)), center.a);
}
fn closestDepth(p: vec2i) -> f32 {
    let q = coord(p); var value = textureLoad(depth, q, 0); let pick = textureLoad(picking, q, 0);
    if (pick.x > 0u) { value = min(value, bitcast<f32>(pick.w)); }
    return value;
}
fn circleOfConfusion(p: vec2f) -> f32 {
    let view = viewPosition(p, closestDepth(vec2i(p)));
    var value = (abs(view.z) - settings.dof.z) / max(settings.dof.w, 0.000001);
    if (settings.dofCenter.w > 0.0) { value = distance(view, settings.dofCenter.xyz) / max(settings.dof.w, 0.000001); }
    return clamp(value, -1.0, 1.0);
}
@fragment fn dof(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let p = position.xy; let original = textureLoad(color, vec2i(p), 0);
    if (!insideViewport(p) || original.a == 0.0) { return original; }
    let coc = smoothstep(0.0, 1.0, abs(circleOfConfusion(p)));
    if (coc == 0.0) { return original; }
    var weighted = vec4f(0.0); var count = 0.0;
    for (var x = 0u; x < u32(settings.dof.x); x++) { for (var y = 0u; y < u32(settings.dof.x); y++) {
        let offset = (vec2f(f32(x), f32(y)) - settings.dof.x * 0.5) * settings.dof.y;
        let samplePoint = clamp(p + offset, settings.viewport.xy + 0.5, settings.viewport.xy + settings.viewport.zw - 0.5);
        let weight = smoothstep(0.0, 1.0, abs(circleOfConfusion(samplePoint)));
        weighted += sampleColor(samplePoint / settings.dimensions.xy) * weight; count += weight;
    } }
    if (count <= 0.0 || weighted.a <= 0.0) { return original; }
    // Preserve alpha while filtering straight color, then premultiply the result.
    let blurred = weighted.rgb / weighted.a * original.a;
    return vec4f(mix(original.rgb, blurred, coc), original.a);
}
@fragment fn bloomCompose(@builtin(position) position: vec4f) -> @location(0) vec4f {
    let p = vec2i(position.xy); let original = textureLoad(color, p, 0);
    let glow = textureLoad(bloomTexture, p, 0).rgb;
    if (all(glow == vec3f(0.0))) { return original; }
    let rgb = original.rgb + glow;
    return vec4f(rgb, min(max(original.a, max(rgb.r, max(rgb.g, rgb.b))), 1.0));
}
`;
