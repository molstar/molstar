/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/**
 * Screen space shadows (ray marched towards each light through the depth), bounded by the
 * \`oSteps\` and \`oLightCount\` override constants. Depth textures are always available on WebGPU.
 */
export const shadows_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tDepth')}

override oSteps: i32 = 1;
override oLightCount: i32 = 1;

// map full [0,1] coords -> bounds [0,1] coords
fn fullCoordsToBoundsCoords(fullCoords: vec2<f32>) -> vec2<f32> {
    let b = uObject.uBounds;
    return (fullCoords - b.xy) / (b.zw - b.xy);
}

// map bounds [0,1] coords -> full [0,1] coords
fn boundsCoordsToFullCoords(boundsCoords: vec2<f32>) -> vec2<f32> {
    let b = uObject.uBounds;
    return b.xy + (b.zw - b.xy) * boundsCoords;
}

fn isBackground(depth: f32) -> bool {
    return depth == 0.0;
}

fn outsideBounds(p: vec2<f32>) -> bool {
    let b = uObject.uBounds;
    return p.x < b.x || p.y < b.y || p.x > b.z || p.y > b.w;
}

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

fn getDepth(coords: vec2<f32>) -> f32 {
    return read_tDepth(coords).r;
}

fn screenFade(coords: vec2<f32>) -> f32 {
    let c = fullCoordsToBoundsCoords(coords);
    let fade = max(12.0 * abs(c - 0.5) - 5.0, vec2<f32>(0.0));
    return clamp(1.0 - dot(fade, fade), 0.0, 1.0);
}

// based on https://panoskarabelas.com/posts/screen_space_shadows/
fn screenSpaceShadow(position: vec3<f32>, lightDirection: vec3<f32>, lightColor: vec3<f32>, stepLength: f32) -> vec3<f32> {
    // ray position and direction (in view-space)
    var rayPos = position;
    let rayStep = -lightDirection * stepLength;

    // ray march towards the light
    for (var i = 0; i < oSteps; i++) {
        rayPos += rayStep;

        var rayCoords = uObject.uProjection * vec4<f32>(rayPos, 1.0);
        let ndc = (rayCoords.xy / rayCoords.w) * vec2<f32>(0.5, -0.5) + 0.5;
        let coords = boundsCoordsToFullCoords(ndc);

        if (outsideBounds(coords)) {
            return lightColor;
        }

        // difference between the ray's and the camera's depth
        let viewZ = getViewZ(getDepth(coords));
        let zDelta = rayPos.z - viewZ;

        if (zDelta < uObject.uTolerance) {
            // fade out as the edges of the screen are approached
            return mix(vec3<f32>(0.0), lightColor, 1.0 - screenFade(coords));
        }
    }

    return lightColor;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let selfCoords = fragmentUv(in.position, uObject.uTexSize);
    let selfDepth = getDepth(selfCoords);

    if (isBackground(selfDepth)) {
        return vec4<f32>(0.0);
    }

    let selfBoundsCoords = fullCoordsToBoundsCoords(selfCoords);
    if (any(selfBoundsCoords < vec2<f32>(0.0)) || any(selfBoundsCoords > vec2<f32>(1.0))) {
        return vec4<f32>(0.0);
    }

    let selfViewPos = screenSpaceToViewSpace(vec3<f32>(selfBoundsCoords, selfDepth), uObject.uInvProjection);
    let stepLength = uObject.uMaxDistance / f32(oSteps);

    var l = length(uObject.uAmbientColor);
    var a = l;
    #if dLightCount != 0
        for (var i = 0; i < oLightCount; i++) {
            let s = screenSpaceShadow(selfViewPos, uObject.uLightDirection[i], uObject.uLightColor[i], stepLength);
            l += length(s);
            a += length(uObject.uLightColor[i]);
        }
    #endif

    // no light: no shadow
    return vec4<f32>(select(1.0, l / a, a > 0.0));
}
`;
