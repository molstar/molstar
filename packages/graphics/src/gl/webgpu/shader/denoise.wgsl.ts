/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { packing_wgsl } from './chunks/packing.wgsl.js';

/**
 * The illumination denoiser (\`smartDeNoise\` of the compose program: a 13 x 13 bilateral filter, weighted by distance,
 * normal and color) as a kernel: each workgroup loads its tile and a border of 6 texels once into workgroup memory,
 * the color and the normal decoded, so that a texel is read and its normal decoded once, not by each of its 169
 * neighbours. The same weights in the same order as the compose program; the background (no normal) keeps its color.
 */
export const denoise_wgsl = /* wgsl */ `
${packing_wgsl}

const Radius = 6;
const TileX = WorkgroupSizeX + 2u * u32(Radius);
const TileY = WorkgroupSizeY + 2u * u32(Radius);
const TileCount = TileX * TileY;
const EPSILON = 1.0e-6;
const INV_PI = 0.31830988618379067154;
const INV_SQRT_OF_2PI = 0.39894228040143267794;

/** exp(-i² / (2 sigma²)) for i in 0..6, sigma 3: the spatial weight of an offset (x, y) is w[|x|] w[|y|] */
const SpatialWeights = array<f32, 7>(1.0, 0.94595947, 0.80073740, 0.60653066, 0.41111229, 0.24935221, 0.13533528);

// normals as three arrays: an array of vec3 takes 16 bytes an element, too many for the 16 KiB of the core limits
var<workgroup> tileColor: array<vec4<f32>, TileCount>;
var<workgroup> tileNormalX: array<f32, TileCount>;
var<workgroup> tileNormalY: array<f32, TileCount>;
var<workgroup> tileNormalZ: array<f32, TileCount>;

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(
    @builtin(workgroup_id) group: vec3<u32>,
    @builtin(local_invocation_id) local: vec3<u32>,
    @builtin(local_invocation_index) index: u32,
) {
    let size = vec2<i32>(textureDimensions(tColor));
    let origin = vec2<i32>(group.xy * vec2<u32>(WorkgroupSizeX, WorkgroupSizeY)) - vec2<i32>(Radius);
    // clamped at the edges, as the compose program's reads
    for (var i = index; i < TileCount; i += WorkgroupSizeX * WorkgroupSizeY) {
        let texel = clamp(origin + vec2<i32>(i32(i % TileX), i32(i / TileX)), vec2<i32>(0), size - 1);
        tileColor[i] = textureLoad(tColor, texel, 0);
        let normal = unpackRGBAToNormal(textureLoad(tNormal, texel, 0));
        tileNormalX[i] = normal.x;
        tileNormalY[i] = normal.y;
        tileNormalZ[i] = normal.z;
    }
    workgroupBarrier();

    let pixel = origin + vec2<i32>(Radius) + vec2<i32>(local.xy);
    if (any(pixel >= size)) {
        return;
    }
    let center = (local.y + u32(Radius)) * TileX + local.x + u32(Radius);
    let centrPx = tileColor[center];
    let normal = vec3<f32>(tileNormalX[center], tileNormalY[center], tileNormalZ[center]);
    let threshold = uParams.uDenoiseThreshold;
    if (threshold == 0.0 || all(normal == vec3<f32>(0.0))) {
        textureStore(outColor, pixel, centrPx);
        return;
    }

    let sigma = 3.0;
    let invSigmaQx2PI = INV_PI * (0.5 / (sigma * sigma));
    let invThresholdSqx2 = 0.5 / (threshold * threshold);
    let invThresholdSqrt2PI = INV_SQRT_OF_2PI / threshold;

    var zBuff = 0.0;
    var colorBuff = vec4<f32>(0.0);
    for (var x = -Radius; x <= Radius; x++) {
        let wx = SpatialWeights[abs(x)] * invSigmaQx2PI;
        for (var y = -Radius; y <= Radius; y++) {
            let s = u32(i32(center) + y * i32(TileX) + x);
            var blurFactor = wx * SpatialWeights[abs(y)];
            // the normal weight to the power of 6
            let n1 = saturate(dot(normal, vec3<f32>(tileNormalX[s], tileNormalY[s], tileNormalZ[s])));
            let n2 = n1 * n1;
            blurFactor *= n2 * n2 * n2;

            let walkPx = tileColor[s];
            let dC = walkPx - centrPx;
            let deltaFactor = exp(-dot(dC, dC) * invThresholdSqx2) * invThresholdSqrt2PI * blurFactor;
            zBuff += deltaFactor;
            colorBuff += deltaFactor * walkPx;
        }
    }
    textureStore(outColor, pixel, select(colorBuff / zBuff, centrPx, zBuff <= EPSILON));
}
`;
