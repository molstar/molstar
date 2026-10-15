/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/**
 * Fast approximate anti-aliasing (FXAA 3.11 by Timothy Lottes, as explained in
 * http://blog.simonrodriguez.fr/articles/2016/07/implementing_fxaa.html), in GL window space; thresholds, iterations
 * and subpixel quality are override constants
 */
export const fxaa_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tColor')}

override oEdgeThresholdMin: f32 = 0.0312;
override oEdgeThresholdMax: f32 = 0.063;
override oIterations: i32 = 12;
override oSubpixelQuality: f32 = 0.3;

fn quality(q: i32) -> f32 {
    if (q < 5) {
        return 1.0;
    } else if (q == 5) {
        return 1.5;
    } else if (q < 10) {
        return 2.0;
    } else if (q < 11) {
        return 4.0;
    }
    return 8.0;
}

fn rgb2luma(rgb: vec3<f32>) -> f32 {
    return sqrt(dot(rgb, vec3<f32>(0.299, 0.587, 0.114)));
}

fn sampleLuma(uv: vec2<f32>) -> f32 {
    return rgb2luma(read_tColor(uv).rgb);
}

fn sampleLumaOffset(uv: vec2<f32>, offsetX: f32, offsetY: f32) -> f32 {
    return sampleLuma(uv + uObject.uTexSizeInv * vec2<f32>(offsetX, offsetY));
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let inverseScreenSize = uObject.uTexSizeInv;
    let coords = glCoords(in.position, 1.0 / inverseScreenSize);

    let colorCenter = read_tColor(coords);

    // luma at the current fragment and its four direct neighbours
    let lumaCenter = rgb2luma(colorCenter.rgb);
    let lumaDown = sampleLumaOffset(coords, 0.0, -1.0);
    let lumaUp = sampleLumaOffset(coords, 0.0, 1.0);
    let lumaLeft = sampleLumaOffset(coords, -1.0, 0.0);
    let lumaRight = sampleLumaOffset(coords, 1.0, 0.0);

    let lumaMin = min(lumaCenter, min(min(lumaDown, lumaUp), min(lumaLeft, lumaRight)));
    let lumaMax = max(lumaCenter, max(max(lumaDown, lumaUp), max(lumaLeft, lumaRight)));
    let lumaRange = lumaMax - lumaMin;

    // not on an edge (or in a really dark area): no anti-aliasing
    if (lumaRange < max(oEdgeThresholdMin, lumaMax * oEdgeThresholdMax)) {
        return colorCenter;
    }

    // the four corners
    let lumaDownLeft = sampleLumaOffset(coords, -1.0, -1.0);
    let lumaUpRight = sampleLumaOffset(coords, 1.0, 1.0);
    let lumaUpLeft = sampleLumaOffset(coords, -1.0, 1.0);
    let lumaDownRight = sampleLumaOffset(coords, 1.0, -1.0);

    let lumaDownUp = lumaDown + lumaUp;
    let lumaLeftRight = lumaLeft + lumaRight;

    let lumaLeftCorners = lumaDownLeft + lumaUpLeft;
    let lumaDownCorners = lumaDownLeft + lumaDownRight;
    let lumaRightCorners = lumaDownRight + lumaUpRight;
    let lumaUpCorners = lumaUpRight + lumaUpLeft;

    // gradient estimation along the horizontal and vertical axis
    let edgeHorizontal = abs(-2.0 * lumaLeft + lumaLeftCorners) + abs(-2.0 * lumaCenter + lumaDownUp) * 2.0 + abs(-2.0 * lumaRight + lumaRightCorners);
    let edgeVertical = abs(-2.0 * lumaUp + lumaUpCorners) + abs(-2.0 * lumaCenter + lumaLeftRight) * 2.0 + abs(-2.0 * lumaDown + lumaDownCorners);

    let isHorizontal = edgeHorizontal >= edgeVertical;

    // one pixel step across the edge
    var stepLength = select(inverseScreenSize.x, inverseScreenSize.y, isHorizontal);

    let luma1 = select(lumaLeft, lumaDown, isHorizontal);
    let luma2 = select(lumaRight, lumaUp, isHorizontal);
    let gradient1 = luma1 - lumaCenter;
    let gradient2 = luma2 - lumaCenter;

    let is1Steepest = abs(gradient1) >= abs(gradient2);
    let gradientScaled = 0.25 * max(abs(gradient1), abs(gradient2));

    var lumaLocalAverage = 0.0;
    if (is1Steepest) {
        stepLength = -stepLength;
        lumaLocalAverage = 0.5 * (luma1 + lumaCenter);
    } else {
        lumaLocalAverage = 0.5 * (luma2 + lumaCenter);
    }

    // shift by half a pixel towards the edge
    var currentUv = coords;
    if (isHorizontal) {
        currentUv.y += stepLength * 0.5;
    } else {
        currentUv.x += stepLength * 0.5;
    }

    // explore along the edge on both sides
    let offset = select(vec2<f32>(0.0, inverseScreenSize.y), vec2<f32>(inverseScreenSize.x, 0.0), isHorizontal);
    var uv1 = currentUv - offset * quality(0);
    var uv2 = currentUv + offset * quality(0);

    var lumaEnd1 = sampleLuma(uv1) - lumaLocalAverage;
    var lumaEnd2 = sampleLuma(uv2) - lumaLocalAverage;

    var reached1 = abs(lumaEnd1) >= gradientScaled;
    var reached2 = abs(lumaEnd2) >= gradientScaled;
    var reachedBoth = reached1 && reached2;

    if (!reached1) {
        uv1 -= offset * quality(1);
    }
    if (!reached2) {
        uv2 += offset * quality(1);
    }

    if (!reachedBoth) {
        for (var i = 2; i < oIterations; i++) {
            if (!reached1) {
                lumaEnd1 = sampleLuma(uv1) - lumaLocalAverage;
            }
            if (!reached2) {
                lumaEnd2 = sampleLuma(uv2) - lumaLocalAverage;
            }
            reached1 = abs(lumaEnd1) >= gradientScaled;
            reached2 = abs(lumaEnd2) >= gradientScaled;
            reachedBoth = reached1 && reached2;

            if (!reached1) {
                uv1 -= offset * quality(i);
            }
            if (!reached2) {
                uv2 += offset * quality(i);
            }
            if (reachedBoth) {
                break;
            }
        }
    }

    // distances to the ends of the edge
    let distance1 = select(coords.y - uv1.y, coords.x - uv1.x, isHorizontal);
    let distance2 = select(uv2.y - coords.y, uv2.x - coords.x, isHorizontal);

    let isDirection1 = distance1 < distance2;
    let distanceFinal = min(distance1, distance2);
    let edgeThickness = distance1 + distance2;

    // the luma deltas at the ends should vary as the center luma does
    let isLumaCenterSmaller = lumaCenter < lumaLocalAverage;
    let correctVariation1 = (lumaEnd1 < 0.0) != isLumaCenterSmaller;
    let correctVariation2 = (lumaEnd2 < 0.0) != isLumaCenterSmaller;
    let correctVariation = select(correctVariation2, correctVariation1, isDirection1);

    let pixelOffset = -distanceFinal / edgeThickness + 0.5;
    var finalOffset = select(0.0, pixelOffset, correctVariation);

    // sub-pixel shifting from the weighted average luma of the 3x3 neighborhood
    let lumaAverage = (1.0 / 12.0) * (2.0 * (lumaDownUp + lumaLeftRight) + lumaLeftCorners + lumaRightCorners);
    let subPixelOffset1 = clamp(abs(lumaAverage - lumaCenter) / lumaRange, 0.0, 1.0);
    let subPixelOffset2 = (-2.0 * subPixelOffset1 + 3.0) * subPixelOffset1 * subPixelOffset1;
    let subPixelOffsetFinal = subPixelOffset2 * subPixelOffset2 * oSubpixelQuality;

    finalOffset = max(finalOffset, subPixelOffsetFinal);

    var finalUv = coords;
    if (isHorizontal) {
        finalUv.y += finalOffset * stepLength;
    } else {
        finalUv.x += finalOffset * stepLength;
    }
    return read_tColor(finalUv);
}
`;
