/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';
import { outline_read_wgsl, pack_2x4_wgsl } from './outlines.wgsl.js';
import type { DefineLookup } from '../../renderable/schema.js';

const TraceTextures = ['tColor', 'tNormal', 'tShaded', 'tThickness', 'tDepth'];

/**
 * Screen space path tracing of the opaque scene, blended into the average of the previous frames (the pass's blend
 * constant, see `TracingPass`) with the frame's weight in alpha; the random numbers are seeded with the GL fragment
 * coordinates, so that the noise matches WebGL
 *
 * parts adapted from
 * - https://blog.demofox.org/2020/05/25/casual-shadertoy-path-tracing-1-basic-camera-diffuse-emissive/
 * - https://github.com/0beqz/realism-effects/blob/v2-debug/src/ssgi/shader/ssgi.frag
 */
export const trace_wgsl = /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${TraceTextures.map(texture_read_wgsl).join('\n')}

override oLightCount: i32 = 0;
override oBounces: i32 = 4;

const TWO_PI = 6.28318530718;

// after a hit, it moves the ray this far along the normal away from a surface.
// Helps prevent incorrect intersections when rays bounce off of objects.
const RayPosNormalNudge = 0.0001;

// https://www.pcg-random.org/
// https://jcgt.org/published/0009/03/02/
fn pcg(seed: ptr<function, u32>) -> u32 {
    *seed = *seed * 747796405u + 2891336453u;
    let word = ((*seed >> ((*seed >> 28u) + 4u)) ^ *seed) * 277803737u;
    return (word >> 22u) ^ word;
}

fn randomFloat(state: ptr<function, u32>) -> f32 {
    return f32(pcg(state)) / 4294967296.0;
}

fn randomUnitVector(state: ptr<function, u32>) -> vec3<f32> {
    let z = randomFloat(state) * 2.0 - 1.0;
    let a = randomFloat(state) * TWO_PI;
    let r = sqrt(1.0 - z * z);
    let x = r * cos(a);
    let y = r * sin(a);
    return vec3<f32>(x, y, z);
}

struct RayHitInfo {
    missed: bool,
    position: vec3<f32>,
    normal: vec3<f32>,
    color: vec3<f32>,
    emissive: vec3<f32>,
};

//

// inset by half a texel to avoid sampling texels outside the viewport
fn clampToBounds(coords: vec2<f32>) -> vec2<f32> {
    let halfTexel = 0.5 / uObject.uTexSize;
    return clamp(coords, uObject.uBounds.xy + halfTexel, uObject.uBounds.zw - halfTexel);
}

fn getDepth(coords: vec2<f32>) -> f32 {
    return read_tDepth(clampToBounds(coords)).r;
}

fn getThickness(coords: vec2<f32>) -> f32 {
    return unpackRGBAToDepth(read_tThickness(clampToBounds(coords)));
}

fn isBackground(depth: f32) -> bool {
    return depth == 0.0;
}

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

fn viewSpaceToScreenSpace(position: vec3<f32>) -> vec2<f32> {
    let projectedCoord = uObject.uProjection * vec4<f32>(position, 1.0);
    // NDC to window uv (top-left origin)
    return (projectedCoord.xy / projectedCoord.w) * vec2<f32>(0.5, -0.5) + 0.5;
}

fn getRayOffset(position: vec3<f32>) -> f32 {
    let projectedCoord = uObject.uProjection * vec4<f32>(position, 1.0);
    // window coordinates (depth already in [0, 1])
    let screenCoord = vec3<f32>((projectedCoord.xy / projectedCoord.w) * vec2<f32>(0.5, -0.5) + 0.5, projectedCoord.z / projectedCoord.w);
    let viewportWidth = uObject.uTexSize.x * (uObject.uBounds.z - uObject.uBounds.x);
    let adjacentPosition = screenSpaceToViewSpace(vec3<f32>(screenCoord.xy + vec2<f32>(1.0 / viewportWidth, 0.0), screenCoord.z), uObject.uInvProjection);
    return max(distance(position, adjacentPosition) * 0.1, 0.001);
}

fn binarySearch(dir: ptr<function, vec3<f32>>, hitPos: ptr<function, vec3<f32>>) -> vec2<f32> {
    *dir *= 0.5;
    *hitPos -= *dir;

    for (var i = 0; i < uObject.uRefineSteps; i++) {
        let coords = viewSpaceToScreenSpace(*hitPos);
        let z = getViewZ(getDepth(coords));
        let rayHitDepthDifference = z - (*hitPos).z;

        *dir *= 0.5;
        if (rayHitDepthDifference >= 0.0) {
            *hitPos -= *dir;
        } else {
            *hitPos += *dir;
        }
    }

    return clampToBounds(viewSpaceToScreenSpace(*hitPos));
}

fn calculateGrowthFactor(begin: f32, end: f32, steps: f32) -> f32 {
    return pow(end / begin, 1.0 / steps);
}

fn rayMarch(dirIn: vec3<f32>, thickness: f32, hitPos: ptr<function, vec3<f32>>, missed: ptr<function, bool>) -> vec2<f32> {
    var coords = vec2<f32>(0.0);

    let begin = getRayOffset(*hitPos);
    var dir = dirIn * begin;
    *missed = false;
    let gf = calculateGrowthFactor(begin, uObject.uRayDistance, f32(uObject.uSteps));

    for (var i = 1; i < uObject.uSteps; i++) {
        *hitPos += dir;
        let stepZ = abs(dir.z);
        dir *= gf;

        coords = clampToBounds(viewSpaceToScreenSpace(*hitPos));
        let z = getViewZ(getDepth(coords));
        let rayHitDepthDifference = z - (*hitPos).z;

        if (rayHitDepthDifference >= 0.0) {
            var t = thickness;
            if (t == 0.0) {
                #if defined(dThicknessMode_auto)
                    t = max(uObject.uMinThickness, (z - getViewZ(getThickness(coords))) * uObject.uThicknessFactor * read_tColor(coords).a);
                #else
                    t = uObject.uThickness;
                #endif
            }

            if (rayHitDepthDifference < max(t, stepZ)) {
                if (uObject.uRefineSteps == 0) {
                    return coords;
                } else {
                    return binarySearch(&dir, hitPos);
                }
            }
        }
    }

    *missed = true;

    return coords;
}

//

fn trace(rayPos: vec3<f32>, rayDir: vec3<f32>) -> RayHitInfo {
    var hitPos = rayPos;
    var missed = false;
    let coords = rayMarch(rayDir, 0.0, &hitPos, &missed);

    var hitInfo: RayHitInfo;
    hitInfo.missed = missed;
    hitInfo.position = hitPos;
    hitInfo.normal = -unpackRGBAToNormal(read_tNormal(coords));
    hitInfo.color = read_tColor(coords).rgb;
    hitInfo.emissive = read_tColor(coords).rgb * read_tNormal(coords).a * 2.0;

    // the background, whatever color it is cleared to: no emitter or surface normal (as transparent, cleared to zero)
    if (isBackground(getDepth(coords))) {
        hitInfo.emissive = vec3<f32>(0.0);
        hitInfo.normal = vec3<f32>(0.0);
    }
    return hitInfo;
}

var<private> viewPos: vec3<f32>;
/** uv of the fragment (top-left origin) */
var<private> fragCoords: vec2<f32>;

fn colorForRay(startRayPos: vec3<f32>, startRayDir: vec3<f32>, rngState: ptr<function, u32>) -> vec3<f32> {
    var ret = vec3<f32>(0.0);

    var throughput = vec3<f32>(1.0);
    var rayPos = startRayPos;
    var rayDir = startRayDir;

    var hitInfo: RayHitInfo;
    var prevHitInfo: RayHitInfo;

    for (var bounceIndex = 0; bounceIndex <= oBounces; bounceIndex++) {
        // shoot a ray out into the world
        if (bounceIndex == 0) {
            let depth = getDepth(fragCoords);

            hitInfo.missed = false;
            hitInfo.position = screenSpaceToViewSpace(vec3<f32>(fragCoords, depth), uObject.uInvProjection);
            hitInfo.normal = -unpackRGBAToNormal(read_tNormal(fragCoords));
            hitInfo.color = read_tShaded(fragCoords).rgb;
            hitInfo.emissive = read_tColor(fragCoords).rgb * read_tNormal(fragCoords).a;

            // shadow
            #if defined(dShadowEnable)
                #if dLightCount != 0
                    // attenuate by the fraction of unoccluded irradiance, weighting each
                    // light by its actual contribution so back-facing lights don't darken
                    var directLight = uObject.uAmbientColor;
                    var fullLight = uObject.uAmbientColor;
                    for (var i = 0; i < oLightCount; i++) {
                        let lightDirection = uObject.uLightDirection[i];
                        let ndotl = saturate(dot(hitInfo.normal, -lightDirection));
                        let irradiance = ndotl * uObject.uLightColor[i];
                        fullLight += irradiance;
                        if (ndotl > 0.0) {
                            var missed = false;
                            var hitPos = viewPos + hitInfo.normal * RayPosNormalNudge;
                            hitPos += -lightDirection * randomFloat(rngState);
                            rayMarch(-lightDirection + randomUnitVector(rngState) * uObject.uShadowSoftness, uObject.uShadowThickness, &hitPos, &missed);
                            if (missed) {
                                directLight += irradiance;
                            }
                        }
                    }
                    hitInfo.color *= directLight / max(fullLight, vec3<f32>(0.0001));
                #endif
            #endif

            if (all(hitInfo.normal == vec3<f32>(0.0))) {
                hitInfo.missed = true;
            }
        } else {
            prevHitInfo = hitInfo;
            hitInfo = trace(rayPos, rayDir);
        }

        // if the ray missed, we are done
        if (hitInfo.missed) {
            var escapeColor = prevHitInfo.color;
            if (bounceIndex > 1) {
                var accIrradiance = uObject.uAmbientColor;
                for (var i = 0; i < oLightCount; i++) {
                    accIrradiance += saturate(dot(prevHitInfo.normal, -uObject.uLightDirection[i])) * uObject.uLightColor[i];
                }
                escapeColor = min(prevHitInfo.color * accIrradiance, vec3<f32>(0.99)) * uObject.uExposure;
            }
            ret += escapeColor * throughput;
            break;
        }

        // add emissive light
        ret += hitInfo.emissive * throughput;

        // update the ray position
        rayPos = hitInfo.position + hitInfo.normal * RayPosNormalNudge;

        // new ray direction from normal oriented cosine weighted hemisphere sample
        rayDir = normalize(hitInfo.normal + randomUnitVector(rngState));

        if (bounceIndex == 0) {
            continue;
        }

        // update the colorMultiplier.
        throughput *= hitInfo.color;

        // Russian Roulette
        // As the throughput gets smaller, the ray is more likely to get terminated early.
        // Survivors have their value boosted to make up for fewer samples being in the average.
        let p = max(throughput.r, max(throughput.g, throughput.b));
        // \`p\` is 0 after a black hit and \`randomFloat\` may return 0: \`1 / p\` would make the throughput NaN
        if (p <= 0.0 || randomFloat(rngState) >= p) {
            break;
        }

        // Add the energy we 'lose' by randomly terminating paths
        throughput *= 1.0 / p;
    }

    // return pixel color
    return ret;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    // the target is as large as uTexSize
    fragCoords = in.position.xy / uObject.uTexSize;
    let depth = getDepth(fragCoords);

    if (isBackground(depth)) {
        return read_tColor(fragCoords);
    }

    // initialize a random number state based on the pixel and uFrameNo, the pixel with a bottom-left origin as the
    // GLSL program's \`gl_FragCoord\` (the same random sequences)
    let glFragCoord = vec2<f32>(in.position.x, uObject.uTexSize.y - in.position.y);
    var rngState = (u32(glFragCoord.x) * 1973u + u32(glFragCoord.y) * 9277u + u32(uObject.uFrameNo) * 26699u) | 1u;

    let cameraPos = vec3<f32>(0.0);
    viewPos = screenSpaceToViewSpace(vec3<f32>(fragCoords, depth), uObject.uInvProjection);
    let rayDir = normalize(viewPos);

    // raytrace for this pixel
    var color = vec3<f32>(0.0);
    for (var index = 0; index < uObject.uRendersPerFrame; index++) {
        color += colorForRay(cameraPos, rayDir, &rngState) / f32(uObject.uRendersPerFrame);
    }

    // the frame's weight in the running average (blended with it), kept for the denoiser
    return vec4<f32>(color, 1.0 / f32(uObject.uFrameNo + 1));
}
`;

/** The textures a variant of the compose pass reads through `read_<name>` (upsampled occlusion samples directly) */
function getComposeTextures(d: DefineLookup) {
  const blend = d('dBlendTransparency');
  const occlusion = d('dOcclusionEnable');
  const upsample = d('dOcclusionUpsample');
  const transparentDepth = d('dTransparentOutline') || occlusion || d('dBloomEnable');
  return [
    'tColor',
    'tDepthOpaque',
    d('dDenoise') && 'tNormal',
    blend && 'tTransparentColor',
    blend && transparentDepth && 'tDepthTransparent',
    occlusion && !upsample && 'tSsaoDepth',
    blend && occlusion && !upsample && 'tSsaoDepthTransparent',
    d('dOutlineEnable') && 'tOutlines',
    d('dBloomEnable') && 'tBloom',
  ].filter((t): t is string => !!t);
}

/** `getDepthTransparent(coords)`: depth of the blended transparency, for the compose pass with it */
const depth_transparent_wgsl = /* wgsl */ `
fn getDepthTransparent(coords: vec2<f32>) -> f32 {
    #if defined(dTransparentOutline) || defined(dOcclusionEnable) || defined(dBloomEnable)
        return unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords)).x;
    #else
        return 0.0;
    #endif
}
`;

/**
 * Denoises the traced image and composes fog, occlusion, outlines, blended transparency and bloom onto it. Depth
 * textures are always available on WebGPU, so the opaque depth is never packed. The texture reads, the denoiser and
 * the outline helper are included by the variant's defines (`WgslProgram.source`).
 */
export const illumination_compose_wgsl = (d: DefineLookup) => /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${getComposeTextures(d).map(texture_read_wgsl).join('\n')}
${d('dOutlineEnable') ? `${pack_2x4_wgsl}\n${outline_read_wgsl}` : ''}
${d('dBlendTransparency') ? depth_transparent_wgsl : ''}
${d('dDenoise') ? denoise_wgsl : ''}

fn getViewZ(depth: f32) -> f32 {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #else
        return perspectiveDepthToViewZ(depth, uObject.uNear, uObject.uFar);
    #endif
}

fn isBackground(depth: f32) -> bool {
    // the far plane (depth is reversed), packed or raw
    return depth == 0.0;
}

${compose_fragment_wgsl}`;

/** `smartDeNoise(uv)`: the traced color denoised with normal weights, for the compose pass with `dDenoise` */
const denoise_wgsl = /* wgsl */ `
const EPSILON = 1e-6;
const INV_SQRT_OF_2PI = 0.39894228040143267793994605993439; // 1.0/SQRT_OF_2PI
const INV_PI = 0.31830988618379067153776752674503;

// https://github.com/BrutPitt/glslSmartDeNoise, with normal weights; sigma 3, kSigma 2

/** exp(-i² / (2 sigma²)) for i in 0..6: the spatial weight of an offset (x, y) is w[|x|] w[|y|] */
const SpatialWeights = array<f32, 7>(1.0, 0.94595947, 0.80073740, 0.60653066, 0.41111229, 0.24935221, 0.13533528);

fn smartDeNoise(uv: vec2<f32>) -> vec4<f32> {
    let sigma = 3.0;
    let threshold = uObject.uDenoiseThreshold;

    let centrPx = read_tColor(uv);
    if (threshold == 0.0) {
        return centrPx;
    }

    let invSigmaQx2 = 0.5 / (sigma * sigma); // 1.0 / (sigma^2 * 2.0)
    let invSigmaQx2PI = INV_PI * invSigmaQx2; // 1.0 / (sqrt(PI) * sigma)

    let invThresholdSqx2 = 0.5 / (threshold * threshold); // 1.0 / (sigma^2 * 2.0)
    let invThresholdSqrt2PI = INV_SQRT_OF_2PI / threshold; // 1.0 / (sqrt(2*PI) * sigma)

    var zBuff = 0.0;
    var colorBuff = vec4<f32>(0.0);

    let normal = unpackRGBAToNormal(read_tNormal(uv));
    // the background: every sample's normal weight is 0
    if (all(normal == vec3<f32>(0.0))) {
        return centrPx;
    }

    for (var x = -6; x <= 6; x++) {
        let wx = SpatialWeights[abs(x)] * invSigmaQx2PI;
        for (var y = -6; y <= 6; y++) {
            let d = vec2<f32>(f32(x), f32(y));

            var blurFactor = wx * SpatialWeights[abs(y)];
            let uvSample = uv + d / uObject.uTexSize;

            let normalSample = unpackRGBAToNormal(read_tNormal(uvSample));
            // the normal weight to the power of 6
            let n1 = saturate(dot(normal, normalSample));
            let n2 = n1 * n1;
            blurFactor *= n2 * n2 * n2;

            let walkPx = read_tColor(uvSample);

            let dC = walkPx - centrPx;
            let deltaFactor = exp(-dot(dC, dC) * invThresholdSqx2) * invThresholdSqrt2PI * blurFactor;

            zBuff += deltaFactor;
            colorBuff += deltaFactor * walkPx;
        }
    }

    if (zBuff <= EPSILON) {
        return centrPx;
    }

    return colorBuff / zBuff;
}
`;

const compose_fragment_wgsl = /* wgsl */ `
@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, uObject.uTexSize);
    let transparentBackground = uObject.uTransparentBackground != 0u;

    #if defined(dDenoise)
        var color = smartDeNoise(coords);
    #else
        var color = read_tColor(coords);
    #endif

    let opaqueDepth = read_tDepthOpaque(coords).r;
    let backgroundViewZ = 2.0 * uObject.uFar;
    let opaqueSelfViewZ = select(getViewZ(opaqueDepth), backgroundViewZ, isBackground(opaqueDepth));
    let fogFactor = smoothstep(uObject.uFogNear, uObject.uFogFar, abs(opaqueSelfViewZ));
    let fogAlpha = 1.0 - fogFactor;

    #if defined(dBlendTransparency)
        var blendTransparency = true;
        var transparentColor = read_tTransparentColor(coords);
        let transparentDepth = getDepthTransparent(coords);
    #endif

    var alpha = 1.0;
    if (!transparentBackground) {
        // mix opaque objects with background color
        color = vec4<f32>(mix(color.rgb, uObject.uFogColor, fogFactor), color.a);
    } else {
        // pre-multiplied alpha expected for transparent background
        alpha = fogAlpha;
        color = vec4<f32>(color.rgb * fogAlpha, color.a);
    }

    #if defined(dOcclusionEnable)
        if (!isBackground(opaqueDepth)) {
            #if defined(dOcclusionUpsample)
                let occlusionFactor = getSsaoUpsampled(tSsaoDepth, coords, opaqueSelfViewZ, uObject.uNear, uObject.uFar);
            #else
                let occlusionFactor = unpackSsao(read_tSsaoDepth(coords));
            #endif
            if (!transparentBackground) {
                color = vec4<f32>(mix(mix(uObject.uOcclusionColor, uObject.uFogColor, fogFactor), color.rgb, occlusionFactor), color.a);
            } else {
                color = vec4<f32>(mix(uObject.uOcclusionColor * (1.0 - fogFactor), color.rgb, occlusionFactor), color.a);
            }
        }
        #if defined(dBlendTransparency)
            if (!isBackground(transparentDepth)) {
                let transparentFogFactor = smoothstep(uObject.uFogNear, uObject.uFogFar, abs(getViewZ(transparentDepth)));
                #if defined(dOcclusionUpsample)
                    let occlusionFactor = getSsaoUpsampled(tSsaoDepthTransparent, coords, getViewZ(transparentDepth), uObject.uNear, uObject.uFar);
                #else
                    let occlusionFactor = unpackSsao(read_tSsaoDepthTransparent(coords));
                #endif
                transparentColor = vec4<f32>(mix(uObject.uOcclusionColor * (1.0 - transparentFogFactor), transparentColor.rgb, occlusionFactor), transparentColor.a);
            }
        #endif
    #endif

    #if defined(dOutlineEnable)
        let outline = getOutline(coords);

        if (outline.hasOpaque) {
            let outlineFogFactor = smoothstep(uObject.uFogNear, uObject.uFogFar, abs(getViewZ(outline.opaqueDepth)));
            if (!transparentBackground) {
                color = vec4<f32>(mix(uObject.uOutlineColor, uObject.uFogColor, outlineFogFactor), color.a);
            } else {
                alpha = 1.0 - outlineFogFactor;
                color = vec4<f32>(mix(uObject.uOutlineColor, vec3<f32>(0.0), outlineFogFactor), color.a);
            }
        }

        #if defined(dBlendTransparency)
            if (outline.hasTransparent) {
                if (outline.hasOpaque && outline.opaqueDepth > outline.transparentDepth) {
                    blendTransparency = false;
                } else {
                    let finalOutlineAlpha = clamp(outline.alpha * 2.0, 0.0, 1.0);
                    let outlineFogFactor = smoothstep(uObject.uFogNear, uObject.uFogFar, abs(getViewZ(outline.transparentDepth)));
                    let finalAlpha = max(transparentColor.a, finalOutlineAlpha * (1.0 - outlineFogFactor));
                    transparentColor = vec4<f32>(uObject.uOutlineColor * finalAlpha, finalAlpha);
                }
            }
        #endif
    #endif

    #if defined(dBlendTransparency)
        if (blendTransparency && transparentColor.a != 0.0) {
            // blending
            color = transparentColor + color * (1.0 - transparentColor.a);
            alpha = transparentColor.a + alpha * (1.0 - transparentColor.a);
        }
    #endif

    #if defined(dBloomEnable)
        let bloom = read_tBloom(coords);
        // additive glow over geometry and background alike; alpha just covers the premultiplied rgb
        color = vec4<f32>(color.rgb + bloom.rgb, color.a);
        alpha = min(max(alpha, max(color.r, max(color.g, color.b))), 1.0);
    #endif

    return vec4<f32>(color.rgb, alpha);
}
`;
