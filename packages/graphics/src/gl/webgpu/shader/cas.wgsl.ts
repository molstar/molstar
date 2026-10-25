/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';

/**
 * Robust contrast adaptive sharpening (RCAS of AMD FidelityFX FSR 1, MIT License, https://gpuopen.com/fsr), as
 * adapted for WebGL in https://www.shadertoy.com/view/stXSWB
 */
export const cas_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tColor')}

const FsrRcasLimit = 0.25 - (1.0 / 16.0);

/** pixel \`p\` (top-left origin), at its center */
fn fsrRcasLoad(p: vec2<f32>) -> vec4<f32> {
    return read_tColor(p * uObject.uTexSizeInv);
}

fn fsrRcas(sp: vec2<f32>, con: f32) -> vec3<f32> {
    // minimal 3x3 pixel neighborhood
    //    b
    //  d e f
    //    h
    let b = fsrRcasLoad(sp + vec2<f32>(0.0, -1.0)).rgb;
    let d = fsrRcasLoad(sp + vec2<f32>(-1.0, 0.0)).rgb;
    let e = fsrRcasLoad(sp).rgb;
    let f = fsrRcasLoad(sp + vec2<f32>(1.0, 0.0)).rgb;
    let h = fsrRcasLoad(sp + vec2<f32>(0.0, 1.0)).rgb;

    // luma times 2
    let bL = b.g + 0.5 * (b.b + b.r);
    let dL = d.g + 0.5 * (d.b + d.r);
    let eL = e.g + 0.5 * (e.b + e.r);
    let fL = f.g + 0.5 * (f.b + f.r);
    let hL = h.g + 0.5 * (h.b + h.r);

    // noise detection
    #if defined(dDenoise)
        var nz = 0.25 * (bL + dL + fL + hL) - eL;
        // a flat area has no range (0 / 0)
        nz = clamp(abs(nz) / max(max(max(bL, dL), max(eL, max(fL, hL))) - min(min(bL, dL), min(eL, min(fL, hL))), 1e-5), 0.0, 1.0);
        nz = 1.0 - 0.5 * nz;
    #endif

    // min and max of ring
    let mn4 = min(min(b, d), min(f, h));
    let mx4 = max(max(b, d), max(f, h));

    // immediate constants for peak range
    let peakC = vec2<f32>(1.0, -4.0);

    // limiters
    // black (mx4 = 0) and white (mn4 = 1) channels are 0 / 0
    let hitMin = mn4 / max(4.0 * mx4, vec3<f32>(1e-5));
    let hitMax = (peakC.x - mx4) / min(4.0 * mn4 + peakC.y, vec3<f32>(-1e-5));
    let lobeRGB = max(-hitMin, hitMax);
    var lobe = max(-FsrRcasLimit, min(max(lobeRGB.r, max(lobeRGB.g, lobeRGB.b)), 0.0)) * con;

    // noise removal
    #if defined(dDenoise)
        lobe *= nz;
    #endif

    return (lobe * (b + d + h + f) + e) / (4.0 * lobe + 1.0);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    // sharpness in stops (0 is maximum) to linear
    let con = exp2(-uObject.uSharpness);
    let fragCoord = in.position.xy;
    let col = fsrRcas(fragCoord, con);
    return vec4<f32>(col, fsrRcasLoad(fragCoord).a);
}
`;
