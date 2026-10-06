/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

// Spherical harmonic polynomials match shader.frag.ts and spherical-functions.ts.
export const orbitalShader = /* wgsl */ `
struct Parameters { dimensions: vec4u, origin: vec4f, delta: vec4f, batch: vec4u, counts: vec4u };
@group(0) @binding(0) var<uniform> params: Parameters;
@group(0) @binding(1) var<storage, read> centers: array<vec4f>;
@group(0) @binding(2) var<storage, read> infos: array<vec4f>;
@group(0) @binding(3) var<storage, read> coefficients: array<f32>;
@group(0) @binding(4) var<storage, read> alphas: array<f32>;
@group(0) @binding(5) var<storage, read> occupancies: array<f32>;
@group(0) @binding(6) var<storage, read_write> values: array<f32>;
fn L1(p: vec3f, a0: f32, a1: f32, a2: f32) -> f32 {
    return a0 * p.z + a1 * p.x + a2 * p.y;
}

fn L2(p: vec3f, a0: f32, a1: f32, a2: f32, a3: f32, a4: f32) -> f32 {
    let x = p.x; let y = p.y; let z = p.z;
    let xx = x * x; let yy = y * y; let zz = z * z;
    return (
        a0 * (-0.5 * xx - 0.5 * yy + zz) +
        a1 * (1.7320508075688772 * x * z) +
        a2 * (1.7320508075688772 * y * z) +
        a3 * (0.8660254037844386 * xx - 0.8660254037844386 * yy) +
        a4 * (1.7320508075688772 * x * y)
    );
}

fn L3(p: vec3f, a0: f32, a1: f32, a2: f32, a3: f32, a4: f32, a5: f32, a6: f32) -> f32 {
    let x = p.x; let y = p.y; let z = p.z;
    let xx = x * x; let yy = y * y; let zz = z * z;
    let xxx = xx * x; let yyy = yy * y; let zzz = zz * z;
    return (
        a0 * (-1.5 * xx * z - 1.5 * yy * z + zzz) +
        a1 * (-0.6123724356957945 * xxx - 0.6123724356957945 * x * yy + 2.449489742783178 * x * zz) +
        a2 * (-0.6123724356957945 * xx * y - 0.6123724356957945 * yyy + 2.449489742783178 * y * zz) +
        a3 * (1.9364916731037085 * xx * z - 1.9364916731037085 * yy * z) +
        a4 * (3.872983346207417 * x * y * z) +
        a5 * (0.7905694150420949 * xxx - 2.3717082451262845 * x * yy) +
        a6 * (2.3717082451262845 * xx * y - 0.7905694150420949 * yyy)
    );
}

fn L4(p: vec3f, a0: f32, a1: f32, a2: f32, a3: f32, a4: f32, a5: f32, a6: f32, a7: f32, a8: f32) -> f32 {
    let x = p.x; let y = p.y; let z = p.z;
    let xx = x * x; let yy = y * y; let zz = z * z;
    let xxx = xx * x; let yyy = yy * y; let zzz = zz * z;
    let xxxx = xxx * x; let yyyy = yyy * y; let zzzz = zzz * z;
    return (
        a0 * (0.375 * xxxx + 0.75 * xx * yy + 0.375 * yyyy - 3.0 * xx * zz - 3.0 * yy * zz + zzzz) +
        a1 * (-2.3717082451262845 * xxx * z - 2.3717082451262845 * x * yy * z + 3.1622776601683795 * x * zzz) +
        a2 * (-2.3717082451262845 * xx * y * z - 2.3717082451262845 * yyy * z + 3.1622776601683795 * y * zzz) +
        a3 * (-0.5590169943749475 * xxxx + 0.5590169943749475 * yyyy + 3.3541019662496847 * xx * zz - 3.3541019662496847 * yy * zz) +
        a4 * (-1.118033988749895 * xxx * y - 1.118033988749895 * x * yyy + 6.708203932499369 * x * y * zz) +
        a5 * (2.091650066335189 * xxx * z + -6.274950199005566 * x * yy * z) +
        a6 * (6.274950199005566 * xx * y * z + -2.091650066335189 * yyy * z) +
        a7 * (0.739509972887452 * xxxx - 4.437059837324712 * xx * yy + 0.739509972887452 * yyyy) +
        a8 * (2.958039891549808 * xxx * y + -2.958039891549808 * x * yyy)
    );
}

fn spherical(l: u32, p: vec3f, a: u32) -> f32 {
    if (l == 0u) { return alphas[a]; }
    if (l == 1u) { return L1(p, alphas[a + 0u], alphas[a + 1u], alphas[a + 2u]); }
    if (l == 2u) { return L2(p, alphas[a + 0u], alphas[a + 1u], alphas[a + 2u], alphas[a + 3u], alphas[a + 4u]); }
    if (l == 3u) { return L3(p, alphas[a + 0u], alphas[a + 1u], alphas[a + 2u], alphas[a + 3u], alphas[a + 4u], alphas[a + 5u], alphas[a + 6u]); }
    if (l == 4u) { return L4(p, alphas[a + 0u], alphas[a + 1u], alphas[a + 2u], alphas[a + 3u], alphas[a + 4u], alphas[a + 5u], alphas[a + 6u], alphas[a + 7u], alphas[a + 8u]); }
    return 0.0;
}
@compute @workgroup_size(64)
fn compute(@builtin(global_invocation_id) invocation: vec3u) {
    let local = invocation.x + invocation.y * params.batch.z * 64u;
    if (local >= params.batch.y) { return; }
    let index = params.batch.x + local;
    let dim = params.dimensions.xyz;
    let coordinate = vec3u(index / (dim.y * dim.z), (index / dim.z) % dim.y, index % dim.z);
    let point = params.origin.xyz + vec3f(coordinate) * params.delta.xyz;
    var total = 0.0;
    for (var orbital = 0u; orbital < params.counts.y; orbital++) {
        if (params.counts.z != 0u && occupancies[orbital] == 0.0) { continue; }
        var value = 0.0;
        for (var i = 0u; i < params.counts.x; i++) {
            let p = point - centers[i].xyz;
            let r2 = dot(p, p);
            if (r2 > centers[i].w) { continue; }
            let info = infos[i];
            var radial = 0.0;
            for (var c = u32(info.z); c < u32(info.w); c++) {
                radial += coefficients[c * 3u] * exp(-coefficients[c * 3u + 1u] * r2);
            }
            value += radial * spherical(u32(info.x), p, u32(info.y) + orbital * params.counts.w);
        }
        if (params.counts.z == 0u) { total = value; }
        else { total += occupancies[orbital] * value * value; }
    }
    values[local] = total;
}
`;
