/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The instances' data of a GPU culler (`merged-culler.ts`, `sCulledData`), read by index: vec4s in one buffer (one
 * binding: the vertex stage's storage buffers are few) in blocks, for `n` instances of `segmentCount` segments: the
 * transforms' affine rows (three vec4s each, from 0: as the sphere raster's transforms, `affine-rows.ts`), `aInstance`
 * (four a vec4), each instance's segment index (four a vec4), each segment's groups (a vec4: `aSegment.xyz`, the same
 * for its instances); 56 bytes an instance (a buffer of `maxBufferSize`, 2 GiB, holds 38 M)
 */
export const CulledData = {
  /** vec4s of `n` instances' values: four a vec4 */
  quarter: (n: number) => Math.ceil(n / 4),
  vec4s: (n: number, segmentCount: number) => 3 * n + 2 * Math.ceil(n / 4) + segmentCount,
  instanceBase: (n: number) => 3 * n,
  segmentIndexBase: (n: number) => 3 * n + Math.ceil(n / 4),
  segmentGroupsBase: (n: number) => 3 * n + 2 * Math.ceil(n / 4),
  /**
   * WGSL readers of instance `i` (`u32`) of `data` (an `array<vec4<f32>>`) holding `n` instances (an `u32`
   * expression): `culledTransform0..2`, `culledInstance`, `culledSegment` (as `aSegment`: the segment's groups, its
   * index)
   */
  wgsl: (data: string, n: string) => /* wgsl */ `
fn culledTransform0(i: u32) -> vec4<f32> { return ${data}[i * 3u]; }
fn culledTransform1(i: u32) -> vec4<f32> { return ${data}[i * 3u + 1u]; }
fn culledTransform2(i: u32) -> vec4<f32> { return ${data}[i * 3u + 2u]; }
fn culledInstance(i: u32) -> f32 { return ${data}[3u * (${n}) + i / 4u][i % 4u]; }
fn culledSegment(i: u32) -> vec4<f32> {
    let n = ${n};
    let q = (n + 3u) / 4u;
    let s = ${data}[3u * n + q + i / 4u][i % 4u];
    return vec4<f32>(${data}[3u * n + 2u * q + u32(s)].xyz, s);
}
`,
};
