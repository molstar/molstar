/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
} from '../compute.js';
import type { CulledDataLayout, CulledPoseSlots } from '../culled-data.js';
import { UniformSpec } from '../renderable/schema.js';

/**
 * The poses' slots of `CulledData` expanded on the GPU: the particle list's values uploaded as
 * stored (float32 or grid positions, float32 or snorm16 rotations, float32 scales or none: 16 B a particle of a
 * pack where a slot is 32 B), a thread a slot writing what `CulledData.layout` writes on the CPU (`writePoses`): the
 * translation `t = p - s R c` and the scale, the unit quaternion, a part's similarity folded in. A part of a range of
 * elements uploads slices of the list's arrays as they are (no pass over them); one by indices gathers its elements
 */
export const CulledPosesKernel = KernelCode('culled-poses', [64]);

/** vec4s of a part's record (`CulledPosesSchema.parts`) */
export const CulledPoseRecordVec4s = 9;

export const CulledPosesSchema = {
  /** the slots of the batch */
  uThreadCount: UniformSpec('i'),
  uPartCount: UniformSpec('i'),
  /**
   * per part (`CulledPoseRecordVec4s`, the bits of f32s where not counts): its positions' first word and kind (0 float32 xyz, 1
   * grid), its rotations' first word and kind (0 none, 1 float32, 2 snorm16); its scales' first word, whether it has
   * them, the uniform scale, its slots; the grid's origin and step; the center and the similarity's scale; its rotation
   * (a quaternion); its matrix's rows; its first slot (vec4) and first thread
   */
  parts: StorageBufferSpec('vec4<u32>', 'read'),
  /** the parts' values as stored, a part's ids from 0 */
  sources: StorageBufferSpec('u32', 'read'),
  /** the culled data, its poses' slots written */
  data: StorageBufferSpec('vec4<f32>', 'read-write'),
};

export const CulledPosesKernels: readonly KernelRequest[] = [{ code: CulledPosesKernel, schema: CulledPosesSchema }];

/** The bytes a batch's values take at most (the binding limit where lower): a submit each */
const MaxBatchBytes = 1 << 26;

/** A batch of part records and their values, under the binding limit */
interface PoseBatch {
  readonly records: Uint32Array;
  readonly parts: number;
  readonly threads: number;
  readonly words: number;
  /** the values' arrays and their first word in `sources` */
  readonly writes: readonly { readonly view: ArrayBufferView; readonly word: number }[];
}

/** words an id of `p`'s positions, rotations and scales take */
function getIdWords(p: CulledPoseSlots) {
  const { positions, rotations, scales } = p.poses;
  const pw = positions instanceof Float32Array ? 3 : 2;
  const rw = !rotations ? 0 : rotations instanceof Float32Array ? 4 : 2;
  return { pw, rw, sw: scales ? 1 : 0 };
}

/** ids `from` to `from + m` of `p`'s array `a` (`w` values an element, `Float32Array`, `Uint32Array` or `Int16Array`) */
function getValues<T extends Float32Array | Uint32Array | Int16Array>(
  p: CulledPoseSlots,
  a: T,
  w: number,
  from: number,
  m: number,
): T {
  const { indices } = p.poses;
  if (!indices) {
    const e = (p.poses.first ?? 0) + from;
    const view = a.subarray(e * w, (e + m) * w) as T;
    // writes take 4-byte aligned bytes: a snorm16 slice at an odd pair copied
    return view.byteOffset % 4 === 0 ? view : (view.slice() as T);
  }
  const out = new (a.constructor as { new (n: number): T })(m * w);
  for (let i = 0; i < m; ++i) {
    const e = indices[from + i] * w;
    for (let k = 0; k < w; ++k) out[i * w + k] = a[e + k];
  }
  return out;
}

/** The batches of `slots`, a part split where it is larger than a batch */
function* getPoseBatches(slots: readonly CulledPoseSlots[], slotFirst: number, batchBytes: number) {
  const batchWords = Math.floor(batchBytes / 4);
  const maxParts = Math.max(1, Math.floor(batchBytes / (CulledPoseRecordVec4s * 16)));
  let records: number[] = [];
  let writes: { view: ArrayBufferView; word: number }[] = [];
  let parts = 0;
  let threads = 0;
  let words = 0;
  const take = (): PoseBatch => {
    const batch = { records: new Uint32Array(records), parts, threads, words, writes };
    records = [];
    writes = [];
    parts = threads = words = 0;
    return batch;
  };
  const f = new Float32Array(1);
  const u = new Uint32Array(f.buffer);
  /** a float's bits */
  const fb = (v: number) => {
    f[0] = v;
    return u[0];
  };
  for (const p of slots) {
    const { positions, rotations, scales, scale, center } = p.poses;
    const { pw, rw, sw } = getIdWords(p);
    const idWords = pw + rw + sw;
    const grid = positions instanceof Float32Array ? undefined : positions;
    const s = p.similarity;
    const m = s?.matrix;
    const [ax, ay, az, aw] = s ? s.rotation : [0, 0, 0, 1];
    for (let from = 0; from < p.ids; ) {
      if (parts === maxParts || words + idWords > batchWords) yield take();
      const n = Math.min(p.ids - from, Math.max(1, Math.floor((batchWords - words) / idWords)));
      const posWord = words;
      const rotWord = posWord + n * pw;
      const scaleWord = rotWord + n * rw;
      writes.push({ view: getValues(p, grid ? grid.data : (positions as Float32Array), pw, from, n), word: posWord });
      if (rotations) writes.push({ view: getValues(p, rotations, 4, from, n), word: rotWord });
      if (scales) writes.push({ view: getValues(p, scales, 1, from, n), word: scaleWord });
      const rows = m
        ? [m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13], m[2], m[6], m[10], m[14]]
        : [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0];
      records.push(
        posWord,
        grid ? 1 : 0,
        rotWord,
        !rotations ? 0 : rw === 4 ? 1 : 2,
        scaleWord,
        scales ? 1 : 0,
        fb(scale !== undefined && scale > 0 ? scale : 1),
        n,
        ...[grid?.origin[0] ?? 0, grid?.origin[1] ?? 0, grid?.origin[2] ?? 0, grid?.step ?? 0].map(fb),
        ...[center[0], center[1], center[2], s ? s.k : 1, ax, ay, az, aw, ...rows].map(fb),
        slotFirst + p.base + 2 * from,
        threads,
        0,
        0,
      );
      parts += 1;
      threads += n;
      words += n * idWords;
      from += n;
    }
  }
  if (parts > 0) yield take();
}

/** Writes the poses' slots of a `CulledData` layout on the GPU */
export interface CulledPoseExpander {
  /**
   * the poses' slots of `layout` (written without them, `poses: false`) into `data`, submitted batch by batch (later
   * writes land after); `dispatchItems` a kernel over a count of items (in rows past the limit)
   */
  expand: (
    layout: CulledDataLayout,
    data: StorageBuffer,
    dispatchItems: (kernel: ComputeKernel, n: number) => void,
  ) => void;
  destroy: () => void;
}

/**
 * The expander's kernel created at once, not prepared in the background: a culled object's data is complete before
 * its first draw (the culling's own kernels are not: frames draw none meanwhile), so the first submit waits for
 * its compile
 */
export function createCulledPoseExpander(compute: ComputeContext, batchBytes = MaxBatchBytes): CulledPoseExpander {
  const bytes = Math.min(batchBytes, compute.maxStorageBufferSize);
  const placeholder = compute.createBuffer(16, 'culled poses placeholder');
  const values = {
    uThreadCount: ValueCell.create(0),
    uPartCount: ValueCell.create(0),
    parts: ValueCell.create(placeholder),
    sources: ValueCell.create(placeholder),
    data: ValueCell.create(placeholder),
  };
  const kernel = compute.createKernel(CulledPosesKernel, CulledPosesSchema, values);
  /** `submitLast` false: the last batch left to the next submit (its buffers not written again) */
  const run = (
    batches: Iterable<PoseBatch>,
    data: StorageBuffer,
    dispatchItems: (kernel: ComputeKernel, n: number) => void,
  ) => {
    ValueCell.updateIfChanged(values.data, data);
    let parts: StorageBuffer | undefined;
    let sources: StorageBuffer | undefined;
    let pending = false;
    for (const batch of batches) {
      if (batch.threads === 0) continue;
      // the batch before read what this one writes over
      if (pending) compute.submit();
      // a batch's buffers reused by the next: its writes after this one's submit
      if (!parts || parts.byteLength < batch.records.byteLength) {
        parts?.destroy();
        parts = compute.createBuffer(Math.max(16, batch.records.byteLength), 'culled pose parts');
      }
      if (!sources || sources.byteLength < batch.words * 4) {
        sources?.destroy();
        sources = compute.createBuffer(Math.max(16, batch.words * 4), 'culled pose values');
      }
      parts.write(batch.records);
      for (const { view, word } of batch.writes) sources.write(view, word * 4);
      ValueCell.updateIfChanged(values.parts, parts);
      ValueCell.updateIfChanged(values.sources, sources);
      ValueCell.updateIfChanged(values.uThreadCount, batch.threads);
      ValueCell.updateIfChanged(values.uPartCount, batch.parts);
      compute.beginComputePass({ label: 'culled poses' });
      dispatchItems(kernel, batch.threads);
      compute.endComputePass();
      pending = true;
    }
    if (pending) compute.submit();
    // destroyed once the submitted work is done (the device keeps them until then)
    parts?.destroy();
    sources?.destroy();
    ValueCell.updateIfChanged(values.parts, placeholder);
    ValueCell.updateIfChanged(values.sources, placeholder);
  };
  return {
    expand: (layout, data, dispatchItems) => {
      if (layout.poseSlots.length === 0) return;
      run(getPoseBatches(layout.poseSlots, layout.slotFirst, bytes), data, dispatchItems);
    },
    destroy: () => {
      kernel.destroy();
      placeholder.destroy();
    },
  };
}
