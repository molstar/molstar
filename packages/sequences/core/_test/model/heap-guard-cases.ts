/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The synthetic steps of the measurement helpers' self-tests, for heap-guard.heap.ts (`heapGrowth`, `heapScaling`, the
// guard template's control) and its quick twin (heap-guard.quick.ts, `quickCase`): a step that allocates nothing, one
// that allocates 208 bytes (two `subarray` views, BK1's case), one that allocates once, and garbage by the object.

export const sink: { v: unknown } = { v: null };

/** A run of a million additions that allocates nothing (`heapGrowth`'s run). */
export const clean = () => {
  let s = 0;
  for (let i = 0; i < 1_000_000; i++) s += i & 7;
  sink.v = s;
};

/** `count` short-lived objects, none kept: 40 bytes each in Node (three header words and two fields, uncompressed). */
export const garbage = (count: number) => () => {
  for (let i = 0; i < count; i++) sink.v = { a: i, b: i };
};

/** A step that allocates nothing: a small integer a step, kept on an object's field. */
export const flatStep = () => {
  let s = 0;
  return (i: number) => {
    s = (s + (i & 7)) | 0;
    sink.v = s;
  };
};

const f32 = new Float32Array(256);
/**
 * A step of two `subarray` views, 104 bytes each (BK1's `PathBatch.append`, fixed there): 812 KB over 4,000 steps, under
 * `LIMIT` for one run.
 */
export const viewsStep = (i: number) => {
  sink.v = f32.subarray(i & 7);
  sink.v = f32.subarray(0, 8 + (i & 7));
};

/** A step that allocates 512 KB at its first index and nothing after: once per run, not per step. */
export const onceStep = (i: number) => {
  if (i === 0) sink.v = new Array(1 << 16).fill(i); // On the heap, as a typed array's store is not.
};
