/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twins' view of the conformance kit's `heap` hook (rn2-draws.quick.ts, rn2-hot-paths.quick.ts,
// values-lines.quick.ts): the kit hands its frames over one at a time and disposes their viewer when the geometry ends,
// so a twin that cycles through a case's frames in one step needs them all while they live. The kit's order is
// deterministic: a first run counts the linear frames, a second one collects them and, at the last, hands them over.

type Hook = (frame: () => void, where: string) => void | Promise<void>;

/**
 * Runs `run` (a case through the kit, given the hook) twice and calls `measure` with the linear geometry's frames, all
 * live, at the last of them. `measure` may throw (the kit records it as a failed heap check, as a guard's hook does).
 */
export async function withLinearFrames<T>(
  run: (heap: Hook) => Promise<T>,
  measure: (frames: (() => void)[]) => void,
): Promise<T> {
  let count = 0;
  await run((_, where) => {
    if (where.startsWith('linear')) count++;
  });
  const frames: (() => void)[] = [];
  return run((frame, where) => {
    if (!where.startsWith('linear')) return;
    frames.push(frame);
    if (frames.length === count) measure(frames);
  });
}

/** One step cycling through `frames`: step `i` runs frame `i % frames.length`. */
export const cycling = (frames: readonly (() => void)[]) => (i: number) => frames[i % frames.length]();
