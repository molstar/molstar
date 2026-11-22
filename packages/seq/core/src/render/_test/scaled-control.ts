/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A positive control for an allocation case that allocates hundreds of kilobytes a step by design (a layout): the
// harness's control (`quickCase`, `leaking`) reads the case and the case with an 88-byte object apart, each a median of
// its windows, and over a step of 100 to 500 KB (300 to 1,200 collections a window) the windows wander by 20 to 250
// bytes a step: the 88 bytes read from -13 to 105 over the case, against the control's 66 (the 3,000-contig drag, the
// bridge's relayout). This control adds an array of 4 KB a step (garbage at the next step), well above that noise, and
// checks the measurement sees three quarters of it.

import { type Alloc, measureBytesPerStep } from '@molstar/test-harness/quick';
import { expect } from 'vitest';

/** The control's array: 510 small integers, 4 KB with its header. */
const Slots = 510;
export const ScaledControlBytes = 4096;
const garbage: unknown[] = [null];

/** The why of `quickCase`'s `control: false` for a case checked by `scaledControl`. */
export const ScaledControlWhy =
  'the 88-byte control is within a layout’s window noise: checked by scaledControl (4 KB)';

/** Measures `step` with 4 KB more a step, right after its case `r` on the same state, and checks the difference. */
export function scaledControl(id: string, step: (i: number) => void, r: Alloc, steps = 1_000): Alloc {
  const c = measureBytesPerStep(
    (i) => {
      step(i);
      garbage[0] = new Array(Slots).fill(i);
    },
    { warmSteps: steps, steps, windows: 3, maxWarmSteps: steps },
  );
  const over = c.bytesPerStep - r.bytesPerStep;
  console.info(`${id} +4 KB control: ${over.toFixed(0)} B over the case`);
  expect(over, `${id} +4 KB: ${c.note}`).toBeGreaterThan(0.75 * ScaledControlBytes);
  return c;
}
