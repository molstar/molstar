/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7's heap guard: reading a sparse composed path (seqspec:model#paths) allocates nothing per step. A DMD-sized protein
// (79 exons over 2.2 Mb) coded on a lazy 200 Mb chromosome, its path read per element both ways (`Paths.at`,
// `Paths.spanAt`) and by ranges either way (`forEachFrom`: a streamed request's runs; `forEachTo`: an item's image), as
// the stream manager's pass and its draws read it. The positive control is the template's: the same step with a leak.
// Untimed, as every core guard: `heapScaling` warms by count. The cases are paths-st7-cases.ts's, shared with the quick
// twin (paths-st7.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { down, read, sums, up } from './paths-st7-cases.js';

const N = 250_000;

describe("ST7: a sparse path's reads allocate nothing", () => {
  it('holds the CDS, not the chromosome', () => {
    expect([down.length, down.heldTo.length < 12_000, up.heldSpans?.length]).toEqual([200_000_000, true, 3 * 3685]);
  });

  it('Paths.at, spanAt, forEachFrom and forEachTo over a 79-exon path through 200 Mb', () => {
    const r = heapScaling(N, read);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(Number.isFinite(sums())).toBe(true);
  });

  // CM-05: the control kept arrays on a list, whose reallocation in the shorter run can cancel the longer's growth.
  it('the same step with a leak is caught (the positive control)', () => {
    expect(controlCaught(() => heapScaling(N, leaking(read)))).toBe(true);
  });
});
