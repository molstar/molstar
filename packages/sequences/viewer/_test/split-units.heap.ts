/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover frames over split rows on units of their own allocate nothing after warm-up (HR1; a row's own marks): the
// pointer moving over the chains' rows (a unit's residue off its items, an item on them) with a highlight on one unit
// and an accent on the shared sequence, each move an overlay frame on the recording backend: the hit's unit residue,
// the channels' hover rebuilt per own unit, the column split into its rows' bands. The positive control runs the same
// frames allocating an object each (`leaking`) through the same warm-up and measurement, which must catch it; both are
// warmed and settled by count (`cheapScaling`, heap.ts). The case loads the viewer afresh (`vi.resetModules`)
// and is split-units-cases.ts's, shared with the quick twin (split-units.quick.ts).

import { describe, expect, it } from 'vitest';
import { cheapScaling, controlCaught, describeScaling, FRAME_N, FRAME_WARM, isFlat, leaking } from './heap.js';
import { unitsHover } from './split-units-cases.js';

describe('hover frames over split rows on their own units', () => {
  it('allocate nothing; an object allocated per frame is caught', async () => {
    const { frame, hovers } = await unitsHover();
    const r = cheapScaling(frame, FRAME_N);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(hovers()).toBeGreaterThan(FRAME_WARM);
    expect(controlCaught(() => cheapScaling(leaking(frame), FRAME_N))).toBe(true);
  }, 300_000);
});
