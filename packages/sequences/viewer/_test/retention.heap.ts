/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Retention (IN1; seqspec:interaction#channels: nothing of the projector is module state): a disposed viewer frees its
// SeqSet and resolved view once the host drops it, after its projectors (the viewports' edge-rate one, a surface's)
// have projected, closed and focused through them. The positive control: the same, with the host keeping the set,
// collected and counted alike, is caught. The check is once-only and async, so it has no step: its projection core is
// retention-cases.ts's `projecting`, measured by the quick twin (retention.quick.ts).

import { expect, it } from 'vitest';
import { controlCaught } from './heap.js';
import { collected, used } from './retention-cases.js';

it('a disposed viewer frees its SeqSet and resolved view', async () => {
  expect(await collected(await used())).toBe(0);
});

it('the same with the host keeping the set is caught (positive control)', async () => {
  const kept: unknown[] = [];
  const live = await collected(await used(kept));
  expect(controlCaught(() => live, 1)).toBe(true);
  kept.length = 0;
});
