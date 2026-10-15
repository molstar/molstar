/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of hover.heap.ts (`test:perf`): the highlight path's hover changes, (a) to (f) and the source used
// twice, the guard's cases (hover-cases.ts), each with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from '../../core/_test/model/quick.js';
import { type Case, codon, containing, direct, fragments, linked, pointer, twice } from './hover-cases.js';

describe('the highlight path, quick', () => {
  const cases: { name: string; make: () => Promise<Case> }[] = [
    { name: '(a) direct block', make: direct },
    { name: '(b) codon route', make: codon },
    { name: '(c) every manager across the link', make: linked },
    { name: '(d) containing closure', make: containing },
    { name: 'a source used twice', make: twice },
    { name: '(e) hoverAt across positions', make: pointer },
    { name: '(f) fragments as axis runs', make: fragments },
  ];
  it.for(cases)('$name', async ({ name, make }, { task }) => {
    const c = await make();
    quickCase(task, `hover: ${name}`, c.step);
    expect(c.seen()).toBeGreaterThan(0);
  });
});
