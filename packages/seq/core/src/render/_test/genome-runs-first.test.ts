/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Runs of small blocks are drawn the same whatever opened first (_spec/README.md "Render", Blocks): a spec of plain
// blocks in a fresh module registry, no genome ever loaded, draws its runs as one after an assembly would. Runs
// isolated (`vi.resetModules`: a registry where no genome chunk has landed).

import { describe, expect, it, vi } from 'vitest';

describe('runs whatever loaded first', () => {
  it('a row of 40 small blocks over one large one draws its run with no genome ever loaded', async () => {
    vi.resetModules();
    const { settle } = await import('./helpers.js');
    const { mixedSpec } = await import('./genome-fixtures.js');
    const tiny = Array.from({ length: 40 }, (_, k) => ({ seq: 1, range: [k * 5, k * 5 + 5] as const }));
    const s = await settle(mixedSpec([{ seq: 0, range: [0, 400], scale: 10 }, ...tiny]));
    const f = await s.frame({ size: { width: 900, height: 400 }, dpr: 2 });
    const p = f.layout.sections[0].placement as { runN: number; row: { run: Int32Array } };
    expect(p.runN).toBe(1);
    expect(p.row.run[1]).toBe(1);
  });
});
