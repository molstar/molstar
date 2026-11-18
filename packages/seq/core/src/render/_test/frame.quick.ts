/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The allocation target (performance.md, Memory): a warm pan, zoom and hover frame allocates 0 bytes. `quickCase`
// (test-harness) reads the bytes allocated a step exactly; the inputs of a walk are made before it starts, so a step's
// bytes are the renderer's. By hand, never in CI: `npx vitest run -c packages/seq/core/vitest.perf.config.ts`. The
// bytes a step are printed (`--reporter=verbose`) and left on the report; there is no baseline until the numbers are
// accepted, so a case fails only on its positive control.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, it } from 'vitest';
import { hoverInput, hovers, prepare, Scenarios, STEPS, walk } from './scenarios.js';

// a frame of the 3M bars costs 9 ms: the warm-up is capped by steps, so a case takes under a minute
const Options = { warmSteps: 500, steps: 200, windows: 3, maxWarmSteps: 4000 };
// a pan or zoom frame allocates kilobytes to megabytes today (the numbers above): the control's 88 bytes are lost in
// that, so it only runs on the hover, which is within a control's reach
const NoControl = { ...Options, control: false, why: 'a frame allocating kilobytes hides the control' } as const;

describe('warm frames allocate nothing', () => {
  for (const s of Scenarios) {
    for (const kind of ['pan', 'zoom'] as const)
      it(`${s.name}: ${kind}`, async ({ task }) => {
        const settled = await prepare(s);
        const inputs = walk(s, kind);
        const { r } = quickCase(
          task,
          `${s.name} ${kind}`,
          (i) => void settled.renderer.frame(inputs[i % STEPS]),
          NoControl,
        );
        console.info(`${s.name} ${kind}: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
      });

    it(`${s.name}: hover`, async ({ task }) => {
      const settled = await prepare(s);
      const input = hoverInput(s);
      const { emphasis } = hovers(settled, input);
      const inputs = emphasis.map((e) => ({ ...input, emphasis: e }));
      const { r } = quickCase(task, `${s.name} hover`, (i) => void settled.renderer.frame(inputs[i % STEPS]), Options);
      console.info(`${s.name} hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    });
  }
});
