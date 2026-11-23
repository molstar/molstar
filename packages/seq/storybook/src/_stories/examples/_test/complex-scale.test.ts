/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews } from '@molstar/seqx-examples';
import { copiesScale } from '@molstar/seqx-examples/complex/blocks';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, describe, expect, vi } from 'vitest';
import { openExample } from './draw-checks.js';

// The Complex at its scale: 200 copies of CDK2 split into a block each, a pan
// frame (the move, its layout and its draw) under 8 ms; and the PAE over the complex (730 × 730, four strips) a pan
// frame under 8 ms once its tiles are built. Quick: each runs once unbounded; heavy: bounded, the best of three.

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const base = ComplexViews.find((v) => v.id === 'blocks')!;

describe('the Complex at scale', () => {
  itWithin('200 copies of CDK2, a block each: a pan frame under 8 ms', async (within) => {
    const view = { ...base, id: 'copies-scale', build: () => copiesScale(200, 'each') };
    viewer = await openExample(view, { fixtures: true, realClock: true });
    const v = viewer;
    expect(v.frame!.layout.sections[0].section.blocks!.length).toBe(200);
    v.view.setRange('/views/main/copies', 0, 2000);
    v.draw();
    // a warm frame first: the move's lays then kept
    v.view.pan('/views/main/copies', 40);
    v.draw();
    within(8, () => {
      v.view.pan('/views/main/copies', 40);
      v.draw();
    });
    expect(v.frame!.layout.issues).toEqual([]);
  });

  itWithin('the PAE over the complex: a pan frame under 8 ms once its tiles are built', async (within) => {
    viewer = await openExample(ComplexViews.find((v) => v.id === 'complex-grid')!, { fixtures: true, realClock: true });
    const v = viewer;
    v.view.setRange('/views/main/chains', 100, 500);
    for (let i = 0; i < 4; i++) {
      v.view.pan('/views/main/chains', 2);
      v.draw();
    }
    within(8, () => {
      v.view.pan('/views/main/chains', 2);
      v.draw();
    });
    expect(v.frame!.layout.sections[1].placement.strips.n).toBe(4);
  });
});
