/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Complex's 200 copies, split and merged again as the Copies story does (one patch each, its frame): what a split
// allocates follows its blocks, each block's new lays bounded (a ranges lay made as one literal, its arrays made at
// their first need). By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { ComplexViews } from '@molstar/seqx-examples';
import { copiesScale } from '@molstar/seqx-examples/complex/blocks';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import type { EditingBehavior } from '@molstar/seqx-viewer/behaviors/editing';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, it } from 'vitest';
import { composeEditing } from '../../features/editing/compose.js';

const Copies = 200;

/**
 * Bytes a split and a merge (the pair, a step) may allocate a copy: measured 2026-11-24 (load 3) at about 58 KB, the
 * split's new lays and their channels; 78 KB while a ranges lay was a spread of its base (dictionary mode, about 25 KB a
 * lay made).
 */
const BytesACopy = 64_000;

const patch = (copies: 'each' | 'merge') => ({
  ops: [{ op: 'set' as const, at: ['views', 'main', 'copies'], key: 'copies', value: copies }],
});

describe('the Complex’s 200 copies', () => {
  it('a split and a merge, each its patch and its frame', async ({ task }) => {
    const base = ComplexViews.find((v) => v.id === 'blocks')!;
    const build = () => copiesScale(Copies);
    const s = composeEditing(
      { ...base, id: 'copies-scale', build },
      { ports: { assets: () => nodeFixtureAssets() }, delayMs: 0 },
    );
    const v = s.viewer;
    v.resize(1440, 900);
    await v.load(build());
    for (let i = 0; i < 20; i++) {
      await new Promise((r) => setTimeout(r, 30));
      v.draw();
    }
    const editing = v.behaviors.get('editing') as EditingBehavior;
    // a step: the patch previewed and its frame, split on even steps, merged on odd (the preview stands for the edit:
    // the frame is the same, and the host's stack is not the measure)
    const r = quickCase(
      task,
      'copies 200 split and merge',
      (i) => {
        editing.preview(patch(i % 2 ? 'merge' : 'each'));
        v.draw();
      },
      { warmSteps: 20, steps: 20, windows: 3, maxWarmSteps: 100, control: false },
    ).r;
    console.info(`200 copies, a split and a merge: ${(r.bytesPerStep / 1024).toFixed(0)} kB a step`);
    expect(v.frame!.layout.sections[0].section.blocks!.length).toBe(1);
    s.dispose();
    // a pair of steps is one split and one merge
    expect(2 * r.bytesPerStep, r.note).toBeLessThan(Copies * BytesACopy);
  }, 120_000);
});
