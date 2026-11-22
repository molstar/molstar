/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover contract's quick representative (_spec/README.md "Render", Emphasis, The hover contract): the sweep
// (`hover-contract.test.ts`, medium lane, every view) over one view a layout: tracks, text, a ring, blocks, rows, a
// grid and a bridge.

import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { AcrossViews, type ExampleView, ExampleViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';
import { sweepHover } from './hover-sweep.js';

/** A view a layout, by `dataset/id`. */
const ByLayout: readonly (readonly [string, string])[] = [
  ['tracks', 'p53/tracks'],
  ['text', 'p53/text'],
  ['rings', 'plasmid/ring'],
  ['blocks', 'complex/blocks'],
  ['rows', 'globins/family'],
  ['grid', 'complex/pae-monomer'],
  ['bridge', 'across/across-fit'],
];

const all: readonly ExampleView[] = [...ExampleViews, ...AcrossViews];

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
  await import('@molstar/seqx-viewer/bridge-view');
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

describe('the hover contract, one view a layout', () => {
  for (const [layout, id] of ByLayout)
    it(`${layout}: ${id}`, async () => {
      const view = all.find((v) => `${v.dataset}/${v.id}` === id)!;
      viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000, awaitOpens: true });
      const r = sweepHover(viewer, id);
      expect(r.probes).toBeGreaterThan(0);
      expect(r.layouts).toContain(layout);
      expect(r.wrong).toEqual([]);
    });
});
