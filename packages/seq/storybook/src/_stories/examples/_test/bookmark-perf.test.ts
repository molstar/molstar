/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { itWithin } from '@molstar/test-harness/timing';
import { afterEach, describe, expect, vi } from 'vitest';
import { openExample } from './draw-checks.js';

// The performance representative of `selection.bookmark` (_spec/README.md "Commands"): the command on the p53 gene
// view (Gene by region, its transcripts and domains loaded) under 4 ms, the selection's address printed and the ref
// added as one patch. The body runs again for each try of the bound (its setup too: each bookmarks a fresh viewer).

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

describe('selection.bookmark on the gene view', () => {
  itWithin('bookmarks a selection of the gene in under 4 ms', async (within) => {
    // the clock of the bound is the real one: its timers alone are faked
    viewer = await openExample(P53Views.find((v) => v.id === 'gene-by-region')!, { gene: '', realClock: true });
    const v = viewer;
    v.addresses.select('/data/tp53-region/(range :from 2100 :to 2400)');
    expect(v.addresses.selectionAddress()).toBeDefined();
    const before = v.addresses.refs.value.length;
    const message = within(4, () => v.execute('selection.bookmark', {}));
    expect(message).toBe('Bookmarked');
    expect(v.addresses.refs.value.length).toBe(before + 1);
  });
});
