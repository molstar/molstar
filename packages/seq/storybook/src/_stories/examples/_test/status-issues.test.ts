/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The example stories' status line (example-root.ts): the export's issues and the viewer's model's as it says them, so
// `seq.resolve.length.unknown` of a ring read from a file is held back by the viewer while the file opens
// (viewer/src/_test/model-issues.test.ts), and gone once it opened.

import { PlasmidViews } from '@molstar/seqx-examples';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { statusIssues } from '../status-issues.js';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

describe('the status line’s issues', () => {
  it('the crossing ring read from its file, opened: no wait for its length said (the viewer held it back meanwhile)', async () => {
    viewer = await openExample(PlasmidViews.find((v) => v.id === 'crossing')!, { gene: '' });
    const v = viewer;
    const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);
    expect(codes(v.spec.model.value.issues)).toEqual([]);
    expect(codes(statusIssues([], v.modelIssues.value))).toEqual([]);
  });

  it('the export’s first, then what only the model says, each once', () => {
    const a = { code: 'a', severity: 'warning', message: 'a', address: '/x' } as const;
    const b = { code: 'b', severity: 'error', message: 'b', address: '/y' } as const;
    expect(statusIssues([a], [a, b])).toEqual([a, b]);
    expect(statusIssues([a], [])).toEqual([a]);
  });
});
