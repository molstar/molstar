/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Globins stories, headless on a fake clock, composed as the stories compose them (a viewer of the view and its SVG
// export): each view loads and draws its rows (the figure its text), and its export has no issue.

import type { RowsLay } from '@molstar/seqx-core/render';
import { GlobinsViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

function rowsLay(v: SeqViewer): RowsLay | undefined {
  for (const s of v.frame!.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  return undefined;
}

describe('Globins stories, headless', () => {
  for (const id of ['family', 'folded', 'figure', 'scale', 'lazy'])
    it(`${id}: loads, draws its rows and exports them with no issue`, async () => {
      viewer = await openExample(GlobinsViews.find((v) => v.id === id)!, { gene: '' });
      const v = viewer;
      const lay = rowsLay(v);
      expect(lay, id).toBeDefined();
      expect(lay!.display.count).toBeGreaterThan(id === 'folded' ? 2 : 10);
      const p = v.renderSvg({ idPrefix: 'g-' });
      await vi.advanceTimersByTimeAsync(400);
      const r = await p;
      expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
      expect(r.svg.length).toBeGreaterThan(1000);
    });
});
