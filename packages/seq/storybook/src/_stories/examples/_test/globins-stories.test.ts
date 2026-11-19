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
  // `scale`'s 10,000 rows (its tiles, its export) are slow by nature: the medium lane; `lazy`'s are the quick representative
  for (const id of ['family', 'folded', 'figure', 'scale', 'lazy'])
    it(`${id}: loads, draws its rows and exports them with no issue`, {
      tags: id === 'scale' ? ['medium'] : [],
    }, async () => {
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

  it('family and folded: labels read whole beside the tree (its own width), and the wedges stop before them', async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'family')!, { gene: '' });
    const fam = rowsLay(viewer)!;
    expect(fam.treeX1 - fam.treeX0).toBe(120);
    expect(fam.labelText).toContain('HBB_HUMAN/1-147');
    viewer.retire();
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'folded')!, { gene: '' });
    const folded = rowsLay(viewer)!;
    const summaries = folded.labelText.filter((t) => /\+\d+$/.test(t));
    expect(summaries.length).toBe(2);
    for (let k = 0; k + 4 < folded.wedges.length; k += 5) expect(folded.wedges[k + 2]).toBeLessThan(folded.labelX);
  });

  it('scale at the fit: the rows are image tiles with their pixels drawn, and the export carries them as images', {
    tags: ['medium'],
  }, async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'scale')!, { gene: '' });
    // every row in the band: 10,000 rows in 240 px
    const key = rowsLay(viewer)!.key;
    viewer.view.setBand(key.slice(0, key.lastIndexOf('/')), { bandTop: 0, rowScale: 240 / 10_000 });
    await vi.advanceTimersByTimeAsync(50);
    viewer.draw();
    const lay = rowsLay(viewer)!;
    expect(lay.display.count * lay.rowScale).toBeLessThanOrEqual(lay.box.height + 1e-6);
    expect(lay.images.length).toBeGreaterThan(0);
    // every tile has opaque pixels (the base level holds no shares: a sample there is whole)
    for (const img of lay.images) {
      let alpha = 0;
      for (let i = 3; i < img.data.length; i += 4) alpha += img.data[i];
      expect(alpha).toBeGreaterThan(0);
    }
    const p = viewer.renderSvg({ idPrefix: 's-' });
    await vi.advanceTimersByTimeAsync(400);
    const svg = (await p).svg;
    expect(svg.match(/<image [^>]*href="data:image\/png;base64,/g)?.length).toBe(lay.images.length);
  });

  for (const id of ['family', 'scale'])
    it(`${id}: the rows rail at the band's right edge, live (its rows scroll), drawn; a drag on it scrolls them`, {
      tags: id === 'scale' ? ['medium'] : [],
    }, async () => {
      viewer = await openExample(GlobinsViews.find((v) => v.id === id)!, { gene: '' });
      const v = viewer;
      const s = v.frame!.layout.sections.find((x) => x.rowsRails?.length)!;
      expect(s, id).toBeDefined();
      const r = s.rowsRails![0];
      const lay = rowsLay(v)!;
      expect(r.vertical).toBe(true);
      expect(r.active).toBe(true);
      // beside the rows' last numbers, within the canvas
      expect(r.top).toBeGreaterThanOrEqual(lay.box.x + lay.box.width);
      expect(r.bottom).toBeLessThanOrEqual(v.width);
      // drawn in the window pass: a rect within the rail's column
      const rects = v.frame!.lists.passes.window.rects;
      let inside = 0;
      for (let i = 0; i < rects.n; i++) if (rects.x[i] >= r.top && rects.x[i] + rects.w[i] <= r.bottom) inside++;
      expect(inside).toBeGreaterThan(1);
      // a drag on its window scrolls the rows
      const c = v.context.renderConfig.navigator.window;
      const y = r.x + Math.max(c.minPx, ((r.end - r.start) / r.total) * r.width) / 2;
      const x = (r.top + r.bottom) / 2;
      expect(v.input.down(x, y)).toBe(true);
      v.input.move(x, y + r.width / 4);
      v.input.up();
      expect(v.view.band(lay.key.slice(0, lay.key.lastIndexOf('/')))?.bandTop).toBeGreaterThan(r.total / 8);
    });
});
