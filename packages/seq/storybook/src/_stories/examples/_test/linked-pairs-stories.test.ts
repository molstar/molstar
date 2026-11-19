/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `linkedPairs` behavior on the stories (_spec/README.md "Behaviors", Linked pairs), composed as the stories
// compose them: on Complex › Contacts the cyclin box hovered lights exactly the contacts with an end in it (both copies,
// merged), a selected one lights them on the channel and the SVG export of the selection draws them; on the tRNA's dot
// plot a modified residue hovered lights its pairs on the arcs and on the dots alike.

import { rangesLocation } from '@molstar/seqx-core/query';
import { recordDrawLists, type SectionLayout } from '@molstar/seqx-core/render';
import { ComplexViews, TrnaViews } from '@molstar/seqx-examples';
import { constructOffset, CyclinFeatures, Data } from '@molstar/seqx-examples/complex/data';
import { LinkedPairsChannel } from '@molstar/seqx-viewer/behaviors';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const Owner = {};
const lit = (v: SeqViewer) => (v.highlight.get(LinkedPairsChannel)?.elements.keys ?? []).map((k) => k.key);

/** The highlight pass's marks in an ink: arcs (curves) and dots (rects), stroked. */
function marks(lists: Parameters<typeof recordDrawLists>[0], ink: string) {
  const lines = recordDrawLists(lists, ['highlight']).split('\n');
  return {
    curves: lines.filter((l) => l.startsWith('path') && l.includes(' C ') && l.includes(`stroke ${ink}`)).length,
    rects: lines.filter((l) => l.startsWith('rect') && l.includes(`stroke ${ink}`)).length,
  };
}

/** The channel's ink as the lists record it. */
const inkOf = (v: SeqViewer) => {
  const c = v.emphasis.highlights?.find((h) => h.elements)?.color ?? 0;
  return `#${(c >>> 0).toString(16).padStart(8, '0').slice(0, 6)}`;
};

describe('linkedPairs on the stories', () => {
  it('Complex › Contacts: the cyclin box lights exactly the contacts with an end in it, both copies; a selection’s export draws them', async () => {
    viewer = await openExample(ComplexViews.find((x) => x.id === 'contacts')!, { fixtures: true });
    const v = viewer;
    const shift = constructOffset();
    const [from, to] = [CyclinFeatures.start[0] - shift, CyclinFeatures.end[0] - shift];
    const expected = Data.contacts
      .filter((c) => c[2] >= from && c[2] < to)
      .map((c) => `${c[0]}:${c[1] + 1}-${c[2] + 1}`)
      .sort();
    expect(new Set(Data.contacts.filter((c) => c[2] >= from && c[2] < to).map((c) => c[0]))).toEqual(new Set([1, 2]));
    v.hover.set(
      {
        layer: '/views/main/complex/features/cyclin',
        system: '/data/cyclin',
        element: { owner: '/data/cyclin-features', key: 'box-n' },
      },
      Owner,
    );
    expect([...lit(v)].sort()).toEqual(expected);
    v.draw();
    expect(marks(v.frame!.lists, inkOf(v)).curves).toBe(expected.length);
    v.hover.clear(Owner);
    expect(lit(v)).toEqual([]);
    // selected: on the channel, never in the selection; the export of the selection draws the lit pairs
    v.addresses.select('/data/cyclin-features/(item :id box-n)');
    expect([...lit(v)].sort()).toEqual(expected);
    v.draw();
    const ink = inkOf(v);
    const svg = v.renderSvg({ idPrefix: 'lp-' });
    await vi.advanceTimersByTimeAsync(400);
    const r = await svg;
    expect(marks(r.lists, ink).curves).toBe(expected.length);
  });

  it('tRNA dot plot: a modified residue hovered lights its pairs on the arcs and the dots alike', async () => {
    viewer = await openExample(TrnaViews.find((x) => x.id === 'dot-plot')!, { fixtures: true });
    const v = viewer;
    const s: SectionLayout = v.frame!.layout.sections[0];
    const modified = s.tracks.find((t) => t.track.address.endsWith('/modified'))!.layers[0];
    const items = modified.type.items!(modified.lay);
    const pairs = v.spec.model.value.annotations.get('/data/pairs')!.pairs!;
    const paired = (p: number) => {
      const out: (string | number)[] = [];
      for (let k = 0; k < pairs.n; k++)
        if ((p >= pairs.aStart[k] && p < pairs.aEnd[k]) || (p >= pairs.bStart[k] && p < pairs.bEnd[k]))
          out.push(pairs.keys[k]);
      return out;
    };
    // a modified residue at a paired base
    const item = items.find((d) => paired(d.start).length)!;
    expect(item).toBeDefined();
    // as the ranges layer's hit says it: the item and its extent (read from the file: no address names it before)
    v.hover.set(
      {
        layer: modified.key,
        system: '/data/trna',
        element: { owner: '/data/modified', key: item.key },
        extent: rangesLocation('/data/trna', [item.start, item.end]),
      },
      Owner,
    );
    expect(lit(v)).toEqual(paired(item.start));
    v.draw();
    const m = marks(v.frame!.lists, inkOf(v));
    expect(m.curves).toBe(paired(item.start).length);
    // the dots mirrored: each pair twice on the plot
    expect(m.rects).toBeGreaterThanOrEqual(paired(item.start).length);
  });
});
