/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A legend row names its layer and lists what the layer draws (_spec/README.md "Themes", Legends; the viewer's
// Legends; FX29: a "Sequence" row listing ATP alone, its eye hiding the whole track). Every example view, settled as
// its story composes it: each row's title is its layer's name where the layer is a letters layer (its key is a scheme,
// never the track's name alone), and every color a letters layer's cells are drawn in is a swatch of its row (a
// color-less theme draws the missing color, a gap or an unknown letter, which is no category). Medium lane: it opens
// every view, headless on a fake clock; the quick representatives are `legend-layer.test.ts` in core.

// @module-tag medium

import { type LettersLay } from '@molstar/seqx-core/render';
import { ExampleViews } from '@molstar/seqx-examples';
import { legendRows } from '@molstar/seqx-viewer/legends';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

/** The views the sweep leaves out, each with why (as `empty-wells.test.ts`). */
const Skipped: ReadonlyMap<string, string> = new Map([
  [
    'vcf-zoom-in',
    'its variants come after a zoom in (the story’s point): the band says so until then, its regions refused',
  ],
]);

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** What a settled view's rows break: a letters row not titled by its layer, a color a letters layer draws no swatch lists. */
function problems(v: SeqViewer): { rows: number; lines: string[] } {
  const lines: string[] = [];
  const rows = legendRows(v, { maxItems: Number.POSITIVE_INFINITY });
  const layout = v.frame!.layout;
  const missing = layout.sections[0]?.tokens.colors['seq-missing'];
  const sequence = layout.sections[0]?.tokens.colors['seq-sequence'];
  const byLayer = new Map<string, (typeof rows)[number]>();
  for (const r of rows) for (const l of r.layers) byLayer.set(l.address, r);
  for (const s of layout.sections)
    for (const band of [...s.tracks.map((t) => t.layers), ...s.lays])
      for (const l of band ?? []) {
        if (l.type.name !== 'letters') continue;
        const row = byLayer.get(l.key);
        if (!row) continue;
        if (row.layers.length === 1 && row.title !== row.layers[0].name)
          lines.push(`${l.key}: titled "${row.title}", its layer is "${row.layers[0].name}"`);
        // never the track's name alone: a letters layer without a name of its own says its kind
        if (s.tracks.some((t) => t.track.label === row.title))
          lines.push(`${l.key}: titled by its track alone, "${row.title}"`);
        const colors = new Set(row.entries.map((e) => e.color));
        const lay = l.lay as LettersLay;
        // Two exemptions, each measured needed (removing either fails 12 views): a molecule's cell of a layer with no theme
        // of its own is drawn in the sequence's neutral (no category to key), and a row noting 'values in view' (hashed
        // values painted) keys the scale, not each color
        const neutral = lay.labels !== undefined && l.layer.theme.color === undefined;
        const fills = new Set<number>();
        for (let j = 0; j < lay.n; j++) fills.add(lay.fill[j]);
        for (const c of fills)
          if (c !== missing && !(neutral && c === sequence) && !colors.has(c) && !row.note?.includes('values in view'))
            lines.push(
              `${l.key}: a cell drawn in #${c.toString(16)} is no swatch of "${row.title}" (${row.entries.map((e) => `${e.label} #${e.color?.toString(16)}`).join(', ')}; ${lay.labels ? 'molecule' : 'letters'}, plain ${lay.plain?.toString(16)})`,
            );
      }
  // every row of every kind (ramps, sizes, categories, rows too): it has a title, and each of its layers a name its eye
  // says ("Hide <name>"; a row of several layers: "Layers of <title>"). A letters row's title is its layer's name (above);
  // another row is titled by its track's label or its scale's (themes.md, Legends: "Dot plot", "PAE"), so its eye's name
  // is the layer's own and the two may differ; what never differs is a row with no title or a layer with no name.
  for (const r of rows) {
    if (!r.title.trim()) lines.push(`a row of ${r.layers.length} layer(s) with no title (${r.layer})`);
    for (const l of r.layers) if (!l.name.trim()) lines.push(`row "${r.title}": layer ${l.address} has no name`);
  }
  return { rows: rows.length, lines };
}

describe('a legend row names its layer and lists what it draws', () => {
  for (const view of ExampleViews.filter((x) => !Skipped.has(x.id))) {
    it(`${view.dataset}/${view.id}`, { tags: ['medium'] }, async () => {
      viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000, awaitOpens: true });
      const p = problems(viewer);
      if (p.lines.length) throw new Error(p.lines.join('\n'));
    });
  }
});
