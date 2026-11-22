/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Tooltips stories, headless on a fake clock, opened as the stories open them (a default viewer of each view): the
// spec checks clean, the view opens, and the card at the story's probe ends with the tooltip's row (by data, by range
// inline and from a source, over an assembly, Markdown with an in-view link); the SVG story's export titles its items by
// their cards; a canvas's card and the SVG's title say the same.

import { loadGenome, loadTooltips } from '@molstar/seqx-core/resolve';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { describerFor } from '@molstar/seqx-viewer/describe';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { TooltipViews } from '../tooltip-views.js';

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
  // the tooltips' door landed first: its import is real I/O, never on the fake clock
  await loadTooltips();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** The card at a story's probe: its section's track, a position on the section's system (a block's, on a genome). */
function probeCard(v: SeqViewer, id: string) {
  const view = TooltipViews.find((x) => x.id === id)!;
  const { probe } = view;
  const s = v.frame!.layout.sections[0];
  const t =
    typeof probe.track === 'number' ? s.tracks[probe.track] : s.tracks.find((x) => x.track.label === probe.track)!;
  // a section of blocks lays its layers out a strip each (the first strip's: chr1)
  const l = t.layers[0] ?? s.lays.flat().find((x) => t.track.layers.includes(x.layer))!;
  // a genome's section: the block holding the position (chr1 first)
  const block = s.section.blocks ? 0 : undefined;
  const system = s.section.blocks ? s.placement.strips.axes[0].system : s.section.axis.system;
  // `first`: the item the layer drew over the position (else its first)
  let element = probe.element === 'first' ? undefined : probe.element;
  if (probe.element === 'first') {
    const items = l.type.items?.(l.lay) ?? [];
    const it = items.find((i) => i.start <= probe.position && probe.position < i.end) ?? items[0];
    if (it) element = { owner: l.layer.data!.address, key: it.key };
  }
  return describerFor({ ...v, layout: () => v.frame?.layout }).describe({
    section: s.section.address,
    layer: l.key,
    system,
    position: probe.position,
    ...(block !== undefined ? { block } : {}),
    ...(element ? { element } : {}),
  });
}

/** Whether a card has a row of the probe's label whose value matches. */
const holds = (rows: readonly { label: string; value: string }[] | undefined, label: string, value: RegExp) =>
  !!rows?.some((r) => r.label === label && value.test(r.value));

describe('Tooltips stories, headless', () => {
  for (const view of TooltipViews)
    it(`${view.id}: checks clean, opens, and its card holds the tooltip’s row`, async () => {
      expect(checkSeqSpec(view.build()).filter((i) => i.severity === 'error')).toEqual([]);
      viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000, awaitOpens: true });
      if (view.id === 'over-assembly')
        viewer.view.setRange(viewer.frame!.layout.sections[0].section.address, 0, 10_000_000);
      await vi.advanceTimersByTimeAsync(500);
      viewer.draw();
      const card = probeCard(viewer, view.id);
      expect(holds(card?.rows, view.probe.label, view.probe.value), JSON.stringify(card?.rows)).toBe(true);
      expect(card!.text).toContain(view.probe.label);
    });

  it('markdown-link: the template’s Markdown and its in-view link on the row', async () => {
    const view = TooltipViews.find((x) => x.id === 'markdown-link')!;
    viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000 });
    const row = probeCard(viewer, view.id)!.rows.find((r) => r.label === 'Note' && r.value.startsWith('R175H'))!;
    expect(JSON.stringify(row.rich)).toContain('!seq-focus=/refs/hotspots');
    expect(JSON.stringify(row.rich)).toContain('"kind":"strong"');
  });

  it('svg-titles: the export titles each item by its card, as the canvas’s card says it', async () => {
    const view = TooltipViews.find((x) => x.id === 'svg-titles')!;
    viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000 });
    const svg = viewer.renderSvg({ idPrefix: 't-' });
    await vi.advanceTimersByTimeAsync(400);
    const r = await svg;
    expect(r.svg).toContain('Domain: DNA-binding (DBD)');
    const named = viewer.renderSvg({ idPrefix: 'n-', titles: 'name' });
    await vi.advanceTimersByTimeAsync(400);
    expect((await named).svg).not.toContain('Domain: DNA-binding (DBD)');
  });
});
