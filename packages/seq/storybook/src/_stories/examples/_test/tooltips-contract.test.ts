/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover contract with the spec's tooltips (the viewer's _spec "Hover card", Tooltips): every example view with a
// tooltip by range over each of its sequences (a row covering the whole), swept as the hover contract is: the contract
// holds (a hit, a mark here and linked, a card), and every view whose sweep reaches a mark over a sequence shows the
// tooltip's row there; the SVG with `card` titles over every view settles and titles items alone (a bin, an image tile,
// a track, a summary none).

// @module-tag medium

import type { Node } from '@molstar/core/repr/tree';
import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome, loadTooltips } from '@molstar/seqx-core/resolve';
import { type ExampleView, ExampleViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';
import { sweepHover } from './hover-sweep.js';

let viewer: SeqViewer | undefined;
/** Each view's tooltips, cards with the tooltip's row, and cards: what the last case reads. */
const stats: number[][] = [];

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
  // the tooltips' door, landed before the first view (the fake clock waits on no import)
  await loadTooltips();
  await import('@molstar/seqx-viewer/bridge-view');
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** The spec with a tooltip over each of its sequences, whole (`Tip`). */
function withTips(spec: Node): { spec: Node; count: number } {
  const seqs = (spec.children ?? [])
    .filter((c) => c.kind === 'data')
    .flatMap((d) =>
      (d.children ?? []).filter((n) => n.kind === 'sequence' && n.id !== undefined).map((n) => `/data/${n.id}`),
    );
  const tips: Node[] = seqs.map((target, i) => ({
    kind: 'tooltip',
    id: `tip-${i}`,
    params: { target, rows: { start: [0], end: [0x7fffffff], tooltip: ['over the sequence'] } as never, label: 'Tip' },
  }));
  return {
    count: tips.length,
    spec: {
      ...spec,
      children: (spec.children ?? []).map((c) =>
        c.kind === 'views' ? { ...c, children: [...(c.children ?? []), ...tips] } : c,
      ),
    },
  };
}

/** Whether a sweep's cards are of a place over one of the spec's sequences (titled by its name and a position). */
function overSequences(cards: readonly string[], spec: Node): boolean {
  const names = (spec.children ?? [])
    .filter((c) => c.kind === 'data')
    .flatMap((d) => (d.children ?? []).filter((n) => n.kind === 'sequence' && n.id !== undefined))
    .map((n) => (typeof n.params?.name === 'string' ? n.params.name : n.id!));
  return cards.some((c) => names.some((n) => c.startsWith(`{"title":"${n} `)));
}

describe('the hover contract with tooltips over every sequence', () => {
  for (const view of ExampleViews)
    it(`${view.dataset}/${view.id}`, { tags: ['medium'] }, async () => {
      const { spec, count } = withTips(view.build());
      const tipped: ExampleView = { ...view, build: () => spec };
      viewer = await openExample(tipped, { fixtures: true, gene: '', delayMs: 0, wait: 5000, awaitOpens: true });
      const r = sweepHover(viewer, `${view.dataset}/${view.id}`);
      expect(r.wrong).toEqual([]);
      const cards = r.cards.split('\n').filter(Boolean);
      const withTip = cards.filter((c) => c.includes('"label":"Tip"')).length;
      // a view over a sequence answers with its row on a mark over it; one over none never does
      if (count === 0) expect(withTip).toBe(0);
      else if (cards.length && overSequences(cards, spec)) {
        expect(withTip).toBeGreaterThan(0);
      }
      // the SVG with card titles: settles; a title an item, never more items titled than named
      const svg = viewer.renderSvg({ idPrefix: 't-', titles: 'card' });
      await vi.advanceTimersByTimeAsync(400);
      const out = await svg;
      expect(out.svg).toContain('<svg');
      // an item's title is its aria-label and its `<title>` (twice); nothing else (a bin, a tile, a track) carries one
      expect(out.svg.match(/<g data-seq-(?!item)[a-z-]*="" role="[^"]*" aria-label="[^"]*Tip: over/g) ?? []).toEqual(
        [],
      );
      stats.push([count, withTip, cards.length]);
    });

  it('the sequences’ views show the tooltip’s row on their marks', { tags: ['medium'] }, () => {
    const over = stats.filter(([count, , cards]) => count > 0 && cards > 0);
    expect(over.length).toBeGreaterThan(0);
  });
});
