/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Coverage in levels (seqspec:format#placement, seqspec:format#layer-types; seqspec:libraries#the-mol-adapter):
// `ranges` with `display: "union"` packs every item into
// one lane, so overlapping items draw as their union; a track with `summary: true` shows under its group's header only
// while the group is collapsed (the union of what the group holds), and its members only while it is open. And a split
// row's line hands its value to a host's track menu items (seqspec:interaction#the-default-ui).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { trackMenu } from '../src/chrome.js';
import { type CustomTrackAction, EntryKind, SeqViewerConfig } from '../src/index.js';

import type { Features } from '../src/prepare.js';
import { layer, opened, preparedRenderers } from './layout-fixtures.js';

const SECTION = 'v.main';
function doc() {
  const b = createSeqViewSpec({ title: 'Coverage' });
  const s = b.data.sequence({ ref: 's', value: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRM' });
  const cov = b.data.annotation({
    ref: 'cov',
    geometry: 'ranges',
    target: s,
    value: {
      id: ['a/A', 'a/B', 'b/A', 'b/A'],
      start: [2, 10, 20, 40],
      stop: [30, 35, 32, 60],
      entry: ['a', 'a', 'b', 'b'],
      chain: ['A', 'B', 'A', 'A'],
    },
  });
  const sec = b.view({ ref: 'v', label: 'V' }).block({ sequence: s });
  sec.track({ label: '' }).ruler();
  const g = sec.group({ ref: 'g', label: 'Coverage', collapsed: true } as never);
  g.track({ ref: 'all', label: 'All', summary: true } as never).ranges({ data: cov, display: 'union' } as never);
  for (const e of ['a', 'b'])
    g.track({ ref: `t-${e}`, label: e } as never).ranges({
      data: cov,
      filter: { field: 'entry', values: [e] },
      display: 'union',
      split: 'chain',
    } as never);
  return b.getDocument();
}

describe('coverage in levels', () => {
  it('a collapsed group shows its summary track; open, its members and not the summary', async () => {
    const t = await opened(doc(), { renderers: preparedRenderers() });
    const shown = () => t.entries(SECTION).map(([k, n]) => `${k}:${n}`);
    expect(shown()).toEqual([`${EntryKind.track}:v.ruler`, `${EntryKind.group}:g`, `${EntryKind.track}:all`]);
    t.ctx.ui.toggleGroup('g');
    await t.ctx.graph.update();
    const open = shown();
    expect(open).toContain(`${EntryKind.track}:t-a`);
    expect(open).toContain(`${EntryKind.track}:t-b`);
    expect(open).not.toContain(`${EntryKind.track}:all`);
  });

  it('union: every item in one lane (the drawn ranges are the union); lanes otherwise', async () => {
    const t = await opened(doc(), { renderers: preparedRenderers() });
    const all = layer(t.ctx, 'All').state as Features;
    expect([all.lanes, all.count]).toEqual([1, 4]);
    t.ctx.ui.toggleGroup('g');
    await t.ctx.graph.update();
    const a = layer(t.ctx, 'a').state as Features & { split?: { labels: string[]; rows: Features[] } };
    // Collapsed structure row: its chains merged in one lane; split by chain: each chain's own.
    expect(a.lanes).toBe(1);
    expect(a.split?.labels).toEqual(['A', 'B']);
    expect(a.split?.rows.map((r) => r.lanes)).toEqual([1, 1]);
  });

  it("a split row's line hands its value to a host's track menu items; the track's own line, null", async () => {
    const t = await opened(doc(), { renderers: preparedRenderers() });
    t.ctx.ui.toggleGroup('g');
    t.ctx.ui.toggleGroup('t-a');
    await t.ctx.graph.update();
    const seen: (string | null)[] = [];
    const action: CustomTrackAction = {
      id: 'host:apply',
      label: 'Apply',
      labelFor: (_c, _s, _t, split) => (split === null ? 'Apply' : `Apply ${split}`),
      checked: (_c, _s, _t, split) => split === 'B',
      run: (_c, _s, _t, split) => void seen.push(split),
    };
    t.ctx.config.set(SeqViewerConfig.UI.TrackActions, [action]);
    const at = (row: number) => trackMenu(t.ctx, SECTION, 't-a', undefined, row);
    expect(at(-1).map((i) => [i.label, i.checked])).toEqual([['Apply', false]]);
    expect(at(1).map((i) => [i.label, i.checked])).toEqual([['Apply B', true]]);
    at(0)[0].run();
    at(-1)[0].run();
    expect(seen).toEqual(['A', null]);
  });
});
