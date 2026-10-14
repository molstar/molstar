/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Split rows (LY1; seqspec:format#layer-types), the layout's and the hit's half: a split layer expanded is a label
// entry and a row per value in the layer's order, each row `splitHeight` tall, hit per row (an item only on its value's
// row); the toggle is `ui.expanded`, reset by a patch or a snapshot that changes the track's `expanded`. The stub
// renderers prepare each value's rows (`splitRows`), as `Features.split`; LY2's prepare gives each row its own features
//. The track menu and the action are UX's and CM1's.

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { case7, case13 } from '../../core/_test/format/cases.js';
import { EntryKind, Hit, HoverKind } from '../src/index.js';
import type { SplitFeatures } from '../src/layout.js';

import type { Features } from '../src/prepare.js';
import { opened, preparedRenderers } from './layout-fixtures.js';

const NV = 'features.natural-variants';
const ROWS = [...'GAVLISTCMDNEQRKHFYWP', '-', '*'];

async function variants(doc: SeqViewSpecDoc = case13(), prepared = false) {
  const renderers = prepared ? preparedRenderers() : undefined;
  const t = await opened(doc, { size: { 'features.main': [1000, 900] }, renderers });
  const rows = () => t.entries('features.main').filter(([, name]) => name === NV);
  const state = () => {
    const m = t.ctx.graph.assembled.value?.view;
    const k = m?.resolved.tracks.split[t.trackIndex(NV)] ?? -1;
    return m?.prepared[k]?.state as SplitFeatures;
  };
  return { ...t, rows, state };
}

describe('split rows', () => {
  it("case 13: the variants grid opens expanded, a thin row per residue in UniProt's order under the label's entry", async () => {
    const t = await variants();
    expect(t.rows().map((r) => [r[0], r[2], r[4]])).toEqual([
      [EntryKind.track, -2, 22],
      ...ROWS.map((_, r) => [EntryKind.track, r, 8]),
    ]);
    expect(t.state().split?.labels).toEqual(ROWS);
  });

  it("LY2's prepare: each row its value's variants alone (R175H, R273H on H; R248W, R282W on W), the compact form whole", async () => {
    const t = await variants(case13(), true);
    const s = t.state() as unknown as Features;
    expect(s.split?.labels).toEqual(ROWS);
    const starts = (f: Features | undefined) => Array.from(f?.start ?? []);
    expect(starts(s.split?.rows[ROWS.indexOf('H')])).toEqual([174, 272]);
    expect(starts(s.split?.rows[ROWS.indexOf('W')])).toEqual([247, 281]);
    expect(starts(s.split?.rows[ROWS.indexOf('G')])).toEqual([]);
    expect(s.count).toBe(10);
    // Hits per row through the prepared rows, as with the stubs.
    const hit = Hit.create();
    t.ctx.layout.hitTest('features.main', t.xOf('features.main', 0, 174), t.yOf('features.main', NV, 15), hit);
    expect([hit.kind, hit.split]).toEqual([HoverKind.item, 15]);
  });

  it('`empty: "hide"`: only the residues with a variant get a row, in the grid\'s order; hits follow', async () => {
    const doc = case13();
    const walk = (n: { params?: Record<string, unknown>; children?: unknown[] }): void => {
      const sp = n.params?.split as { empty?: string } | null | undefined;
      if (sp) sp.empty = 'hide';
      for (const c of n.children ?? []) walk(c as never);
    };
    walk(doc.root as never);
    const t = await variants(doc);
    const present = t.state().split?.labels ?? [];
    expect(present.length).toBeLessThan(ROWS.length);
    expect(present).toEqual(ROWS.filter((v) => present.includes(v)));
    expect(t.rows().map((r) => r[2])).toEqual([-2, ...present.map((_, r) => r)]);
    const hit = Hit.create();
    t.ctx.layout.hitTest(
      'features.main',
      t.xOf('features.main', 0, 174),
      t.yOf('features.main', NV, present.indexOf('H')),
      hit,
    );
    expect([hit.kind, hit.split]).toEqual([HoverKind.item, present.indexOf('H')]);
  });

  it("hits per row: an item only on its value's row, the position elsewhere; the label entry hits no variant", async () => {
    const t = await variants();
    const hit = Hit.create();
    const x = t.xOf('features.main', 0, 174);
    const at = (row: number) => {
      t.ctx.layout.hitTest('features.main', x, t.yOf('features.main', NV, row), hit);
      return [hit.kind, hit.kind === HoverKind.item ? hit.item : -1, hit.split];
    };
    const ann = t.ctx.graph.assembled.value?.set.annotations.find((a) => a?.count === 10);
    expect(at(ROWS.indexOf('H'))).toEqual([HoverKind.item, ann?.item[2], ROWS.indexOf('H')]);
    expect(at(ROWS.indexOf('W'))).toEqual([HoverKind.seq, -1, ROWS.indexOf('W')]);
    expect(at(-2)).toEqual([HoverKind.seq, -1, -2]);
  });

  it('the toggle: `ui.expanded` collapses to one compact row and back; ranges open collapsed', async () => {
    const t = await variants();
    t.ctx.ui.toggleGroup(NV);
    expect(t.rows().map((r) => r[2])).toEqual([-1]);
    t.ctx.ui.toggleGroup(NV);
    expect(t.rows().length).toBe(ROWS.length + 1);
    // Case 7's UniProt track: collapsed at open (a ranges layer's default); expanded, a row per feature type.
    const c = await opened(case7());
    const uni = () => c.entries('complex-one.main').filter(([, n]) => n === 'complex-one.uniprot');
    expect(uni().map((r) => r[2])).toEqual([-1]);
    c.ctx.ui.toggleGroup('complex-one.uniprot');
    expect(uni().length).toBeGreaterThan(2);
  });

  it('a patch that sets `expanded` opens or closes the rows', async () => {
    const t = await variants();
    expect(t.rows().length).toBe(ROWS.length + 1);
    await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: NV, params: { expanded: false } }],
    });
    expect(t.rows().map((r) => r[2])).toEqual([-1]);
    await t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: NV, params: { expanded: true } }],
    });
    expect(t.rows().length).toBe(ROWS.length + 1);
  });
});
