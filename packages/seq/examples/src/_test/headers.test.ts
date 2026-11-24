/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headers drawn against the data in every example view (_spec/README.md "Render", Track headers): every shown
// track of a section laid out as tracks has its header line; with a host's entries and the actions shown, each slot it
// takes is drawn in the tip pass and hit as its entry, and every state is drawn and hit as its entry. The p53 views are
// the quick representatives; the rest are medium.

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import {
  actionColumns,
  builtInHeaderStates,
  type HeaderPartsOf,
  hitView,
  openModel,
  recordDrawLists,
  Renderer,
  renderConfig,
} from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets } from '../node-assets.js';
import { type ExampleView, ExampleViews } from '../views.js';

const Cross = { commands: [0, 1, 0, 1], points: [4, 4, 20, 20, 4, 20, 20, 4] };

/** A host's entries: a slot with an icon and one without, the menu, the built-in states and a state of its own. */
const Host: HeaderPartsOf = (f) => {
  const own = builtInHeaderStates(f)?.states ?? [];
  if (!f.track) return { states: own, menu: true };
  return {
    slots: [
      { name: 'seq:hide', label: 'Hide track', icon: Cross, enabled: true },
      { name: 'host:color-by', label: 'Color by', enabled: true },
    ],
    states: [...own, { name: 'host:bound', text: 'bound' }],
    menu: true,
  };
};

const Size = { width: 960, height: 480 };

/** A number as the record writes it: two decimals at most, no `-0`. */
const fmt = (n: number) => {
  const r = Math.round(n * 100) / 100;
  return Object.is(r, -0) ? '0' : String(r);
};

async function frameOf(view: ExampleView) {
  const ctx = createSeqContext({
    registry: [{ kind: 'format-runtime', runtime: SeqRuntime }, ...(view.registry?.({ delayMs: 0 }) ?? [])],
    ports: { assets: () => nodeFixtureAssets() },
  });
  const store = new DataStore();
  const spec = view.build();
  const opened = await openModel(spec, store, {
    runtimes: ctx.runtimes,
    config: ctx.resolveConfig,
    assets: ctx.assets,
  });
  const renderer = new Renderer({
    model: opened.model,
    view: opened.view!,
    data: store,
    tokens: ctx.tokens,
    text: ctx.text,
    types: ctx.types,
    layouts: ctx.layouts,
    sizings: ctx.sizings,
    config: renderConfig({ header: { actions: { show: 'always' } } }),
  });
  const input = { size: Size, settled: true, headers: Host };
  let f = renderer.frame(input);
  for (let i = 0; i < 6; i++) {
    store.want(f.wants);
    await store.whenSettled();
    const next = renderer.frame(input);
    if (!next.rebuilt) break;
    f = next;
  }
  return { ctx, f: renderer.frame(input) };
}

describe.each(ExampleViews.map((v) => [`${v.dataset}/${v.id}`, v] as const))('%s', (_, view) => {
  it(
    'every shown track has its header line; each slot and state drawn and hit as its entry',
    view.dataset === 'p53' ? {} : { tags: ['medium'] },
    async () => {
      const { ctx, f } = await frameOf(view);
      const layout = f.layout;
      const tip = recordDrawLists(f.lists, ['tip']);
      const labels = recordDrawLists(f.lists, ['labels']);
      const cols = actionColumns(layout.config);
      let icons = 0;
      for (const s of layout.sections) {
        if (s.layout.name !== 'tracks') {
          for (const t of s.tracks) expect(t.header).toBeUndefined();
          continue;
        }
        for (let i = 0; i < s.tracks.length; i++) {
          const t = s.tracks[i];
          if (t.well.y >= layout.height) continue;
          // a folded group's row: one header line, its first band's, the group's
          if (s.groups?.some((g) => g.folded && i > g.first && i <= g.last)) {
            expect(t.header).toBeUndefined();
            continue;
          }
          const h = t.header!;
          expect(h, t.track.address).toBeDefined();
          const rows = t.track.layers.some((l) => l.type === 'rows');
          // the slots a track takes: none on a short track or a rows track (its title row's), else the columns'
          const want = h.short || rows ? 0 : Math.min(cols.n, 2);
          expect(h.slots.length, t.track.address).toBe(want);
          for (const slot of h.slots) {
            if (slot.part.icon) icons++;
            if (slot.y + slot.height > layout.height) continue;
            expect(hitView(layout, slot.x + slot.width / 2, slot.y + slot.height / 2)).toMatchObject({
              target: 'header-slot',
              state: slot.part.name,
              layer: h.group ?? t.track.address,
            });
          }
          // a rows track draws no states (its header is its title row's, not built); the rest every state its entries
          // give, each drawn where it was laid out (a chip its text, a mark a squircle) and hit as its entry
          if (rows) expect(h.states).toEqual([]);
          else expect(h.states.map((x) => x.name)).toContain('host:bound');
          for (const st of h.states) {
            const drawn = st.mark
              ? `rect ${fmt(st.x)} ${fmt(st.y)} ${fmt(st.width)} ${fmt(st.width)} `
              : `text "${st.text}" ${fmt(st.x)} ${fmt(st.y)} start`;
            expect(labels, `${t.track.address} ${st.name}`).toContain(drawn);
            if (st.box.y + st.box.height > layout.height) continue;
            expect(hitView(layout, st.box.x + st.box.width / 2, st.box.y + st.box.height / 2)).toMatchObject({
              target: 'header-state',
              state: st.name,
              layer: h.group ?? t.track.address,
            });
          }
        }
      }
      // every slot with an icon drawn once in the tip pass (`always`), the ones without as their label's first letter
      expect(tip.match(/^path /gm)?.length ?? 0).toBe(icons);
      ctx.retire();
    },
  );
});
