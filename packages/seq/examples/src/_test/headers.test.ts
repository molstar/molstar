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
        for (const t of s.tracks) {
          if (t.well.y >= layout.height) continue;
          const h = t.header!;
          expect(h, t.track.address).toBeDefined();
          const rows = t.track.layers.some((l) => l.type === 'rows');
          // the slots a track takes: none on a short track or a rows track (its title row's), else the columns'
          const want = h.short || rows ? 0 : Math.min(cols.n, 2);
          expect(h.slots.length, t.track.address).toBe(want);
          for (const slot of h.slots) {
            if (slot.icon) icons++;
            if (slot.y + slot.height > layout.height) continue;
            expect(hitView(layout, slot.x + slot.width / 2, slot.y + slot.height / 2)).toMatchObject({
              target: 'header-slot',
              state: slot.name,
              layer: t.track.address,
            });
          }
          // every state its entries give: drawn (a chip its text, a mark a squircle) and hit as its entry
          expect(h.states.map((x) => x.name)).toContain('host:bound');
          for (const st of h.states) {
            if (!st.mark) expect(labels).toContain(`text "${st.text}"`);
            if (st.box.y + st.box.height > layout.height) continue;
            expect(hitView(layout, st.box.x + st.box.width / 2, st.box.y + st.box.height / 2)).toMatchObject({
              target: 'header-state',
              state: st.name,
              layer: t.track.address,
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
