// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The card composed by providers (_spec/README.md "Hover card", Providers): the registry's `card-provider` entries by
// order, then name; the built-ins in the default spec; a host's added, filtered out, replaced; the floor said; a
// provider that throws dropped and said once (EC-262).

import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type CardProviderEntry, composeCardProviders } from '../card-provider.js';
import { BuiltInCardProviders, CardProviderNames } from '../card-providers.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { describerFor } from '../describe.js';
import type { SeqViewerRegistryEntry } from '../header-entries.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames;
const viewers: SeqViewer[] = [];

beforeEach(() => {
  frames = installFrames();
  vi.spyOn(DomHost, 'now').mockImplementation(() => 1000);
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
});

/** A protein with a domain and its letters. */
function spec() {
  const b = createSeqSpec({ title: 'providers' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
  const d = b.data.ranges('d', { target: p, rows: { start: [2], end: [8], id: ['A'], label: ['Domain'] } });
  b.views.view('v', (vw) =>
    vw.section('s', { of: p }, (s) => {
      s.track('letters').layer('l', { type: 'letters' });
      s.track('domains', { label: 'Domains' }).layer('d', { type: 'ranges', data: d });
    }),
  );
  return b.build();
}

async function make(registry: (base: readonly SeqViewerRegistryEntry[]) => readonly SeqViewerRegistryEntry[]) {
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({ ...base, registry: registry(base.registry ?? []), ui: undefined });
  viewers.push(v);
  await v.load(spec());
  v.resize(900);
  await v.settle();
  frames.flush();
  v.draw();
  return v;
}

/** The viewer's card at a letter (position 4, under the domain) and on the domain's item. */
function cards(v: SeqViewer) {
  const d = describerFor({ ...v, layout: () => v.frame?.layout });
  const s = v.frame!.layout.sections[0];
  const letters = s.tracks[0].layers[0];
  const domains = s.tracks[1].layers[0];
  const place = { section: s.section.address, system: '/data/p', position: 4 };
  return {
    letter: d.describe({ ...place, layer: letters.key }),
    item: d.describe({ ...place, layer: domains.key, element: { owner: '/data/d', key: 'A' } }),
  };
}

const note = (name: string, order: number | undefined, label = name): CardProviderEntry => ({
  kind: 'card-provider',
  name,
  ...(order === undefined ? {} : { order }),
  describe(at, out) {
    if (Number.isNaN(at.position)) return;
    out.row(label, `at ${at.position}`);
    out.words.push(`${label} at ${at.position}`);
  },
});

describe('card providers', () => {
  it('orders by order, then name; an entry with no order after the ordered ones', () => {
    const names = composeCardProviders([
      note('lab:b', 5),
      note('lab:none', undefined),
      note('lab:a', 5),
      ...BuiltInCardProviders,
    ]).map((p) => p.name);
    // a layer type's own words first (`seq:layer`, 1)
    expect(names.slice(0, 4)).toEqual([CardProviderNames.layer, 'lab:a', 'lab:b', CardProviderNames.region]);
    expect(names.at(-1)).toBe('lab:none');
    expect(names.indexOf(CardProviderNames.place)).toBeLessThan(names.indexOf(CardProviderNames.block));
  });

  it('the default spec lists every built-in provider: the card as the describer’s own', async () => {
    const v = await make((r) => r);
    const listed = (v.registry ?? []).filter((e) => e.kind === 'card-provider').map((e) => e.name);
    expect(listed).toEqual(BuiltInCardProviders.map((p) => p.name));
    const { letter, item } = cards(v);
    expect(letter?.rows.map((r) => r.label)).toEqual(['index', 'Letter', 'Domain']);
    expect(item?.title).toBe('Domain');
    expect(v.behaviorIssues.value).toEqual([]);
  });

  it('a host’s entry sits by its order: after the place’s rows, before the block', async () => {
    const v = await make((r) => [...r, note('lab:assay', 95, 'Assay')]);
    const { letter } = cards(v);
    expect(letter?.rows.map((r) => r.label)).toEqual(['index', 'Letter', 'Domain', 'Assay']);
    expect(letter?.text).toMatch(/, Assay at 4$/);
  });

  it('an entry under a built-in’s name replaces it; one filtered out says nothing', async () => {
    const own: CardProviderEntry = {
      kind: 'card-provider',
      name: CardProviderNames.item,
      order: 40,
      describe(at, out) {
        if (!at.place.element) return;
        out.title = `item ${String(at.place.element.key)}`;
        out.words.push(out.title);
        out.whole = true;
      },
    };
    const v = await make((r) => [...r.filter((e) => e.kind !== 'card-provider' || e.name !== 'seq:block'), own]);
    expect(cards(v).item).toEqual({ title: 'item A', rows: [], text: 'item A' });
  });

  it('a composition without seq:item or seq:place is said: the hover contract’s floor (EC-262)', async () => {
    const v = await make((r) => r.filter((e) => e.kind !== 'card-provider' || e.name !== CardProviderNames.place));
    expect(v.behaviorIssues.value.map((i) => i.code)).toEqual(['seq.viewer.card.floor']);
    // the letter has no place's rows: the card says nothing there
    expect(cards(v).letter).toBeUndefined();
  });

  it('a provider that throws: its rows dropped, the rest shown, said once (EC-262)', async () => {
    const bad: CardProviderEntry = {
      kind: 'card-provider',
      name: 'lab:bad',
      order: 95,
      describe(_, out) {
        out.row('half', 'written');
        out.title = 'broken';
        throw new Error('no assay');
      },
    };
    const v = await make((r) => [...r, bad]);
    const first = cards(v).letter;
    const again = cards(v).letter;
    expect(first?.rows.map((r) => r.label)).toEqual(['index', 'Letter', 'Domain']);
    expect(first?.title).not.toBe('broken');
    expect(again).toEqual(first);
    const said = v.behaviorIssues.value.filter((i) => i.code === 'seq.viewer.card.provider');
    expect(said).toHaveLength(1);
    expect(said[0].message).toContain('lab:bad');
  });

  it('with no layer hook the built-in cards are byte-equal with and without `seq:layer`', async () => {
    const withIt = cards(await make((r) => r));
    const without = cards(
      await make((r) => r.filter((e) => e.kind !== 'card-provider' || e.name !== CardProviderNames.layer)),
    );
    expect(JSON.stringify(withIt)).toBe(JSON.stringify(without));
  });

  it('a host provider’s title by its order: an earlier one under the item’s title, a later one over it', async () => {
    const titled = (name: string, order: number): CardProviderEntry => ({
      kind: 'card-provider',
      name,
      order,
      describe(at, out) {
        if (at.place.element) out.title = `${name} title`;
      },
    });
    expect(cards(await make((r) => [...r, titled('lab:early', 5)])).item?.title).toBe('Domain');
    expect(cards(await make((r) => [...r, titled('lab:late', 95)])).item?.title).toBe('lab:late title');
  });

  it('a genome’s base has no index row (its number is its position); a protein’s residue keeps it', async () => {
    const b = createSeqSpec({ title: 'genome' });
    const g = b.data.assembly('toy', assembly('Toy', [{ name: 'chr1', length: 10_000 }]));
    b.views.view('v', (vw) => {
      genome(vw, 'g', g);
    });
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    await v.load(b.build());
    v.resize(900);
    await v.settle();
    frames.flush();
    v.draw();
    const d = describerFor({ ...v, layout: () => v.frame?.layout });
    const s = v.frame!.layout.sections[0];
    const card = d.describe({ section: s.section.address, system: '/data/toy/(sequence :name chr1)', position: 41 })!;
    expect(card.title).toBe('chr1 42');
    expect(card.rows.map((r) => r.label)).not.toContain('index');
    expect(cards(await make((r) => r)).letter?.rows.map((r) => r.label)).toContain('index');
  });
});
