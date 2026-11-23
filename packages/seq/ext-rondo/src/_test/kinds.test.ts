/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Kinds of links (_spec/README.md, Kinds of links): `rondoRing`'s `kinds`, a `rondo:ribbons` layer a kind over one
// pairs node, each its data filter on `type`, in its own track (the last drawn on top), its width by its own rule
// (count, or a field summed), its legend entry; each neck the sum of its layer's ribbons, stubs and fan; a layer hidden
// draws nothing of its kind.

import { createSeqContext } from '@molstar/seqx-core/context/default';
import { viewLegends } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { DefaultRondoConfig } from '../config.js';
import { rondo } from '../index.js';
import type { RondoLay } from '../lay.js';
import { type Node, lays, ringSpec, Size, settle, toyChromosomes, toyLinks } from './helpers.js';

const chroms = toyChromosomes();
const links = toyLinks(chroms, 2000, 5);
const Red = '#d62728';
const kinds = [
  { type: 'contact', label: 'Contacts', opacity: 0.5 },
  { type: 'other', label: 'Others', color: Red },
] as const;

async function frame(spec: Node) {
  const s = await settle(spec);
  return s.frame({ size: Size });
}

/** The links of a kind (held rows), and their q (the helpers' `q`: 0.001 + (i·37 mod 100) / 1000). */
const ofKind = (type: string) =>
  links.map((l, i) => ({ ...l, q: 0.001 + ((i * 37) % 100) / 1000 })).filter((l) => l.type === type);

/** Everything a layer read, summed: its ribbons, its stubs and the links within its stems. */
function total(l: RondoLay): number {
  let t = 0;
  for (let r = 0; r < l.nr; r++) t += l.ribCount[r];
  for (let b = 0; b < l.nb; b++) t += l.stubCount[b];
  for (let s = 0; s < l.ns; s++) t += l.stemWithin[s];
  return t;
}

/** Each placed stem's neck: its slices (ribbons, stubs) and its fan fill it exactly (k × its count). */
function necksAdd(l: RondoLay) {
  const f = l.frame;
  const perPx = f.perPx(DefaultRondoConfig.stemDepth * f.room);
  for (let s = 0; s < l.ns; s++) {
    if (!l.stemPlaced[s] || l.stemCount[s] <= 0) continue;
    if (l.k * l.stemCount[s] < 1 / l.dpr) continue; // a neck floored to a device px: its parts are thinner
    const mid = (l.stemT0[s] + l.stemT1[s]) / 2;
    const end = l.stemCut[s] + l.k * l.stemFan[s] * perPx;
    expect(end, `stem ${s}`).toBeCloseTo(mid + l.stemNeck[s], 9);
  }
}

describe('kinds of links', () => {
  it('a layer a kind over one pairs node, in draw order, each reading its kind alone; totals add up', async () => {
    const f = await frame(ringSpec({ held: true, links, kinds }));
    const [contacts, others] = lays(f);
    expect([contacts.label, others.label]).toEqual(['Contacts', 'Others']);
    // the last listed drawn last: on top
    expect(contacts.key.includes('links-contact')).toBe(true);
    expect(others.key.includes('links-other')).toBe(true);
    for (const [l, type] of [
      [contacts, 'contact'],
      [others, 'other'],
    ] as const) {
      const d = l.read!.data!;
      for (let k = 0; k < d.n; k++) if (l.read!.shown[k]) expect(d.kinds[d.kind[k]]).toBe(type);
      expect(total(l), type).toBe(ofKind(type).length);
      necksAdd(l);
    }
    expect(total(contacts) + total(others)).toBe(links.length);
    // one rim for both: the hole's
    expect(contacts.frame.r).toBe(others.frame.r);
    // each its own k (the widest stem of its own kind)
    const widest = (l: RondoLay) => Math.max(...Array.from(l.stemCount.subarray(0, l.ns)));
    expect(contacts.k).toBeCloseTo(DefaultRondoConfig.maxWidthPx / widest(contacts), 9);
    expect(others.k).toBeCloseTo(DefaultRondoConfig.maxWidthPx / widest(others), 9);
  });

  it('a width by a field: q summed, at the item level and at the bins level', async () => {
    const f = await frame(ringSpec({ held: true, links, kinds: kinds.map((k) => ({ ...k, width: 'q' })) }));
    const [contacts, others] = lays(f);
    expect(contacts.items).toBe(false);
    expect(others.items).toBe(true);
    for (const [l, type] of [
      [contacts, 'contact'],
      [others, 'other'],
    ] as const) {
      const q = ofKind(type).reduce((a, x) => a + x.q, 0);
      expect(total(l), type).toBeCloseTo(q, 3);
      necksAdd(l);
    }
    expect(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.width')).toEqual([]);
  });

  it('a width by a field the links lack: by count, said', async () => {
    const f = await frame(ringSpec({ held: true, links, kinds: [{ type: 'other', width: 'score' }] }));
    expect(total(lays(f)[0])).toBe(ofKind('other').length);
    expect(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.width').map((i) => i.message)).toEqual([
      'Ribbons wide by count: the links read carry no numeric score',
    ]);
  });

  it('a legend entry each: the hued kind says so, the coloured one its colour', async () => {
    const f = await frame(ringSpec({ held: true, links, kinds }));
    const { legends } = viewLegends(f.layout);
    const mine = legends.filter((l) => l.layers.some((k) => k.includes('/links-')));
    expect(
      mine.map((l) => [l.title, l.entries.kind === 'swatches' ? l.entries.items.map((i) => i.label) : []]),
    ).toEqual([
      ['Contacts', ['Contacts: by chromosome']],
      ['Others', ['Others']],
    ]);
    // the hued kind's entry is its words alone (no one hue's swatch: it draws in every chromosome's)
    expect(mine[0].entries.kind === 'swatches' && mine[0].entries.items[0].color & 0xff).toBe(0);
    const red = mine[1].entries.kind === 'swatches' ? mine[1].entries.items[0].color : 0;
    expect((red >>> 8).toString(16)).toBe('d62728');
  });

  it('the spec checks clean; a reader without the extension draws each kind as `pairs`, its filter kept', async () => {
    const spec = ringSpec({ held: true, links, kinds });
    expect(createSeqContext({ extensions: [rondo] }).check(spec)).toEqual([]);
    const s = await settle(spec, []);
    const f = await s.frame({ size: Size });
    const pairs = f.layout.sections.flatMap((x) => x.sectionLays).filter((l) => l.type.name === 'pairs');
    expect(pairs.length).toBe(2);
    expect(pairs.map((l) => (l.lay as { prep?: { filtered: number } }).prep?.filtered)).toEqual([0, 0]);
  });

  it('kinds over links read by region: their items each its kind; their bins drawn once, said', async () => {
    // items (under the gate): core's region reads keep each layer's base, so each layer reads its kind alone
    const few = toyLinks(chroms, 300, 5);
    const byKind = (type: string) => few.filter((l) => l.type === type).length;
    const items = await frame(ringSpec({ links: few, kinds }));
    const [c, o] = lays(items);
    expect([c.items, o.items]).toEqual([true, true]);
    expect([total(c), total(o)]).toEqual([byKind('contact'), byKind('other')]);
    expect(items.layout.issues.filter((i) => i.code === 'rondo.ribbons.filter')).toEqual([]);
    // bins: no fields to filter by, so the first layer draws them, the next draws nothing, both say why
    const f = await frame(ringSpec({ links, kinds }));
    const [first, next] = lays(f);
    expect(first.items).toBe(false);
    expect(total(first)).toBe(links.length);
    expect([next.empty, next.polys.n]).toEqual(['shared', 0]);
    expect(new Set(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.filter').map((i) => i.address)).size).toBe(
      2,
    );
  });

  it('a kind hidden (`layer.hidden`): its layer draws nothing, the other unchanged', async () => {
    const spec = ringSpec({ held: true, links, kinds });
    const both = lays(await frame(spec));
    const hide = (n: Node, inContacts = false): Node => {
      const here = inContacts || (n.kind === 'track' && n.id === 'links-contact');
      if (here && n.kind === 'layer') return { ...n, params: { ...n.params, hidden: true } };
      return n.children ? { ...n, children: n.children.map((c) => hide(c, here)) } : n;
    };
    const one = lays(await frame(hide(spec)));
    expect(one.map((l) => l.label)).toEqual(['Others']);
    expect(total(one[0])).toBe(total(both[1]));
  });
});
