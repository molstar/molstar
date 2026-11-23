/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ribbons layer's colour (_spec/README.md, Colour): its `color` channel as other layer types read it. Absent or by
// `system`, the blocks' hues (unchanged); a fixed colour on every polygon (stems and stubs too); a field read link by
// link: at the item level each ribbon its links' colour (a stem pair split by it), at the bins level the layer's one
// colour, the ranges' ink and `rondo.ribbons.color` said when the links' values are mixed; `opacity` scales the fill.

import { theme } from '@molstar/seqx-core/lang';
import { LightTokens, viewLegends } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { ColorMode, Kind, type RondoLay } from '../lay.js';
import { lay, ringSpec, Size, settle, toyChromosomes, toyLinks } from './helpers.js';

const rgb = (c: number) => (c >>> 8).toString(16).padStart(6, '0');
const Red = 'd62728';
const Blue = '1f77b4';
const byType = theme.field('type', {
  palette: theme.categorical({ colors: [`#${Blue}`, `#${Red}`], values: ['contact', 'other'] }),
});

/** The colours of the polygons of a kind (the hole left out). */
function colours(l: RondoLay, kind: number): Set<string> {
  const out = new Set<string>();
  for (let i = 0; i < l.polys.n; i++) if (l.polyKind[i] === kind) out.add(rgb(l.polyColor[i]));
  return out;
}

async function frame(o: Parameters<typeof ringSpec>[0]) {
  const s = await settle(ringSpec({ held: true, ...o }));
  const f = await s.frame({ size: Size });
  return { f, l: lay(f)! };
}

const chroms = toyChromosomes();

describe('a ribbons layer’s colour', () => {
  it('a fixed colour: every ribbon, stem, stub and fan in it, no issue', async () => {
    const { f, l } = await frame({ include: ['chr1', 'chr3'], layer: { color: `#${Red}` } });
    expect(l.mode).toBe(ColorMode.fixed);
    expect(l.nr).toBeGreaterThan(0);
    expect(l.nb).toBeGreaterThan(0);
    for (const k of [Kind.ribbon, Kind.trumpet, Kind.stub]) expect([...colours(l, k)], `kind ${k}`).toEqual([Red]);
    expect(f.layout.issues.filter((i) => i.code.startsWith('rondo.'))).toEqual([]);
  });

  it('by `system` (as `rondoRing` writes it): the blocks’ hues, as without a colour', async () => {
    const { l } = await frame({});
    expect(l.mode).toBe(ColorMode.system);
    expect(l.one).toBe(-1);
    expect(colours(l, Kind.ribbon).size).toBeGreaterThan(1);
  });

  it('`rondo:hues` (a scheme over `system` by default): the blocks’ hues, never one colour', async () => {
    const { l } = await frame({ hues: true });
    expect(l.mode).toBe(ColorMode.system);
    expect(colours(l, Kind.trumpet).size).toBeGreaterThan(1);
    expect(colours(l, Kind.ribbon).size).toBeGreaterThan(1);
  });

  it('a field at the item level: each ribbon its links’ colour, a stem pair split by it; stems the ranges’ ink', async () => {
    const { f, l } = await frame({ links: toyLinks(chroms, 300, 3), layer: { color: byType } });
    expect(l.items).toBe(true);
    expect(l.mode).toBe(ColorMode.field);
    expect(l.byRow).toBe(true);
    expect(colours(l, Kind.ribbon)).toEqual(new Set([Blue, Red]));
    // every ribbon's links are of its colour's type
    const d = l.read!.data!;
    const want = (k: number) => (d.kinds[d.kind[k]] === 'contact' ? Blue : Red);
    const pair = new Map<string, Set<string>>();
    for (let r = 0; r < l.nr; r++) {
      expect(rgb(l.ribColor[r])).toBe(want(l.ribRow[r]));
      const key = `${Math.min(l.ribA[r], l.ribB[r])}|${Math.max(l.ribA[r], l.ribB[r])}|${l.ribStrand[r]}`;
      if (!pair.has(key)) pair.set(key, new Set());
      pair.get(key)!.add(rgb(l.ribColor[r]));
    }
    // a stem pair holding both kinds is two ribbons, never one of mixed links
    expect([...pair.values()].some((c) => c.size === 2)).toBe(true);
    expect([...colours(l, Kind.trumpet)]).toEqual([rgb(LightTokens.colors['seq-range'])]);
    // the ribbons' hits told apart by colour
    expect(new Set(l.ribKey.slice(0, l.nr)).size).toBe(l.nr);
    expect(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.color')).toEqual([]);
  });

  it('held bins coloured by `type`: split by it, each ribbon its type’s colour, nothing said', async () => {
    const { f, l } = await frame({ links: toyLinks(chroms, 3000, 3), layer: { color: byType } });
    expect(l.items).toBe(false);
    expect(l.byRow).toBe(true);
    expect(colours(l, Kind.ribbon)).toEqual(new Set([Blue, Red]));
    const d = l.read!.data!;
    for (let r = 0; r < l.nr; r++)
      expect(rgb(l.ribColor[r])).toBe(d.kinds[d.kind[l.ribRow[r]]] === 'contact' ? Blue : Red);
    expect(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.color')).toEqual([]);
  });

  it('bins that carry not the colour’s field (q): one colour (the ink), said once, the legend saying so', async () => {
    const byQ = theme.field('q', { palette: theme.continuous({ colors: [`#${Blue}`, `#${Red}`] }) });
    const s = await settle(ringSpec({ held: true, links: toyLinks(chroms, 3000, 3), layer: { color: byQ } }));
    const f = await s.frame({ size: Size });
    const l = lay(f)!;
    expect(l.items).toBe(false);
    expect([...colours(l, Kind.ribbon)]).toEqual([rgb(LightTokens.colors['seq-range'])]);
    expect(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.color').map((i) => [i.severity, i.message])).toEqual(
      [['info', 'A colour by q is read link by link: these bins carry no q, drawn in one colour']],
    );
    const legend = viewLegends(f.layout).legends.find((x) => x.layers.includes(l.key))!;
    expect(legend.entries.kind === 'swatches' && legend.entries.items[0].label).toBe(
      'Links (bins in one colour: q read link by link)',
    );
  });

  it('more colours than a stem pair’s key holds: said, never merged silently', async () => {
    // 100 links between one pair of stems, each a kind of its own with its own colour, at the item level
    const one = Array.from({ length: 100 }, (_, i) => ({
      a: 0,
      b: 1,
      as: 1_000_000 + i,
      bs: 1_000_000 + i,
      type: `t${i}`,
    }));
    const values = one.map((l) => l.type);
    const colors = values.map((_, i) => `#${(0x100000 + i * 0x2a3b1).toString(16).slice(-6)}`);
    const byKind = theme.field('type', { palette: theme.categorical({ colors, values }) });
    const s = await settle(ringSpec({ held: true, links: one, layer: { color: byKind } }));
    const f = await s.frame({ size: Size });
    const l = lay(f)!;
    expect(l.items).toBe(true);
    expect(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.colors').length).toBe(1);
  });

  it('the colour’s scan kept by the read: a warm layout scans none', async () => {
    const { l } = await frame({ links: toyLinks(chroms, 300, 3), layer: { color: byType } });
    const kept = l.rowColor;
    const read = l.colorRead;
    expect(read).toBe(l.read);
    expect(l.rowColor).toBe(kept);
  });

  it('a field at the bins level over one kind (a filtered layer): its value’s colour, nothing said', async () => {
    const { f, l } = await frame({
      links: toyLinks(chroms, 8000, 3),
      query: '/(item :type other)',
      layer: { color: byType },
    });
    expect(l.items).toBe(false);
    expect([...colours(l, Kind.ribbon)]).toEqual([Red]);
    expect([...colours(l, Kind.trumpet)]).toEqual([Red]);
    expect(f.layout.issues.filter((i) => i.code === 'rondo.ribbons.color')).toEqual([]);
  });

  it('`opacity` scales the fill', async () => {
    const full = (await frame({ layer: { color: `#${Red}` } })).l;
    const faint = (await frame({ layer: { color: `#${Red}`, opacity: 0.5 } })).l;
    const a = (l: RondoLay) => l.polyColor[l.ribPoly[0]] & 0xff;
    expect(a(faint)).toBe(Math.round(a(full) / 2));
  });
});
