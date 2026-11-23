/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The edge cases the design takes (_spec/README.md): bins drawn until items land, nothing pending left blank
// (EC-169); a hue per chromosome stable across sections (EC-218); a minimum by the data filter, counted (EC-225);
// two assemblies on one ring, sided (EC-210); the words a reader hears (the card, ARIA, the cursor).

import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { rondo } from '../index.js';
import { Kind, type RondoLay } from '../lay.js';
import { RibbonsType } from '../ribbons.js';
import { lay, ringSpec, Size, settle, toyChromosomes } from './helpers.js';

const hueOfSystem = (l: RondoLay, name: string) => l.hues[l.names.indexOf(name) % l.hues.length];

describe('edge cases', () => {
  it('EC-169: before the bins land the layer is pending (said, nothing drawn); once they land, ribbons', async () => {
    const s = await settle(ringSpec());
    const first = lay(s.renderer.frame({ size: Size }))!;
    expect(first.empty).toBe('pending');
    expect(RibbonsType.summary!(first)).toBe('links loading');
    const l = lay(await s.frame({ size: Size }))!;
    expect([l.empty, l.nr > 0]).toEqual(['', true]);
  });

  it('EC-169: a frame drawn settled while its bins are wanted is still pending, never “no links”', async () => {
    // a host settles for a frame between a source's open and its tiles' wants: the reader is waiting there
    const s = await settle(ringSpec());
    const first = lay(s.renderer.frame({ size: Size, settled: true }))!;
    expect(first.reader.waiting).toBe(true);
    expect([first.empty, RibbonsType.summary!(first)]).toEqual(['pending', 'links loading']);
  });

  it('EC-218: a chromosome keeps its hue on the genome’s ring, on a pair’s and on a stub pointing at it', async () => {
    const all = lay(await (await settle(ringSpec())).frame({ size: Size }))!;
    const pair = lay(await (await settle(ringSpec({ include: ['chr5', 'chr3'] }))).frame({ size: Size }))!;
    const one = lay(await (await settle(ringSpec({ include: ['chr1'] }))).frame({ size: Size }))!;
    for (const name of ['chr3', 'chr5']) expect(hueOfSystem(pair, name)).toBe(hueOfSystem(all, name));
    // the domain is the assembly's order, before `include`: chr3 the third hue everywhere
    expect(pair.names.indexOf('chr3')).toBe(2);
    const stubTo3 = [...Array(one.nb).keys()].find(
      (b) => one.names[one.stubPartner[b]] === 'chr3' && one.stubPoly[b] >= 0,
    )!;
    expect(one.polyColor[one.stubPoly[stubTo3]] >>> 8).toBe(hueOfSystem(all, 'chr3') >>> 8);
  });

  it('EC-225: the data filter leaves links out of every count, and says how many', async () => {
    const all = lay(await (await settle(ringSpec({ held: true }))).frame({ size: Size }))!;
    const some = lay(
      await (await settle(ringSpec({ held: true, query: '/(except :of (item :type other))' }))).frame({ size: Size }),
    )!;
    const sum = (l: RondoLay) => {
      let n = 0;
      for (let s = 0; s < l.ns; s++) n += l.stemCount[s] + 2 * l.stemWithin[s];
      return n;
    };
    expect(RibbonsType.filtered!(some)).toBeGreaterThan(0);
    expect(sum(some)).toBeLessThan(sum(all));
    expect(RibbonsType.filtered!(all)).toBe(0);
  });

  it('EC-210: two assemblies on one ring: synteny ribbons between them, each end on its own assembly', async () => {
    const chroms = toyChromosomes(3);
    const b = createSeqSpec({ title: 'two genomes', extensions: [rondo] });
    const peach = b.data.assembly('peach', assembly('peach', chroms));
    const grape = b.data.assembly('grape', assembly('grape', chroms));
    const n = 60;
    const links = b.data.pairs('synteny', {
      a: peach,
      b: grape,
      rows: {
        a_name: Array.from({ length: n }, (_, i) => chroms[i % 3].name),
        b_name: Array.from({ length: n }, (_, i) => chroms[(i + 1) % 3].name),
        a_start: Int32Array.from({ length: n }, (_, i) => i * 50_000),
        a_end: Int32Array.from({ length: n }, (_, i) => i * 50_000 + 40_000),
        b_start: Int32Array.from({ length: n }, (_, i) => i * 60_000),
        b_end: Int32Array.from({ length: n }, (_, i) => i * 60_000 + 40_000),
      },
    } as never);
    b.views.view('v', (v) =>
      v
        .section('ring', { of: [peach, grape], layout: 'rings', center: 'none' } as never)
        .track('links')
        .layer('ribbons', { type: 'rondo:ribbons', fallback: 'pairs', data: links } as never),
    );
    const l = lay(await (await settle(b.build())).frame({ size: Size }))!;
    expect(l.nr).toBeGreaterThan(0);
    // every ribbon joins a peach stem to a grape stem: its two ends' systems differ
    for (let r = 0; r < l.nr; r++)
      expect(l.systems[l.stemSystem[l.ribA[r]]]).not.toBe(l.systems[l.stemSystem[l.ribB[r]]]);
  });
});

describe('the words', () => {
  it('a ribbon’s card names its two stems and its count; a stem’s its links and partners; the ARIA its counts', async () => {
    const l = lay(await (await settle(ringSpec())).frame({ size: Size }))!;
    const r = l.strongest;
    const card = RibbonsType.describe!({ element: { owner: l.owner, key: l.ribKey[r] }, position: 0 }, l)!;
    expect(card[0]).toMatch(/^chr\d:[\d.]+–[\d.]+ Mb ↔ chr\d:[\d.]+–[\d.]+ Mb$/);
    expect(card[1]).toMatch(new RegExp(`^Links: ${Math.round(l.ribCount[r]).toLocaleString('en-US')} links`));
    expect(card).toContain(
      `Bins: ${l.level >= 1e3 ? Math.round(l.level / 1024) : l.level}`.slice(0, 6) +
        card.find((x) => x.startsWith('Bins: '))!.slice(6),
    );
    const stem = l.ribA[r];
    const s = RibbonsType.describe!({ element: { owner: l.owner, key: l.stemKey[stem] }, position: 0 }, l)!;
    expect(s[1]).toMatch(/^Links: [\d,]+ links to \d+ chromosomes?$/);
    expect(RibbonsType.summary!(l)).toMatch(/^[\d,]+ ribbons, 6 chromosomes; strongest chr\d ↔ chr\d, [\d,]+$/);
    // the cursor at a stem's start: its links
    expect(RibbonsType.announce!(l, l.systems[l.stemSystem[stem]], l.stemStart[stem])).toBe(
      `stem of ${Math.round(l.stemCount[stem]).toLocaleString('en-US')} links`,
    );
  });

  it('a stub says its partner, not shown; a ring with no links says so', async () => {
    const one = lay(await (await settle(ringSpec({ include: ['chr2'] }))).frame({ size: Size }))!;
    const i = [...Array(one.polys.n).keys()].find((k) => one.polyKind[k] === Kind.stub)!;
    const lines = RibbonsType.describe!({ ...one.hits[i], position: one.hits[i].position ?? Number.NaN }, one)!;
    expect(lines[0]).toMatch(/^chr2:.+ → chr\d$/);
    expect(lines[1]).toMatch(/, not shown$/);
    const none = lay(await (await settle(ringSpec({ links: [] }))).frame({ size: Size }))!;
    expect(RibbonsType.summary!(none)).toBe('no links');
  });
});
