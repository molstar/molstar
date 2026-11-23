/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The card names the kind (_spec/README.md, Kinds of links, The card): a link's card lists its fields (its type, its
// ends' strands, `q`, gene names where the data has them), through the `describe` hook; a bin's and a stem's card their
// count on each kind's layer of the section.

import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { rondo } from '../index.js';
import { Kind, type RondoLay } from '../lay.js';
import { rondoRing } from '../ring.js';
import { describeRondo } from '../words.js';
import { lays, Size, settle, toyChromosomes, toyLinks } from './helpers.js';

const chroms = toyChromosomes();

/** 3,000 contacts (bins at the fit) and three translocations (items), two of them fusions with gene names. */
function spec() {
  const contacts = toyLinks(chroms, 3000, 9).map((l) => ({ ...l, type: 'contact' }));
  const tx = [
    { a: 0, as: 2_000_000, b: 3, bs: 4_000_000, sa: '+', sb: '-', genes: ['XBCR', 'XABL1'] },
    { a: 1, as: 1_000_000, b: 4, bs: 3_000_000, sa: '-', sb: '-', genes: ['XETV6', 'XRUNX'] },
    { a: 2, as: 5_000_000, b: 5, bs: 1_000_000, sa: '+', sb: '+', genes: [null, null] },
  ];
  const rows = {
    a_name: [...contacts.map((l) => chroms[l.a].name), ...tx.map((t) => chroms[t.a].name)],
    b_name: [...contacts.map((l) => chroms[l.b].name), ...tx.map((t) => chroms[t.b].name)],
    a_start: [...contacts.map((l) => l.as), ...tx.map((t) => t.as)],
    a_end: [...contacts.map((l) => l.as + 500), ...tx.map((t) => t.as + 500)],
    b_start: [...contacts.map((l) => l.bs), ...tx.map((t) => t.bs)],
    b_end: [...contacts.map((l) => l.bs + 500), ...tx.map((t) => t.bs + 500)],
    a_strand: [...contacts.map(() => '+'), ...tx.map((t) => t.sa)],
    b_strand: [...contacts.map(() => '+'), ...tx.map((t) => t.sb)],
    type: [...contacts.map(() => 'contact'), ...tx.map(() => 'translocation')],
    q: [...contacts.map(() => 0.01), 0.0012, 0.004, 0.03],
    gene_a: [...contacts.map(() => null), ...tx.map((t) => t.genes[0])],
    gene_b: [...contacts.map(() => null), ...tx.map((t) => t.genes[1])],
  };
  const b = createSeqSpec({ title: 'kinds', extensions: [rondo] });
  const g = b.data.assembly('toy', assembly('Toy', chroms));
  const links = b.data.pairs('links', { a: g, rows } as never);
  b.views.view('v', (v) => {
    rondoRing(v, {
      assembly: g,
      links,
      square: false,
      kinds: [
        { type: 'contact', label: 'Contacts' },
        { type: 'translocation', label: 'Translocations', color: '#d62728' },
      ],
    });
  });
  return b.build();
}

/** The card of a polygon of a kind (the first such, or the one `pick` chooses). */
function card(l: RondoLay, kind: number, pick: (ref: number) => boolean = () => true) {
  for (let i = 0; i < l.polys.n; i++)
    if (l.polyKind[i] === kind && pick(l.polyRef[i]))
      return { ref: l.polyRef[i], lines: describeRondo(l.hits[i] as never, l)! };
  throw new Error('no such polygon');
}

describe('the card names the kind', () => {
  it('a link: its type, strands, q and gene names (none where the data has none)', async () => {
    const f = await (await settle(spec())).frame({ size: Size });
    const [, tx] = lays(f);
    expect(tx.items).toBe(true);
    expect(tx.nr).toBe(3);
    const fusion = card(tx, Kind.ribbon, (r) => tx.read!.data!.aStart[tx.ribRow[r]] === 2_000_000).lines;
    expect(fusion.slice(1)).toEqual([
      'Links: 1 link (translocation)',
      'Strand: reversed',
      'Type: translocation',
      'Strands: + −',
      'q: 0.0012',
      'gene_a: XBCR',
      'gene_b: XABL1',
    ]);
    const plain = card(tx, Kind.ribbon, (r) => tx.read!.data!.aStart[tx.ribRow[r]] === 5_000_000).lines;
    expect(plain.slice(1)).toEqual(['Links: 1 link (translocation)', 'Type: translocation', 'Strands: + +', 'q: 0.03']);
  });

  it('a bin and a stem: their count on each kind’s layer', async () => {
    const f = await (await settle(spec())).frame({ size: Size });
    const [contacts] = lays(f);
    expect(contacts.items).toBe(false);
    const [, tx] = lays(f);
    const n = (k: number) => `${k} link${k === 1 ? '' : 's'}`;
    // every drawn bin ribbon: its contacts (its own links: one strand) and the translocations between its stems
    let checked = 0;
    for (let r = 0; r < contacts.nr; r++) {
      if (!contacts.ribDrawn[r]) continue;
      const i = contacts.ribPoly[r];
      const lines = describeRondo(contacts.hits[i] as never, contacts)!;
      const [a, b] = [contacts.ribA[r], contacts.ribB[r]];
      const d = tx.read!.data!;
      let t = 0;
      for (let k = 0; k < d.n; k++) {
        if (!tx.read!.shown[k]) continue;
        const ea = tx.names[d.aUnit![k]];
        const eb = tx.names[d.bUnit![k]];
        const inS = (s: number, name: string, x: number) =>
          contacts.names[contacts.stemSystem[s]] === name && contacts.stemStart[s] <= x && x < contacts.stemEnd[s];
        if (
          (inS(a, ea, d.aStart[k]) && inS(b, eb, d.bStart[k])) ||
          (inS(b, ea, d.aStart[k]) && inS(a, eb, d.bStart[k]))
        )
          t++;
      }
      expect(lines).toContain(`Contacts: ${n(contacts.ribLinks[r])}`);
      expect(lines).toContain(`Translocations: ${n(t)}`);
      checked++;
    }
    expect(checked).toBeGreaterThan(10);
    // a stem: its links on each layer
    const at = (s: number, name: string, x: number) =>
      contacts.names[contacts.stemSystem[s]] === name && contacts.stemStart[s] <= x && x < contacts.stemEnd[s];
    const stem = card(contacts, Kind.trumpet, (s) => at(s, 'chr1', 2_000_000)).lines;
    expect(stem.some((x) => /^Contacts: \d+ links$/.test(x))).toBe(true);
    expect(stem).toContain('Translocations: 1 link');
  });

  it('one layer alone: no count per layer', async () => {
    const f = await (await settle(spec())).frame({ size: Size });
    const [contacts] = lays(f);
    const only = { ...contacts, peers: [contacts] } as RondoLay;
    const { lines } = card(only, Kind.trumpet);
    expect(lines.some((x) => x.startsWith('Contacts:'))).toBe(false);
  });
});
