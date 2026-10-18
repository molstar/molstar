/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Cyclic coordinate systems in the format (CY1; seqspec:format#positions, seqspec:format#data-nodes,
// seqspec:format#blocks, seqspec:format#focus): `sequence.cyclic`, a block's range across the origin and its `origin`,
// as pieces with two-segment layout maps; a positions filter and a focus across the origin (an annotation's as its
// shortest arc); and what is reported where a system is not cyclic.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { Coords, OrderedSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  compile,
  type ResolvedView,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
  validate,
} from '../../src/format/index.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
const all = supportedAll(reg);
async function resolve(doc: SeqViewSpecDoc, view: string) {
  const { set, issues } = await compile(doc, reg, fixtureLoaders());
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  return { set, ...resolveView(doc, set, view, reg, all) };
}
/** Block b's indices in shown order, through its layout map alone. */
function shownOrder(v: ResolvedView, b: number) {
  const out: number[] = [];
  const map = v.blocks.layout[b];
  for (let k = 0; k < map.toStop; k++) Coords.forEachFrom(map, k, k + 1, (s) => out.push(s), null);
  return out;
}
const span = (a: number, b: number) => Array.from({ length: b - a }, (_, k) => a + k);

/** A 40-nucleotide cyclic `ring` and a linear `line`; `ori` covers [36, 40) and [0, 4), written as [36, 4). */
function plasmid(view: (b: ReturnType<typeof createSeqViewSpec>) => void) {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'ring', value: 'ACGTACGTTA'.repeat(4), alphabet: 'dna', cyclic: true });
  b.data.sequence({ ref: 'line', value: 'ACGTACGTTA'.repeat(4), alphabet: 'dna' });
  b.data.annotation({ ref: 'ori', geometry: 'ranges', target: 'ring', value: { start: [36], stop: [4] } });
  view(b);
  return b.getDocument();
}

describe('cyclic coordinate systems in the format (CY1)', () => {
  it('a range across the origin is one piece of two segments; origin starts a whole piece elsewhere, either way', async () => {
    const doc = plasmid((b) => {
      const v = b.view({ ref: 'v', label: 'V' });
      v.block({ sequence: 'ring', range: [36, 4] });
      v.block({ sequence: 'ring', origin: 30 });
      v.block({ sequence: 'ring', origin: 30, direction: -1 });
      v.block({ sequence: 'ring' });
      v.track().letters();
    });
    const { view, set } = await resolve(doc, 'v');
    expect(set.molecules[0]?.cyclic).toBe(true);
    const b = view.blocks;
    expect(Array.from({ length: 4 }, (_, k) => [b.start[k], b.stop[k], b.first[k], b.layout[k].count])).toEqual([
      [36, 4, 36, 2],
      [30, 30, 30, 2],
      [30, 30, 29, 2],
      [0, 40, 0, 1],
    ]);
    expect(shownOrder(view, 0)).toEqual([...span(36, 40), ...span(0, 4)]);
    expect(shownOrder(view, 1)).toEqual([...span(30, 40), ...span(0, 30)]);
    expect(shownOrder(view, 2)).toEqual([...span(0, 30).reverse(), ...span(30, 40).reverse()]);
    // Copies of one unit merge only with the same piece: four pieces, four blocks.
    expect(b.start.length).toBe(4);
  });

  it("a positions filter and a focus across the origin; an annotation's focus is its shortest arc", async () => {
    const doc = plasmid((b) => {
      const v = b.view({ ref: 'v', label: 'V', positions: 'ori' });
      v.block({ sequence: 'ring', origin: 30 });
      v.track().letters();
      const f = b.view({ ref: 'f', label: 'F' });
      f.block({ sequence: 'ring' });
      f.track().letters();
      f.focus({ annotation: 'ori' });
      const t = b.view({ ref: 't', label: 'T' });
      t.block({ sequence: 'ring', range: [30, 10] });
      t.track().letters();
      t.focus({ target: 'ring', start: 38, stop: 2 });
    });
    const v = (await resolve(doc, 'v')).view;
    const shown = v.blocks.positions[0];
    expect(shown && Array.from({ length: OrderedSet.size(shown) }, (_, k) => OrderedSet.getAt(shown, k))).toEqual([
      ...span(0, 4),
      ...span(36, 40),
    ]);
    expect((await resolve(doc, 'f')).view.sections.focus[0]).toEqual({ block: 0, start: 36, stop: 4 });
    expect((await resolve(doc, 't')).view.sections.focus[0]).toEqual({ block: 0, start: 38, stop: 2 });
  });

  it('start > stop or an origin on a linear system is svs.position.range, saying what wraps', async () => {
    const doc = plasmid((b) => {
      const v = b.view({ ref: 'v', label: 'V' });
      v.block({ sequence: 'line', range: [36, 4] });
      v.block({ sequence: 'line', origin: 5 });
      v.block({ sequence: 'ring', origin: 45 });
      v.track().letters();
      v.focus({ target: 'line', start: 38, stop: 2 });
    });
    const { view, issues } = await resolve(doc, 'v');
    expect(issues.map((i) => [i.code, i.message.replace(/.*(wraps only|origin \d+ is [^:]*).*/, '$1')])).toEqual([
      ['svs.position.range', 'wraps only'],
      ['svs.position.range', "origin 5 is on a cyclic sequence's whole unit only"],
      ['svs.position.range', 'origin 45 is within [0, 40) only'],
      ['svs.focus.range', expect.any(String)],
    ]);
    // The line's range left out; the origins ignored.
    expect(Array.from(view.blocks.start)).toEqual([0, 0]);
  });

  it('the plasmid fixture: a TSV row across the origin is two runs of one item', async () => {
    const b = createSeqViewSpec();
    const p = b.data.sequence({
      ref: 'p',
      alphabet: 'dna',
      cyclic: true,
      source: { host_file: 'pDemo.fasta', format: 'fasta' },
    });
    const f = b.data.annotation({
      ref: 'f',
      geometry: 'ranges',
      target: p,
      source: { host_file: 'pDemo.tsv', format: 'tsv' },
    });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: p, origin: 200 }).track().ranges({ data: f });
    const { set } = await resolve(b.getDocument(), 'v');
    const ann = set.annotations[0];
    expect([set.molecules[0]?.length, ann?.count, ann?.itemCount]).toEqual([240, 5, 4]);
    expect([ann?.start[0], ann?.stop[0], ann?.start[1], ann?.stop[1], ann?.item[0], ann?.item[1]]).toEqual([
      200, 240, 0, 30, 0, 0,
    ]);
  });

  it('the validator: origin with a range, or on a complex or an alignment, or negative', () => {
    const doc = plasmid((b) => {
      const aln = b.data.alignment({ ref: 'aln', value: ['AC', 'AC'] });
      b.view({ ref: 'v', label: 'V' })
        .block({ sequence: 'ring', range: [36, 4], origin: 3 })
        .track()
        .letters();
      b.view({ ref: 'w', label: 'W' }).block({ alignment: aln, origin: 1 }).track().rows();
      b.view({ ref: 'x', label: 'X' }).block({ sequence: 'ring', origin: -2 }).track().letters();
    });
    expect(validate(doc, reg).map((i) => [i.code, i.path.replace(/.*params\//, ''), i.message])).toEqual([
      ['svs.params.exclusive', 'origin', expect.stringMatching(/a range starts where it starts/)],
      ['svs.params.exclusive', 'origin', expect.stringMatching(/an alignment's columns have none/)],
      ['svs.position.range', 'origin', expect.stringMatching(/negative/)],
    ]);
  });
});
