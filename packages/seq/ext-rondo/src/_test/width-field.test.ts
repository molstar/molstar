/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A width by a field (_spec/README.md, Kinds of links): the ribbons, stubs and necks are wide by the field's sum, but
// folding and the cards count links: a stub folds by its links, never by its `q`, and every card's count is links.

import { describe, expect, it } from 'vitest';
import { Kind, type RondoLay } from '../lay.js';
import { describeRondo } from '../words.js';
import { lay, ringSpec, Size, settle, toyChromosomes, toyLinks } from './helpers.js';

const chroms = toyChromosomes();
const links = toyLinks(chroms, 1500, 4);

/** Rows of chr1 toward each partner chromosome, by stem: what the stubs count. */
function stubRows(l: RondoLay) {
  const d = l.read!.data!;
  const out = new Map<string, number>();
  for (let k = 0; k < d.n; k++) {
    if (!l.read!.shown[k]) continue;
    const a = l.names[d.aUnit![k]];
    const b = l.names[d.bUnit![k]];
    if ((a === 'chr1') === (b === 'chr1')) continue;
    const [pos, partner] = a === 'chr1' ? [d.aStart[k], b] : [d.bStart[k], a];
    const key = `${Math.floor(pos / l.stemLevel)}|${partner}`;
    out.set(key, (out.get(key) ?? 0) + 1);
  }
  return out;
}

describe('a width by a field: folding and the cards count links', () => {
  it('chr1 alone, wide by q: stubs fold by their links, cards and labels say links', async () => {
    const s = await settle(ringSpec({ held: true, links, include: ['chr1'], layer: { width: 'q' } }));
    const l = lay(await s.frame({ size: Size }))!;
    expect(l.width).toBe('q');
    expect(l.nb).toBeGreaterThan(0);
    const rows = stubRows(l);
    // a stub a partner holds its links (folded by links, never for a small q)
    for (let b = 0; b < l.nb; b++) {
      const partner = l.stubPartner[b];
      if (partner < 0) continue;
      const want = rows.get(`${l.stemBin[l.stubStem[b]]}|${l.names[partner]}`);
      expect(l.stubLinks[b]).toBe(want);
      // a small stub stays only alone on its stem (two or more fold)
      if (want! < 2)
        for (let c = 0; c < l.nb; c++)
          if (c !== b && l.stubStem[c] === l.stubStem[b] && l.stubPartner[c] >= 0)
            expect(l.stubLinks[c]).toBeGreaterThanOrEqual(2);
    }
    for (const [key, n] of rows)
      if (n >= 2) {
        const [bin, partner] = key.split('|');
        const kept = Array.from({ length: l.nb }, (_, b) => b).some(
          (b) => l.stemBin[l.stubStem[b]] === Number(bin) && l.names[l.stubPartner[b]] === partner,
        );
        expect(kept, key).toBe(true);
      }
    // every card's and label's count is a whole number of links, at least one
    for (let i = 0; i < l.polys.n; i++) {
      if (l.polyKind[i] !== Kind.stub && l.polyKind[i] !== Kind.trumpet) continue;
      const lines = describeRondo(l.hits[i] as never, l)!;
      const m = lines.join('\n').match(/Links: (\d+) links?/);
      expect(Number(m![1]), lines.join(' | ')).toBeGreaterThanOrEqual(1);
      for (const p of lines.join('\n').matchAll(/chr\d+ \((\d+)\)/g)) expect(Number(p[1])).toBeGreaterThanOrEqual(1);
    }
    for (let b = 0; b < l.nb; b++) if (l.stubLabel[b]) expect(l.stubLabel[b]).toMatch(/ · [1-9]\d*$/);
  });
});
