/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories' own search (lib/seqstar/storybook/src/_stories/support/seqstar-story-search.ts) over case 18's Tracks view: item
// names match case-insensitively anywhere in the name, letters match from one letter up, names first.
// The suggested searches each find what the view shows.

import { compile } from '@molstar/seq-core/format';
import type { SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { fixtureLoaders, P53 } from '../../core/_test/format/load-fixtures.js';
import { compiledCase, reg } from '../src/_stories/support/seqstar-cases.js';
import { type SearchContext, storySearch } from '../src/_stories/support/seqstar-story-search.js';

async function search(text: string) {
  const { set } = await compiledCase('case18');
  return searchIn(set, text);
}
async function searchIn(set: SeqSet, text: string) {
  const ctx: SearchContext = {
    graph: { assembled: { value: { set } } },
    refs: { toLoci: (name, range) => ({ kind: 'test', name, ...range }) as never },
  };
  const r = await storySearch(ctx)(text, new AbortController().signal);
  return r.labels ?? [];
}

describe("the story's search", () => {
  it('matches item names anywhere in the name, in any case', async () => {
    // Case 18's query results and its bookmarks are items with names too.
    const binding = ['DNA binding (Domains)', 'DNA binding (Cysteine-rich domains)', 'MDM2 binding (Bookmarks)'];
    expect(await search('binding')).toEqual(binding);
    expect(await search('BIND')).toEqual(binding);
    expect(await search('zinc')).toHaveLength(4);
    expect(await search('tetra')).toEqual(['Tetramerization (Domains)']);
    expect((await search('contact')).length).toBe(8);
    expect(await search('dna bind')).toEqual(binding.slice(0, 2));
  });

  it('matches letters from two letters up, and from one', async () => {
    // p53 reads D, A at 352 to 353 (1-based), and E, L, K from 349.
    expect(await search('DA')).toContain('DA at 352');
    expect(await search('ELK')).toContain('ELK at 349');
    const one = await search('W');
    expect(one.length).toBeGreaterThan(0);
    expect(one.every((l) => l.startsWith('W at '))).toBe(true);
  });

  it('lists name hits first, then letter hits, and finds what the suggestions promise', async () => {
    const dna = await search('DNA');
    // "DNA binding" (a domain, and a rich one) and eight "DNA contact" residues; no p53 letters D-N-A in a row.
    expect(dna.filter((l) => l.includes('(')).length).toBe(10);
    const ti = await search('TI');
    const firstLetter = ti.findIndex((l) => l.startsWith('TI at '));
    expect(firstLetter).toBeGreaterThan(0);
    expect(ti.slice(0, firstLetter).every((l) => l.includes('('))).toBe(true);
    expect(ti.slice(firstLetter).every((l) => l.startsWith('TI at '))).toBe(true);
    expect(await search('RHSVV')).toEqual(['RHSVV at 213']);
    expect(await search('')).toEqual([]);
  });

  it("reads a text with a regular expression's characters as a motif over the letters", async () => {
    // PROSITE's bipartite nuclear localization signal: the Motif feature's 305 to 321.
    expect(await search('[KR]{2}.{10,12}[KR]{3,}')).toEqual(['KRALPNNTSSSPQPKKK at 305']);
    expect(await search('N[^P][ST]')).toEqual(['NSS at 239', 'NNT at 310', 'NTS at 311']);
    // Not a regular expression: no hits, no error.
    expect(await search('[KR')).toEqual([]);
  });

  it('matches names from two characters; one letter searches the letters only, once per residue', async () => {
    const a = await search('A');
    expect(a.every((l) => l.startsWith('A at '))).toBe(true);
    expect(a.length).toBe([...P53.sequence].filter((c) => c === 'A').length);
    expect(new Set(a).size).toBe(a.length);
    expect((await search('a')).length).toBe(a.length);
  });

  it('does not search the sets it keeps, so a kept query does not multiply the next one', async () => {
    const { doc } = await compiledCase('case18');
    const withKept = structuredClone(doc);
    const data = withKept.root.children?.find((c) => c.kind === 'data');
    // What `results.save` writes for the query "A": a ranges set under the section's name, one item per hit.
    data?.children?.push({
      kind: 'annotation',
      ref: 'tracks.main.keep.a',
      params: {
        geometry: 'ranges',
        label: 'A',
        value: [118, 128, 138].map((start, k) => ({
          target: 'p53/unit',
          start,
          stop: start + 1,
          id: `r${k + 1}`,
          label: `A at ${start + 1}`,
        })),
      },
    });
    const { set, issues } = await compile(withKept, reg, fixtureLoaders());
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    expect((await searchIn(set, 'A')).length).toBe((await search('A')).length);
    expect(await searchIn(set, 'at 119')).toEqual([]);
  });
});
