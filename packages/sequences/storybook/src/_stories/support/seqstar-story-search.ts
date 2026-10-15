/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Seq*/Query and search/Search stories' search (an application's, not the viewer's): the names of the items of p53's
// annotations through a `TextIndex`, matched case-insensitively anywhere in the name from two characters, and then the
// letters of the sequence through an exact matcher, from one letter up, or, when the text has a regular expression's
// characters (`[KR]{2}.{10,12}[KR]{3,}`), through the regular-expression matcher. Name hits come first, then letter
// hits in sequence order. The answer is capped at `MAX_HITS` (the status says `n+`); the kept result sets the Search part
// makes (`<section>.keep.<query>`) are not searched.

import { exactMatcher, regexMatcher, type TextIndex, textIndex, toCodes } from '@molstar/seq-core/algorithms/search';
import { Annotations, type ItemLoci, KeyKind, OrderedSet, type SeqSet } from '@molstar/seq-core/model';
import type { Result, SearchFn, SeqViewerContext } from '@molstar/seq-viewer';

/** What the search reads of a viewer context: the assembled set, and a name and range as a loci. */
export interface SearchContext {
  readonly graph: { readonly assembled: { readonly value: { readonly set: SeqSet } | null | undefined } };
  readonly refs: { toLoci: SeqViewerContext['refs']['toLoci'] };
}

/** Every item of the set's annotations by its label: the index, and per id its annotation, item and name. */
function labelIndex(set: SeqSet) {
  const index: TextIndex = textIndex();
  const items: { annotation: number; item: number; name: string }[] = [];
  const texts: string[] = [];
  set.annotations.forEach((a, annotation) => {
    if (!a || kept(set.strings.get(a.key))) return;
    const rows = Annotations.itemRows(a);
    for (let item = 0; item < a.itemCount; item++) {
      if (rows.start[item] === rows.start[item + 1]) continue;
      const label = Annotations.rowValue(a, set, rows.rows[rows.start[item]], 'label');
      if (typeof label !== 'string') continue;
      items.push({ annotation, item, name: `${label} (${set.strings.get(a.label) || set.strings.get(a.key)})` });
      texts.push(label);
    }
  });
  index.add(
    Int32Array.from(texts, (_, k) => k),
    texts,
    0,
    texts.length,
  );
  return { index, items };
}

/** What the search returns at most: the viewer lights thousands as readily, but a person reads none of the tail. */
export const MAX_HITS = 2000;
/** The shortest text that searches names ("A" would match every name with an a); letters search from one. */
export const MIN_NAME_CHARS = 2;
/** The sets the Search part keeps are results, not data to search (`results.save`'s refs, `<section>.keep.<query>`). */
const kept = (key: string) => key.includes('.keep.');

/** A text with these is read as a regular expression over the letters; the longest hit it finds is `MOTIF_LENGTH`. */
const REGEX_CHARS = /^[A-Za-z0-9[\]^{},.()|?*+-]+$/;
const MOTIF_LENGTH = 64;

const escaped = (s: string) => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

/** The story's search: item names containing the text (any case, anywhere; two characters or more), then the letters of p53. */
export function storySearch(ctx: SearchContext): SearchFn {
  let built: { set: SeqSet; index: TextIndex; items: ReturnType<typeof labelIndex>['items'] } | null = null;
  return async (given) => {
    const text = given.trim();
    const set = ctx.graph.assembled.value?.set;
    if (!set || text.length === 0) return { results: [] };
    if (built?.set !== set) built = { set, ...labelIndex(set) };
    const results: Result[] = [];
    const labels: string[] = [];
    const { items } = built;
    // The index normalizes to lower case: the expression is written for that.
    const full = () => results.length >= MAX_HITS;
    if (text.length >= MIN_NAME_CHARS)
      built.index.matches(new RegExp(escaped(text.toLowerCase())), (id) => {
        if (full()) return;
        const it = items[id];
        const elements = [{ annotation: it.annotation as never, items: OrderedSet.ofSingleton(it.item) as never }];
        results.push({ kind: 'item-loci', set, elements } as ItemLoci);
        labels.push(it.name);
      });
    const t = { kind: -1, id: -1, row: -1 };
    const codes = set.byName('p53', t) && t.kind === KeyKind.molecule ? set.molecules[t.id]?.code : null;
    if (codes && /^[A-Za-z]+$/.test(text)) {
      exactMatcher(toCodes(text.toUpperCase())).find(codes, 0, codes.length, false, (start, stop) => {
        if (full()) return;
        const loci = ctx.refs.toLoci('p53', { start, stop });
        if (!loci) return;
        results.push(loci);
        labels.push(`${text.toUpperCase()} at ${start + 1}`);
      });
    } else if (codes && REGEX_CHARS.test(text)) {
      // A motif as a regular expression over the letters (PROSITE's patterns are written this way once translated); a
      // text that is not one searches names only.
      let motif: ReturnType<typeof regexMatcher> | null = null;
      try {
        motif = regexMatcher(new RegExp(text, 'i'), MOTIF_LENGTH);
      } catch {
        motif = null;
      }
      motif?.find(codes, 0, codes.length, false, (start, stop) => {
        if (full()) return;
        const loci = ctx.refs.toLoci('p53', { start, stop });
        if (!loci) return;
        results.push(loci);
        labels.push(`${String.fromCharCode(...codes.subarray(start, stop))} at ${start + 1}`);
      });
    }
    return { results, labels, more: full() };
  };
}
