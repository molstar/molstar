/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's `highlights` (seqspec:format#sections, LS1): the param checked with the document (its keys, its data
// refs, its style), and resolved per section as item loci on the set, a query node's as any annotation's.

import { OrderedSet } from '@molstar/seq-core/model';
import { Q } from '@molstar/seq-core/query';
import { describe, expect, it } from 'vitest';
import {
  compile,
  resolveView,
  type SeqViewSpecDoc,
  sectionHighlights,
  supportedAll,
  validate,
} from '../../src/format/index.js';
import { CASES } from './cases.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
type J = { kind: string; ref?: string; params?: Record<string, unknown>; children?: J[] };
/** Case 1 with a query node and the Tracks view's section given `highlights`. */
function withHighlights(highlights: unknown): SeqViewSpecDoc {
  const doc = CASES.case1();
  const kids = (doc.root.children ?? []) as J[];
  const data = kids.find((c) => c.kind === 'data') as J;
  data.children?.push({
    kind: 'computed',
    ref: 'cysteines',
    params: { provider: 'query', query: Q.read('(and p53 (letters C))') },
  });
  const section = kids.find((c) => c.ref === 'tracks')?.children?.[0] as J;
  section.params = { ...section.params, highlights };
  return doc;
}
const said = (doc: SeqViewSpecDoc) => validate(doc, reg).map((i) => `${i.code} ${i.path.replace(/.*\/params/, '')}`);

describe('section.highlights', () => {
  it('is checked with the document: its keys, its refs and its style', () => {
    expect(said(withHighlights([{ data: 'regions' }, { data: 'cysteines', items: null, style: 'accent' }]))).toEqual(
      [],
    );
    expect(
      said(
        withHighlights([
          { data: 'nope' },
          { data: 'regions', style: 'glow' },
          { data: 'p53' },
          { data: 'regions', color: 'red' },
        ]),
      ),
    ).toEqual([
      'svs.params.type /highlights/1/style',
      'svs.params.unknown /highlights/3/color',
      'svs.ref.unresolved /highlights/0/data',
      'svs.ref.wrong_kind /highlights/2/data',
    ]);
  });

  it("resolves per section into item loci: every item, the ones named, a query node's, none of an unavailable node", async () => {
    const doc = withHighlights([
      { data: 'regions' },
      { data: 'cysteines', items: ['nope', 'q0'], style: 'highlight' },
      { data: 'cysteines', style: 'accent' },
      { data: 'regions', items: ['nope'] },
    ]);
    const { set, issues } = await compile(doc, reg, fixtureLoaders());
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    const { view } = resolveView(doc, set, 'tracks', reg, supportedAll(reg));
    const hs = sectionHighlights(view, set, 0);
    expect(hs.map((h) => [h.data, h.items, h.style])).toEqual([
      ['regions', null, 'band'],
      ['cysteines', ['nope', 'q0'], 'highlight'],
      ['cysteines', null, 'accent'],
      ['regions', ['nope'], 'band'],
    ]);
    const items = (k: number) => {
      const e = hs[k].loci?.elements[0];
      return e ? Array.from({ length: OrderedSet.size(e.items) }, (_, i) => OrderedSet.getAt(e.items, i)) : null;
    };
    expect(items(0)).toEqual([0, 1, 2]);
    // A query node's one group is its item `q0` (seqspec:format#placement).
    expect([items(1), items(2)]).toEqual([[0], [0]]);
    expect(hs[1].loci?.elements[0].annotation).not.toBe(hs[0].loci?.elements[0].annotation);
    expect(hs[3].loci).toBeNull();
    // The other view's section has none.
    const text = resolveView(doc, set, 'text', reg, supportedAll(reg)).view;
    expect(sectionHighlights(text, set, 0)).toEqual([]);
  });
});
