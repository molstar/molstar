/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `nodeRefs` (seqspec:format#nodes-and-refs): a node's refs as validation reads them, the one source of a driver's
// edges (VS-24), over `indexDoc`'s nodes (implied refs filled, VS-23).
import { Q, type SeqQLQuery } from '@molstar/seq-core/query';
import { describe, expect, it } from 'vitest';
import { indexDoc, nodeRefs, type SeqViewSpecDoc, validate } from '../../src/format/index.js';
import { TEMPLATE_NODE } from '../../src/format/schema/field-schema.js';
import { CASES } from './cases.js';
import { registries } from './load-fixtures.js';

const reg = registries();
const heads = (doc: SeqViewSpecDoc, ref: string) => {
  const n = indexDoc(doc).nodes.get(ref);
  if (!n) throw new Error(`no node ${ref}`);
  return [...new Set(nodeRefs(n, reg).map((u) => u.ref.split('/')[0]))].sort();
};

function withQueries(queries: Record<string, SeqQLQuery>): SeqViewSpecDoc {
  const doc = CASES.case1();
  const data = (doc.root.children ?? []).find((c) => c.kind === 'data') as { children: unknown[] };
  for (const [ref, query] of Object.entries(queries))
    data.children.push({ kind: 'computed', ref, params: { provider: 'query', query } });
  return doc;
}

describe('nodeRefs', () => {
  it("gives a query's names, sel's expanded, and a theme's data: what the hand-kept lists missed", () => {
    const doc = withQueries({ cys: Q.read('(and p53 (letters C))'), big: Q.read('(or (sel "p53[1-5]") regions)') });
    expect(heads(doc, 'cys')).toEqual(['p53']);
    expect(heads(doc, 'big')).toEqual(['p53', 'regions']);
    const members = indexDoc(CASES.case3()).order.find((n) => n.kind === 'layer' && n.params.type === 'rows');
    expect(members && nodeRefs(members, reg).map((u) => u.path.slice(u.path.indexOf('/params')))).toContain(
      '/params/color/params/data',
    );
  });

  it('over every case, names nodes of the document, the implied root and data included in the index', () => {
    let count = 0;
    for (const [name, make] of Object.entries(CASES)) {
      const doc = make();
      expect(
        validate(doc, reg).filter((i) => i.severity === 'error'),
        name,
      ).toEqual([]);
      const index = indexDoc(doc);
      expect([index.nodes.has('root'), index.nodes.has('data')], name).toEqual([true, true]);
      for (const n of index.order) {
        for (const u of nodeRefs(n, reg)) {
          count++;
          if (u.ref !== TEMPLATE_NODE) expect(index.nodes.has(u.ref.split('/')[0]), `${name} ${u.path}`).toBe(true);
          expect(u.from, name).toBe(n.ref);
        }
      }
    }
    expect(count).toBeGreaterThan(100);
  });

  it("answers from validation's memo: the same params object, the same list", () => {
    const doc = CASES.case3();
    const n = indexDoc(doc).order.find((x) => x.kind === 'layer' && x.params.type === 'rows');
    if (!n) throw new Error('no rows layer');
    expect(nodeRefs(n, reg)).toBe(nodeRefs(n, reg));
  });
});
