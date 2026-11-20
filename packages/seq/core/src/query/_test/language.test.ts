/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The language's laws, seeded at small fixed counts on the 3SGJ glycan fixture: items as params, columns and records
// one table (SoA ≡ AoS), the long names read as the short, printed short and read back equal as a value; `elements`
// equal to its `filter` sugar for random items; monadic composition (a form over a union of prefixes is the union of
// the form over each); types checked when read; and the index walk's representative: rows plus matches, never rows ×
// units.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { type Atom, type Node, parseAddress, printAddress } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import {
  evaluateQuery,
  location,
  locationsEqual,
  type QueryModel,
  type QueryUnit,
  unionLocations,
} from '@molstar/seqx-core/query';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { elementsForm, type Items, itemsEqual, itemsOf, itemsPredicate, readItems } from '../items.js';
import { BuiltInPropertyTable } from '../properties.js';
import { codes, count, ex, glycanSpec, loc, prop, specRun } from './structures.js';

const T = BuiltInPropertyTable;
const g = specRun(glycanSpec());
const D = '/data/c/D';

const formOf = (text: string): Node => {
  const s = parseAddress(`/${text}`).pattern?.segments[0];
  if (!s || typeof s === 'string') throw new Error(`not a form: ${text}`);
  return s;
};
const items = (text: string) => readItems(T, formOf(text));
const printed = (i: Items) => printAddress([elementsForm(i)]).slice(1);

describe('items: one table however written', () => {
  it('params, columns and records read alike; JSON’s object and array of objects too (SoA ≡ AoS)', () => {
    const cols = items('(elements :chain [A C] :auth-seq 87)');
    const recs = items('(elements :items [(fields :chain A :auth-seq 87) (fields :chain C :auth-seq 87)])');
    expect(itemsEqual(cols, recs)).toBe(true);
    expect(itemsEqual(cols, itemsOf(T, { chain: ['A', 'C'], 'auth-seq': 87 }))).toBe(true);
    expect(
      itemsEqual(
        cols,
        itemsOf(T, [
          { chain: 'A', auth_seq_id: 87 },
          { label_asym_id: 'C', 'auth-seq': 87 },
        ]),
      ),
    ).toBe(true);
    expect(cols.rows).toBe(2);
    // a row given twice is one; nil is no test on that field
    expect(items('(elements :chain [A A])').rows).toBe(1);
    expect(items('(elements :chain [A nil] :seq [1 2])').columns.chain).toEqual(['A', null]);
  });

  it('long names read, short printed; the shorter layout; parse(print(a)) equal as a value', () => {
    expect(printed(items('(elements :label_asym_id A :label_seq_id 33)'))).toBe('(elements :chain A :seq 33)');
    expect(printed(items('(elements :chain [A C] :auth-seq 87)'))).toBe('(elements :auth-seq 87 :chain [A C])');
    const sparse =
      '(elements :items [(fields :atom CA :chain A :ins B :seq 3) (fields :auth-chain C :auth-seq 7 :entity 2 :operator X)])';
    // an entity id is a word: printed quoted when it reads as a number
    expect(printed(items(sparse))).toBe(sparse.replace(':entity 2', ':entity "2"'));
    const mixed = items('(elements :items [(fields :chain A :seq 33) (fields :chain B :beg 140 :end 146)])');
    expect(printed(mixed)).toBe('(elements :beg [nil 140] :chain [A B] :end [nil 146] :seq [33 nil])');
    const rng = new PCG(41);
    const pick = <T>(xs: readonly T[]) => xs[rng.int() % xs.length];
    for (let trial = 0; trial < 40; trial++) {
      const rows = 1 + (rng.int() % 4);
      const recs = Array.from({ length: rows }, () => {
        const r: Record<string, Atom> = {};
        if (rng.float() < 0.7) r[pick(['chain', 'label_asym_id'])] = pick(['A', 'B', '1']);
        if (rng.float() < 0.5) r.seq = 1 + (rng.int() % 9);
        else if (rng.float() < 0.5) [r.beg, r.end] = [1 + (rng.int() % 3), 4 + (rng.int() % 3)];
        if (rng.float() < 0.5) r[pick(['atom', 'label_atom_id'])] = pick(['C1', 'CA', 'FE']);
        if (rng.float() < 0.3) r.standard = rng.float() < 0.5;
        return r;
      });
      const a = itemsOf(T, recs);
      const back = items(printed(a));
      expect(itemsEqual(a, back), printed(a)).toBe(true);
    }
  });

  it('bounds: unequal columns, a record that is not one, seq with beg, beg past end, a field twice', () => {
    const err = (text: string) => {
      try {
        items(text);
        return 'none';
      } catch (e) {
        return (e as { issues: { code: string }[] }).issues[0].code;
      }
    };
    expect(err('(elements :chain [A B] :seq [1 2 3])')).toBe('seq.query.form.param');
    expect(err('(elements :items [(range :from 1)])')).toBe('seq.query.form.param');
    expect(err('(elements :seq 3 :beg 1)')).toBe('seq.query.form.param');
    expect(err('(elements :beg 9 :end 3)')).toBe('seq.query.form.param');
    expect(err('(elements :items [(fields :chain A :label_asym_id A)])')).toBe('seq.query.form.param');
    expect(err('(elements :seq A)')).toBe('seq.query.form.type');
    expect(err('(elements :"af:plddt" 3)')).toBe('seq.query.extension.missing');
  });
});

describe('items: an empty selection, and the rows bounded first', () => {
  it('no row selects nothing and prints so; print∘parse = id over random selections, empty ones included', () => {
    const none = itemsOf(T, { chain: [] });
    expect(none.rows).toBe(0);
    expect(printed(none)).toBe('(elements :items [])');
    expect(count(g.run(`/data/c/${printed(none)}`))).toBe(0);
    expect(count(g.run('/data/c/(elements)'))).toBe(17);
    const rng = new PCG(53);
    for (let trial = 0; trial < 60; trial++) {
      const n = rng.int() % 4;
      const cols: Record<string, Atom[]> = {};
      if (rng.float() < 0.7) cols.chain = Array.from({ length: n }, () => (rng.float() < 0.5 ? 'A' : 'D'));
      if (rng.float() < 0.5) cols.seq = Array.from({ length: n }, () => 1 + (rng.int() % 9));
      if (rng.float() < 0.4) cols.atom = Array.from({ length: n }, () => (rng.float() < 0.3 ? null : 'C1'));
      const a = itemsOf(T, cols.chain || cols.seq || cols.atom ? cols : { chain: [] });
      const text = printed(a);
      const back = items(text);
      expect(itemsEqual(back, a), text).toBe(true);
      expect(printed(back)).toBe(text);
      // and what it selects is what the table selects
      if (!a.rows) expect(count(g.run(`/data/c/${text}`)), text).toBe(0);
    }
  });

  it('more rows than the limit: refused before they are read, a column’s length or the records', () => {
    const over = (f: () => unknown) => {
      try {
        f();
        return 'none';
      } catch (e) {
        return (e as { issues: { code: string }[] }).issues[0].code;
      }
    };
    expect(over(() => items('(elements :seq [1 2 3 4 5])'))).toBe('none');
    expect(over(() => readItems(T, formOf('(elements :seq [1 2 3 4 5])'), 3))).toBe('repr.tree.address.limit_exceeded');
    expect(over(() => readItems(T, formOf('(elements :items [(fields :seq 1) (fields :seq 2)])'), 1))).toBe(
      'repr.tree.address.limit_exceeded',
    );
    expect(over(() => itemsOf(T, { seq: Array.from({ length: 10 }, (_, i) => i + 1) }, 9))).toBe(
      'repr.tree.address.limit_exceeded',
    );
    const seqs = Array.from({ length: 50 }, (_, i) => i + 1).join(' ');
    const r = evaluateQuery(g.spec, g.query, `/data/c/(elements :seq [${seqs}])`, { limit: 20 });
    expect(codes(r)).toEqual(['repr.tree.address.limit_exceeded']);
  });
});

describe('elements is its filter (a property)', () => {
  it('for random items, the walk equals the filter over the same elements', () => {
    const rng = new PCG(7);
    const pick = <T>(xs: readonly T[]) => xs[rng.int() % xs.length];
    for (let trial = 0; trial < 30; trial++) {
      const atomLevel = rng.float() < 0.4;
      const recs = Array.from({ length: 1 + (rng.int() % 3) }, () => {
        const r: Record<string, Atom> = {};
        if (rng.float() < 0.6) r.chain = pick(['A', 'D']);
        if (rng.float() < 0.4) r.seq = 1 + (rng.int() % 9);
        else if (rng.float() < 0.4) [r.beg, r.end] = [1 + (rng.int() % 4), 3 + (rng.int() % 5)];
        if (rng.float() < 0.4) r.component = pick(['NAG', 'MAN', 'ASN', 'FUC']);
        if (atomLevel) r.atom = pick(['C1', 'C6', 'O5']);
        return r;
      });
      const it = itemsOf(T, recs);
      const walk = g.run(`/data/c/${printed(it)}`);
      // the same elements: every residue of the complex's units, or every atom
      const prefix = atomLevel ? `${D}/(atoms)/(all)` : '/data/c/(elements)';
      const filtered = g.form(prefix, { kind: 'filter', params: { test: itemsPredicate(it) } });
      const a = walk.result?.type === 'location' ? walk.result.location : undefined;
      const b = filtered.result?.type === 'location' ? filtered.result.location : undefined;
      expect(a && b && locationsEqual(a, b), printed(it)).toBe(true);
    }
  });
});

describe('monadic composition', () => {
  it('a form over a union of prefixes is the union of the form over each', () => {
    const rng = new PCG(19);
    for (let trial = 0; trial < 20; trial++) {
      const r = () => {
        const a = 1 + (rng.int() % 8);
        return `${D}/(range :from ${a} :to ${a + (rng.int() % 3)})`;
      };
      const [p, q] = [r(), r()];
      for (const form of ['(elements :atom C1)', '(elements :component NAG)', '(expand :to unit)']) {
        const both = g.run(`/(merge :of [${p} ${q}])/${form}`);
        const each = [g.run(`${p}/${form}`), g.run(`${q}/${form}`)];
        const loc2 = (e: (typeof each)[0]) => (e.result?.type === 'location' ? e.result.location : location({}));
        expect(locationsEqual(loc2(both), unionLocations(loc2(each[0]), loc2(each[1]))), `${p} ${q} ${form}`).toBe(
          true,
        );
      }
    }
  });

  it('an element lacking a property yields nothing, the rest stays: one warning', () => {
    const r = g.run('/data/c/(elements :atom-index 0)');
    expect([count(r), codes(r)]).toEqual([1, ['seq.query.property.absent']]);
    // a polymer's atoms read by its residues' standard components: none has a C1, nothing lacks
    expect(codes(g.run('/data/c/(elements :atom C1)'))).toEqual([]);
  });

  it('combinators over one type: except of a location, one operand its complement within the prefix', () => {
    expect(loc(g.run(`${D}/(except :of (elements :component NAG))`), D)).toEqual([2, 4, 5, 6, 7, 8]);
    expect(loc(g.run(`${D}/(except :of [(elements) (elements :component NAG)])`), D)).toEqual([2, 4, 5, 6, 7, 8]);
    expect(codes(g.run(`${D}/(merge :of [(elements) (atoms)])`))).toEqual(['seq.query.form.type']);
  });
});

describe('types checked when read', () => {
  const f = (test: Node) => codes(g.form(`${D}/(elements)`, { kind: 'filter', params: { test } }));
  it('an operator over a word and a number, a test that is not boolean, logic over a value', () => {
    expect(f(ex('<', prop('chain'), 3))).toEqual(['seq.query.form.type']);
    expect(f(prop('seq'))).toEqual(['seq.query.form.type']);
    expect(f(ex('and', prop('seq')))).toEqual(['seq.query.form.type']);
    expect(f(ex('=', prop('seq'), prop('chain')))).toEqual(['seq.query.form.type']);
    // a literal reads as the property's type: an entity id written as a number is a word
    expect(f(ex('=', prop('entity'), 3))).toEqual(['seq.query.property.absent']);
    expect(f(ex('in-range', prop('seq'), 2, 4))).toEqual([]);
  });

  it('a property unknown: bare a param error with the names, namespaced an extension missing', () => {
    const bare = g.form(`${D}/(elements)`, { kind: 'filter', params: { test: ex('=', prop('label_comp'), 'NAG') } });
    expect(codes(bare)).toEqual(['seq.query.form.param']);
    expect(bare.issues[0].hint).toMatch(/properties: entity/);
    expect(f(ex('=', prop('af:plddt'), 1))).toEqual(['seq.query.extension.missing']);
    expect(codes(g.run(`${D}/(af:confident)`))).toEqual(['seq.query.extension.missing']);
  });
});

/** A complex of `n` units, each of one 10-residue sequence, chain `U<i>`; the units' `chain` reads counted. */
function manyUnits(n: number): { model: QueryModel; spec: Node; reads: () => number } {
  const b = createSeqSpec();
  const s = b.data.sequence('s', { alphabet: 'protein', letters: 'ACDEFGHIKL' });
  b.data.complex('c', {}, (c) => {
    for (let i = 0; i < n; i++) c.unit(`U${i}`, { of: s, chain: `U${i}` });
  });
  const spec = b.build();
  const model = resolveSpec(spec);
  let reads = 0;
  const counted = (u: QueryUnit): QueryUnit =>
    Object.defineProperty({ ...u }, 'chain', {
      get: () => {
        reads++;
        return u.chain;
      },
    });
  const units = new Map<string, readonly QueryUnit[] | undefined>();
  const proxy = new Proxy(model, {
    get(target, key) {
      if (key === 'units')
        return (c: string) => {
          if (!units.has(c)) units.set(c, target.units?.(c)?.map(counted));
          return units.get(c);
        };
      const v = Reflect.get(target, key, target);
      return typeof v === 'function' ? v.bind(target) : v;
    },
  });
  return { model: proxy, spec, reads: () => reads };
}

describe('the index walk (the quick representative)', () => {
  it('rows grouped by their unit fields, each unit looked up once: rows plus matches, never rows × units', async () => {
    const n = 200;
    const { model, spec, reads } = manyUnits(n);
    const chains = Array.from({ length: n }, (_, i) => `U${i}`).join(' ');
    const seqs = Array.from({ length: n }, (_, i) => 1 + (i % 10)).join(' ');
    const r = evaluateQuery(spec, model, `/data/c/(elements :chain [${chains}] :seq [${seqs}])`, { limit: 10000 });
    expect(r.issues).toEqual([]);
    expect(count(r)).toBe(n);
    // each unit's chain read once (a scan of rows × units would read it n × n times)
    expect(reads()).toBeLessThanOrEqual(2 * n);
  });
});

describe('the designs’ examples as text (locations-and-queries.md, properties.md, addressing.md)', () => {
  const examples = [
    '/data/4hhb/(elements :component HEM :atom FE)',
    '/data/4hhb/(elements :auth-chain [A C] :auth-seq 87)',
    '/data/1fin/(elements :chain A :beg 10 :end 20 :atom CA)',
    '/data/1fin/(elements :chain A :atom CA)/(filter :test (> (prop b-factor) 60))',
    '/data/1fin/(elements :chain A)/(filter :test (> (prop b-factor) 60))',
    '/data/af-p53/(elements)/(filter :test (>= (prop af:plddt) 70))',
    '/data/af-p53/(elements)/(filter :test (>= (value /data/af-p53/plddt :field score) 70))',
    '/data/glycan/(elements :atom C1)',
    '/data/1fin/(elements :chain A :seq 33)/(linked :via /data/1fin/links)',
    '/data/1ehz/(groups :residue-test (= (prop modified) true))',
    '/data/1ehz/(elements :component MSE :element SE)',
    '/data/x/(elements :items [(fields :chain A :seq 33) (fields :chain B :beg 140 :end 146)])',
    '/data/x/(filter :test (= (prop component-class) saccharide))',
    '/data/x/(molecule :entity 3)',
    '/data/x/(elements :atom CA)/(expand :to residue)',
    '/data/x/(filter :test (and (in-range (prop seq) 10 20) (not (in (prop component) [GLY PRO]))))',
    '/data/x/(except :of [(elements :chain A) (elements :chain A :seq 33)])',
  ];
  it('each parses, prints canonically and reads back equal (keys sorted)', () => {
    for (const t of examples) {
      const p = parseAddress(t);
      expect(p.issues, t).toEqual([]);
      const printed = printAddress(p.pattern!.segments);
      expect(parseAddress(printed).pattern, t).toEqual(p.pattern);
      expect(printAddress(parseAddress(printed).pattern!.segments)).toBe(printed);
    }
    // one written in canonical order prints as written
    expect(printAddress(parseAddress(examples[4]).pattern!.segments)).toBe(examples[4]);
  });
});

/** A protein of `n` residues, its components' reads counted (`at` and `standard`). */
function longProtein(n: number): { model: QueryModel; spec: Node; reads: () => number } {
  const b = createSeqSpec();
  const amino = 'ACDEFGHIKLMNPQRSTVWY';
  b.data.sequence('p', {
    alphabet: 'protein',
    letters: Array.from({ length: n }, (_, i) => amino[(i * 7) % 20]).join(''),
  });
  const spec = b.build();
  const model = resolveSpec(spec);
  let reads = 0;
  const proxy = new Proxy(model, {
    get(target, key) {
      if (key === 'components')
        return (s: string) => {
          const c = target.components?.(s);
          if (!c) return c;
          return {
            ...c,
            at: (p: number) => {
              reads++;
              return c.at(p);
            },
            standard: (p: number) => {
              reads++;
              return c.standard?.(p);
            },
          };
        };
      const v = Reflect.get(target, key, target);
      return typeof v === 'function' ? v.bind(target) : v;
    },
  });
  return { model: proxy, spec, reads: () => reads };
}

const fieldRows = (n: number, rows: number, field: string) => {
  const amino = ['ALA', 'CYS', 'TRP', 'MET', 'GLY'];
  const values =
    field === 'index'
      ? Array.from({ length: rows }, (_, i) => (i * 97) % n)
      : field === 'atom'
        ? Array.from({ length: rows }, (_, i) => ['SD', 'SG', 'NE1', 'CA'][i % 4])
        : Array.from({ length: rows }, (_, i) => amino[i % 5]);
  return `/data/p/(elements :${field} [${values.join(' ')}])`;
};

describe('the fine fields indexed (a row naming no position)', () => {
  it('`index` named by index; a component or an atom field from one map a frame and field: rows plus matches', () => {
    const n = 2000;
    for (const field of ['index', 'component', 'atom']) {
      const { model, spec, reads } = longProtein(n);
      const r = evaluateQuery(spec, model, fieldRows(n, 300, field), { limit: 10000 });
      expect(codes(r), field).toEqual([]);
      const matches = count(r) as number;
      expect(matches, field).toBeGreaterThan(0);
      // a scan per row would read 300 × 2000 components
      expect(reads(), field).toBeLessThanOrEqual(n + 2 * matches);
    }
    // the index field reads none, and names one position each
    const { model, spec } = longProtein(n);
    expect(count(evaluateQuery(spec, model, '/data/p/(elements :index [0 5 1999 2000])', { limit: 10 }))).toBe(3);
  });

  itWithin('1,000 rows over 100,000 residues: index, component and atom fields', (within) => {
    const n = 100_000;
    for (const field of ['index', 'component', 'atom']) {
      const { model, spec } = longProtein(n);
      const query = fieldRows(n, 1000, field);
      evaluateQuery(spec, model, query, { limit: 1e6 }); // warm
      // measured: tens of ms (an M-series laptop); a scan per row took seconds
      const r = within(400, () => evaluateQuery(spec, model, query, { limit: 1e6 }));
      expect(codes(r), field).toEqual([]);
    }
  });
});
