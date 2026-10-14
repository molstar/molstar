/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL2a: the binder (seqspec:seqql#binding-and-checking): short names by operand kind and slot chain, positionals and
// keywords, coercion, macros, and its errors with path, position, signature and a suggestion.
import { describe, expect, it } from 'vitest';
import { expandMacros } from '../../src/query/bind.js';
import {
  type Apply,
  bindQuery,
  createLanguage,
  defineSymbol,
  type Expression,
  type Language,
  Q,
  type QueryIssue,
  type QueryNode,
  readQuery,
  renderIssue,
} from '../../src/query/index.js';

const BASE = createLanguage().language;

function bound(text: string, language: Language = BASE) {
  const r = readQuery(text);
  expect(r.issues).toEqual([]);
  return bindQuery(r.node as QueryNode, language);
}
const expr = (text: string, language?: Language) => {
  const b = bound(text, language);
  expect(b.issues.map(renderIssue)).toEqual([]);
  return b.expr as Expression;
};
const issues = (text: string, language?: Language) => {
  const b = bound(text, language);
  expect(b.expr).toBeNull();
  return b.issues;
};
const one = (text: string, language?: Language): QueryIssue => {
  const list = issues(text, language);
  expect(list.map(renderIssue)).toHaveLength(1);
  return list[0] as QueryIssue;
};
/** The full name at a JSON pointer of the bound expression. */
function headAt(e: Expression, pointer: string): string {
  let x: unknown = e;
  for (const part of pointer.split('/').slice(1)) x = (x as Record<string, unknown>)[part];
  return (x as Apply).head;
}

describe('resolving names', () => {
  it("binds the spec's stored example exactly", () => {
    expect(JSON.stringify(expr('(pick (items :of domains) :where (matches (label) kinase))'))).toBe(
      '{"head":"seq.filter.pick","args":{"0":{"head":"seq.generator.items","args":{"of":"domains"}},"where":{"head":"core.str.matches","args":{"0":{"head":"seq.item.label"},"1":"kinase"}}}}',
    );
  });

  it('takes a full name as written, wherever it is', () => {
    expect(expr('(seq.generator.all)')).toEqual({ head: 'seq.generator.all' });
    // a full name wins over the slot chain: the checker (QL2b) judges whether the slot fits
    expect(headAt(expr('(items :where (= (seq.system.label) x))'), '/args/where/args/0')).toBe('seq.system.label');
    expect(headAt(expr('(core.logic.and 1brs/A 1brs/B)'), '')).toBe('core.logic.and');
  });

  it('step 1: an unknown name gets did-you-mean, nearest first, with its signature', () => {
    const i = one('(pick (item :of domains) :where (matches (label) kinase))');
    expect(renderIssue(i)).toBe(
      '1:8 seqql.symbol.unknown: unknown symbol "item"; did you mean "items"?\n  (items :of A :where T)',
    );
    expect(i.path).toBe('/args/0');
    expect(i.known?.[0]).toBe('items');
    expect(one('(seq.generator.item)').message).toBe(
      'unknown symbol "seq.generator.item"; did you mean "seq.generator.items"?',
    );
    const far = one('(frobnicate)');
    expect([far.message, far.signature, far.known]).toEqual([
      'unknown symbol "frobnicate": no symbol has a similar name',
      null,
      [],
    ]);
  });

  it('step 2: and, or, not by operand kind; a bare word is a selection', () => {
    expect(headAt(expr('(and 1brs/A domains)'), '')).toBe('seq.combinator.intersect');
    expect(headAt(expr('(or 1brs/A (sel "1brs/B[1-5]"))'), '')).toBe('seq.combinator.merge');
    expect(headAt(expr('(not domains :in 1brs/A)'), '')).toBe('seq.modifier.complement');
    const logic = expr('(items :where (and (> (length) 3) (or (= (id) a) (not (= (id) b)))))');
    expect(headAt(logic, '/args/where')).toBe('core.logic.and');
    expect(headAt(logic, '/args/where/args/1/0')).toBe('core.logic.or');
    expect(headAt(logic, '/args/where/args/1/0/args/1/0')).toBe('core.logic.not');
    // nested: the inner and is a selection, so the outer one is too
    expect(headAt(expr('(and (and 1brs/A domains) (letters C))'), '')).toBe('seq.combinator.intersect');
    // no operand kind known: a selection outside a :where, a test inside one
    expect(headAt(expr('(and)'), '')).toBe('seq.combinator.intersect');
    expect(headAt(expr('(items :where (and))'), '/args/where')).toBe('core.logic.and');
  });

  it('step 2: a selection mixed with a test is an instructive type error', () => {
    const i = one('(and 1brs/A (> (length) 3))');
    expect(renderIssue(i)).toBe(
      '1:2 seqql.type: and mixes a selection and a test: filter with (pick Q :where T)\n  (and T T...) or (and Q Q...)',
    );
    expect(one('(pick (all) :where (or (> (length) 3) domains))').path).toBe('/args/where');
  });

  it("step 3: properties along the innermost binder's slot chain, the first fitting signature winning", () => {
    const pick = (t: string) => headAt(expr(`(pick (items :of domains) :where ${t})`), '/args/where/args/0');
    expect(pick('(matches (label) kinase)')).toBe('seq.item.label'); // group → item
    expect(pick('(> (length) 50)')).toBe('seq.group.length');
    expect(pick('(= (kind) unit)')).toBe('seq.system.kind');
    expect(pick('(< (field evalue) 0.1)')).toBe('seq.item.field');
    // the spec's example: the item's field has no :of, so the system's is taken
    expect(pick('(>= (field identity :of stats) 0.9)')).toBe('seq.system.field');
    const at = (t: string, p: string) => headAt(expr(t), p);
    expect(at('(positions :in p53 :where (= (label) 45A))', '/args/where/args/0')).toBe('seq.position.label');
    expect(at('(positions :where (> (length) 100))', '/args/where/args/0')).toBe('seq.system.length');
    expect(at('(rows :of msa :where (matches (label) HUMAN))', '/args/where/args/0')).toBe('seq.system.label');
    expect(at('(partners 1brs/B :via contacts :where (= (category) hbond))', '/args/where/args/0')).toBe(
      'seq.item.category',
    );
    // the innermost binder rules: a generator inside a test brings its own chain
    expect(
      at('(pick (all) :where (overlaps (items :where (> (length) 3))))', '/args/where/args/0/args/where/args/0'),
    ).toBe('seq.item.length');
  });

  it('step 3: with no enclosing :where, the error names where the name works', () => {
    expect(renderIssue(one('(label)'))).toBe(
      '1:2 seqql.slot: label reads a position, a system or an item: use it inside a :where\n  (label :numbering N) or (label)',
    );
    expect(one('(letter)').message).toBe('letter reads a position: use it in (positions :where …)');
    // a selection inside a test reads no slot (seqspec:seqql#elements-and-groups): computed once, not per element
    const i = one('(pick (all) :where (> (count (letters (letter))) 2))');
    expect([i.message, i.path]).toEqual([
      'letter reads a position: use it in (positions :where …)',
      '/args/where/args/0/args/0/args/0',
    ]);
  });

  it('step 3: a chain with no fitting candidate names the binders where the name works', () => {
    expect(renderIssue(one('(items :where (= (letter) C))'))).toBe(
      '1:19 seqql.slot: letter reads a position: use it in (positions :where …)\n  (letter)',
    );
    expect(one('(items :where (> (score) 1))').message).toBe(
      'score reads a system or a group: use it in (systems :where …) or (rows :where …) or (pick Q :where …)',
    );
    expect(one('(items :where (= (label :numbering auth) x))').message).toBe(
      'label reads a position: use it in (positions :where …)',
    );
  });

  it("an extension's property shares a short name across slots and resolves like the base's", () => {
    const prop = (name: string, reads: 'item' | 'system') =>
      defineSymbol({
        name,
        alias: 'depth',
        kind: 'property',
        reads,
        args: [],
        result: 'num',
        signature: '(depth)',
        doc: 'a depth',
        examples: [],
        hints: [],
      });
    const { language } = createLanguage([
      { symbols: [prop('x.item.depth', 'item'), prop('x.system.depth', 'system')] },
    ]);
    expect(headAt(expr('(pick (all) :where (> (depth) 3))', language), '/args/where/args/0')).toBe('x.item.depth');
    expect(headAt(expr('(systems :where (> (depth) 3))', language), '/args/where/args/0')).toBe('x.system.depth');
  });
});

describe('arguments and coercion', () => {
  it('names positionals, gathers a rest positional into an array, and keeps words and strings as strings', () => {
    expect(expr('(and A B C)')).toEqual({ head: 'seq.combinator.intersect', args: { '0': 'A', '1': ['B', 'C'] } });
    expect(expr('(slice "1brs/A" 9 20 30 40)')).toEqual(Q.seq.modifier.slice('1brs/A', 9, 20, 30, 40));
    expect(expr('(range 1brs/A 10-20 45 :numbering uniprot)')).toEqual({
      head: 'seq.modifier.range',
      args: { '0': '1brs/A', '1': ['10-20', 45], numbering: 'uniprot' },
    });
    expect(expr('(input x)')).toEqual({ head: 'seq.generator.input', args: { '0': 'x' } });
  });

  it("gives the builder's JSON: keywords in signature order, defaults omitted", () => {
    expect(JSON.stringify(expr('(items :where (= (id) a) :of domains)'))).toBe(
      JSON.stringify(Q.seq.generator.items({ of: 'domains', where: Q.core.rel.eq(Q.seq.item.id(), 'a') })),
    );
    expect(expr('(first (all) 1)')).toEqual({ head: 'seq.filter.first', args: { '0': { head: 'seq.generator.all' } } });
    expect(expr('(systems :kind unit :of 1brs)')).toEqual({ head: 'seq.generator.systems', args: { of: '1brs' } });
    expect(expr('(items :where (matches (id) x :case false))')).toEqual(
      Q.seq.generator.items({ where: Q.core.str.matches(Q.seq.item.id(), 'x', { case: false }) }),
    );
  });

  it('binds (field label) as (label), and keeps a full-named field as written', () => {
    expect(headAt(expr('(items :where (= (field label) x))'), '/args/where/args/0')).toBe('seq.item.label');
    expect(headAt(expr('(systems :where (= (field kind) row))'), '/args/where/args/0')).toBe('seq.system.kind');
    expect(expr('(systems :where (= (field molecule) p53))')).toEqual(
      Q.seq.generator.systems({ where: Q.core.rel.eq(Q.seq.system.field('molecule'), 'p53') }),
    );
    expect(headAt(expr('(items :where (= (seq.item.field label) x))'), '/args/where/args/0')).toBe('seq.item.field');
  });

  it('keyword errors: unknown (did-you-mean or the list), repeated, missing a value', () => {
    expect(renderIssue(one('(items :in domains)'))).toBe(
      '1:8 seqql.arg.unknown: items has no :in; did you mean :of?\n  (items :of A :where T)',
    );
    expect(one('(rows :frobs x)').message).toBe('rows has no :frobs; it takes :of, :where');
    // the symbol's hint, when it has one, says what it takes
    expect(one('(items :frobs x)').message).toBe(
      'items has no :frobs; it takes :of and :where; items touching a selection: (pick (items :of A) :where (overlaps Q))',
    );
    expect(one('(all :of x)').message).toBe('all has no :of; it takes no keywords');
    expect(renderIssue(one('(items :of a :of b)'))).toBe(
      '1:14 seqql.arg.repeated: :of is given twice\n  (items :of A :where T)',
    );
    expect(one('(items :of :where (= (id) a))').message).toBe(':of needs a value: (items :of A :where T)');
    expect(one('(items :of)').code).toBe('seqql.arg.missing');
  });

  it('positional errors: after a keyword, and too many', () => {
    expect(renderIssue(one('(range 1brs/A :numbering uniprot 10-20)'))).toBe(
      '1:34 seqql.arg.order: positional arguments come before keywords: (range 1brs/A 10-20 :numbering uniprot)\n  (range Q part... :numbering N)',
    );
    expect(one('(expand Q 5 6)').message).toBe(
      'expand takes one positional after the selection, the width: (expand Q 5)',
    );
    expect(one('(all x)').message).toBe('all takes 0 positional arguments: (all)');
    expect(one('(onto Q R S)').message).toBe('onto takes 2 positional arguments: (onto Q R :kind K :keep b)');
  });

  it('reports every problem, each with its JSON path, position and signature', () => {
    const list = issues('(or (item :of a)\n    (items :in b))');
    expect(list.map((i) => [i.code, i.path, i.at, i.signature])).toEqual([
      ['seqql.symbol.unknown', '/args/0', { line: 1, column: 6 }, '(items :of A :where T)'],
      ['seqql.arg.unknown', '/args/1/0', { line: 2, column: 12 }, '(items :of A :where T)'],
    ]);
  });

  it('maps every JSON path to the node it came from', () => {
    const b = bound('(and 1brs/A\n  (sel "p53@0:50"))');
    expect([...b.source].map(([p, n]) => [p, n.at.line, n.at.column])).toEqual([
      ['', 1, 1],
      ['/args/0', 1, 6],
      ['/args/1/0', 2, 3],
      ['/args/1/0/args/0', 2, 8],
    ]);
  });
});

describe('macros', () => {
  // An app's macro, as Stories' chain: (chain A) is (systems :where (= (field auth_asym_id :ns mvs) A)); it reports a
  // lower-case chain at the letter's index in its string.
  const chain = defineSymbol({
    name: 'stories.chain',
    alias: 'chain',
    kind: 'macro',
    args: [{ name: 'C', type: 'str', positional: true, doc: 'the chain' }],
    result: 'query',
    signature: '(chain C)',
    doc: 'units of chain C',
    examples: [],
    hints: [],
    expand: (call, report) => {
      const c = String(call.args?.['0']);
      const lower = c.search(/[a-z]/);
      if (lower >= 0) report({ code: 'x.chain.case', message: 'chain ids are upper-case', index: lower });
      const field = Q.seq.system.field('auth_asym_id', { ns: 'mvs' });
      return Q.seq.generator.systems({ where: Q.core.rel.eq(field, c) });
    },
  });
  const { language } = createLanguage([{ symbols: [chain] }]);

  it('keeps a macro as written in the JSON, and expands it for checking and compiling', () => {
    const e = expr('(and (chain A) (letters C))', language);
    expect(e).toEqual(
      Q.seq.combinator.intersect({ head: 'stories.chain', args: { '0': 'A' } }, Q.seq.generator.letters('C')),
    );
    const x = expandMacros(e, language);
    expect(x.issues).toEqual([]);
    expect(x.expr).toEqual(
      Q.seq.combinator.intersect(
        Q.seq.generator.systems({ where: Q.core.rel.eq(Q.seq.system.field('auth_asym_id', { ns: 'mvs' }), 'A') }),
        Q.seq.generator.letters('C'),
      ),
    );
  });

  it("places a macro's report inside its string, escapes counted as written", () => {
    const b = bound('(or (chain AB)\n  (chain "\\u0041\\"Bc"))', language);
    const x = expandMacros(b.expr as Expression, language, b.source);
    expect(x.issues.map((i) => [i.path, i.at, renderIssue(i)])).toEqual([
      ['/args/1/0', { line: 2, column: 20 }, '2:20 x.chain.case: chain ids are upper-case\n  (chain C)'],
    ]);
  });
});
