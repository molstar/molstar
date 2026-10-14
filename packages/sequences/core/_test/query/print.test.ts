/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL2a: the canonical printer (seqspec:seqql#the-printer): short names where they bind back, the layout, the round trip
// over every example of the table and over random expressions built with Q; and a fuzz: no input makes the pipeline
// throw.
import { describe, expect, it } from 'vitest';
import {
  type ArgDef,
  BASE_SYMBOLS,
  bindQuery,
  createLanguage,
  defineSymbol,
  type Expression,
  type Keywords,
  printQuery,
  Q,
  type QueryNode,
  readQuery,
  type SymbolDef,
  type Type,
} from '../../src/query/index.js';

const BASE = createLanguage().language;

/** read, then bind; null when either reports an issue. */
function parse(text: string, language = BASE): Expression | null {
  const r = readQuery(text);
  return r.node ? bindQuery(r.node, language).expr : null;
}

describe('printing', () => {
  it('prints every example of the table as written: each is canonical and binds back', () => {
    let count = 0;
    for (const s of BASE_SYMBOLS) {
      for (const ex of s.examples) {
        const e = parse(ex.text);
        expect(e, ex.text).not.toBeNull();
        expect(printQuery(e as Expression)).toBe(ex.text);
        count++;
      }
    }
    expect(count).toBeGreaterThan(60);
  });

  it('writes strings as words when they read back as the same word, else as JSON strings', () => {
    const p = (...v: Expression[]) => printQuery(Q.seq.filter.ids('x', ...v));
    expect(p('PF00069', '1brs/entity:1', '45A', '-5-10', '007', 'é')).toBe(
      '(ids x PF00069 1brs/entity:1 45A -5-10 007 é)',
    );
    expect(p('123', '1.5e3', 'true', 'false', ':of', '$s', 'a b', '')).toBe(
      '(ids x "123" "1.5e3" "true" "false" ":of" "$s" "a b" "")',
    );
    expect(p('x;y', '[1]', 'a,b', 'q"', '(', 't\tx')).toBe('(ids x "x;y" "[1]" "a,b" "q\\"" "(" "t\\tx")');
    expect(printQuery(Q.seq.modifier.range('1brs/A', 45, 1.5, -5, 1e21))).toBe('(range 1brs/A 45 1.5 -5 1e+21)');
    // sel's string is always quoted
    expect(printQuery(Q.seq.macro.sel('p53@0:50'))).toBe('(sel "p53@0:50")');
  });

  it('writes positionals, then keywords in signature order, defaults omitted', () => {
    const e = { head: 'seq.generator.search', args: { regex: true, '0': 'K..K', in: 'p53', edits: 0 } };
    expect(printQuery(e)).toBe('(search K..K :in p53 :regex true)');
    expect(printQuery({ head: 'seq.filter.first', args: { '0': 'x', '1': 1 } })).toBe('(first x)');
  });

  it('writes the full name where the short name would bind to another symbol', () => {
    const pick = (t: Expression) =>
      printQuery(Q.seq.filter.pick(Q.seq.generator.all(), { where: Q.core.rel.eq(t, 'x') }));
    expect(pick(Q.seq.item.label())).toBe('(pick (all) :where (= (label) x))');
    expect(pick(Q.seq.system.label())).toBe('(pick (all) :where (= (seq.system.label) x))');
    expect(pick(Q.seq.system.field('identity', { of: 'stats' }))).toBe(
      '(pick (all) :where (= (field identity :of stats) x))',
    );
    expect(pick(Q.seq.system.field('identity'))).toBe('(pick (all) :where (= (seq.system.field identity) x))');
    // (field label) would bind as (label)
    expect(pick(Q.seq.item.field('label'))).toBe('(pick (all) :where (= (seq.item.field label) x))');
    expect(printQuery(Q.seq.item.label())).toBe('(seq.item.label)');
    expect(printQuery(Q.core.logic.and('1brs/A', '1brs/B'))).toBe('(core.logic.and 1brs/A 1brs/B)');
    expect(printQuery(Q.seq.combinator.intersect('1brs/A', '1brs/B'))).toBe('(and 1brs/A 1brs/B)');
    // an extension's symbol without a short name, and a head the language lacks, by full name
    const bare = defineSymbol({
      name: 'x.thing',
      alias: null,
      kind: 'generator',
      args: [],
      result: 'query',
      signature: '(x.thing)',
      doc: 'a thing',
      examples: [],
      hints: [],
    });
    const { language } = createLanguage([{ symbols: [bare] }]);
    expect(printQuery(Q.seq.combinator.merge(bare(), 'a'), language)).toBe('(or (x.thing) a)');
    expect(printQuery({ head: 'y.other', args: { '0': 'a', k: 1 } })).toBe('(y.other a :k 1)');
  });

  it('breaks a long form: the head and the positionals that fit, then one argument a line, two spaces in', () => {
    const a =
      '(onto (and 1brs/A (partners 1brs/B :via contacts :where (= (category) hbond)) (onto (values :of cons :above 0.8) 1brs/A)) msa :keep true)';
    expect(printQuery(parse(a) as Expression)).toBe(
      [
        '(onto',
        '  (and 1brs/A (partners 1brs/B :via contacts :where (= (category) hbond))',
        '    (onto (values :of cons :above 0.8) 1brs/A))',
        '  msa',
        '  :keep true)',
      ].join('\n'),
    );
    const b =
      '(pick (items :of variants :where (and (matches (label) "pathogenic|likely pathogenic") (> (length) 1))) :where (overlaps (sel "domains[#PF00069,#PF07714]")))';
    const printed = printQuery(parse(b) as Expression);
    expect(printed).toBe(
      [
        '(pick',
        '  (items',
        '    :of variants',
        '    :where (and (matches (label) "pathogenic|likely pathogenic") (> (length) 1)))',
        '  :where (overlaps (sel "domains[#PF00069,#PF07714]")))',
      ].join('\n'),
    );
    expect(printQuery(parse(printed) as Expression)).toBe(printed);
  });

  it('Q.query stores the version, the JSON and its canonical text', () => {
    const kinases = Q.seq.filter.pick(Q.seq.generator.items({ of: 'domains' }), {
      where: Q.core.str.matches(Q.seq.item.label(), 'kinase'),
    });
    expect(Q.query(kinases)).toEqual({
      seqql: '1.0',
      expr: kinases,
      text: '(pick (items :of domains) :where (matches (label) kinase))',
    });
  });
});

/** mulberry32: a seeded generator, so a failure repeats. */
function random(seed: number) {
  let a = seed;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

describe('random expressions built with Q round-trip', () => {
  const STRINGS = [
    '1brs/A',
    'msa/3',
    'domains',
    'a b',
    '123',
    'true',
    '',
    ':x',
    '$s',
    'x;y',
    'q"',
    'é',
    '[1-2]',
    '-5-10',
    '45A',
    ' ',
    'back\\slash',
  ];
  const queries = BASE_SYMBOLS.filter((s) => s.result === 'query');
  const tests = BASE_SYMBOLS.filter((s) => s.result === 'bool');
  const properties = BASE_SYMBOLS.filter((s) => s.kind === 'property');

  function generate(next: () => number, depth: number) {
    const pickOf = <T>(xs: readonly T[]) => xs[Math.floor(next() * xs.length)] as T;
    const value = (t: Type): Expression => {
      // now and then a value of the wrong kind: the printer must fall back to full names where a short one would
      // bind elsewhere
      const wrong = next() < 0.08;
      if (typeof t === 'object') return pickOf(t.enum);
      if ((t === 'query' || wrong) && depth > 0 && next() < 0.6) return call(pickOf(queries));
      if ((t === 'bool' || t === 'test' || wrong) && depth > 0 && next() < 0.7) return call(pickOf(tests));
      if ((t === 'json' || t === 'num' || t === 'int' || wrong) && depth > 0 && next() < 0.4)
        return call(pickOf(properties));
      if (t === 'int') return Math.floor(next() * 100);
      if (t === 'num' || (t === 'json' && next() < 0.3)) return Math.round((next() - 0.3) * 1e4) / 100 || 1;
      if (t === 'bool') return next() < 0.5;
      return pickOf(STRINGS);
    };
    function call(s: SymbolDef): Expression {
      depth--;
      const positional: Expression[] = [];
      const keywords: Record<string, Expression> = {};
      for (const a of s.args as readonly ArgDef[]) {
        if (a.positional && a.rest) for (let k = Math.floor(next() * 3); k > 0; k--) positional.push(value(a.type));
        else if (a.positional && (!a.optional || next() < 0.5)) positional.push(value(a.type));
        else if (!a.positional && (!a.optional || next() < 0.4)) keywords[a.name] = value(a.type);
      }
      depth++;
      return s(...positional, keywords as Keywords);
    }
    return call(pickOf(BASE_SYMBOLS));
  }

  it('binds each print back to the expression, and printing is a fixed point', () => {
    const next = random(20261012);
    let multiline = 0;
    let fullNames = 0;
    for (let k = 0; k < 600; k++) {
      const e = generate(next, 4);
      const text = printQuery(e);
      const r = readQuery(text);
      expect(r.issues, text).toEqual([]);
      const b = bindQuery(r.node as QueryNode, BASE);
      expect(b.issues, text).toEqual([]);
      expect(b.expr, text).toEqual(e);
      expect(printQuery(b.expr as Expression)).toBe(text);
      if (text.includes('\n')) multiline++;
      if (/\((core|seq)\./.test(text)) fullNames++;
    }
    expect(multiline).toBeGreaterThan(30);
    expect(fullNames).toBeGreaterThan(30);
  });
});

describe('fuzz: random text never throws, only reports issues', () => {
  const TOKENS = [
    '(',
    ')',
    '"',
    '\\',
    ';',
    '\n',
    ' ',
    '[',
    ',',
    '$',
    ':',
    ':of',
    ':where',
    '1brs/A',
    '45',
    '-',
    'true',
    'and',
    'or',
    'not',
    'label',
    'pick',
    'items',
    'field',
    'sel',
    'seq.item.label',
    'é',
    '\u0000',
    '(range',
    '(and',
    '"x"',
  ];

  it('reads, binds and prints whatever it is given', () => {
    const next = random(7);
    let bound = 0;
    for (let k = 0; k < 3000; k++) {
      let text = '';
      const n = Math.floor(next() * 24);
      for (let j = 0; j < n; j++) text += TOKENS[Math.floor(next() * TOKENS.length)];
      if (next() < 0.3) text = `(pick (items :of domains) :where ${text})`;
      const r = readQuery(text);
      for (const i of r.issues) expect(i.code).toMatch(/^seqql\.read\./);
      if (!r.node) continue;
      const b = bindQuery(r.node, BASE);
      for (const i of b.issues) expect(i.code).toMatch(/^seqql\./);
      if (!b.expr) continue;
      bound++;
      const printed = printQuery(b.expr);
      expect(parse(printed), printed).toEqual(b.expr);
    }
    expect(bound).toBeGreaterThan(100);
  });
});
