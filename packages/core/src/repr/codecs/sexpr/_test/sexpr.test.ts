/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CodecIssues, ConformanceLanguage, checkCodec } from '../../../codec/index.js';
import { Mini } from '../../../schema/_test/fixtures.js';
import { categoryColumn, hashNode, type Node, pathToPointer, TreeIssues } from '../../../tree/index.js';
import { parseSexpr, printSexpr, printSexprLine, SexprCodec, sexprCodec } from '../index.js';

describe('repr/codecs/sexpr', () => {
  it('conforms to the codec interface', () => {
    expect(checkCodec(SexprCodec)).toEqual([]);
    expect(checkCodec(sexprCodec({ language: ConformanceLanguage }))).toEqual([]);
  });

  it('prints a spec for people, and reads the same spelling back', () => {
    const node: Node = {
      kind: 'spec',
      params: { version: '0.0.0', title: 'p53' },
      children: [
        {
          kind: 'data',
          id: 'data',
          children: [
            { kind: 'sequence', id: 'p53', params: { alphabet: 'protein', cyclic: false } },
            {
              kind: 'ranges',
              id: 'domains',
              params: {
                target: '/data/p53',
                rows: { start: new Int32Array([93, 291]), label: ['DNA-binding', 'tetramer'] },
              },
              custom: { 'ext:color': '#ff0000' },
            },
          ],
        },
        {
          kind: 'views',
          id: 'views',
          children: [{ kind: 'layer', params: { data: '/data/p53/(range :from 94 :to 292)', label: '94A' } }],
        },
      ],
    };
    const text = printSexpr(node, { language: Mini });
    expect(text).toBe(
      [
        '(spec :version 0.0.0 :title p53',
        '  (data data',
        '    (sequence p53 :alphabet protein :cyclic false)',
        '    (ranges domains',
        '      :target /data/p53',
        '      :rows (:start [93 291] :label [DNA-binding tetramer])',
        '      :$custom (:"ext:color" #ff0000)))',
        '  (views views (layer :data /data/p53/(range :from 94 :to 292) :label 94A)))',
      ].join('\n'),
    );
    const back = parseSexpr(text);
    expect(back.issues).toEqual([]);
    expect(printSexpr(back.node!, { language: Mini })).toBe(text);
  });

  it('elides a column past its head, and its reader refuses it', () => {
    const node: Node = { kind: 'values', params: { v: new Float32Array(3_000_000).fill(0.5) } };
    const text = printSexprLine(node);
    expect(text).toBe('(values :v {float32 × 3,000,000: 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 …})');
    const { node: back, issues } = parseSexpr(text);
    expect(back).toBeUndefined();
    expect(issues.map((i) => [i.code, i.offset])).toEqual([['repr.codec.sexpr.elided', 11]]);
    expect(printSexprLine({ kind: 'x', params: { c: categoryColumn(['a', null, 'a']) } }, { head: 2 })).toBe(
      '(x :c {category × 3: a nil …})',
    );
    // a longer head keeps it whole; a list of the same values is a list, never an elided column
    expect(sexprCodec({ head: 4 }).encode({ kind: 'x', params: { c: new Int8Array([1, 2, 3, 4]) } })).toBe(
      '(x :c (:$column int8 :values [1 2 3 4]))',
    );
    expect(parseSexpr('(x :c [float32 × 3])').node?.params?.c).toEqual(['float32', '×', 3]);
    expect(parseSexpr('(x :c {int8})').issues.map((i) => i.code)).toEqual(['repr.codec.sexpr.syntax']);
  });

  it('spells a column outside a declared param tagged, read back typed', () => {
    const node: Node = {
      kind: 'x',
      params: { a: categoryColumn(['b', null, 'a']) },
      custom: { f: new Float64Array([1.5, Number.NaN]) },
    };
    const text = printSexprLine(node);
    expect(text).toBe(
      '(x :a (:$column category :dict [b a] :codes [0 -1 1]) :$custom (:f (:$column float64 :values [1.5 nan])))',
    );
    expect(hashNode(parseSexpr(text).node!)).toBe(hashNode(node));
    for (const bad of ['(:$column int8 :values [300])', '(:$column category :dict [a] :codes [2])', '(:$column int8)'])
      expect(
        parseSexpr(`(x :c ${bad})`).issues.map((i) => i.code),
        bad,
      ).toEqual(['repr.tree.value.tag']);
  });

  it('reads a word opening with / as an address, canonically, and a bad one as an issue', () => {
    const { node } = parseSexpr('(layer :data /data/p53/( range  :from 94 ) :at "/x")');
    expect(node?.params).toEqual({ data: '/data/p53/(range :from 94)', at: '/x' });
    expect(parseSexpr('(layer :data /a//b)').issues.map((i) => i.code)).toEqual(['repr.tree.address.syntax']);
  });

  it('reads atoms as tree/atoms.ts spells them', () => {
    const { node } = parseSexpr('(x :a 94A :b -1 :c .5 :d 0x10 :e nan :f 7n :g nil :h 2nm)');
    expect(node?.params).toEqual({ a: '94A', b: -1, c: 0.5, d: '0x10', e: Number.NaN, f: 7n, g: null, h: '2nm' });
  });

  it('reports what is not a tree, by code, offset and path', () => {
    const cases: [string, string, number][] = [
      ['(a', 'repr.codec.sexpr.syntax', 0],
      ['(a) (b)', 'repr.codec.sexpr.shape', 4],
      ['(a {)', 'repr.codec.sexpr.syntax', 3],
      ['(a :x 1 :x 2)', 'repr.codec.sexpr.duplicate_key', 8],
      ['(a :x)', 'repr.codec.sexpr.shape', 3],
      ['(A)', 'repr.tree.node.kind', 1],
      ['(a "b c")', 'repr.tree.node.id', 3],
      ['(a :kind 1)', 'repr.tree.value.reserved_key', 3],
      ['(a :$custom 1)', 'repr.codec.sexpr.shape', 3],
      ['(a b 1)', 'repr.codec.sexpr.shape', 5],
    ];
    for (const [text, code, offset] of cases) {
      const { node, issues } = parseSexpr(text);
      expect(node, text).toBeUndefined();
      expect(
        issues.map((i) => [i.code, i.offset]),
        text,
      ).toEqual([[code, offset]]);
      const catalog = code.startsWith('repr.codec.') ? CodecIssues.codes : TreeIssues.codes;
      expect(catalog).toHaveProperty([code]);
    }
    expect(pathToPointer(parseSexpr('(a (b :x 1 :x 2))').issues[0].path ?? [])).toBe('/children/0/params/x');
  });

  it('reads data forms when told which kinds they are', () => {
    const text = '(expr /data/p53 2 :n 1 (expr [1 2]))';
    const { node } = parseSexpr(text, { dataForms: new Set(['expr']) });
    expect(node).toEqual({ kind: 'expr', params: { n: 1 }, args: ['/data/p53', 2, { kind: 'expr', args: [[1, 2]] }] });
    expect(printSexprLine(node!)).toBe(text);
  });
});
