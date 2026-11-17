/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The codecs' conformance suite (_spec/README.md "Codecs"): trees covering the value model, and `checkCodec`, which
// says how a codec breaks the interface's guarantees (codec.ts). Data and a function, no test framework: a codec's
// own test runs it and expects no failure, so a codec outside this package is held to the same suite. The trees are
// normalized by `ConformanceLanguage`, which declares their columns; their columns are short (8 values at most), so a
// printing that elides past a head of 8 keeps them whole.

import { normalize } from '../schema/normalize.js';
import { defineLanguage } from '../schema/check.js';
import { dataFormsOf } from '../schema/language.js';
import { categoryColumn, hashNode, type Node, pathToPointer, type Value } from '../tree/index.js';
import type { Codec } from './codec.js';

const any = { type: 'any', doc: 'Anything.' } as const;
const column = <C extends string>(of: C) => ({ type: 'column', of, doc: `A ${of} column.` }) as const;

/** The language that declares the conformance trees' columns. */
export const ConformanceLanguage = defineLanguage({
  id: 'conformance',
  version: '1.0.0',
  root: 'doc',
  doc: "The codecs' conformance trees.",
  kinds: {
    doc: {
      doc: 'The root.',
      parents: [],
      params: {
        any,
        slot: { type: 'node', language: 'conformance', doc: 'A slot.' },
        list: { type: 'list', of: { type: 'node', language: 'conformance', doc: 'A slot.' }, doc: 'Slots.' },
      },
    },
    item: {
      doc: 'An item: every column type.',
      parents: ['doc', 'item'],
      params: {
        any,
        i8: column('int8'),
        u8: column('uint8'),
        i16: column('int16'),
        u16: column('uint16'),
        i32: column('int32'),
        u32: column('uint32'),
        f32: column('float32'),
        f64: column('float64'),
        cat: { type: 'column', of: 'category', missing: true, doc: 'A category column.' },
        rows: {
          type: 'table',
          columns: {
            start: { type: 'column', of: 'int32', required: true, doc: 'Starts.' },
            score: { type: 'column', of: 'float64', doc: 'Scores.' },
            label: { type: 'list', of: { type: 'string', doc: 'A label.' }, doc: 'Labels.' },
          },
          doc: 'A table.',
        },
      },
    },
    expr: {
      doc: 'A data form.',
      parents: [],
      data: { args: [{ name: 'x', ...any }], rest: true, max: 8 },
      params: { any },
    },
  },
});

const atoms: Record<string, Value> = {
  string: 'text',
  empty: '',
  unicode: 'ä α ☃',
  quote: 'a "quoted" \\ line\nbreak\ttab',
  numberLike: '01',
  wordLike: 'true',
  label: '94A',
  slash: 'HBA_HUMAN/1-141',
  path: '/not//an/address',
  address: '/data/p53/(range :from 94 :to 292)',
  zero: 0,
  negZero: -0,
  float: 1.5,
  big: 1e21,
  small: -1e-7,
  nan: Number.NaN,
  inf: Number.POSITIVE_INFINITY,
  negInf: Number.NEGATIVE_INFINITY,
  bigint: 12345678901234567890n,
  negBigint: -7n,
  yes: true,
  no: false,
  nil: null,
};

/** The conformance trees, normalized. */
export function conformanceTrees(): readonly { readonly name: string; readonly node: Node }[] {
  const trees: { name: string; node: Node }[] = [
    { name: 'the root alone', node: { kind: 'doc' } },
    { name: 'every atom', node: { kind: 'doc', params: { any: atoms } } },
    {
      name: 'lists, records, slots and odd keys',
      node: {
        kind: 'doc',
        params: {
          any: {
            list: [1, 'a', [], [[null]], { k: [true] }],
            record: { nested: { deeper: {} } },
            'p.ALB': 1,
            'a b': 2,
            ä: 3,
          },
          slot: { kind: 'item', params: { any: 'in a slot' } },
          list: [{ kind: 'doc' }, { kind: 'expr', args: [1] }],
        },
      },
    },
    {
      name: 'ids, custom data and extension kinds',
      node: {
        kind: 'doc',
        custom: { note: 'kept', n: [1, 2] },
        children: [
          { kind: 'item', id: 'α.1+b-c_9', custom: { 'ext:color': '#ff0000' } },
          { kind: 'item', id: '0' },
          { kind: 'ext:thing', id: 'e', params: { any: [1, { x: 'y' }] }, children: [{ kind: 'ext:part' }] },
          { kind: 'item' },
        ],
      },
    },
    {
      name: 'data forms',
      node: {
        kind: 'doc',
        params: {
          slot: { kind: 'expr', args: ['/data/p53', 2, { kind: 'expr', args: [[1, 2]] }], params: { any: 1 } },
        },
      },
    },
    {
      name: 'typed columns of every element type',
      node: {
        kind: 'doc',
        children: [
          {
            kind: 'item',
            id: 'columns',
            params: {
              i8: new Int8Array([-128, 0, 127]),
              u8: new Uint8Array([0, 255]),
              i16: new Int16Array([-32768, 32767]),
              u16: new Uint16Array([0, 65535]),
              i32: new Int32Array([-2147483648, 0, 2147483647]),
              u32: new Uint32Array([0, 4294967295]),
              f32: new Float32Array([0.1, Number.NaN, -0, Number.POSITIVE_INFINITY, -1.5e-30]),
              f64: new Float64Array([Number.NaN, -0, 0, Number.NEGATIVE_INFINITY, Math.PI, 1e300, 5e-324]),
              cat: categoryColumn(['b', 'a', null, 'b']),
            },
          },
          {
            kind: 'item',
            id: 'empty',
            params: { i32: new Int32Array(0), f64: new Float64Array(0), cat: categoryColumn([]) },
          },
        ],
      },
    },
    {
      name: 'a table with an extra column',
      node: {
        kind: 'doc',
        children: [
          {
            kind: 'item',
            id: 't',
            params: { rows: { start: [0, 10], score: [0.5, null], label: ['x', 'y'], extra: [true, false] } },
          },
        ],
      },
    },
    { name: 'a deep tree', node: deep(64) },
  ];
  return trees.map(({ name, node }) => ({ name, node: normalize(ConformanceLanguage, node) }));
}

function deep(depth: number): Node {
  let node: Node = { kind: 'item', id: 'leaf', params: { i32: [depth] } };
  for (let i = depth - 1; i >= 0; i--) node = { kind: 'item', id: `n${i}`, children: [node] };
  return { kind: 'doc', children: [node] };
}

/**
 * How a codec breaks the guarantees on the conformance trees (none: it conforms): a decode issue, a tree that does
 * not come back (by hash, after normalization), an encoding that is not deterministic.
 */
export function checkCodec<T extends string | Uint8Array>(codec: Codec<T>): string[] {
  const failures: string[] = [];
  const dataForms = dataFormsOf([ConformanceLanguage]);
  for (const { name, node } of conformanceTrees()) {
    let encoded: T;
    try {
      encoded = codec.encode(node);
    } catch (e) {
      failures.push(`${name}: encode threw (${String(e)})`);
      continue;
    }
    const decoded = codec.decode(encoded, { dataForms });
    if (!decoded.node) {
      const why = decoded.issues.map((i) => `${i.code} at ${i.offset ?? pathToPointer(i.path ?? [])}`).join(', ');
      failures.push(`${name}: decode gave issues (${why})`);
      continue;
    }
    const back = normalize(ConformanceLanguage, decoded.node);
    if (hashNode(back) !== hashNode(node)) failures.push(`${name}: the decoded tree is not the tree`);
    if (!sameOutput(codec.encode(back), encoded)) failures.push(`${name}: encoding is not deterministic`);
  }
  return failures;
}

function sameOutput(a: string | Uint8Array, b: string | Uint8Array): boolean {
  if (typeof a === 'string' || typeof b === 'string') return a === b;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}
