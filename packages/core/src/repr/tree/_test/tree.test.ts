/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  CategoryColumn,
  categoryColumn,
  categoryFromCodes,
  columnsEqual,
  hashBytes,
  hashColumn,
  hashNode,
  issueOf,
  KindPattern,
  kindNamespace,
  type Node,
  nodesEqual,
  printAtom,
  readWord,
  toColumn,
  TreeIssues,
  valuesEqual,
} from '../index.js';

describe('repr/tree: nodes and kinds', () => {
  it('takes a kind in a namespace', () => {
    for (const k of ['spec', 'a-b_c1', 'ext:track', 'my-ext:layer_2']) expect(KindPattern.test(k)).toBe(true);
    for (const k of ['Spec', '1a', 'ext:', ':a', 'a:b:c', 'ext:Track']) expect(KindPattern.test(k)).toBe(false);
    expect(kindNamespace('ext:track')).toBe('ext');
    expect(kindNamespace('track')).toBeUndefined();
  });
});

describe('repr/tree: hashNode', () => {
  it("is meso*'s for a tree without the delta: the canonical text and the hash function are a contract", () => {
    const node: Node = { kind: 'population', id: 'p1', params: { name: 'ä', count: 3 } };
    expect(hashNode(node)).toBe('066e920268ab0cde469132d9e0a89ea2');
  });

  it('is stable across key order, and absent params, custom or children are empty ones', () => {
    const a: Node = { kind: 'p', params: { x: 1, y: { b: 2, a: 1 } } };
    const b: Node = { params: { y: { a: 1, b: 2 }, x: 1 }, kind: 'p' };
    expect(hashNode(a)).toBe(hashNode(b));
    expect(hashNode({ kind: 'p' })).toBe(hashNode({ kind: 'p', params: {}, custom: {}, children: [] }));
  });

  it('sees every change: custom data, a column type, a column byte, NaN against null', () => {
    const base: Node = { kind: 'p', params: { c: new Int32Array([1, 2]) }, custom: { x: 1 } };
    const variants: Node[] = [
      { ...base, custom: { x: 2 } },
      { ...base, custom: undefined },
      { ...base, params: { c: new Int32Array([1, 3]) } },
      { ...base, params: { c: new Uint32Array([1, 2]) } },
      { ...base, params: { c: new Float64Array([1, 2]) } },
      { ...base, params: { c: [1, 2] } },
      { ...base, params: { c: categoryColumn(['1', '2']) } },
      { ...base, params: { c: new Float32Array([Number.NaN]) } },
      { ...base, params: { c: new Float32Array([0]) } },
      { ...base, params: { c: new Float32Array([-0]) } },
    ];
    const hashes = new Set([base, ...variants].map(hashNode));
    expect(hashes.size).toBe(variants.length + 1);
  });

  it('equal columns hash and compare equal, whatever their buffers', () => {
    const big = new Int32Array([9, 1, 2, 9]);
    const a: Node = { kind: 'p', params: { c: big.subarray(1, 3) } };
    const b: Node = { kind: 'p', params: { c: new Int32Array([1, 2]) } };
    expect(nodesEqual(a, b)).toBe(true);
    expect(valuesEqual(a.params!.c, b.params!.c)).toBe(true);
    expect(valuesEqual(new Int32Array([1]), [1])).toBe(false);
  });

  it('hashes bytes as content', () => {
    const a = new Uint8Array([1, 2, 3, 4, 5]);
    expect(hashBytes(a)).toMatch(/^[0-9a-f]{32}$/);
    expect(hashBytes(new Uint8Array([9, 1, 2, 3, 4, 5]).subarray(1))).toBe(hashBytes(a));
  });
});

describe('repr/tree: columns', () => {
  it('a category is canonical: its dictionary in the order of first use', () => {
    const a = categoryColumn(['b', 'a', null, 'b']);
    expect(a.dict).toEqual(['b', 'a']);
    expect(Array.from(a.codes)).toEqual([0, 1, -1, 0]);
    const b = categoryFromCodes(['a', 'b', 'unused'], [1, 0, -1, 1]);
    expect(columnsEqual(a, b)).toBe(true);
    expect(hashColumn(a)).toBe(hashColumn(b));
    expect(b.get(2)).toBeNull();
    expect(() => categoryFromCodes(['a'], [1])).toThrow();
  });

  it('every NaN of a float column is the canonical NaN, made so in place (the tree owns the column)', () => {
    const odd = new Float32Array(2);
    new DataView(odd.buffer).setUint32(0, 0x7fc00001, true); // a payload
    new DataView(odd.buffer).setUint32(4, 0xffc00000, true); // the sign bit
    const plain = new Float32Array([Number.NaN, Number.NaN]);
    expect(hashColumn(odd)).toBe(hashColumn(plain));
    expect(columnsEqual(odd, plain)).toBe(true);
    expect(Array.from(new Uint32Array(odd.buffer))).toEqual([0x7fc00000, 0x7fc00000]);
    const odd64 = new Float64Array(1);
    new DataView(odd64.buffer).setUint32(4, 0xfff80000, true);
    expect(toColumn('float64', odd64)).toBe(odd64);
    expect(columnsEqual(odd64, new Float64Array([Number.NaN]))).toBe(true);
  });

  it('takes ownership: no copy, and the hash is computed once', () => {
    const c = new Int32Array([1, 2]);
    expect(toColumn('int32', c)).toBe(c);
    const h = hashColumn(c);
    c[0] = 9; // never done to an owned column: the cached hash shows it is not looked at again
    expect(hashColumn(c)).toBe(h);
  });

  it('identity is by bytes: NaN is NaN, -0 is not 0', () => {
    expect(columnsEqual(new Float64Array([Number.NaN]), new Float64Array([Number.NaN]))).toBe(true);
    expect(columnsEqual(new Float64Array([-0]), new Float64Array([0]))).toBe(false);
    expect(columnsEqual(new Int8Array([1]), new Uint8Array([1]))).toBe(false);
  });

  it('types a plain array by its declared element type, or refuses it', () => {
    expect(toColumn('int32', [1, -2])).toEqual(new Int32Array([1, -2]));
    expect(toColumn('uint8', [256])).toBeUndefined();
    expect(toColumn('int16', [1.5])).toBeUndefined();
    const f = toColumn('float32', [0.5, null]) as Float32Array;
    expect(f[0]).toBe(0.5);
    expect(Number.isNaN(f[1])).toBe(true);
    expect(toColumn('category', ['a', null])).toBeUndefined();
    expect(toColumn('category', ['a', null], true)).toBeInstanceOf(CategoryColumn);
    const same = new Int32Array([1]);
    expect(toColumn('int32', same)).toBe(same);
    expect(toColumn('int32', new Uint32Array([1]))).toBeUndefined();
  });
});

describe('repr/tree: atoms', () => {
  it('a word that parses as a JavaScript number is a number; any other is a symbol', () => {
    expect(readWord('-1')).toBe(-1);
    expect(readWord('1e5')).toBe(1e5);
    expect(readWord('.5')).toBe(0.5);
    expect(Object.is(readWord('-0'), -0)).toBe(true);
    expect(readWord('94A')).toBe('94A');
    expect(readWord('0x10')).toBe('0x10');
    expect(readWord('2nm')).toBe('2nm');
    expect(readWord('nil')).toBeNull();
    expect(readWord('nan')).toBeNaN();
    expect(readWord('12n')).toBe(12n);
    expect(readWord('a b')).toBeUndefined();
    expect(readWord(':a')).toBeUndefined();
  });

  it('prints an atom so that it reads back as itself', () => {
    const atoms = ['94A', '01', 'true', 'nil', 'a b', '', '/x', -0, 1e21, Number.NaN, -Infinity, 7n, null, false, 'ä'];
    for (const a of atoms) {
      const text = printAtom(a);
      const back = text.startsWith('"') ? JSON.parse(text) : readWord(text);
      expect(Object.is(back, a), text).toBe(true);
    }
    expect(printAtom('94A')).toBe('94A');
    expect(printAtom('01')).toBe('"01"');
  });
});

describe('repr/tree: issues', () => {
  it('carry their catalog severity', () => {
    const i = issueOf(TreeIssues, 'repr.tree.node.kind', 'Not a kind', { path: ['kind'] });
    expect(i).toEqual({ code: 'repr.tree.node.kind', severity: 'error', message: 'Not a kind', path: ['kind'] });
  });
});
