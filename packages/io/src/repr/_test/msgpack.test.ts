/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { checkCodec, ConformanceLanguage } from '@molstar/core/repr/codec/index';
import { JsonCodec } from '@molstar/core/repr/codecs/json/index';
import { normalize } from '@molstar/core/repr/schema/normalize';
import { CategoryColumn, categoryColumn, hashNode, type Node, pathToPointer } from '@molstar/core/repr/tree/index';
import { describe, expect, it } from 'vitest';
import { encodeMsgPack } from '../../common/msgpack/encode.js';
import { MsgPackCodec, MsgPackExt } from '../msgpack.js';

/** A spec of the shape a viewer's spec has: a sequence, ranges with columns, a track of values, custom data. */
function spec(): Node {
  return {
    kind: 'spec',
    params: { title: 'p53', version: '1.0.0' },
    children: [
      {
        kind: 'data',
        children: [
          {
            kind: 'sequence',
            id: 'p53',
            params: { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' },
          },
          {
            kind: 'ranges',
            id: 'domains',
            params: {
              target: '/data/p53',
              start: new Int32Array([0, 17]),
              end: new Int32Array([14, 30]),
              label: categoryColumn(['Disordered', 'TAD']),
            },
          },
          {
            kind: 'values',
            id: 'conservation',
            params: { target: '/data/p53', values: new Float32Array([0.5, Number.NaN, 1, Number.NEGATIVE_INFINITY]) },
            custom: { note: 'ä α ☃', big: 12345678901234567890n, zero: -0, flags: [true, false, null] },
          },
        ],
      },
    ],
  };
}

/** Bytes of `{ kind: 'x', params: { a: <an extension> } }`, the extension written by hand. */
function withExtension(type: number, data: number[]): Uint8Array {
  const head = [0x82, 0xa4, ...Buffer.from('kind'), 0xa1, 0x78, 0xa6, ...Buffer.from('params'), 0x81, 0xa1, 0x61];
  const body = data.length === 1 ? [0xd4, type, ...data] : [0xc7, data.length, type, ...data];
  return Uint8Array.from([...head, ...body]);
}

describe('io/repr/msgpack', () => {
  it('conforms to the codec interface', () => {
    expect(checkCodec(MsgPackCodec)).toEqual([]);
  });

  it('round-trips a spec through MessagePack and JSON by hash (json(msgpack(x)) is json(x))', () => {
    const node = spec();
    const bytes = MsgPackCodec.encode(node);
    const viaMsgPack = MsgPackCodec.decode(bytes).node!;
    expect(hashNode(viaMsgPack)).toBe(hashNode(node));
    const text = JsonCodec.encode(viaMsgPack);
    const viaJson = JsonCodec.decode(text).node!;
    expect(hashNode(viaJson)).toBe(hashNode(node));
    expect(JsonCodec.encode(MsgPackCodec.decode(MsgPackCodec.encode(viaJson)).node!)).toBe(text);
    // the first byte tells the encodings apart: a map, never `{`
    expect(bytes[0] >> 4).toBe(0x8);
  });

  it('gives every column type back typed, and the atoms JSON lacks', () => {
    const node: Node = {
      kind: 'x',
      params: {
        i8: new Int8Array([-1, 2]),
        u16: new Uint16Array([65535]),
        f64: new Float64Array([Number.NaN, -0, 0.1]),
        c: categoryColumn(['a', null, 'a', 'b']),
        n: Number.NaN,
        nz: -0,
        inf: Number.POSITIVE_INFINITY,
        b: -7n,
      },
    };
    const back = MsgPackCodec.decode(MsgPackCodec.encode(node)).node!;
    expect(back.params?.i8).toEqual(new Int8Array([-1, 2]));
    expect(back.params?.u16).toBeInstanceOf(Uint16Array);
    expect(Object.is((back.params?.f64 as Float64Array)[1], -0)).toBe(true);
    expect(back.params?.c).toBeInstanceOf(CategoryColumn);
    expect(back.params?.n).toBeNaN();
    expect(Object.is(back.params?.nz, -0)).toBe(true);
    expect(back.params?.inf).toBe(Number.POSITIVE_INFINITY);
    expect(back.params?.b).toBe(-7n);
  });

  it('writes the typed columns as extensions 0x11-0x18, little-endian', () => {
    const bytes = MsgPackCodec.encode({ kind: 'x', params: { c: new Uint16Array([0x0102]) } });
    // fixext 2, type 0x14, then 02 01
    const at = bytes.findIndex((b, i) => b === 0xd5 && bytes[i + 1] === MsgPackExt.uint16);
    expect([...bytes.subarray(at, at + 4)]).toEqual([0xd5, 0x14, 0x02, 0x01]);
  });

  it('decodes a column in place where it is aligned, with one copy where it is not', () => {
    const column = new Float64Array(100).map((_, i) => i + 0.5);
    const encoded = MsgPackCodec.encode({ kind: 'x', params: { c: column } });
    // the payload ends the document: shifting the document by (8 - offset mod 8) aligns it
    const payload = encoded.length - column.byteLength;
    const pad = (8 - (payload % 8)) % 8;
    for (const [shift, aligned] of [
      [pad, true],
      [pad + 1, false],
    ] as const) {
      const host = new Uint8Array(encoded.length + 16);
      host.set(encoded, shift);
      const input = host.subarray(shift, shift + encoded.length);
      const out = MsgPackCodec.decode(input).node!.params!.c as Float64Array;
      expect([...out]).toEqual([...column]);
      expect(out.buffer === host.buffer, `in place at shift ${shift}`).toBe(aligned);
      expect(out.byteOffset % 8).toBe(0);
    }
  });

  it('refuses what is not a tree, with issues and never a throw', () => {
    const issue = (bytes: Uint8Array) => MsgPackCodec.decode(bytes).issues.map((i) => i.code);
    expect(issue(new Uint8Array([0xc1]))).toEqual(['repr.codec.msgpack.syntax']);
    expect(issue(new Uint8Array([]))).toEqual(['repr.codec.msgpack.syntax']);
    const good = MsgPackCodec.encode({ kind: 'x' });
    expect(issue(good.subarray(0, good.length - 1))).toEqual(['repr.codec.msgpack.syntax']);
    expect(issue(Uint8Array.of(...good, 0))).toEqual(['repr.codec.msgpack.syntax']);
    // not a node, a bin taken for a column, a reserved key, a bad kind, a key a node lacks
    expect(issue(encodeMsgPack(5))).toEqual(['repr.tree.node.shape']);
    expect(issue(encodeMsgPack({ kind: 'x', params: { a: new Uint8Array([1]) } }))).toEqual(['repr.tree.value.shape']);
    expect(issue(encodeMsgPack({ kind: 'x', params: { $a: 1 } }))).toEqual(['repr.tree.value.reserved_key']);
    expect(issue(encodeMsgPack({ kind: 'X' }))).toEqual(['repr.tree.node.kind']);
    expect(issue(encodeMsgPack({ kind: 'x', other: 1 }))).toEqual(['repr.tree.node.shape']);
  });

  it('reports a malformed extension with its path', () => {
    const where = [['repr.tree.value.tag', '/params/a']];
    const bad = (type: number, data: number[]) =>
      MsgPackCodec.decode(withExtension(type, data)).issues.map((i) => [i.code, pathToPointer(i.path ?? [])]);
    expect(bad(MsgPackExt.number, [9])).toEqual(where);
    expect(bad(MsgPackExt.int32, [1, 2, 3])).toEqual(where);
    expect(bad(MsgPackExt.category, [0, 0, 0])).toEqual(where);
    // a code out of its dictionary
    expect(bad(MsgPackExt.category, [0, 0, 0, 0, 1, 0, 0, 0])).toEqual(where);
    expect(bad(MsgPackExt.bigint, [0x78])).toEqual(where);
    expect(bad(0x42, [1])).toEqual(where);
    // and one that reads
    expect(MsgPackCodec.decode(withExtension(MsgPackExt.number, [3])).node?.params?.a).toBe(-0);
  });

  it('bounds the depth', () => {
    let node: Node = { kind: 'x' };
    for (let i = 0; i < 600; i++) node = { kind: 'x', children: [node] };
    const codes = MsgPackCodec.decode(MsgPackCodec.encode(node)).issues.map((i) => i.code);
    expect(codes).toEqual(['repr.tree.depth_exceeded']);
  });

  it('is typed by the language like any other codec', () => {
    const node: Node = { kind: 'item', params: { i8: new Int8Array([1]), cat: categoryColumn(['a', null]) } };
    const back = normalize(ConformanceLanguage, MsgPackCodec.decode(MsgPackCodec.encode(node)).node!);
    expect(hashNode(back)).toBe(hashNode(node));
  });
});
