/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { checkCodec } from '../../../codec/index.js';
import { categoryColumn, pathToPointer, TreeIssues } from '../../../tree/index.js';
import { fromJson, JsonCodec, jsonCodec, parseJson, printJson } from '../index.js';

describe('repr/codecs/json', () => {
  it('conforms to the codec interface', () => {
    expect(checkCodec(JsonCodec)).toEqual([]);
    expect(checkCodec(jsonCodec({ indent: 2 }))).toEqual([]);
  });

  it('writes NaN as null inside a float column, and tags the atoms JSON lacks elsewhere', () => {
    const text = printJson({
      kind: 'x',
      params: {
        f: new Float32Array([0.5, Number.NaN, Number.POSITIVE_INFINITY, -0]),
        i: new Int8Array([-1, 2]),
        c: categoryColumn(['a', null]),
        n: Number.NaN,
        z: -0,
        b: 7n,
      },
      custom: { note: 1 },
    });
    expect(text).toBe(
      '{"kind":"x","params":{"f":[0.5,null,{"$number":"Infinity"},{"$number":"-0"}],"i":[-1,2],"c":["a",null],' +
        '"n":{"$number":"NaN"},"z":{"$number":"-0"},"b":{"$bigint":"7"}},"custom":{"note":1}}',
    );
    // schema-free: columns come back as plain arrays (normalization types them)
    expect(parseJson(text).node?.params?.f).toEqual([0.5, null, Number.POSITIVE_INFINITY, -0]);
  });

  it('reports every issue, by code, with its path', () => {
    const { node, issues } = fromJson({
      kind: 'Scene',
      extra: 1,
      id: '',
      params: { $x: 1, ok: { $number: 'nope' }, f: undefined },
      custom: [],
      children: [{ params: {} }, 'x'],
    });
    expect(node).toBeUndefined();
    expect(issues.map((i) => [i.code, pathToPointer(i.path ?? [])])).toEqual([
      ['repr.tree.node.shape', '/extra'],
      ['repr.tree.node.kind', '/kind'],
      ['repr.tree.node.id', '/id'],
      ['repr.tree.value.reserved_key', '/params/$x'],
      ['repr.tree.value.tag', '/params/ok'],
      ['repr.tree.value.shape', '/params/f'],
      ['repr.tree.node.shape', '/custom'],
      ['repr.tree.node.shape', '/children/0'],
      ['repr.tree.node.shape', '/children/1'],
    ]);
    for (const i of issues) {
      expect(TreeIssues.codes).toHaveProperty([i.code]);
      expect(i.severity).toBe('error');
    }
  });

  it('reports text that is not JSON, and skips a byte order mark', () => {
    const { issues } = parseJson('{"kind": "x",');
    expect(issues.map((i) => [i.code, i.severity])).toEqual([['repr.codec.json.syntax', 'error']]);
    expect(parseJson('﻿{"kind":"ext:x"}').node).toEqual({ kind: 'ext:x' });
  });

  it('bounds the depth it builds', () => {
    let v: unknown = 1;
    for (let i = 0; i < 600; i++) v = [v];
    expect(fromJson({ kind: 'x', params: { v } }).issues.map((i) => i.code)).toEqual(['repr.tree.depth_exceeded']);
  });
});
