/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ConformanceLanguage, checkCodec } from '../../../codec/index.js';
import { Mini } from '../../../schema/_test/fixtures.js';
import { CategoryColumn, categoryColumn, type Node, pathToPointer, TreeIssues } from '../../../tree/index.js';
import { fromJson, JsonCodec, jsonCodec, parseJson, printJson } from '../index.js';

describe('repr/codecs/json', () => {
  it('conforms to the codec interface, with and without the language', () => {
    expect(checkCodec(JsonCodec)).toEqual([]);
    expect(checkCodec(jsonCodec({ indent: 2 }))).toEqual([]);
    expect(checkCodec(jsonCodec({ language: ConformanceLanguage }))).toEqual([]);
  });

  const values: Node = {
    kind: 'values',
    id: 'v',
    params: {
      target: '/data/p53',
      values: new Float32Array([0.5, Number.NaN, Number.POSITIVE_INFINITY, -0]),
      scale: Number.NaN,
    },
    custom: { c: categoryColumn(['a', null]), i: new Int8Array([-1, 2]), z: -0, b: 7n },
  };

  it('writes a declared column as a plain array (NaN null), tags any other column and the atoms JSON lacks', () => {
    const text = printJson(values, { language: Mini });
    expect(text).toBe(
      '{"kind":"values","id":"v","params":{"target":"/data/p53","values":[0.5,null,{"$number":"Infinity"},{"$number":"-0"}],' +
        '"scale":{"$number":"NaN"}},"custom":{"c":{"$column":"category","dict":["a"],"codes":[0,-1]},' +
        '"i":{"$column":"int8","values":[-1,2]},"z":{"$number":"-0"},"b":{"$bigint":"7"}}}',
    );
    // schema-free: a declared column comes back plain (normalization types it), a tagged one typed
    const back = parseJson(text).node!;
    expect(back.params?.values).toEqual([0.5, null, Number.POSITIVE_INFINITY, -0]);
    expect(back.custom?.i).toEqual(new Int8Array([-1, 2]));
    expect(back.custom?.c).toBeInstanceOf(CategoryColumn);
    // without the language every column is tagged
    expect(printJson(values)).toContain('"values":{"$column":"float32","values":[0.5,null,');
  });

  it('refuses a malformed tagged column and a key given twice', () => {
    for (const bad of [
      '{"$column":"int8","values":[300]}',
      '{"$column":"int9","values":[]}',
      '{"$column":"category","dict":["a"],"codes":[1]}',
      '{"$column":"int8","values":[],"x":1}',
    ])
      expect(
        parseJson(`{"kind":"x","params":{"c":${bad}}}`).issues.map((i) => i.code),
        bad,
      ).toEqual(['repr.tree.value.tag']);
    const dup = parseJson('{"kind":"x","params":{"a":1,"b":[{"c":1,"\\u0063":2}],"a":3}}');
    expect(dup.node).toBeUndefined();
    expect(dup.issues.map((i) => [i.code, pathToPointer(i.path ?? [])])).toEqual([
      ['repr.codec.json.duplicate_key', '/params/b/0/c'],
      ['repr.codec.json.duplicate_key', '/params/a'],
    ]);
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
