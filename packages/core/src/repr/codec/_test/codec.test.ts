/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { JsonCodec, parseJson, printJson } from '../../codecs/json/index.js';
import { parseSexpr, printSexpr } from '../../codecs/sexpr/index.js';
import { LanguageIssues, normalize, SchemaIssues } from '../../schema/index.js';
import { hashNode, type Node, TreeIssues } from '../../tree/index.js';
import { type Codec, CodecIssues, ConformanceLanguage, checkCodec, conformanceTrees } from '../index.js';

describe('repr/codec: the conformance suite', () => {
  it('finds a codec that loses data, or does not decode what it encodes', () => {
    const dropsCustom: Codec<string> = {
      ...JsonCodec,
      encode: (node) => printJson(strip(node)),
    };
    expect(checkCodec(dropsCustom)).toEqual([
      'ids, custom data and extension kinds: the decoded tree is not the tree',
      'columns outside declared params: custom data, an any param, a NaN with a payload and a sign bit: the decoded tree is not the tree',
    ]);
    const refuses: Codec<string> = { ...JsonCodec, decode: () => parseJson('{') };
    expect(checkCodec(refuses)).toHaveLength(conformanceTrees().length);
    const throws: Codec<string> = {
      ...JsonCodec,
      encode: () => {
        throw new Error('no');
      },
    };
    expect(checkCodec(throws)[0]).toMatch(/encode threw/);
  });

  it('JSON and the printing carry the same tree', () => {
    for (const { node } of conformanceTrees()) {
      const json = parseJson(printJson(node)).node!;
      const sexpr = parseSexpr(printSexpr(node), { dataForms: new Set(['expr']) }).node!;
      expect(hashNode(normalize(ConformanceLanguage, json))).toBe(hashNode(normalize(ConformanceLanguage, sexpr)));
    }
  });
});

function strip(node: Node): Node {
  const { custom: _, ...rest } = node;
  return rest.children ? { ...rest, children: rest.children.map(strip) } : rest;
}

describe('repr: issue codes', () => {
  it('every code a source uses is cataloged, with its severity', () => {
    const catalogs = [TreeIssues, SchemaIssues, LanguageIssues, CodecIssues];
    const known = new Set(catalogs.flatMap((c) => Object.keys(c.codes)));
    for (const c of catalogs)
      for (const info of Object.values(c.codes)) expect(info.severity).toMatch(/^(error|warning|info)$/);
    const root = fileURLToPath(new URL('../../', import.meta.url));
    const used = new Set<string>();
    const walk = (dir: string) => {
      for (const name of readdirSync(dir)) {
        const path = join(dir, name);
        if (statSync(path).isDirectory()) {
          if (name !== '_test' && name !== '_spec') walk(path);
        } else if (name.endsWith('.ts')) {
          for (const m of readFileSync(path, 'utf8').matchAll(/'(repr\.[a-z_.]*[a-z_])'/g)) used.add(m[1]);
        }
      }
    };
    walk(root);
    const namespaces = new Set(catalogs.map((c) => c.namespace));
    const missing = [...used].filter((code) => !known.has(code) && !namespaces.has(code));
    expect(missing).toEqual([]);
    expect(used.size).toBeGreaterThan(40);
  });
});
