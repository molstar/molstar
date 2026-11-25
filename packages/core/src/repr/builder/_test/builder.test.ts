/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder's engine on the schema's test language: the normalized tree it builds, the id and body rules thrown
// where written, the description's issues by `check`, misuse as type errors (`@ts-expect-error`, checked when the
// tests are type-checked), and that its files import only the tree and the schema.

import { readdirSync, readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { parseJson, printJson } from '../../codecs/json/index.js';
import { normalize } from '../../schema/index.js';
import { Mini } from '../../schema/_test/fixtures.js';
import { CategoryColumn, hashNode } from '../../tree/index.js';
import { createBuilder } from '../index.js';

const mini = createBuilder(Mini);

describe('repr/builder: the engine', () => {
  it('builds the normalized tree, leaving out what is absent', () => {
    const node = mini.root({ version: '0.0.0', title: undefined }, (s) => {
      s.data('data', {}, (d) => {
        d.sequence('p53', { alphabet: 'protein', letters: 'MEEPQ' });
        d.ranges('domains', {
          target: '/data/p53',
          rows: { start: [0, 10], end: new Int32Array([5, 20]), id: ['a', 'b'], label: ['x', 'y'] },
          color: mini.forms.field(['id']),
        });
      });
    });
    expect(node.params).toEqual({ version: '0.0.0' });
    const ranges = node.children![0].children![1];
    const rows = ranges.params!.rows as Record<string, unknown>;
    expect(rows.start).toEqual(new Int32Array([0, 10]));
    expect(rows.id).toBeInstanceOf(CategoryColumn);
    expect(ranges.params!.color).toEqual({ kind: 'field', args: ['id'] });
    expect(mini.check(node)).toEqual([]);
    // the tree a codec gives back, once normalized, is the builder's
    expect(hashNode(normalize(Mini, parseJson(printJson(node)).node!))).toBe(hashNode(node));
  });

  it('returns each new node’s address as its handle', () => {
    const handles: (string | undefined)[] = [];
    mini.root({ version: '0' }, (s) => {
      const data: string = s.data('data', {}, (d) => {
        const p53: string = d.sequence('p53', { alphabet: 'protein' });
        handles.push(p53, d.ranges('r', { target: p53, rows: { start: [0], end: [1] } }));
      });
      handles.push(data);
    });
    expect(handles).toEqual(['/data/p53', '/data/r', '/data']);
    // a slot's node has no address, nor what is built under it
    mini.forms.data('d', {}, (d) => {
      const none: undefined = d.sequence('x', { alphabet: 'dna' });
      handles.push(none);
    });
    expect(handles[3]).toBeUndefined();
  });

  it('takes ownership of a typed array it is given: no copy, hashed once', () => {
    const start = new Int32Array([0, 10]);
    const node = mini.root({ version: '0' }, (s) =>
      s.data('data', {}, (d) => d.ranges('r', { target: '/data/p53', rows: { start, end: [5, 20] } })),
    );
    const rows = node.children![0].children![0].params!.rows as Record<string, unknown>;
    expect(rows.start).toBe(start);
    const h = hashNode(node);
    start[0] = 99; // never done to an owned column: the cached hashes show it is not looked at again
    expect(hashNode(node)).toBe(h);
  });

  it('checks a tree against the description', () => {
    const bad = mini.root({ version: '0.0.0' }, (s) =>
      s.data('data', {}, (d) => (d as never as Record<string, (...a: unknown[]) => void>).sequence('x', {})),
    );
    expect(mini.check(bad).map((i) => i.code)).toEqual(['repr.schema.param.missing']);
  });

  it('refuses what is written wrong where it is written', () => {
    const body = (f: (d: Record<string, (...a: unknown[]) => void>) => void) =>
      mini.root({ version: '0' }, (s) => s.data('data', {}, (d) => f(d as never)));
    expect(() =>
      body((d) => {
        d.sequence('a', { alphabet: 'dna' });
        d.sequence('a', { alphabet: 'dna' });
      }),
    ).toThrow(/used twice/);
    expect(() => body((d) => d.sequence('a b', { alphabet: 'dna' }))).toThrow(/is not an id/);
    expect(() => body((d) => d.sequence())).toThrow(/needs an id/);
    expect(() => body((d) => d.data('x'))).toThrow(/does not sit under data/);
    expect(() => body((d) => d.sequence('a', 1))).toThrow(/are an object/);
    expect(() => mini.forms.field('x' as never)).toThrow(/as a list/);
  });

  it('types misuse as errors where it is written', () => {
    const typed = () =>
      mini.root({ version: '0' }, (s) => {
        s.data('data', {}, (d) => {
          // @ts-expect-error: an alphabet the enum lacks
          d.sequence('a', { alphabet: 'klingon' });
          // @ts-expect-error: the required alphabet is missing
          d.sequence('b', {});
          // @ts-expect-error: an int32 column takes numbers
          d.ranges('r', { target: '/data/a', rows: { start: ['x'], end: [1] } });
          // @ts-expect-error: a sequence does not sit under a spec
          s.sequence('c', { alphabet: 'dna' });
        });
      });
    expect(typed).toThrow();
  });

  it('imports only the tree and the schema', () => {
    const dir = new URL('../', import.meta.url);
    for (const f of readdirSync(dir).filter((x) => x.endsWith('.ts'))) {
      const specs = [...readFileSync(new URL(f, dir), 'utf8').matchAll(/from '([^']+)'/g)].map((m) => m[1]);
      for (const spec of specs) expect(spec, f).toMatch(/^\.\.\/(tree|schema)\/index\.js$|^\.\/builder\.js$/);
    }
  });
});
