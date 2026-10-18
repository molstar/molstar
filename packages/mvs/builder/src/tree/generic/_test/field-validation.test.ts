/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import {
  dict,
  float,
  int,
  isCertainlyValid,
  list,
  nullable,
  object,
  partial,
  str,
  tuple,
  union,
} from '../field-schema.js';
import { MVSTreeSchema } from '../../mvs/mvs-tree.js';
import { MVSAnimationSchema } from '../../animation/animation-tree.js';
import { ColorT, Vector3 } from '../../mvs/param-types.js';
import type { Tree } from '../tree-schema.js';

const FixtureDir = path.join(import.meta.dirname, '../../../../../runtime/src/_test/fixtures');

/** Sound: whenever the fast check says valid, io-ts decode succeeds. */
function expectSound(type: any, value: unknown) {
  if (isCertainlyValid(type, value))
    expect(type.decode(value)._tag, JSON.stringify(value)?.slice(0, 200)).toBe('Right');
}

describe('isCertainlyValid', () => {
  it('is sound on every field of the fixtures, and certain for most of them', () => {
    let fields = 0;
    let certain = 0;
    for (const name of ['features', 'motm1', 'exosome']) {
      const data = JSON.parse(fs.readFileSync(path.join(FixtureDir, `${name}.mvsj`), 'utf8'));
      for (const snapshot of data.snapshots) {
        const visit = (node: Tree, schema: any) => {
          const nodeSchema = schema.nodes[node.kind].params;
          const params: any = node.params ?? {};
          const fieldSchemas =
            nodeSchema.type === 'simple'
              ? nodeSchema.fields
              : nodeSchema.cases[params[nodeSchema.discriminator]].fields;
          for (const key in fieldSchemas) {
            if (params[key] === undefined) continue;
            fields++;
            if (isCertainlyValid(fieldSchemas[key].type, params[key])) certain++;
            expectSound(fieldSchemas[key].type, params[key]);
          }
          for (const child of node.children ?? []) visit(child, schema);
        };
        visit(snapshot.root, MVSTreeSchema);
        if (snapshot.animation) visit(snapshot.animation, MVSAnimationSchema);
      }
    }
    // (custom object types, e.g. selectors and MolQL, go to `decode`)
    expect(certain / fields).toBeGreaterThan(0.8);
  });

  it('is sound on invalid and edge values', () => {
    const types: any[] = [
      str,
      int,
      float,
      ColorT,
      Vector3,
      nullable(ColorT),
      list(float),
      tuple([float, float]),
      dict(union(int, str), ColorT),
      object({ a: int }, { b: str }),
      partial({ a: int, b: list(str) }),
      union(nullable(ColorT), dict(union(int, str), ColorT)),
    ];
    const values: unknown[] = [
      'red',
      'notacolor',
      '#ff0000',
      '#ff00',
      'red/blue',
      'red/nope',
      1,
      1.5,
      NaN,
      null,
      undefined,
      true,
      [1, 2],
      [1, 2, 3],
      [1, '2', 3],
      [],
      {},
      { 1: 'red', 2: 'blue' },
      { 1: 'red', 2: 'nope' },
      { a: 1 },
      { a: 1.5 },
      { a: 1, b: 'x' },
      { a: 1, b: 2 },
      { b: ['x'] },
      { b: [1] },
      new Date(0),
    ];
    for (const type of types) for (const value of values) expectSound(type, value);
  });
});
