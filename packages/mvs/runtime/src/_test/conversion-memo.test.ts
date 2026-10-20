/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { generateStateTransition } from '@molstar/mvs/helpers/animation';
import {
  convertMvsToMolstar,
  convertMvsToMolstarMemo,
  createMolstarConversionMemo,
} from '@molstar/mvs/tree/molstar/conversion';

// `exosome`, the largest, is medium (AGENTS.md, Tests)
const Fixtures = ['features', 'motm1', 'exosome'];
const lane = (name: string) => (name === 'exosome' ? ['medium'] : []);

describe('convertMvsToMolstarMemo', () => {
  for (const name of Fixtures) {
    it(`converts every frame of ${name} as convertMvsToMolstar does`, {
      timeout: 60_000,
      tags: lane(name),
    }, async () => {
      const json = fs.readFileSync(path.join(import.meta.dirname, 'fixtures', `${name}.mvsj`), 'utf8');
      const data = MVSData.fromMVSJ(json);
      const snapshots = data.kind === 'multiple' ? data.snapshots : MVSData.stateToStates(data).snapshots;
      for (const sourceUrl of [undefined, 'https://example.org/data/index.mvsj']) {
        let units = 0;
        let reused = 0;
        for (let i = 0; i < snapshots.length; i++) {
          const transition = await generateStateTransition(SyncRuntimeContext, snapshots[i], i, snapshots.length);
          const trees = [snapshots[i].root, ...(transition?.frames.map((f) => f[0]) ?? [])];
          const memo = createMolstarConversionMemo(sourceUrl);
          let previous = new Set<object>();
          for (const tree of trees) {
            const memoized = convertMvsToMolstarMemo(tree, memo);
            expect(JSON.stringify(memoized.tree)).toBe(JSON.stringify(convertMvsToMolstar(tree, sourceUrl)));
            for (const u of memoized.units) {
              units++;
              if (previous.has(u)) reused++;
            }
            previous = memoized.units;
          }
        }
        if (name !== 'features') expect(reused / units).toBeGreaterThan(0.5);
      }
    });
  }

  describe('downloads that do not convert to one Molstar download', () => {
    const structure = (format: string) => ({
      kind: 'parse',
      params: { format },
      children: [
        {
          kind: 'structure',
          params: { type: 'model' },
          children: [
            {
              kind: 'component',
              params: { selector: 'polymer' },
              children: [{ kind: 'representation', params: { type: 'cartoon' } }],
            },
          ],
        },
      ],
    });
    const canvas = { kind: 'canvas', params: { background_color: 'white' } };
    const cases: [string, any[]][] = [
      ['a download without a parse', [{ kind: 'download', params: { url: 'a.bcif' } }, canvas]],
      [
        'two parses of different binary-ness',
        [{ kind: 'download', params: { url: 'a.cif' }, children: [structure('bcif'), structure('mmcif')] }, canvas],
      ],
      [
        'twin downloads',
        [
          { kind: 'download', params: { url: 'a.bcif' }, children: [structure('bcif')] },
          canvas,
          { kind: 'download', params: { url: 'a.bcif' }, children: [structure('bcif')] },
        ],
      ],
    ];
    for (const [name, children] of cases) {
      it(name, () => {
        const tree = { kind: 'root', children } as any;
        const memo = createMolstarConversionMemo(undefined);
        const expected = JSON.stringify(convertMvsToMolstar(tree, undefined));
        expect(JSON.stringify(convertMvsToMolstarMemo(tree, memo).tree)).toBe(expected);
        // again with the memo filled
        expect(JSON.stringify(convertMvsToMolstarMemo(tree, memo).tree)).toBe(expected);
      });
    }
  });
});
