/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import { Structure } from '@molstar/model/model/structure';
import { CustomStructureProperty } from '../custom-structure-property.js';

const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };

async function load1crn() {
  const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/1crn.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  return Structure.ofModel((await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative);
}

/** A provider whose computations wait until released; `fail` makes the next one throw */
function provider() {
  const defaultParams = { k: PD.Numeric(1) };
  const runs: { k: number; release: () => void; fail: (e: Error) => void }[] = [];
  const p = CustomStructureProperty.createProvider({
    label: 'Test',
    descriptor: CustomPropertyDescriptor({ name: 'test-in-flight' }),
    type: 'root',
    defaultParams,
    getParams: () => defaultParams,
    isApplicable: () => true,
    obtain: (_ctx, _data, props) =>
      new Promise((resolve, reject) => {
        runs.push({ k: props.k, release: () => resolve({ value: props.k * 10 }), fail: reject });
      }),
  });
  return { p, runs };
}

const tick = () => new Promise((r) => setTimeout(r, 0));

// Attaches that overlap with the same props share one computation (perf-audit-2 06/B9).
describe('CustomStructureProperty.attach in flight', () => {
  it('computes once for overlapping attaches with the same props', async () => {
    const structure = await load1crn();
    const { p, runs } = provider();
    const a = p.attach(ctx, structure);
    const b = p.attach(ctx, structure, { k: 1 });
    await tick();
    expect(runs.length).toBe(1);
    runs[0].release();
    await Promise.all([a, b]);
    expect(p.get(structure).value).toBe(10);
  });

  it('computes again for other props, and after a failed computation', async () => {
    const structure = await load1crn();
    const { p, runs } = provider();
    const a = p.attach(ctx, structure);
    const b = p.attach(ctx, structure, { k: 2 });
    await tick();
    expect(runs.map((r) => r.k)).toEqual([1, 2]);
    runs[0].release();
    runs[1].release();
    await Promise.all([a, b]);
    expect(p.get(structure).value).toBe(20);

    // the first requester's computation fails (an abort): the second computes itself
    const c = p.attach(ctx, structure, { k: 3 });
    const d = p.attach(ctx, structure, { k: 3 });
    await tick();
    expect(runs.length).toBe(3);
    runs[2].fail(new Error('aborted'));
    await expect(c).rejects.toThrow('aborted');
    await tick();
    expect(runs.length).toBe(4);
    runs[3].release();
    await d;
    expect(p.get(structure).value).toBe(30);
  });
});
