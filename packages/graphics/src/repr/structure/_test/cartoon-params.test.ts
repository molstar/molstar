/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import { CartoonParams, getCartoonParams } from '../representation/cartoon.js';

describe('cartoon params', () => {
  it('are made once per structure (perf-audit-2 12/P8)', async () => {
    const file = path.resolve(import.meta.dirname, '../../../../../../data/examples/1cbs_updated.cif');
    const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const structure = Structure.ofModel(
      (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative,
    );
    const ctx = {} as any;
    const params = getCartoonParams(ctx, structure);
    expect(getCartoonParams(ctx, structure)).toBe(params);
    expect(params.visuals.defaultValue).toEqual(['polymer-trace']);
    // the shared definitions are not changed
    expect(CartoonParams.visuals.defaultValue).not.toBe(params.visuals.defaultValue);
    // another structure gets its own
    expect(getCartoonParams(ctx, Structure.create(structure.units))).not.toBe(params);
  });
});
