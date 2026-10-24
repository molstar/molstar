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
import { Structure } from '../../structure.js';
import { StructureSubsetBuilder } from '../subset-builder.js';

const examples = path.resolve(import.meta.dirname, '../../../../../../../../data/examples');

describe('StructureSubsetBuilder', () => {
  it('getStructureDeduplicate deduplicates each unit by its own elements (perf-audit-2 05/B5)', async () => {
    const parsed = await CIF.parse(fs.readFileSync(path.join(examples, '1cbs_updated.cif'), 'utf8')).run();
    if (parsed.isError) throw new Error(parsed.message);
    const model = (await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run()).representative;
    const root = Structure.ofModel(model);
    const [u0, u1] = root.units;
    expect(u1).toBeDefined();

    const b = new StructureSubsetBuilder(root, false);
    for (const [u, picks] of [
      [u0, [5, 2, 5, 9, 2]],
      [u1, [1, 0, 1]],
    ] as const) {
      b.beginUnit(u.id);
      for (const i of picks) b.addElement(u.elements[i]);
      b.commitUnit();
    }
    const s = b.getStructureDeduplicate();
    expect(s.units.map((u) => Array.from(u.elements))).toEqual([
      [2, 5, 9].map((i) => u0.elements[i]),
      [0, 1].map((i) => u1.elements[i]),
    ]);
  });
});
