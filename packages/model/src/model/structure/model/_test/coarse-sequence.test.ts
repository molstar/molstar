/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '../../../../formats/structure/mmcif.js';
import { Sequence } from '../../../sequence.js';

describe('coarse sequences from sphere ranges', () => {
  it('Sequence.ofSequenceRanges numbers from the first seq id', () => {
    const seq = Sequence.ofSequenceRanges(Column.ofIntArray([3, 8]), Column.ofIntArray([7, 12]));
    expect(seq.length).toBe(10);
    expect(seq.seqId.value(0)).toBe(3);
    expect(seq.seqId.value(seq.length - 1)).toBe(12);
  });

  // 8zzh without `_entity_poly_seq`: the sequences come from the spheres' seq id ranges (perf-audit-2 04/B4)
  it('covers every sphere of a chain, the last included', async () => {
    const file = path.resolve(import.meta.dirname, '../../../../../../../data/examples/ihm-coarse-8zzh.cif');
    const text = fs.readFileSync(file, 'utf8').replaceAll('_entity_poly_seq.', '_entity_poly_seq_hidden.');
    const parsed = await CIF.parse(text).run();
    if (parsed.isError) throw new Error(parsed.message);
    const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
    const model = trajectory.representative;
    const { spheres } = model.coarseHierarchy;
    expect(spheres.count).toBeGreaterThan(0);

    const { sequences } = model.sequence;
    expect(sequences.length).toBeGreaterThan(0);
    for (const { entityId, sequence } of sequences) {
      // the seq id range of the entity's first chain
      let first = -1;
      for (let i = 0; i < spheres.count; i++) {
        if (spheres.entity_id.value(i) === entityId) {
          first = i;
          break;
        }
      }
      const asym = spheres.asym_id.value(first);
      let min = Infinity;
      let max = -Infinity;
      for (let i = first; i < spheres.count && spheres.asym_id.value(i) === asym; i++) {
        min = Math.min(min, spheres.seq_id_begin.value(i));
        max = Math.max(max, spheres.seq_id_end.value(i));
      }
      expect(sequence.seqId.value(0)).toBe(min);
      expect(sequence.seqId.value(sequence.length - 1)).toBe(max);
      expect(sequence.length).toBe(max - min + 1);
    }
  });
});
