/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { CifWriter } from '@molstar/io/writer/cif';
import { parseCifBinary } from '../parser.js';
import { readCifField } from '../../data-model.js';

function createTestBcif() {
  const rows = [0, 1, 2, 3, 4];
  const encoder = CifWriter.createEncoder({ binary: true, encoderName: 'test-encoder' });
  encoder.startDataBlock('TEST');
  encoder.writeCategory({
    name: 'pdbx_struct_oper_list',
    instance: () =>
      CifWriter.categoryInstance(
        [CifWriter.Field.int<number>('id', (i) => i + 1), CifWriter.Field.float<number>('vector[1]', (i) => i * 1.5)],
        { data: rows, rowCount: rows.length },
      ),
  });
  return encoder.getData() as Uint8Array;
}

describe('binary-cif category', () => {
  it('reads a field once without keeping it (readField), the cached one where read before', async () => {
    const parsed = await parseCifBinary(createTestBcif()).run();
    if (parsed.isError) throw new Error(parsed.message);
    const category = parsed.result.blocks[0].categories.pdbx_struct_oper_list;

    // a new field each read, with the values of the cached one; the cache stays empty
    const read = readCifField(category, 'vector[1]')!;
    expect(Array.from(read.toFloatArray())).toEqual([0, 1.5, 3, 4.5, 6]);
    expect(readCifField(category, 'vector[1]')).not.toBe(read);
    const cached = category.getField('vector[1]')!;
    expect(cached).not.toBe(read);
    expect(cached).toBe(category.getField('vector[1]'));
    // read after: the cached field, no second decode
    expect(readCifField(category, 'vector[1]')).toBe(cached);
    expect(readCifField(category, 'missing')).toBeUndefined();
  });
});
