/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Test helpers: CIF files re-encoded as text or BinaryCIF, and their models.

import { ArrayEncoding, classifyFloatArray, classifyIntArray } from '@molstar/io/common/binary-cif';
import { CIF, type CifFile, getCifFieldType } from '@molstar/io/reader/cif';
import { CifWriter } from '@molstar/io/writer/cif';
import type { Model } from '../../../../model/model/model.js';
import { trajectoryFromMmCIF } from '../../mmcif.js';

export async function models(file: CifFile) {
  const trajectory = await trajectoryFromMmCIF(file.blocks[0], file).run();
  const ret: Model[] = [];
  for (let i = 0; i < trajectory.frameCount; i++) ret.push((await trajectory.getFrameAtIndex(i)) as Model);
  return ret;
}

export async function parseBinary(data: Uint8Array) {
  const parsed = await CIF.parseBinary(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}

/** As `cif2bcif` encodes a file; with `rawFloats`, float columns as Float64 byte arrays. */
export function toBinary(file: CifFile, rawFloats = false) {
  return encode(file, true, rawFloats) as Uint8Array;
}

/** As `cif2bcif` writes text. */
export function toText(file: CifFile) {
  return encode(file, false, false) as string;
}

function encode(file: CifFile, binary: boolean, rawFloats: boolean) {
  const encoder = CifWriter.createEncoder({ binary, binaryAutoClassifyEncoding: true });
  for (const b of file.blocks) {
    encoder.startDataBlock(b.header);
    for (const c of b.categoryNames) {
      const cat = b.categories[c];
      const fields = cat.fieldNames.map((name): CifWriter.Field => {
        const field = cat.getField(name)!;
        const type = getCifFieldType(field)['@type'];
        if (type === 'str')
          return { name, type: CifWriter.Field.Type.Str, value: field.str, valueKind: field.valueKind };
        if (type === 'float') {
          const encoder = rawFloats
            ? ArrayEncoding.by(ArrayEncoding.byteArray)
            : classifyFloatArray(field.toFloatArray({ array: Float64Array }));
          return CifWriter.Field.float(name, field.float, {
            valueKind: field.valueKind,
            encoder,
            typedArray: Float64Array,
          });
        }
        const encoder = classifyIntArray(field.toIntArray({ array: Int32Array }));
        return CifWriter.Field.int(name, field.int, { valueKind: field.valueKind, encoder, typedArray: Int32Array });
      });
      encoder.writeCategory({
        name: cat.name,
        instance: () => CifWriter.categoryInstance(fields, { data: cat, rowCount: cat.rowCount }),
      });
    }
  }
  return encoder.getData();
}

export async function parseText(data: string) {
  const parsed = await CIF.parseText(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}
