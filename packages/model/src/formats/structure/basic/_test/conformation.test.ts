/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { classifyFloatArray, classifyIntArray } from '@molstar/io/common/binary-cif';
import { CIF, type CifFile, getCifFieldType } from '@molstar/io/reader/cif';
import { CifWriter } from '@molstar/io/writer/cif';
import type { Model } from '../../../../model/model/model.js';
import { trajectoryFromMmCIF } from '../../mmcif.js';

// The conformation's coordinates are Float32 from text and BinaryCIF alike: a BinaryCIF `Cartn_*` decodes straight to
// Float32 (binary/field.ts), and a sorted or windowed `atom_site` gathers them from that array (sort.ts). Each value
// equals the source value rounded once to Float32.

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

async function models(file: CifFile) {
  const trajectory = await trajectoryFromMmCIF(file.blocks[0], file).run();
  const ret: Model[] = [];
  for (let i = 0; i < trajectory.frameCount; i++) ret.push((await trajectory.getFrameAtIndex(i)) as Model);
  return ret;
}

async function parseBinary(data: Uint8Array) {
  const parsed = await CIF.parseBinary(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}

/** As `cif2bcif` encodes a file. */
function toBinary(file: CifFile) {
  const encoder = CifWriter.createEncoder({ binary: true, binaryAutoClassifyEncoding: true });
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
          const encoder = classifyFloatArray(field.toFloatArray({ array: Float64Array }));
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
  return encoder.getData() as Uint8Array;
}

function expectFloat32Coordinates(model: Model, file: CifFile) {
  const { x, y, z } = model.atomicConformation;
  const atom_site = file.blocks[0].categories['atom_site'];
  const source = model.atomicHierarchy.atomSourceIndex;
  expect(x).toBeInstanceOf(Float32Array);
  for (const [axis, values] of [
    ['x', x],
    ['y', y],
    ['z', z],
  ] as const) {
    const field = atom_site.getField(`Cartn_${axis}`)!;
    for (let i = 0; i < values.length; i++) {
      const expected = Math.fround(field.float(source.value(i)));
      if (values[i] !== expected) expect([axis, i, values[i]]).toEqual([axis, i, expected]);
    }
  }
}

describe('BinaryCIF conformation', () => {
  it('reads a model in order and a sorted model as Float32', async () => {
    for (const name of ['1tup.bcif', '21nr.bcif']) {
      const file = await parseBinary(new Uint8Array(fs.readFileSync(path.join(examples, name))));
      const [model] = await models(file);
      expectFloat32Coordinates(model, file);
    }
    // 21nr's rows move
    const file = await parseBinary(new Uint8Array(fs.readFileSync(path.join(examples, '21nr.bcif'))));
    const [model] = await models(file);
    const source = model.atomicHierarchy.atomSourceIndex;
    expect(Array.from(source.toArray()).some((s, i) => s !== i)).toBe(true);
  });

  it('reads the models of a multi-model file as Float32, text and binary alike', async () => {
    const text = fs.readFileSync(path.join(examples, '1crn.cif'), 'utf8');
    const block = text.slice(text.indexOf('_atom_site.group_PDB') - 'loop_\n'.length);
    const header = text.slice(0, text.indexOf('_atom_site.group_PDB') - 'loop_\n'.length);
    const lines = block.split('\n');
    const rows = lines.filter((l) => l.startsWith('ATOM') || l.startsWith('HETATM'));
    const names = lines.filter((l) => l.startsWith('_atom_site.')).map((l) => l.trim());
    const modelColumn = names.indexOf('_atom_site.pdbx_PDB_model_num');
    expect(modelColumn).toBeGreaterThan(0);
    // a second model: the same rows with model number 2
    const second = rows.map((r) => {
      const cols = r.trim().split(/\s+/);
      cols[modelColumn] = '2';
      return cols.join(' ');
    });
    const twoModels = `${header}loop_\n${names.join('\n')}\n${[...rows, ...second].join('\n')}\n#\n`;
    const parsed = await CIF.parseText(twoModels).run();
    if (parsed.isError) throw new Error(parsed.message);

    const textModels = await models(parsed.result);
    expect(textModels.length).toBe(2);
    for (const m of textModels) expectFloat32Coordinates(m, parsed.result);

    const binary = await parseBinary(toBinary(parsed.result));
    const binaryModels = await models(binary);
    expect(binaryModels.length).toBe(2);
    for (const m of binaryModels) expectFloat32Coordinates(m, binary);
    expect(Array.from(binaryModels[1].atomicConformation.x)).toEqual(Array.from(textModels[1].atomicConformation.x));
  });
});
