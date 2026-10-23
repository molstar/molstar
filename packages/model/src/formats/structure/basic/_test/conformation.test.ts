/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { CIF, type CifFile } from '@molstar/io/reader/cif';
import type { Model } from '../../../../model/model/model.js';
import { models, parseBinary, toBinary } from './encode.js';

// The conformation's coordinates are Float32 from text and quantized BinaryCIF alike: a quantized `Cartn_*` decodes
// straight to Float32 (binary/field.ts), and a sorted or windowed `atom_site` gathers them from that array (sort.ts).
// Each value equals the source value rounded once to Float32. A raw Float64 column stays Float64 in every path.

const examples = path.resolve(import.meta.dirname, '../../../../../../../data/examples');

function expectFloat32Coordinates(
  model: Model,
  file: CifFile,
  type: Float32ArrayConstructor | Float64ArrayConstructor = Float32Array,
) {
  const { x, y, z } = model.atomicConformation;
  const atom_site = file.blocks[0].categories['atom_site'];
  const source = model.atomicHierarchy.atomSourceIndex;
  expect(x).toBeInstanceOf(type);
  const round = type === Float32Array ? Math.fround : (v: number) => v;
  for (const [axis, values] of [
    ['x', x],
    ['y', y],
    ['z', z],
  ] as const) {
    const field = atom_site.getField(`Cartn_${axis}`)!;
    for (let i = 0; i < values.length; i++) {
      const expected = round(field.float(source.value(i)));
      if (values[i] !== expected) expect([axis, i, values[i]]).toEqual([axis, i, expected]);
    }
  }
}

describe('BinaryCIF conformation', () => {
  it('reads a model in order and a sorted model as Float32 (or a raw column as it is)', async () => {
    for (const name of ['1tup.bcif', '21nr.bcif']) {
      const file = await parseBinary(new Uint8Array(fs.readFileSync(path.join(examples, name))));
      const [model] = await models(file);
      // 1tup's coordinates are FixedPoint, 21nr's raw Float64, which stays so
      const raw = file.blocks[0].categories['atom_site'].getField('Cartn_x')!.binaryEncoding![0].kind === 'ByteArray';
      expectFloat32Coordinates(model, file, raw ? Float64Array : Float32Array);
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

    // a raw Float64 column stays Float64 in every model, as in an in-order one
    const raw = await parseBinary(toBinary(parsed.result, true));
    for (const m of await models(raw)) expectFloat32Coordinates(m, raw, Float64Array);
  });

  it('reads the sorted models of a multi-model file from their own rows', async () => {
    const text = fs.readFileSync(path.join(examples, '1crn.cif'), 'utf8');
    const at = text.indexOf('_atom_site.group_PDB') - 'loop_\n'.length;
    const lines = text.slice(at).split('\n');
    const rows = lines.filter((l) => l.startsWith('ATOM') || l.startsWith('HETATM'));
    const names = lines.filter((l) => l.startsWith('_atom_site.')).map((l) => l.trim());
    const modelColumn = names.indexOf('_atom_site.pdbx_PDB_model_num');
    const second = rows.map((r) => {
      const cols = r.trim().split(/\s+/);
      cols[modelColumn] = '2';
      return cols.join(' ');
    });
    // each model in reverse, so that it sorts
    const reversed = [...rows].reverse().concat([...second].reverse());
    const parsed = await CIF.parseText(
      `${text.slice(0, at)}loop_\n${names.join('\n')}\n${reversed.join('\n')}\n#\n`,
    ).run();
    if (parsed.isError) throw new Error(parsed.message);
    for (const file of [parsed.result, await parseBinary(toBinary(parsed.result))]) {
      const field = file.blocks[0].categories['atom_site'].getField('Cartn_x')!;
      const read = vi.spyOn(field, 'toFloatArray');
      const ms = await models(file);
      expect(ms.length).toBe(2);
      for (const m of ms) {
        expect(Array.from(m.atomicHierarchy.atomSourceIndex.toArray()).some((s, i) => s !== i)).toBe(true);
        expectFloat32Coordinates(m, file);
      }
      // the gather of a model reads its own rows (a text field parses what is read)
      for (const [params] of read.mock.calls) {
        if (params?.array === Float32Array) expect(params.end! - params.start!).toBe(rows.length);
      }
      read.mockRestore();
    }
  });
});
