/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import { utf8Encode } from '@molstar/core/util/utf8';
import { parsePrmtop } from '../prmtop/parser.js';
import { plainResult } from './plain-result.js';

const parm7 = fs.readFileSync(
  path.resolve(import.meta.dirname, '../../../../../data/examples/trajectory/protein.parm7'),
  'utf8',
);

async function parse(data: string | Uint8Array) {
  const parsed = await parsePrmtop(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}

describe('prmtop reader', () => {
  it('reads the pointers and the fixed-width sections', async () => {
    const prmtop = await parse(parm7);
    expect(prmtop.title).toEqual(['alanine-dipeptide.solvated.pdb']);
    const { NATOM, NRES } = prmtop.pointers;
    expect(prmtop.atomName.rowCount).toBe(NATOM);
    expect(prmtop.charge.rowCount).toBe(NATOM);
    expect(prmtop.residueLabel.rowCount).toBe(NRES);
    expect(prmtop.residuePointer.value(0)).toBe(1);
    expect(prmtop.bondsIncHydrogen.rowCount).toBe(prmtop.pointers.NBONH * 3);
  });

  it('parses bytes as it parses the text', async () => {
    expect(plainResult(await parse(utf8Encode(parm7)))).toEqual(plainResult(await parse(parm7)));
  });

  it('reads a non-ASCII title', async () => {
    const text = parm7.replace('alanine-dipeptide', 'alanin-dipeptid é 中');
    const prmtop = await parse(utf8Encode(text));
    expect(prmtop.title).toEqual(['alanin-dipeptid é 中.solvated.pdb']);
    expect(plainResult(prmtop)).toEqual(plainResult(await parse(text)));
  });
});
