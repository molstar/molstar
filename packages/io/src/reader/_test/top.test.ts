/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { utf8Encode } from '@molstar/core/util/utf8';
import { parseTop } from '../top/parser.js';
import { plainResult } from './plain-result.js';

const topString = `; a water and an ion
[ defaults ]
; nbfunc comb-rule gen-pairs
1 2 yes 0.5 0.8333

[ moleculetype ]
; molname nrexcl
Wätér 2

[ atoms ]
; nr type resnr residu atom cgnr charge mass
1 OW 1 SOL OW  1 -0.834 15.9994
2 HW 1 SOL HW1 1  0.417  1.008
3 HW 1 SOL HW2 1  0.417  1.008

[ bonds ]
1 2
1 3

[ moleculetype ]
ION 1

[ atoms ]
1 NA 1 NA NA 1 1.0 22.99

[ system ]
Wätér in a box 中

[ molecules ]
Wätér 2
ION 1
`;

async function parse(data: string | Uint8Array) {
  const parsed = await parseTop(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result;
}

describe('top reader', () => {
  it('reads compounds, the system and the molecules, non-ASCII names included', async () => {
    const top = await parse(topString);
    expect(top.system).toBe('Wätér in a box 中');
    expect(Object.keys(top.compounds)).toEqual(['Wätér', 'ION']);
    const water = top.compounds['Wätér'];
    expect(water.atoms._rowCount).toBe(3);
    expect(water.atoms.atom.value(2)).toBe('HW2');
    expect(water.atoms.charge.value(0)).toBeCloseTo(-0.834, 6);
    expect(water.atoms.mass.value(0)).toBeCloseTo(15.9994, 6);
    expect(water.bonds?.aj.toArray()).toEqual([2, 3]);
    expect(top.compounds['ION'].atoms.atom.value(0)).toBe('NA');
    expect(top.molecules.compound.toArray()).toEqual(['Wätér', 'ION']);
    expect(top.molecules.molCount.toArray()).toEqual([2, 1]);
  });

  it('parses bytes as it parses the text', async () => {
    expect(plainResult(await parse(utf8Encode(topString)))).toEqual(plainResult(await parse(topString)));
  });
});
