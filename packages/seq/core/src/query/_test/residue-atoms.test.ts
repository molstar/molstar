/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A polymer's atoms read by its residues' standard components (the author, 2026-11-20: "Residues, as a structure
// viewer reads them"): the table's codes (PYL among them), the element sets made once a code.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { standardAtoms, standardElements } from '../residue-atoms.js';
import { loc, specRun } from './structures.js';

describe('standard components’ atoms', () => {
  it('PYL and SEC in the table; element sets made once a code', () => {
    expect(standardAtoms('PYL')).toContain('NZ');
    expect(standardAtoms('SEC')).toContain('SE');
    expect(standardElements('MET')).toBe(standardElements('MET'));
    expect([...(standardElements('MET') ?? [])].sort()).toEqual(['C', 'N', 'O', 'S']);
    expect(standardElements('XYZ')).toBeUndefined();
  });

  it('an atom field on a protein with pyrrolysine: the residues whose component has it', () => {
    const b = createSeqSpec();
    b.data.sequence('p', { alphabet: 'protein', letters: 'AOKOU' });
    const t = specRun(b.build());
    expect(loc(t.run('/data/p/(elements :atom CA2)'), '/data/p')).toEqual([1, 2, 3, 4]);
    expect(loc(t.run('/data/p/(elements :element SE)'), '/data/p')).toEqual([4, 5]);
  });
});
