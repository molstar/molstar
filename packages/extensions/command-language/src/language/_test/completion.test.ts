/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { applyCompletion, complete } from '../completion.js';
import { EmptyEnv, type ScriptEnv } from '../env.js';
import type { Dialect } from '../ir.js';

// Completion (_spec §8): the command word, arguments by kind, selections' words and the scene's names.

const env: ScriptEnv = {
  ...EmptyEnv,
  objects: [{ name: '1cbs', molecule: 0, chains: ['A'], ligands: ['REA'] }],
  names: new Map([['lig', { text: 'resn REA', dialect: 'pymol', parts: [] }]]),
  top: '1cbs',
  nextMolecule: 1,
};

const labels = (text: string, dialect: Dialect = 'pymol') => complete(text, dialect, env).items.map((i) => i.label);

describe('completion', () => {
  it('offers PyMOL commands at the start, and replaces the word typed', () => {
    expect(labels('sh')).toContain('show');
    const c = complete('color red; sh', 'pymol', env);
    const show = c.items.find((i) => i.label === 'show')!;
    expect(applyCompletion('color red; sh', c, show)).toBe('color red; show ');
  });

  it('offers by the argument’s kind', () => {
    expect(labels('show st')).toEqual(['sticks']);
    expect(labels('color sla')).toContain('slate');
    expect(labels('set cartoon_tr')).toEqual(['cartoon_transparency']);
    expect(labels('delete ')).toEqual(['1cbs', 'lig']);
  });

  it('offers a selection’s words, names, chains and ligands', () => {
    expect(labels('show sticks, li')).toContain('lig');
    expect(labels('show sticks, chain ')).toEqual(['A']);
    expect(labels('show sticks, resn ')).toEqual(['REA']);
    expect(labels('show sticks, poly')).toContain('polymer');
  });

  it('offers VMD’s commands, subcommands and arguments', () => {
    expect(labels('mo', 'vmd')).toContain('mol');
    expect(labels('mol add', 'vmd')).toEqual(['addrep']);
    expect(labels('mol addrep ', 'vmd')).toEqual(['top', '0', 'all']);
    expect(labels('mol representation New', 'vmd')).toEqual(['NewCartoon', 'NewRibbons', 'NewSurf']);
    expect(labels('mol color Chai', 'vmd')).toEqual(['Chain']);
    expect(labels('mol selection {prot', 'vmd')).toContain('protein');
  });
});
