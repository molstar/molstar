/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import { compileQuery } from '../compile.js';
import { queryLanguages, registerQueryLanguages } from '../languages.js';
import { PyMolLanguage } from '../languages/pymol.js';
import { VmdLanguage } from '../languages/vmd.js';
import { parseTarget } from '../target.js';

// The query languages a viewer reads are its spec's (targets/languages.ts): with none, keywords and address lists still
// read and text says so; with PyMOL only, another language is refused by name. In order: the registry is per module.

describe('query languages', () => {
  it('reads keywords and address lists with none, and says text needs a language', () => {
    expect(queryLanguages()).toEqual([]);
    expect(parseTarget('ligand')).toMatchObject({ ok: true, parsed: { kind: 'keyword' } });
    expect(parseTarget('A:120-130')).toMatchObject({ ok: true, parsed: { kind: 'addresses' } });
    expect(compileQuery('chain A and resn HIS')).toMatchObject({
      ok: false,
      message: 'It is not a keyword or an address list, and this viewer reads no query language.',
    });
  });

  it('reads the languages registered, and refuses another by name', () => {
    registerQueryLanguages([PyMolLanguage]);
    registerQueryLanguages([PyMolLanguage]);
    expect(queryLanguages().map((l) => l.id)).toEqual(['pymol']);
    expect(compileQuery('chain A and resn HIS').ok).toBe(true);
    const vmd = compileQuery('resname HIS', 'vmd');
    expect(vmd).toMatchObject({ ok: false, message: 'This viewer does not read VMD.' });
    expect(!vmd.ok && vmd.hint).toMatch(/^Write it in PyMOL, or as a keyword or an address list/);
  });

  it("a viewer registers its spec's languages when it is made", async () => {
    const v = await createPrismViewerHeadless({ spec: { ...DefaultPrismViewerSpec(), languages: [VmdLanguage] } });
    expect(queryLanguages().map((l) => l.id)).toEqual(['pymol', 'vmd']);
    expect(compileQuery('resname HIS', 'vmd').ok).toBe(true);
    await v.dispose();
  });
});
