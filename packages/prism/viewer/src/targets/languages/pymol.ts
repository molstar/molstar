/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { transpiler } from '@molstar/query-language/transpilers/pymol/parser';
import type { QueryLanguageEntry } from '../languages.js';

export const PyMolLanguage: QueryLanguageEntry = {
  id: 'pymol',
  name: 'PyMOL',
  hint: 'PyMOL selection syntax, e.g. "chain A and resi 100-120", "resn HIS", "byres polymer within 5 of resn HEM"; or an address list ("A:100-120").',
  transpile: (source) => transpiler(source),
};
