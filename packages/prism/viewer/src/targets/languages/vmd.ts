/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { transpiler } from '@molstar/query-language/transpilers/vmd/parser';
import type { QueryLanguageEntry } from '../languages.js';

export const VmdLanguage: QueryLanguageEntry = {
  id: 'vmd',
  name: 'VMD',
  hint: 'VMD selection syntax, e.g. "chain A and resid 100 to 120", "resname HIS"; or an address list ("A:100-120").',
  transpile: (source) => transpiler(source),
};
