/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { transpiler } from '@molstar/query-language/transpilers/jmol/parser';
import type { QueryLanguageEntry } from '../languages.js';

export const JmolLanguage: QueryLanguageEntry = {
  id: 'jmol',
  name: 'Jmol',
  hint: 'Jmol selection syntax, e.g. ":A and 100-120", "[HIS]"; or an address list ("A:100-120").',
  transpile: (source) => transpiler(source),
};
