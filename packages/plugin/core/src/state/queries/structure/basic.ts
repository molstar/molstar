/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as MS } from '@molstar/model/script/language/builder';
import { StructureSelectionQuery } from './query.js';

export const all = StructureSelectionQuery('All', MS.struct.generator.all(), { category: '', priority: 1000 });
export const current = StructureSelectionQuery('Current Selection', MS.internal.generator.current(), {
  category: '',
  referencesCurrent: true,
});

export const BasicSelectionQueries = [all, current];
