/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Tokens } from '../common/text/tokenizer.js';

export interface PdbFile {
    lines: Tokens
    id?: string,
    variant: 'pdb' | 'pdbqt' | 'pqr'
}