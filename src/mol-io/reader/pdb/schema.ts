/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Tokens } from '../common/text/tokenizer';

export interface PdbFile {
    id?: string,
    lines: Tokens
}