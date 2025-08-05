/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Kim Juho <juho_kim@outlook.com>
 */

import { Tokens } from '../common/text/tokenizer';

export interface PdbFile {
    lines: Tokens
    id?: string,
    isPdbqt?: boolean,
    is4LetterResidueName?: boolean
}