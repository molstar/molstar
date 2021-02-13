/*
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Expression } from './expression';

export interface Container {
    source?: string,
    version: string,
    expression: Expression
}
