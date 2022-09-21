

import { loadCheckpoint } from '../../mol-util/debug';
loadCheckpoint(`mol-script/language/container.ts::start`);
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

loadCheckpoint(`mol-script/language/container.ts::end`);
