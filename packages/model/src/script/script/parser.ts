/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexanderose@weirdbyte.de>
 */

import type { Expression } from '../language/expression.js';

export type Parser = (source: string) => Expression
