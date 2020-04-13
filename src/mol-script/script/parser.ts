/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexanderose@weirdbyte.de>
 */

import Expression from '../language/expression';

type Parser = (source: string) => Expression

export default Parser;