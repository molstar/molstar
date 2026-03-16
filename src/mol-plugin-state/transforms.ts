/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './transforms/data';
import * as Misc from './transforms/misc';
import * as Model from './transforms/model';
import * as Volume from './transforms/volume';
import * as Representation from './transforms/representation';
import * as Shape from './transforms/shape';

// Use lazy getters so that namespace imports are resolved at access time
// rather than at object-construction time. This makes the code safe
// regardless of bundler module evaluation order (circular dependency).
// @see https://github.com/molstar/molstar/issues/1791
export const StateTransforms = {
    get Data() { return Data; },
    get Misc() { return Misc; },
    get Model() { return Model; },
    get Volume() { return Volume; },
    get Representation() { return Representation; },
    get Shape() { return Shape; },
};

export type StateTransforms = typeof StateTransforms