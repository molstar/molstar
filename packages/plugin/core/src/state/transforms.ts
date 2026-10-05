/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Data from './transforms/data.js';
import * as Misc from './transforms/misc.js';
import * as Model from './transforms/model.js';
import * as Particles from './transforms/particles.js';
import * as Volume from './transforms/volume.js';
import * as Representation from './transforms/representation.js';
import * as Shape from './transforms/shape.js';

// Use lazy getters so that namespace imports are resolved at access time
// rather than at object-construction time. This makes the code safe
// regardless of bundler module evaluation order (circular dependency).
// @see https://github.com/molstar/molstar/issues/1791
export const StateTransforms = {
    get Data() { return Data; },
    get Misc() { return Misc; },
    get Model() { return Model; },
    get Particles() { return Particles; },
    get Volume() { return Volume; },
    get Representation() { return Representation; },
    get Shape() { return Shape; },
};

export type StateTransforms = typeof StateTransforms