/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { create } from 'mutative/dist/index.js';

/** Apply changes to an immutable-like object */
export function produce<T>(base: T, recipe: (draft: T) => T | void): T {
    if (typeof base === 'object' && !('prototype' in (base as any))) {
        return create({ ...base }, recipe as any) as T;
    }
    return create(base, recipe as any) as T;
}