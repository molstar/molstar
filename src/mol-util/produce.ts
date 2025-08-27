/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { create, rawReturn } from 'mutative';

let currentRecipe: any = undefined;
function recipeWrapper(draft: any) {
    const r = currentRecipe(draft);
    if (r !== undefined && r !== draft) return rawReturn(r);
    return r;
}

/** Apply changes to an immutable-like object */
export function produce<T>(base: T, recipe: (draft: T) => T | void): T {
    currentRecipe = recipe;
    if (typeof base === 'object' && !('prototype' in (base as any))) {
        return create({ ...base }, recipeWrapper) as T;
    }
    return create(base, recipeWrapper) as T;
}