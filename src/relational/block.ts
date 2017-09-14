/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Category, UndefinedCategory } from './category'

export abstract class Block<T> {
    abstract getCategory(name: string): T|undefined
    abstract addCategory(category: T): void

    getCategoriesFromSchema<T extends object> (schema: T) {
        return BlockCategories(this, schema)
    }
}

export type BlockCategories<Categories extends string> = { readonly [name in Categories]: Category }
export function BlockCategories<T extends object>(block: Block<any> | undefined, categories: T): BlockCategories<keyof T> {
    const ret = Object.create(null);
    if (!block) for (const c of Object.keys(categories)) ret[c] = UndefinedCategory;
    else for (const c of Object.keys(categories)) ret[c] = block.getCategory(c);
    return ret;
}
