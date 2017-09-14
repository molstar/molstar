/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Block } from './block'
import { TextCategory } from './text-category'

export class TextBlock extends Block<TextCategory> {
    private categoryMap: Map<string, TextCategory>;
    private categoryList: TextCategory[];

    data: string;

    /**
     * Gets a category by its name.
     */
    getCategory(name: string) {
        return this.categoryMap.get(name);
    }

    /**
     * Adds a category.
     */
    addCategory(category: TextCategory) {
        this.categoryList[this.categoryList.length] = category;
        this.categoryMap.set(category.name, category);
    }

    constructor(data: string) {
        super()

        this.data = data;

        this.categoryMap = new Map()
        this.categoryList = []
    }
}