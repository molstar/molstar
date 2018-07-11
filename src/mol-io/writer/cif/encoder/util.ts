/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Field, Category } from '../encoder';

export function getFieldDigitCount(field: Field) {
    if (field.defaultFormat && typeof field.defaultFormat.digitCount !== 'undefined') return Math.max(0, Math.min(field.defaultFormat.digitCount, 16));
    return 6;
}

export function getIncludedFields(category: Category.Instance) {
    return category.fields.some(f => !!f.shouldInclude)
        ? category.fields.filter(f => !f.shouldInclude || f.shouldInclude(category.data))
        : category.fields;
}