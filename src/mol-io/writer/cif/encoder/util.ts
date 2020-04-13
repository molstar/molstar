/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Iterator } from '../../../../mol-data';
import { Field, Category } from '../encoder';

export function getFieldDigitCount(field: Field) {
    if (field.defaultFormat && typeof field.defaultFormat.digitCount !== 'undefined') return Math.max(0, Math.min(field.defaultFormat.digitCount, 16));
    return 6;
}

export function getIncludedFields(category: Category.Instance) {
    return category.fields.some(f => !!f.shouldInclude)
        ? category.fields.filter(f => !f.shouldInclude || category.source.some(src => f.shouldInclude!(src.data)))
        : category.fields;
}

export interface CategoryInstanceData<Ctx = any> {
    instance: Category.Instance<Ctx>,
    rowCount: number,
    source: { data: any, keys: () => Iterator<any>, rowCount: number }[]
}

export function getCategoryInstanceData<Ctx>(category: Category<Ctx>, ctx?: Ctx): CategoryInstanceData<Ctx> {
    const instance = category.instance(ctx as any);
    let sources = instance.source.filter(s => s.rowCount > 0);
    if (!sources.length) return { instance, rowCount: 0, source: [] };

    const rowCount = sources.reduce((a, c) => a + c.rowCount, 0);
    const source = sources.map(c => ({
        data: c.data,
        keys: () => c.keys ? c.keys() : Iterator.Range(0, c.rowCount - 1),
        rowCount: c.rowCount
    }));

    return { instance, rowCount, source };
}