/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifCategory, CifField, CifFrame, getCifFieldType } from '../../../mol-io/reader/cif';
import { CifWriter } from '../../../mol-io/writer/cif';
import { Task } from '../../../mol-task';
// import { showProgress } from './util';

function getCategoryInstanceProvider(cat: CifCategory, fields: CifWriter.Field[]): CifWriter.Category {
    return {
        name: cat.name,
        instance: () => CifWriter.categoryInstance(fields, { data: cat, rowCount: cat.rowCount })
    };
}

function classify(name: string, field: CifField): CifWriter.Field {
    const type = getCifFieldType(field);
    if (type['@type'] === 'str') {
        return { name, type: CifWriter.Field.Type.Str, value: field.str, valueKind: field.valueKind };
    } else if (type['@type'] === 'float') {
        return CifWriter.Field.float(name, field.float, { valueKind: field.valueKind, typedArray: Float64Array });
    } else {
        return CifWriter.Field.int(name, field.int, { valueKind: field.valueKind, typedArray: Int32Array });
    }
}

export function classifyCif(frame: CifFrame) {
    return Task.create('Classify CIF Data', async ctx => {
        let maxProgress = 0;
        for (const c of frame.categoryNames) maxProgress += frame.categories[c].fieldNames.length;

        const ret: CifWriter.Category[] = [];

        let current = 0;
        for (const c of frame.categoryNames) {
            const cat = frame.categories[c];
            const fields: CifWriter.Field[] = [];
            for (const f of cat.fieldNames) {
                const cifField = classify(f, cat.getField(f)!);
                fields.push(cifField);
                current++;
                if (ctx.shouldUpdate) await ctx.update({ message: 'Classifying...', current, max: maxProgress });
            }
            ret.push(getCategoryInstanceProvider(cat, fields));
        }
        return ret;
    }).run();
}