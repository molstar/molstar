/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import CIF, { CifCategory } from 'mol-io/reader/cif'
import { CifWriter } from 'mol-io/writer/cif'
import * as fs from 'fs'
import classify from './field-classifier'
import { Progress, Task, RuntimeContext } from 'mol-task';

function showProgress(p: Progress) {
    process.stdout.write(`\r${new Array(80).join(' ')}`);
    process.stdout.write(`\r${Progress.format(p)}`);
}

async function getCIF(ctx: RuntimeContext, path: string) {
    const str = fs.readFileSync(path, 'utf8');
    const parsed = await CIF.parseText(str).runInContext(ctx);
    if (parsed.isError) {
        throw new Error(parsed.toString());
    }
    return parsed.result;
}

function getCategoryInstanceProvider(cat: CifCategory, fields: CifWriter.Field[]): CifWriter.Category.Provider {
    return function (ctx: any) {
        return {
            data: cat,
            name: cat.name,
            fields,
            rowCount: cat.rowCount
        };
    }
}

export default function convert(path: string, asText = false) {
    return Task.create<Uint8Array>('BinaryCIF', async ctx => {
        const cif = await getCIF(ctx, path);

        const encoder = CifWriter.createEncoder({ binary: !asText, encoderName: 'mol* cif2bcif' });

        let maxProgress = 0;
        for (const b of cif.blocks) {
            maxProgress += b.categoryNames.length;
            for (const c of b.categoryNames) maxProgress += b.categories[c].fieldNames.length;
        }

        let current = 0;
        for (const b of cif.blocks) {
            encoder.startDataBlock(b.header);
            for (const c of b.categoryNames) {
                const cat = b.categories[c];
                const fields: CifWriter.Field[] = [];
                for (const f of cat.fieldNames) {
                    fields.push(classify(f, cat.getField(f)!))
                    current++;
                    if (ctx.shouldUpdate) await ctx.update({ message: 'Encoding...', current, max: maxProgress });
                }

                encoder.writeCategory(getCategoryInstanceProvider(b.categories[c], fields));
                current++;
                if (ctx.shouldUpdate) await ctx.update({ message: 'Encoding...', current, max: maxProgress });
            }
        }
        await ctx.update('Exporting...');
        const ret = encoder.getData() as Uint8Array;
        await ctx.update('Done.');
        return ret;
    }).run(showProgress, 250);
}