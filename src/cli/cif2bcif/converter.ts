/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { CIF, CifCategory, getCifFieldType, CifField, CifFile } from '../../mol-io/reader/cif';
import { CifWriter, EncodingStrategyHint } from '../../mol-io/writer/cif';
import * as util from 'util';
import * as fs from 'fs';
import * as zlib from 'zlib';
import { Progress, Task, RuntimeContext } from '../../mol-task';
import { classifyFloatArray, classifyIntArray } from '../../mol-io/common/binary-cif';
import { BinaryEncodingProvider } from '../../mol-io/writer/cif/encoder/binary';
import { Category } from '../../mol-io/writer/cif/encoder';
import { ReaderResult } from '../../mol-io/reader/result';

function showProgress(p: Progress) {
    process.stdout.write(`\r${new Array(80).join(' ')}`);
    process.stdout.write(`\r${Progress.format(p)}`);
}

const readFileAsync = util.promisify(fs.readFile);
const unzipAsync = util.promisify<zlib.InputType, Buffer>(zlib.unzip);

async function readFile(ctx: RuntimeContext, filename: string): Promise<ReaderResult<CifFile>> {
    const isGz = /\.gz$/i.test(filename);
    if (filename.match(/\.bcif/)) {
        let input = await readFileAsync(filename);
        if (isGz) input = await unzipAsync(input);
        return await CIF.parseBinary(new Uint8Array(input)).runInContext(ctx);
    } else {
        let str: string;
        if (isGz) {
            const data = await unzipAsync(await readFileAsync(filename));
            str = data.toString('utf8');
        } else {
            str = await readFileAsync(filename, 'utf8');
        }
        return await CIF.parseText(str).runInContext(ctx);
    }
}

async function getCIF(ctx: RuntimeContext, filename: string) {
    const parsed = await readFile(ctx, filename);
    if (parsed.isError) {
        throw new Error(parsed.toString());
    }
    return parsed.result;
}

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
        const encoder = classifyFloatArray(field.toFloatArray({ array: Float64Array }));
        return CifWriter.Field.float(name, field.float, { valueKind: field.valueKind, encoder, typedArray: Float64Array });
    } else {
        const encoder = classifyIntArray(field.toIntArray({ array: Int32Array }));
        return CifWriter.Field.int(name, field.int, { valueKind: field.valueKind, encoder, typedArray: Int32Array });
    }
}

export default function convert(path: string, asText = false, hints?: EncodingStrategyHint[], filter?: string) {
    return Task.create<Uint8Array>('BinaryCIF', async ctx => {
        const encodingProvider: BinaryEncodingProvider = hints
            ? CifWriter.createEncodingProviderFromJsonConfig(hints)
            : { get: (c, f) => void 0 };
        const cif = await getCIF(ctx, path);

        const encoder = CifWriter.createEncoder({
            binary: !asText,
            encoderName: 'mol*/ciftools cif2bcif',
            binaryAutoClassifyEncoding: true,
            binaryEncodingPovider: encodingProvider
        });

        if (filter) {
            encoder.setFilter(Category.filterOf(filter));
        }

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
                    fields.push(classify(f, cat.getField(f)!));
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
        await ctx.update('Done.\n');
        return ret;
    }).run(showProgress, 250);
}