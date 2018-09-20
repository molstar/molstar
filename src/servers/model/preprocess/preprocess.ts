/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readStructure } from '../server/structure-wrapper';
import { classifyCif } from './converter';
// import { ConsoleLogger } from 'mol-util/console-logger';
import { Structure } from 'mol-model/structure';
import { CifWriter } from 'mol-io/writer/cif';
import Writer from 'mol-io/writer/writer';
import { wrapFileToWriter } from '../server/api-local';
import { Task/*, now*/ } from 'mol-task';
import { /*showProgress, clearLine */ } from './util';
import { encode_mmCIF_categories, CifExportContext } from 'mol-model/structure/export/mmcif';

// TODO: error handling
// let linearId = 0;

export async function preprocessFile(filename: string, outputCif?: string, outputBcif?: string) {
    // linearId++;

    //const started = now();
    //ConsoleLogger.log(`${linearId}`, `Reading '${filename}'...`);
    // TODO: support the custom prop provider list here.
    const input = await readStructure('entry', '_local_', filename, void 0);
    //ConsoleLogger.log(`${linearId}`, `Classifying CIF categories...`);
    const categories = await classifyCif(input.cifFrame);
    //clearLine();

    const exportCtx = CifExportContext.create(input.structure, input.structure.models[0]);

    if (outputCif) {
        //ConsoleLogger.log(`${linearId}`, `Encoding CIF...`);
        const writer = wrapFileToWriter(outputCif);
        const encoder = CifWriter.createEncoder({ binary: false });
        await encode(input.structure, input.cifFrame.header, categories, encoder, exportCtx, writer);
        // clearLine();
        writer.end();
    }

    if (outputBcif) {
        // ConsoleLogger.log(`${linearId}`, `Encoding BinaryCIF...`);
        const writer = wrapFileToWriter(outputBcif);
        const encoder = CifWriter.createEncoder({ binary: true, binaryAutoClassifyEncoding: true });
        await encode(input.structure, input.cifFrame.header, categories, encoder, exportCtx, writer);
        //clearLine();
        writer.end();
    }
    // ConsoleLogger.log(`${linearId}`, `Finished '${filename}' in ${Math.round(now() - started)}ms`);
}

function encode(structure: Structure, header: string, categories: CifWriter.Category[], encoder: CifWriter.Encoder, exportCtx: CifExportContext, writer: Writer) {
    return Task.create('Encode', async ctx => {
        const skipCategoryNames = new Set<string>(categories.map(c => c.name));
        encoder.startDataBlock(header);
        // let current = 0;
        for (const cat of categories){
            encoder.writeCategory(cat);
            // current++;
            // if (ctx.shouldUpdate) await ctx.update({ message: 'Encoding...', current, max: categories.length });
        }
        encode_mmCIF_categories(encoder, structure, { skipCategoryNames, exportCtx });
        encoder.encode();
        encoder.writeTo(writer);
    }).run();
}