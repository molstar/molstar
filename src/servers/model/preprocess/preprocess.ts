/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readStructureWrapper, resolveStructures, readDataAndFrame } from '../server/structure-wrapper';
import { classifyCif } from './converter';
import { Structure } from '../../../mol-model/structure';
import { CifWriter } from '../../../mol-io/writer/cif';
import Writer from '../../../mol-io/writer/writer';
import { encode_mmCIF_categories, CifExportContext } from '../../../mol-model/structure/export/mmcif';
import { ModelPropertiesProvider } from '../property-provider';
import { FileResultWriter } from '../utils/writer';

// TODO: error handling

export function preprocessFile(filename: string, propertyProvider?: ModelPropertiesProvider, outputCif?: string, outputBcif?: string) {
    return propertyProvider
        ? preprocess(filename, propertyProvider, outputCif, outputBcif)
        : convert(filename, outputCif, outputBcif);
}

async function preprocess(filename: string, propertyProvider?: ModelPropertiesProvider, outputCif?: string, outputBcif?: string) {
    const input = await readStructureWrapper('entry', '_local_', filename, void 0, propertyProvider);
    const categories = await classifyCif(input.cifFrame);
    const inputStructures = (await resolveStructures(input))!;
    const exportCtx = CifExportContext.create(inputStructures);

    if (outputCif) {
        const writer = new FileResultWriter(outputCif);
        const encoder = CifWriter.createEncoder({ binary: false });
        encode(inputStructures[0], input.cifFrame.header, categories, encoder, exportCtx, writer);
        writer.end();
    }

    if (outputBcif) {
        const writer = new FileResultWriter(outputBcif);
        const encoder = CifWriter.createEncoder({ binary: true, binaryAutoClassifyEncoding: true });
        encode(inputStructures[0], input.cifFrame.header, categories, encoder, exportCtx, writer);
        writer.end();
    }
}

async function convert(filename: string, outputCif?: string, outputBcif?: string) {
    const { frame } = await readDataAndFrame(filename);
    const categories = await classifyCif(frame);

    if (outputCif) {
        const writer = new FileResultWriter(outputCif);
        const encoder = CifWriter.createEncoder({ binary: false });
        encodeConvert(frame.header, categories, encoder, writer);
        writer.end();
    }

    if (outputBcif) {
        const writer = new FileResultWriter(outputBcif);
        const encoder = CifWriter.createEncoder({ binary: true, binaryAutoClassifyEncoding: true });
        encodeConvert(frame.header, categories, encoder, writer);
        writer.end();
    }
}

function encodeConvert(header: string, categories: CifWriter.Category[], encoder: CifWriter.Encoder, writer: Writer) {
    encoder.startDataBlock(header);
    for (const cat of categories) {
        encoder.writeCategory(cat);
    }
    encoder.encode();
    encoder.writeTo(writer);
}

function encode(structure: Structure, header: string, categories: CifWriter.Category[], encoder: CifWriter.Encoder, exportCtx: CifExportContext, writer: Writer) {
    const skipCategoryNames = new Set<string>(categories.map(c => c.name));
    encoder.startDataBlock(header);
    for (const cat of categories) {
        encoder.writeCategory(cat);
    }
    encode_mmCIF_categories(encoder, structure, { skipCategoryNames, exportCtx });
    encoder.encode();
    encoder.writeTo(writer);
}