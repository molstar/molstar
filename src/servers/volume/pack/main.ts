/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Format from './format'
import * as File from '../common/file'
import * as Data from './data-model'
import * as Sampling from './sampling'
import * as DataFormat from '../common/data-format'
import { FileHandle } from 'mol-io/common/file-handle';

export default async function pack(input: { name: string, filename: string }[], blockSize: number, isPeriodic: boolean, outputFilename: string) {
    try {
        await create(outputFilename, input, blockSize, isPeriodic);
    } catch (e) {
        console.error('[Error] ' + e);
    }
}

function getTime() {
    let t = process.hrtime();
    return t[0] * 1000 + t[1] / 1000000;
}

function updateAllocationProgress(progress: Data.Progress, progressDone: number) {
    let old = (100 * progress.current / progress.max).toFixed(0);
    progress.current += progressDone;
    let $new = (100 * progress.current / progress.max).toFixed(0);
    if (old !== $new) {
        process.stdout.write(`\rAllocating...      ${$new}%`);
    }
}

/**
 * Pre allocate the disk space to be able to do "random" writes into the entire file.
 */
async function allocateFile(ctx: Data.Context) {
    const { totalByteSize, file } = ctx;
    const buffer = new Buffer(Math.min(totalByteSize, 8 * 1024 * 1024));
    const progress: Data.Progress = { current: 0, max: Math.ceil(totalByteSize / buffer.byteLength) };
    let written = 0;
    while (written < totalByteSize) {
        written += file.writeBufferSync(written, buffer, Math.min(totalByteSize - written, buffer.byteLength));
        updateAllocationProgress(progress, 1);
    }
}

function determineBlockSize(data: Format.Data, blockSize: number) {
    const { extent } = data.header;
    const maxLayerSize = 1024 * 1024 * 1024;
    const valueCount = extent[0] * extent[1];
    if (valueCount * blockSize <= maxLayerSize) return blockSize;

    while (blockSize > 0) {
        blockSize -= 4;
        if (valueCount * blockSize <= maxLayerSize) return blockSize;
    }

    throw new Error('Could not determine a valid block size.');
}

async function writeHeader(ctx: Data.Context) {
    const header = DataFormat.encodeHeader(Data.createHeader(ctx));
    await File.writeInt(ctx.file, header.byteLength, 0);
    await ctx.file.writeBuffer(4, header);
}

async function create(filename: string, sourceDensities: { name: string, filename: string }[], sourceBlockSize: number, isPeriodic: boolean) {
    const startedTime = getTime();

    if (sourceBlockSize % 4 !== 0 || sourceBlockSize < 4) {
        throw Error('Block size must be a positive number divisible by 4.');
    }

    if (!sourceDensities.length) {
        throw Error('Specify at least one source density.');
    }

    process.stdout.write('Initializing... ');
    const files: FileHandle[] = [];
    try {
        // Step 1a: Read the Format headers
        const channels: Format.Context[] = [];
        for (const s of sourceDensities) {
            channels.push(await Format.open(s.name, s.filename));
        }
        // Step 1b: Check if the Format headers are compatible.
        const isOk = channels.reduce((ok, s) => ok && Format.compareHeaders(channels[0].data.header, s.data.header), true);
        if (!isOk) {
            throw new Error('Input file headers are not compatible (different grid, etc.).');
        }
        const blockSize = determineBlockSize(channels[0].data, sourceBlockSize);
        for (const ch of channels) Format.assignSliceBuffer(ch.data, blockSize);

        // Step 1c: Create data context.
        const context = await Sampling.createContext(filename, channels, blockSize, isPeriodic);
        for (const s of channels) files.push(s.data.file);
        files.push(context.file);
        process.stdout.write('   done.\n');

        console.log(`Block size: ${blockSize}`);

        // Step 2: Allocate disk space.
        process.stdout.write('Allocating...      0%');
        await allocateFile(context);
        process.stdout.write('\rAllocating...      done.\n');

        // Step 3: Process and write the data
        process.stdout.write('Writing data...    0%');
        await Sampling.processData(context);
        process.stdout.write('\rWriting data...    done.\n');

        // Step 4: Write the header at the start of the file.
        // The header is written last because the sigma/min/max values are computed
        // during step 3.
        process.stdout.write('Writing header...  ');
        await writeHeader(context);
        process.stdout.write('done.\n');

        // Step 5: Report the time, d'ph.
        const time = getTime() - startedTime;
        console.log(`[Done] ${time.toFixed(0)}ms.`);
    } finally {
        for (let f of files) f.close();

        // const ff = await File.openRead(filename);
        // const hh = await DataFormat.readHeader(ff);
        // File.close(ff);
        // console.log(hh.header);
    }
}