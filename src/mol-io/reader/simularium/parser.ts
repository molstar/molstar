/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Parser for Simularium trajectory files (JSON and binary variants).
 * The binary layout mirrors the reference implementation in
 * https://github.com/simularium/simularium-viewer (BinaryFileReader.ts).
 */

import { Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';
import { SimulariumFile, SimulariumFrame, SimulariumTrajectoryInfo } from './schema';

/** 16-byte ASCII signature at the start of a binary Simularium file. */
const SIGNATURE = 'SIMULARIUMBINARY';
/** Size in bytes of a block header (`type` int + `size` int). */
const BLOCK_HEADER_SIZE = 8;
/** Number of 32-bit ints in the file header before the block table-of-contents. */
const OFFSET_TO_TABLE_OF_CONTENTS = 3;

const enum BlockType {
    SpatialDataJson = 0,
    TrajectoryInfoJson = 1,
    PlotDataJson = 2,
    SpatialDataBinary = 3,
}

function hasBinarySignature(data: Uint8Array): boolean {
    if (data.length < SIGNATURE.length) return false;
    for (let i = 0; i < SIGNATURE.length; ++i) {
        if (data[i] !== SIGNATURE.charCodeAt(i)) return false;
    }
    return true;
}

function decodeUtf8(data: Uint8Array, byteOffset: number, byteLength: number): string {
    const view = new Uint8Array(data.buffer, data.byteOffset + byteOffset, byteLength);
    // Trim any trailing null padding bytes before decoding.
    return new TextDecoder('utf-8').decode(view).replace(/\0+$/, '');
}

function framesFromBundleData(bundleData: any): SimulariumFrame[] {
    if (!Array.isArray(bundleData)) {
        throw new Error('Simularium spatial data has no bundleData array.');
    }
    const frames: SimulariumFrame[] = [];
    for (let i = 0; i < bundleData.length; ++i) {
        const bundle = bundleData[i];
        const values: ArrayLike<number> = Array.isArray(bundle?.data) ? bundle.data : [];
        frames.push({
            frameNumber: typeof bundle?.frameNumber === 'number' ? bundle.frameNumber : i,
            time: typeof bundle?.time === 'number' ? bundle.time : 0,
            data: Float32Array.from(values),
        });
    }
    return frames;
}

function parseJson(data: Uint8Array): Result<SimulariumFile> {
    const text = decodeUtf8(data, 0, data.byteLength);
    const json = JSON.parse(text);
    console.log('Simularium JSON', json);

    const trajectoryInfo = json?.trajectoryInfo as SimulariumTrajectoryInfo | undefined;
    if (!trajectoryInfo) {
        return Result.error<SimulariumFile>('Simularium JSON is missing the trajectoryInfo object.');
    }
    const frames = framesFromBundleData(json?.spatialData?.bundleData);
    return Result.success<SimulariumFile>({ trajectoryInfo, frames });
}

interface BlockInfo { offset: number, type: number, size: number }

function readHeaderBlocks(view: DataView): BlockInfo[] {
    let p = SIGNATURE.length;
    const headerLength = view.getUint32(p + 0, true);
    const nBlocks = view.getUint32(p + 8, true);
    if (nBlocks < 1) throw new Error('Simularium binary file has no blocks.');

    const blocks: BlockInfo[] = [];
    p += OFFSET_TO_TABLE_OF_CONTENTS * 4;
    for (let i = 0; i < nBlocks; ++i) {
        blocks.push({
            offset: view.getUint32(p + 0, true),
            type: view.getUint32(p + 4, true),
            size: view.getUint32(p + 8, true),
        });
        p += 12;
    }
    if (blocks[0].offset !== headerLength) {
        throw new Error('Simularium binary header length does not match the first block offset.');
    }
    return blocks;
}

function readSpatialDataBinary(view: DataView, block: BlockInfo): SimulariumFrame[] {
    // Block content begins after the block header; frame offsets are stored relative
    // to the start of the block (i.e. including the block header).
    let p = block.offset + BLOCK_HEADER_SIZE;
    p += 4; // skip spatial data version
    const nFrames = view.getUint32(p, true);
    p += 4;

    const frameOffsets = new Int32Array(nFrames);
    const frameLengths = new Int32Array(nFrames);
    for (let i = 0; i < nFrames; ++i) {
        frameOffsets[i] = view.getUint32(p + 0, true);
        frameLengths[i] = view.getUint32(p + 4, true);
        p += 8;
    }

    const frames: SimulariumFrame[] = [];
    for (let i = 0; i < nFrames; ++i) {
        const frameStart = block.offset + frameOffsets[i];
        const frameNumber = view.getUint32(frameStart + 0, true);
        const time = view.getFloat32(frameStart + 4, true);
        // frameStart + 8 holds nAgents; the remaining bytes are the agent buffer.
        const agentByteStart = frameStart + 12;
        const agentFloatCount = Math.max(0, (frameLengths[i] - 12) >> 2);
        const agentData = new Float32Array(agentFloatCount);
        for (let k = 0; k < agentFloatCount; ++k) {
            agentData[k] = view.getFloat32(agentByteStart + k * 4, true);
        }
        frames.push({ frameNumber, time, data: agentData });
    }
    return frames;
}

function parseBinary(data: Uint8Array): Result<SimulariumFile> {
    const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
    const blocks = readHeaderBlocks(view);

    let trajectoryInfo: SimulariumTrajectoryInfo | undefined;
    let frames: SimulariumFrame[] | undefined;

    for (const block of blocks) {
        const blockType = view.getUint32(block.offset, true);
        if (blockType !== block.type) {
            return Result.error<SimulariumFile>(`Simularium block type mismatch: header ${block.type}, block ${blockType}.`);
        }
        switch (blockType) {
            case BlockType.TrajectoryInfoJson: {
                if (trajectoryInfo) break;
                const text = decodeUtf8(data, block.offset + BLOCK_HEADER_SIZE, block.size - BLOCK_HEADER_SIZE);
                trajectoryInfo = JSON.parse(text);
                break;
            }
            case BlockType.SpatialDataBinary: {
                if (!frames) frames = readSpatialDataBinary(view, block);
                break;
            }
            case BlockType.SpatialDataJson: {
                if (!frames) {
                    const text = decodeUtf8(data, block.offset + BLOCK_HEADER_SIZE, block.size - BLOCK_HEADER_SIZE);
                    frames = framesFromBundleData(JSON.parse(text)?.bundleData);
                }
                break;
            }
            case BlockType.PlotDataJson:
                // TODO: implement plot data parsing if needed
                break;
        }
    }

    if (!trajectoryInfo) return Result.error<SimulariumFile>('No Simularium trajectory info block found.');
    if (!frames) return Result.error<SimulariumFile>('No Simularium spatial data block found.');
    console.log('Simularium Binary', { trajectoryInfo, frames });
    return Result.success<SimulariumFile>({ trajectoryInfo, frames });
}

export function parseSimularium(data: Uint8Array) {
    return Task.create<Result<SimulariumFile>>('Parse Simularium', async () => {
        try {
            return hasBinarySignature(data) ? parseBinary(data) : parseJson(data);
        } catch (e) {
            return Result.error<SimulariumFile>(`Failed to parse Simularium file: ${e}`);
        }
    });
}
