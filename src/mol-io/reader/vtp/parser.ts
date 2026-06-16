/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext } from '../../../mol-task';
import { inflate } from '../../../mol-util/zip/zip';
import { VtpFile, VtpDataArrayDescriptor, VtpScalarArray } from './schema';

// --- Locate AppendedData section ---

interface AppendedInfo {
    /** Raw base64 (or raw binary) text/bytes after the '_' marker */
    rawB64: string;
    /** Encoding: "base64" or "raw" */
    encoding: string;
    /** Byte position of '_' + 1 in the original Uint8Array (for raw mode) */
    rawByteStart: number;
    /** Original Uint8Array (for raw mode) */
    rawData: Uint8Array;
}

function findBytesPosition(data: Uint8Array, tag: string, from = 0): number {
    const tagBytes = tag.split('').map(c => c.charCodeAt(0));
    const tlen = tagBytes.length;
    for (let i = from, il = data.length - tlen; i < il; i++) {
        let match = true;
        for (let j = 0; j < tlen; j++) {
            if (data[i + j] !== tagBytes[j]) { match = false; break; }
        }
        if (match) return i;
    }
    return -1;
}

function extractAppendedInfo(data: Uint8Array, headerText: string): AppendedInfo {
    const encMatch = headerText.match(/<AppendedData\s+encoding="([^"]+)"/i);
    const encoding = encMatch ? encMatch[1].toLowerCase() : 'raw';

    const tagPos = findBytesPosition(data, '<AppendedData');
    if (tagPos === -1) throw new Error('AppendedData section not found');

    let pos = tagPos;
    while (pos < data.length && data[pos] !== 0x3E) pos++;
    pos++; // skip '>'
    while (pos < data.length && (data[pos] === 0x20 || data[pos] === 0x09 || data[pos] === 0x0A || data[pos] === 0x0D)) pos++;
    if (data[pos] !== 0x5F) {
        throw new Error(`Expected '_' before appended data, got 0x${data[pos].toString(16)}`);
    }
    const dataStart = pos + 1;

    if (encoding === 'base64') {
        const endTagPos = findBytesPosition(data, '</AppendedData>', dataStart);
        if (endTagPos < 0) throw new Error('</AppendedData> not found');
        const rawB64 = new TextDecoder('ascii').decode(
            new Uint8Array(data.buffer, data.byteOffset + dataStart, endTagPos - dataStart)
        );
        return { rawB64, encoding, rawByteStart: dataStart, rawData: data };
    } else {
        // raw: offsets are byte offsets
        return { rawB64: '', encoding: 'raw', rawByteStart: dataStart, rawData: data };
    }
}

// --- XML header parsing ---

function attrValue(attrStr: string, name: string): string | undefined {
    const m = attrStr.match(new RegExp(name + '=["\']([^"\']*)["\']'));
    return m ? m[1] : undefined;
}

function parseDataArrayAttrs(attrStr: string): VtpDataArrayDescriptor | null {
    const offsetStr = attrValue(attrStr, 'offset');
    if (offsetStr === undefined) return null;
    return {
        name: attrValue(attrStr, 'Name') ?? '',
        type: attrValue(attrStr, 'type') ?? 'Float32',
        numberOfComponents: parseInt(attrValue(attrStr, 'NumberOfComponents') ?? '1'),
        offset: parseInt(offsetStr),
        rangeMin: parseFloat(attrValue(attrStr, 'RangeMin') ?? '0'),
        rangeMax: parseFloat(attrValue(attrStr, 'RangeMax') ?? '1'),
    };
}

function extractSectionArrays(header: string, sectionName: string): VtpDataArrayDescriptor[] {
    const re = new RegExp(`<${sectionName}(?:\\s[^>]*)?>([\\s\\S]*?)</${sectionName}>`, 'i');
    const sm = header.match(re);
    if (!sm) return [];
    const results: VtpDataArrayDescriptor[] = [];
    const daRe = /<DataArray\s+([^>]*?)\/?>/g;
    let m;
    while ((m = daRe.exec(sm[1])) !== null) {
        const desc = parseDataArrayAttrs(m[1]);
        if (desc) results.push(desc);
    }
    return results;
}

interface ParsedHeader {
    nPoints: number;
    nCells: number;
    compressor: string;
    headerType: string;
    pointDataArrays: VtpDataArrayDescriptor[];
    cellDataArrays: VtpDataArrayDescriptor[];
    pointsArrays: VtpDataArrayDescriptor[];
    polysArrays: VtpDataArrayDescriptor[];
}

function parseHeader(header: string): ParsedHeader {
    const pieceMatch = header.match(/NumberOfPoints="(\d+)"[^>]*NumberOfPolys="(\d+)"/);
    if (!pieceMatch) throw new Error('Cannot parse NumberOfPoints/NumberOfPolys from VTP header');
    const vtkFileMatch = header.match(/<VTKFile\s+([^>]*)>/);
    const vtkAttrs = vtkFileMatch ? vtkFileMatch[1] : '';
    return {
        nPoints: parseInt(pieceMatch[1]),
        nCells: parseInt(pieceMatch[2]),
        compressor: attrValue(vtkAttrs, 'compressor') ?? '',
        headerType: attrValue(vtkAttrs, 'header_type') ?? 'UInt32',
        pointDataArrays: extractSectionArrays(header, 'PointData'),
        cellDataArrays: extractSectionArrays(header, 'CellData'),
        pointsArrays: extractSectionArrays(header, 'Points'),
        polysArrays: extractSectionArrays(header, 'Polys'),
    };
}

// --- Base64 helpers for VTK appended format ---
//
// VTK base64 appended: each DataArray's binary block is split into two base64-encoded chunks:
//   1. The compressed-block header  (nblocks, blockSize, lastSize, compSize[i])
//   2. All compressed blocks concatenated
//
// The DataArray XML `offset` attribute is the CHARACTER position of the first chunk in
// the raw base64 text (after '_').  Knowing the number of bytes to decode lets us compute
// the exact character span (ceil(nBytes/3)*4) without scanning for '='.

function decodeBase64Slice(rawB64: string, charPos: number, nBytes: number): Uint8Array {
    if (nBytes === 0) return new Uint8Array(0);
    const nChars = Math.ceil(nBytes / 3) * 4;
    const slice = rawB64.slice(charPos, charPos + nChars);
    const binary = atob(slice);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
    return bytes;
}

// --- Zlib block decompression (base64 path) ---

async function decompressVtkBlockB64(
    ctx: RuntimeContext,
    rawB64: string,
    charOffset: number
): Promise<Uint8Array> {
    // 1. Read first 12 bytes to get nblocks, blockSize, lastSize
    const minHdr = decodeBase64Slice(rawB64, charOffset, 12);
    const dvMin = new DataView(minHdr.buffer, minHdr.byteOffset, 12);
    const nblocks = dvMin.getUint32(0, true);
    const blockSize = dvMin.getUint32(4, true);
    const lastSize = dvMin.getUint32(8, true);

    // 2. Read full header (adds nblocks * 4 bytes for compressed sizes)
    const totalHdrBytes = 12 + nblocks * 4;
    const hdrRaw = decodeBase64Slice(rawB64, charOffset, totalHdrBytes);
    const hdrChars = Math.ceil(totalHdrBytes / 3) * 4;

    if (nblocks === 0) return new Uint8Array(0);

    const dvHdr = new DataView(hdrRaw.buffer, hdrRaw.byteOffset, hdrRaw.byteLength);
    const compSizes = new Uint32Array(nblocks);
    let totalCompressed = 0;
    for (let i = 0; i < nblocks; i++) {
        compSizes[i] = dvHdr.getUint32(12 + i * 4, true);
        totalCompressed += compSizes[i];
    }

    // 3. Decode all compressed blocks as one contiguous base64 slice
    const dataCharOffset = charOffset + hdrChars;
    const compData = decodeBase64Slice(rawB64, dataCharOffset, totalCompressed);

    // 4. Decompress each zlib block
    const totalUncompressed = lastSize > 0
        ? (nblocks - 1) * blockSize + lastSize
        : nblocks * blockSize;

    const result = new Uint8Array(totalUncompressed);
    let compOffset = 0;
    let outOffset = 0;

    for (let i = 0; i < nblocks; i++) {
        const csize = compSizes[i];
        const blockView = new Uint8Array(
            compData.buffer, compData.byteOffset + compOffset, csize
        ) as Uint8Array<ArrayBuffer>;
        const inflated = await inflate(ctx, blockView);
        result.set(inflated, outOffset);
        outOffset += inflated.length;
        compOffset += csize;
    }

    return result;
}

// --- Zlib block decompression (raw path) ---

async function decompressVtkBlockRaw(
    ctx: RuntimeContext,
    rawData: Uint8Array,
    rawByteStart: number,
    byteOffset: number
): Promise<Uint8Array> {
    const absBase = rawByteStart + byteOffset;
    const dv = new DataView(rawData.buffer, rawData.byteOffset, rawData.byteLength);

    const nblocks = dv.getUint32(absBase, true);
    const blockSize = dv.getUint32(absBase + 4, true);
    const lastSize = dv.getUint32(absBase + 8, true);

    if (nblocks === 0) return new Uint8Array(0);

    const hdrBytes = 12 + nblocks * 4;
    const compSizes = new Uint32Array(nblocks);
    for (let i = 0; i < nblocks; i++) {
        compSizes[i] = dv.getUint32(absBase + 12 + i * 4, true);
    }

    const totalUncompressed = lastSize > 0
        ? (nblocks - 1) * blockSize + lastSize
        : nblocks * blockSize;

    const result = new Uint8Array(totalUncompressed);
    let compOff = absBase + hdrBytes;
    let outOff = 0;

    for (let i = 0; i < nblocks; i++) {
        const csize = compSizes[i];
        const blockView = new Uint8Array(
            rawData.buffer, rawData.byteOffset + compOff, csize
        ) as Uint8Array<ArrayBuffer>;
        const inflated = await inflate(ctx, blockView);
        result.set(inflated, outOff);
        outOff += inflated.length;
        compOff += csize;
    }
    return result;
}

// --- Dispatch decompression by encoding ---

async function decompressBlock(
    ctx: RuntimeContext,
    info: AppendedInfo,
    desc: VtpDataArrayDescriptor
): Promise<Uint8Array> {
    if (info.encoding === 'base64') {
        return decompressVtkBlockB64(ctx, info.rawB64, desc.offset);
    } else {
        return decompressVtkBlockRaw(ctx, info.rawData, info.rawByteStart, desc.offset);
    }
}

// --- Typed array decoding ---

function decodeFloat32(raw: Uint8Array, count: number): Float32Array {
    return new Float32Array(raw.buffer, raw.byteOffset, count);
}

function decodeInt64AsInt32(raw: Uint8Array, count: number): Int32Array {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const out = new Int32Array(count);
    for (let i = 0; i < count; i++) out[i] = dv.getInt32(i * 8, true);
    return out;
}

function decodeToFloat64(raw: Uint8Array, type: string, count: number): Float64Array {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const out = new Float64Array(count);
    switch (type) {
        case 'Float32': { const f32 = decodeFloat32(raw, count); for (let i = 0; i < count; i++) out[i] = f32[i]; break; }
        case 'Float64': for (let i = 0; i < count; i++) out[i] = dv.getFloat64(i * 8, true); break;
        case 'Int32': for (let i = 0; i < count; i++) out[i] = dv.getInt32(i * 4, true); break;
        case 'UInt32': for (let i = 0; i < count; i++) out[i] = dv.getUint32(i * 4, true); break;
        case 'Int64': case 'UInt64':
            for (let i = 0; i < count; i++) {
                const lo = dv.getUint32(i * 8, true);
                const hi = dv.getInt32(i * 8 + 4, true);
                out[i] = hi * 0x100000000 + lo;
            }
            break;
        default: for (let i = 0; i < count; i++) out[i] = dv.getInt32(i * 4, true); break;
    }
    return out;
}

// --- Fan-triangulate connectivity from VTP Polys offsets ---
// offsets[i] = cumulative vertex count through cell i; cell i uses
// rawConn[offsets[i-1]..offsets[i]-1] (with offsets[-1] = 0).

function buildTriangles(rawConn: Int32Array, offsets: Int32Array, nCells: number): Int32Array {
    let nTris = 0;
    for (let i = 0; i < nCells; i++) {
        const start = i > 0 ? offsets[i - 1] : 0;
        nTris += offsets[i] - start - 2;
    }
    const tris = new Int32Array(3 * nTris);
    let ti = 0;
    for (let i = 0; i < nCells; i++) {
        const start = i > 0 ? offsets[i - 1] : 0;
        const end = offsets[i];
        const v0 = rawConn[start];
        for (let j = start + 1; j < end - 1; j++) {
            tris[ti++] = v0;
            tris[ti++] = rawConn[j];
            tris[ti++] = rawConn[j + 1];
        }
    }
    return tris;
}

// --- Main parser ---

async function parseInternal(data: Uint8Array, ctx: RuntimeContext): Promise<Result<VtpFile>> {
    ctx.update({ message: 'Parsing VTP header...', current: 0, max: 100 });

    const previewLen = Math.min(data.length, 65536);
    const headerText = new TextDecoder('ascii').decode(data.subarray(0, previewLen));
    const hdr = parseHeader(headerText);

    if (hdr.compressor && !hdr.compressor.includes('ZLib') && !hdr.compressor.includes('zlib')) {
        throw new Error(`Unsupported VTP compressor: "${hdr.compressor}". Only vtkZLibDataCompressor is supported.`);
    }
    if (hdr.headerType !== 'UInt32') {
        throw new Error(`Unsupported VTP header_type: "${hdr.headerType}". Only UInt32 is supported.`);
    }

    ctx.update({ message: 'Locating binary section...', current: 5, max: 100 });
    const info = extractAppendedInfo(data, headerText);

    // Points positions
    const pointsDesc = hdr.pointsArrays.find(d => d.name === 'Points');
    if (!pointsDesc) throw new Error('VTP file missing Points DataArray');

    ctx.update({ message: 'Decompressing positions...', current: 10, max: 100 });
    const posRaw = await decompressBlock(ctx, info, pointsDesc);
    const positions = decodeFloat32(posRaw, hdr.nPoints * 3);

    // Polys connectivity + offsets
    const connDesc = hdr.polysArrays.find(d => d.name === 'connectivity');
    const offsetsDesc = hdr.polysArrays.find(d => d.name === 'offsets');
    if (!connDesc || !offsetsDesc) throw new Error('VTP file missing Polys connectivity/offsets DataArrays');

    ctx.update({ message: 'Decompressing topology...', current: 20, max: 100 });
    const connRaw = await decompressBlock(ctx, info, connDesc);
    const offsetsRaw = await decompressBlock(ctx, info, offsetsDesc);

    const rawConn = decodeInt64AsInt32(connRaw, connRaw.byteLength / 8);
    const rawOffsets = decodeInt64AsInt32(offsetsRaw, hdr.nCells);

    const connectivity = buildTriangles(rawConn, rawOffsets, hdr.nCells);
    const numberOfTriangles = connectivity.length / 3;

    // PointData scalar arrays
    const pointData = new Map<string, VtpScalarArray>();
    const scalarPD = hdr.pointDataArrays.filter(d => d.numberOfComponents === 1);
    for (let i = 0; i < scalarPD.length; i++) {
        const desc = scalarPD[i];
        ctx.update({ message: `Point data: ${desc.name}...`, current: 30 + Math.floor(i / Math.max(1, scalarPD.length) * 30), max: 100 });
        const raw = await decompressBlock(ctx, info, desc);
        const values = decodeToFloat64(raw, desc.type, hdr.nPoints);
        pointData.set(desc.name, { desc, values });
    }

    // CellData scalar arrays
    const cellData = new Map<string, VtpScalarArray>();
    const scalarCD = hdr.cellDataArrays.filter(d => d.numberOfComponents === 1);
    for (let i = 0; i < scalarCD.length; i++) {
        const desc = scalarCD[i];
        ctx.update({ message: `Cell data: ${desc.name}...`, current: 60 + Math.floor(i / Math.max(1, scalarCD.length) * 35), max: 100 });
        const raw = await decompressBlock(ctx, info, desc);
        const values = decodeToFloat64(raw, desc.type, hdr.nCells);
        cellData.set(desc.name, { desc, values });
    }

    return Result.success({
        numberOfPoints: hdr.nPoints,
        numberOfCells: hdr.nCells,
        positions,
        connectivity,
        numberOfTriangles,
        pointData,
        cellData,
    });
}

export function parseVtp(data: Uint8Array) {
    return Task.create<Result<VtpFile>>('Parse VTP', async ctx => {
        return await parseInternal(data, ctx);
    });
}
