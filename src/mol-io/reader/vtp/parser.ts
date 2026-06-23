/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext } from '../../../mol-task';
import { inflate } from '../../../mol-util/zip/zip';
import { VtpFile, VtpDataArrayDescriptor, VtpScalarArray } from './schema';
import { parseInt as fastParseInt, parseFloat as fastParseFloat } from '../common/text/number-parser';
import { StringLike } from '../../common/string-like';
import { decodeB64Bytes, decodeB64Str } from '../common/decode';
import { TypedArray } from '../../../mol-util/type-helpers';

// --- Locate AppendedData section ---

interface AppendedInfo {
    encoding: string;
    /** Byte position of '_' + 1 in the original Uint8Array */
    rawByteStart: number;
    rawData: Uint8Array;
}


// --- Typed-array decoders ---

function decodeFloat32(raw: Uint8Array, count: number): Float32Array {
    if (raw.byteOffset % 4 === 0) return new Float32Array(raw.buffer, raw.byteOffset, count);
    const aligned = raw.slice(0, count * 4);
    return new Float32Array(aligned.buffer, 0, count);
}

function decodeInt64AsInt32(raw: Uint8Array, count: number): Int32Array {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const out = new Int32Array(count);
    for (let i = 0; i < count; i++) out[i] = dv.getUint32(i * 8, true);
    return out;
}

function decodeTyped(raw: Uint8Array, type: string, count: number): TypedArray {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    switch (type) {
        case 'Float32': return decodeFloat32(raw, count);
        case 'Float64': { const out = new Float64Array(count); for (let i = 0; i < count; i++) out[i] = dv.getFloat64(i * 8, true); return out; }
        case 'Int8': { const out = new Int8Array(count); for (let i = 0; i < count; i++) out[i] = dv.getInt8(i); return out; }
        case 'UInt8': { const out = new Uint8Array(count); for (let i = 0; i < count; i++) out[i] = dv.getUint8(i); return out; }
        case 'Int16': { const out = new Int16Array(count); for (let i = 0; i < count; i++) out[i] = dv.getInt16(i * 2, true); return out; }
        case 'UInt16': { const out = new Uint16Array(count); for (let i = 0; i < count; i++) out[i] = dv.getUint16(i * 2, true); return out; }
        case 'Int32': { const out = new Int32Array(count); for (let i = 0; i < count; i++) out[i] = dv.getInt32(i * 4, true); return out; }
        case 'UInt32': { const out = new Uint32Array(count); for (let i = 0; i < count; i++) out[i] = dv.getUint32(i * 4, true); return out; }
        case 'Int64': case 'UInt64': {
            const out = new Float64Array(count);
            for (let i = 0; i < count; i++) {
                const lo = dv.getUint32(i * 8, true);
                const hi = dv.getInt32(i * 8 + 4, true);
                out[i] = hi * 0x100000000 + lo;
            }
            return out;
        }
        default: throw new Error(`Unsupported scalar type: "${type}".`);
    }
}

function decodePositions(raw: Uint8Array, type: string, count: number): Float32Array {
    if (type === 'Float32') return decodeFloat32(raw, count);
    return new Float32Array(decodeTyped(raw, type, count));
}

function decodeConnectivity(raw: Uint8Array, type: string): Int32Array {
    if (type === 'Int64' || type === 'UInt64') {
        return decodeInt64AsInt32(raw, raw.byteLength / 8);
    }
    if (type === 'Int32' || type === 'UInt32') {
        if (raw.byteOffset % 4 === 0) return new Int32Array(raw.buffer, raw.byteOffset, raw.byteLength / 4);
        const aligned = raw.slice();
        return new Int32Array(aligned.buffer, 0, aligned.byteLength / 4);
    }
    throw new Error(`Unsupported connectivity type: "${type}". Expected Int32, UInt32, Int64, or UInt64.`);
}

function findBytesPosition(data: Uint8Array, tag: string, from = 0): number {
    const tagBytes = tag.split('').map(c => c.charCodeAt(0));
    const tlen = tagBytes.length;
    for (let i = from, il = data.length - tlen; i <= il; i++) {
        let match = true;
        for (let j = 0; j < tlen; j++) {
            if (data[i + j] !== tagBytes[j]) { match = false; break; }
        }
        if (match) return i;
    }
    return -1;
}

function extractAppendedInfo(data: Uint8Array): AppendedInfo {
    const tagPos = findBytesPosition(data, '<AppendedData');
    if (tagPos === -1) throw new Error('AppendedData section not found');

    let pos = tagPos;
    while (pos < data.length && data[pos] !== 0x3E) pos++;
    const tagText = new TextDecoder('ascii').decode(data.subarray(tagPos, pos + 1));
    const { attrs: appendedAttrs } = readAttrs(tagText, '<AppendedData'.length);
    const encoding = (appendedAttrs.get('encoding') ?? 'raw').toLowerCase();
    pos++; // skip '>'
    while (pos < data.length && (data[pos] === 0x20 || data[pos] === 0x09 || data[pos] === 0x0A || data[pos] === 0x0D)) pos++;
    if (pos >= data.length || data[pos] !== 0x5F) {
        throw new Error(`Expected '_' before appended data, got 0x${data[pos]?.toString(16) ?? 'EOF'}`);
    }
    const dataStart = pos + 1;

    return { encoding, rawByteStart: dataStart, rawData: data };
}

// --- Tokenizer helpers (replace regex-based XML attribute parsing) ---

function isWS(c: number): boolean {
    return c === 32 || c === 9 || c === 10 || c === 13;
}

/**
 * Find '<tagName' followed by whitespace, '>', or '/' starting from `from`.
 * The char-after guard prevents 'Points' from matching '<PointData>'.
 */
function findTagStart(src: string, tagName: string, from: number): number {
    const needle = '<' + tagName;
    const nl = needle.length;
    outer: for (let i = from, il = src.length - nl; i <= il; i++) {
        if (src.charCodeAt(i) !== 60) continue;
        for (let j = 1; j < nl; j++) {
            if (src.charCodeAt(i + j) !== needle.charCodeAt(j)) continue outer;
        }
        const after = src.charCodeAt(i + nl);
        if (isWS(after) || after === 62 || after === 47) return i;
    }
    return -1;
}

/**
 * Parse attribute name="value" pairs from the body of an opening tag (the part after
 * the tag name).  Returns attrs and the position of the closing '>' or '/'.
 */
function readAttrs(src: string, pos: number): { attrs: Map<string, string>, tagEnd: number } {
    const attrs = new Map<string, string>();
    const len = src.length;
    while (pos < len) {
        while (pos < len && isWS(src.charCodeAt(pos))) pos++;
        const c = src.charCodeAt(pos);
        if (c === 62 || c === 47) break; // '>' or '/'
        const ns = pos;
        while (pos < len && src.charCodeAt(pos) !== 61 && !isWS(src.charCodeAt(pos)) && src.charCodeAt(pos) !== 62 && src.charCodeAt(pos) !== 47) pos++;
        const name = src.slice(ns, pos);
        while (pos < len && isWS(src.charCodeAt(pos))) pos++;
        if (pos < len && src.charCodeAt(pos) === 61) pos++; // '='
        while (pos < len && isWS(src.charCodeAt(pos))) pos++;
        if (pos < len && (src.charCodeAt(pos) === 34 || src.charCodeAt(pos) === 39)) {
            const q = src.charCodeAt(pos++);
            const vs = pos;
            while (pos < len && src.charCodeAt(pos) !== q) pos++;
            if (name) attrs.set(name, src.slice(vs, pos));
            if (pos < len) pos++;
        }
    }
    return { attrs, tagEnd: pos };
}

// --- XML header parsing ---

function parseDataArraysInSection(header: string, sectionName: string): VtpDataArrayDescriptor[] {
    const sStart = findTagStart(header, sectionName, 0);
    if (sStart === -1) return [];
    let pos = sStart + 1 + sectionName.length;
    while (pos < header.length && header.charCodeAt(pos) !== 62) pos++;
    pos++; // skip '>'
    const sEnd = header.indexOf('</' + sectionName + '>', pos);
    if (sEnd === -1) return [];

    const results: VtpDataArrayDescriptor[] = [];
    while (pos < sEnd) {
        const daStart = findTagStart(header, 'DataArray', pos);
        if (daStart === -1 || daStart >= sEnd) break;
        const { attrs, tagEnd } = readAttrs(header, daStart + '<DataArray'.length);
        let content = '';
        let nextPos: number;
        if (header.charCodeAt(tagEnd) === 47) {
            nextPos = tagEnd + 2; // '/>'
        } else {
            nextPos = tagEnd + 1;
            const closeDA = header.indexOf('</DataArray>', nextPos);
            if (closeDA !== -1) {
                content = header.slice(nextPos, closeDA);
                nextPos = closeDA + '</DataArray>'.length;
            } else {
                nextPos = sEnd;
            }
        }
        pos = nextPos;

        const fmt = (attrs.get('format') ?? 'binary') as 'ascii' | 'binary' | 'appended';
        const offsetStr = attrs.get('offset');
        const b64 = fmt === 'binary' ? content.replace(/\s/g, '') : '';
        const asciiText = fmt === 'ascii' ? content : undefined;

        if (fmt === 'appended' && offsetStr === undefined) continue;
        if (fmt === 'binary' && b64.length === 0) continue;

        const rangeMinStr = attrs.get('RangeMin');
        const rangeMaxStr = attrs.get('RangeMax');
        results.push({
            name: attrs.get('Name') ?? '',
            type: attrs.get('type') ?? 'Float32',
            numberOfComponents: parseInt(attrs.get('NumberOfComponents') ?? '1', 10),
            format: fmt,
            offset: offsetStr !== undefined ? parseInt(offsetStr, 10) : -1,
            hasRange: rangeMinStr !== undefined && rangeMaxStr !== undefined,
            rangeMin: parseFloat(rangeMinStr ?? '0'),
            rangeMax: parseFloat(rangeMaxStr ?? '1'),
            ...(b64.length > 0 ? { inlineBase64: b64 } : {}),
            ...(asciiText !== undefined ? { asciiText } : {}),
        });
    }
    return results;
}

interface ParsedHeader {
    nPoints: number;
    nCells: number;
    nVerts: number;
    nLines: number;
    nStrips: number;
    byteOrder: string;
    compressor: string;
    headerType: 'UInt32' | 'UInt64';
    pointDataArrays: VtpDataArrayDescriptor[];
    cellDataArrays: VtpDataArrayDescriptor[];
    pointsArrays: VtpDataArrayDescriptor[];
    polysArrays: VtpDataArrayDescriptor[];
}

function parseHeader(header: string): ParsedHeader {
    const vtkStart = findTagStart(header, 'VTKFile', 0);
    if (vtkStart === -1) throw new Error('Cannot find VTKFile element in VTP header');
    const { attrs: vtkAttrs } = readAttrs(header, vtkStart + '<VTKFile'.length);

    const pieceStart = findTagStart(header, 'Piece', 0);
    if (pieceStart === -1) throw new Error('Cannot find Piece element in VTP header');
    const { attrs: pieceAttrs } = readAttrs(header, pieceStart + '<Piece'.length);

    const nPointsStr = pieceAttrs.get('NumberOfPoints');
    const nPolysStr = pieceAttrs.get('NumberOfPolys');
    if (!nPointsStr || !nPolysStr) {
        throw new Error('Cannot parse NumberOfPoints/NumberOfPolys from VTP Piece element');
    }
    return {
        nPoints: parseInt(nPointsStr, 10),
        nCells: parseInt(nPolysStr, 10),
        nVerts: parseInt(pieceAttrs.get('NumberOfVerts') ?? '0', 10),
        nLines: parseInt(pieceAttrs.get('NumberOfLines') ?? '0', 10),
        nStrips: parseInt(pieceAttrs.get('NumberOfStrips') ?? '0', 10),
        byteOrder: vtkAttrs.get('byte_order') ?? 'LittleEndian',
        compressor: vtkAttrs.get('compressor') ?? '',
        headerType: vtkAttrs.get('header_type') === 'UInt64' ? 'UInt64' : 'UInt32',
        pointDataArrays: parseDataArraysInSection(header, 'PointData'),
        cellDataArrays: parseDataArraysInSection(header, 'CellData'),
        pointsArrays: parseDataArraysInSection(header, 'Points'),
        polysArrays: parseDataArraysInSection(header, 'Polys'),
    };
}

/**
 * --- Base64 helpers for VTK appended format ---
 *
 * VTK base64 appended: each DataArray's binary block is split into two base64-encoded chunks:
 *    1. The compressed-block header  (nblocks, blockSize, lastSize, compSize[i])
 *    2. All compressed blocks concatenated
 * The DataArray XML `offset` attribute is a byte offset in the base64 ASCII stream after '_'.
 * Since base64 is ASCII (1 char == 1 byte), the character offset equals the byte offset.
 *
 * --- Zlib block decompression (base64 paths) ---
 * Two variants: bytes (appended — no string) and str (inline — avoids atob).
 */
async function decompressVtkBlockB64Bytes(
    ctx: RuntimeContext,
    rawData: Uint8Array,
    rawByteStart: number,
    byteOffset: number,
    headerType: 'UInt32' | 'UInt64' = 'UInt32',
    hasCompressor = true
): Promise<Uint8Array> {
    const hdrItemBytes = headerType === 'UInt64' ? 8 : 4;
    const base = rawByteStart + byteOffset;

    if (!hasCompressor) {
        const hdr = decodeB64Bytes(rawData, base, hdrItemBytes);
        const nbytes = new DataView(hdr.buffer, hdr.byteOffset, hdrItemBytes).getUint32(0, true);
        return decodeB64Bytes(rawData, base + Math.ceil(hdrItemBytes / 3) * 4, nbytes);
    }

    const minHdr = decodeB64Bytes(rawData, base, 3 * hdrItemBytes);
    const dvMin = new DataView(minHdr.buffer, minHdr.byteOffset, 3 * hdrItemBytes);
    const nblocks = dvMin.getUint32(0, true);
    const blockSize = dvMin.getUint32(hdrItemBytes, true);
    const lastSize = dvMin.getUint32(2 * hdrItemBytes, true);
    if (nblocks === 0) return new Uint8Array(0);

    const totalHdrBytes = (3 + nblocks) * hdrItemBytes;
    const hdrRaw = decodeB64Bytes(rawData, base, totalHdrBytes);
    const dvHdr = new DataView(hdrRaw.buffer, hdrRaw.byteOffset, hdrRaw.byteLength);
    const compSizes = new Uint32Array(nblocks);
    let totalCompressed = 0;
    for (let i = 0; i < nblocks; i++) {
        compSizes[i] = dvHdr.getUint32((3 + i) * hdrItemBytes, true);
        totalCompressed += compSizes[i];
    }

    const compData = decodeB64Bytes(rawData, base + Math.ceil(totalHdrBytes / 3) * 4, totalCompressed);
    const totalUncompressed = lastSize > 0 ? (nblocks - 1) * blockSize + lastSize : nblocks * blockSize;
    const result = new Uint8Array(totalUncompressed);
    let compOff = 0, outOff = 0;
    for (let i = 0; i < nblocks; i++) {
        const csize = compSizes[i];
        const inflated = await inflate(ctx, new Uint8Array(compData.buffer, compData.byteOffset + compOff, csize) as Uint8Array<ArrayBuffer>);
        result.set(inflated, outOff);
        outOff += inflated.length;
        compOff += csize;
    }
    return result;
}

async function decompressVtkBlockB64Str(
    ctx: RuntimeContext,
    src: StringLike,
    charOffset: number,
    headerType: 'UInt32' | 'UInt64' = 'UInt32',
    hasCompressor = true
): Promise<Uint8Array> {
    const hdrItemBytes = headerType === 'UInt64' ? 8 : 4;

    if (!hasCompressor) {
        const hdr = decodeB64Str(src, charOffset, hdrItemBytes);
        const nbytes = new DataView(hdr.buffer, hdr.byteOffset, hdrItemBytes).getUint32(0, true);
        return decodeB64Str(src, charOffset + Math.ceil(hdrItemBytes / 3) * 4, nbytes);
    }

    const minHdr = decodeB64Str(src, charOffset, 3 * hdrItemBytes);
    const dvMin = new DataView(minHdr.buffer, minHdr.byteOffset, 3 * hdrItemBytes);
    const nblocks = dvMin.getUint32(0, true);
    const blockSize = dvMin.getUint32(hdrItemBytes, true);
    const lastSize = dvMin.getUint32(2 * hdrItemBytes, true);
    if (nblocks === 0) return new Uint8Array(0);

    const totalHdrBytes = (3 + nblocks) * hdrItemBytes;
    const hdrRaw = decodeB64Str(src, charOffset, totalHdrBytes);
    const dvHdr = new DataView(hdrRaw.buffer, hdrRaw.byteOffset, hdrRaw.byteLength);
    const compSizes = new Uint32Array(nblocks);
    let totalCompressed = 0;
    for (let i = 0; i < nblocks; i++) {
        compSizes[i] = dvHdr.getUint32((3 + i) * hdrItemBytes, true);
        totalCompressed += compSizes[i];
    }

    const compData = decodeB64Str(src, charOffset + Math.ceil(totalHdrBytes / 3) * 4, totalCompressed);
    const totalUncompressed = lastSize > 0 ? (nblocks - 1) * blockSize + lastSize : nblocks * blockSize;
    const result = new Uint8Array(totalUncompressed);
    let compOff = 0, outOff = 0;
    for (let i = 0; i < nblocks; i++) {
        const csize = compSizes[i];
        const inflated = await inflate(ctx, new Uint8Array(compData.buffer, compData.byteOffset + compOff, csize) as Uint8Array<ArrayBuffer>);
        result.set(inflated, outOff);
        outOff += inflated.length;
        compOff += csize;
    }
    return result;
}

// --- Zlib block decompression (raw path) ---

async function decompressVtkBlockRaw(
    ctx: RuntimeContext,
    rawData: Uint8Array,
    rawByteStart: number,
    byteOffset: number,
    headerType: 'UInt32' | 'UInt64' = 'UInt32',
    hasCompressor = true
): Promise<Uint8Array> {
    const hdrItemBytes = headerType === 'UInt64' ? 8 : 4;
    const absBase = rawByteStart + byteOffset;
    const dv = new DataView(rawData.buffer, rawData.byteOffset, rawData.byteLength);

    if (!hasCompressor) {
        // Uncompressed raw AppendedData: [nbytes (hdrItemBytes)][raw data bytes]
        const nbytes = dv.getUint32(absBase, true); // low 32 bits safe for practical sizes
        return new Uint8Array(rawData.buffer, rawData.byteOffset + absBase + hdrItemBytes, nbytes).slice();
    }

    const nblocks = dv.getUint32(absBase, true);
    const blockSize = dv.getUint32(absBase + hdrItemBytes, true);
    const lastSize = dv.getUint32(absBase + 2 * hdrItemBytes, true);

    if (nblocks === 0) return new Uint8Array(0);

    const hdrBytes = (3 + nblocks) * hdrItemBytes;
    const compSizes = new Uint32Array(nblocks);
    for (let i = 0; i < nblocks; i++) {
        compSizes[i] = dv.getUint32(absBase + (3 + i) * hdrItemBytes, true);
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
    info: AppendedInfo | null,
    desc: VtpDataArrayDescriptor,
    headerType: 'UInt32' | 'UInt64',
    hasCompressor: boolean
): Promise<Uint8Array> {
    if (desc.inlineBase64) {
        return decompressVtkBlockB64Str(ctx, desc.inlineBase64, 0, headerType, hasCompressor);
    }
    if (!info) throw new Error(`No data source for array "${desc.name}"`);
    if (info.encoding === 'base64') {
        return decompressVtkBlockB64Bytes(ctx, info.rawData, info.rawByteStart, desc.offset, headerType, hasCompressor);
    } else {
        return decompressVtkBlockRaw(ctx, info.rawData, info.rawByteStart, desc.offset, headerType, hasCompressor);
    }
}


/**
 * --- ASCII format helpers ---
 * Scan-based tokenizer: avoids the intermediate array that split() would create.
 * Each token is sliced only when parseFloat/parseInt is called; no other copies.
 */
function parseAsciiNumbers(text: StringLike, out: Float32Array | Float64Array | Int32Array, count: number, asInt: boolean): void {
    let pos = 0, n = 0;
    const len = text.length;
    while (pos < len && n < count) {
        while (pos < len && isWS(text.charCodeAt(pos))) pos++;
        if (pos >= len) break;
        const start = pos;
        while (pos < len && !isWS(text.charCodeAt(pos))) pos++;
        if (pos > start) out[n++] = asInt ? fastParseInt(text, start, pos) : fastParseFloat(text, start, pos);
    }
}

function parseAsciiFloat32(text: StringLike, count: number): Float32Array {
    const out = new Float32Array(count);
    parseAsciiNumbers(text, out, count, false);
    return out;
}

function parseAsciiFloat64(text: StringLike, count: number): Float64Array {
    const out = new Float64Array(count);
    parseAsciiNumbers(text, out, count, false);
    return out;
}

function parseAsciiInts(text: StringLike, count: number): Int32Array {
    const out = new Int32Array(count);
    parseAsciiNumbers(text, out, count, true);
    return out;
}

// For ASCII connectivity: count is not known upfront, so two-pass to avoid over-allocation.
function parseAsciiIntsAll(text: StringLike): Int32Array {
    let count = 0, pos = 0;
    const len = text.length;
    while (pos < len) {
        while (pos < len && isWS(text.charCodeAt(pos))) pos++;
        if (pos >= len) break;
        count++;
        while (pos < len && !isWS(text.charCodeAt(pos))) pos++;
    }
    const out = new Int32Array(count);
    parseAsciiNumbers(text, out, count, true);
    return out;
}

/**
 * Fan-triangulate connectivity from VTP Polys offsets.
 * offsets[i] = cumulative vertex count through cell i; cell i uses
 * rawConn[offsets[i-1]..offsets[i]-1] (with offsets[-1] = 0).
 */
function buildTriangles(
    rawConn: Int32Array, offsets: Int32Array, nCells: number
): { connectivity: Int32Array; triangleCellIndex: Int32Array } {
    let nTris = 0;
    for (let i = 0; i < nCells; i++) {
        const start = i > 0 ? offsets[i - 1] : 0;
        nTris += Math.max(0, offsets[i] - start - 2);
    }
    const tris = new Int32Array(3 * nTris);
    const cellIdx = new Int32Array(nTris);
    let ti = 0;
    let triIdx = 0;
    for (let i = 0; i < nCells; i++) {
        const start = i > 0 ? offsets[i - 1] : 0;
        const end = offsets[i];
        const v0 = rawConn[start];
        for (let j = start + 1; j < end - 1; j++) {
            tris[ti++] = v0;
            tris[ti++] = rawConn[j];
            tris[ti++] = rawConn[j + 1];
            cellIdx[triIdx++] = i;
        }
    }
    return { connectivity: tris, triangleCellIndex: cellIdx };
}

// --- Main parser ---

async function parseInternal(data: Uint8Array, ctx: RuntimeContext): Promise<Result<VtpFile>> {
    ctx.update({ message: 'Parsing VTP header...', current: 0, max: 100 });

    // Locate <AppendedData bytewise so header size is unbounded.
    // For inline VTP (no AppendedData) decode the whole file; it is pure ASCII base64-in-XML.
    const appendedTagPos = findBytesPosition(data, '<AppendedData');
    const isInline = appendedTagPos === -1;
    const headerText = new TextDecoder('ascii').decode(
        isInline ? data : data.subarray(0, appendedTagPos)
    );

    const hdr = parseHeader(headerText);

    if (hdr.byteOrder === 'BigEndian') {
        throw new Error('BigEndian VTP files are not supported.');
    }
    if (hdr.compressor && hdr.compressor !== 'vtkZLibDataCompressor') {
        throw new Error(`Unsupported VTP compressor: "${hdr.compressor}". Only vtkZLibDataCompressor is supported.`);
    }
    if (hdr.nVerts + hdr.nLines + hdr.nStrips > 0) {
        console.warn(`VTP: file contains ${hdr.nVerts} verts, ${hdr.nLines} lines, ${hdr.nStrips} strips — only Polys are rendered.`);
    }
    const headerType = hdr.headerType;
    const hasCompressor = !!hdr.compressor;

    ctx.update({ message: 'Locating binary section...', current: 5, max: 100 });
    const info = isInline ? null : extractAppendedInfo(data);

    // Points positions
    const pointsDesc = hdr.pointsArrays.find(d => d.name === 'Points') ?? hdr.pointsArrays[0];
    if (!pointsDesc) throw new Error('VTP file missing Points DataArray');

    ctx.update({ message: 'Decoding positions...', current: 10, max: 100 });
    let positions: Float32Array;
    if (pointsDesc.format === 'ascii') {
        positions = parseAsciiFloat32(pointsDesc.asciiText ?? '', hdr.nPoints * 3);
    } else {
        const posRaw = await decompressBlock(ctx, info, pointsDesc, headerType, hasCompressor);
        positions = decodePositions(posRaw, pointsDesc.type, hdr.nPoints * 3);
    }

    // Polys connectivity + offsets
    const connDesc = hdr.polysArrays.find(d => d.name === 'connectivity');
    const offsetsDesc = hdr.polysArrays.find(d => d.name === 'offsets');
    if (!connDesc || !offsetsDesc) throw new Error('VTP file missing Polys connectivity/offsets DataArrays');

    ctx.update({ message: 'Decoding topology...', current: 20, max: 100 });
    let rawConn: Int32Array;
    let rawOffsets: Int32Array;
    if (connDesc.format === 'ascii') {
        // For ASCII connectivity, count = offsets[last] which we don't know yet — parse all tokens.
        rawConn = parseAsciiIntsAll(connDesc.asciiText ?? '');
        rawOffsets = parseAsciiInts(offsetsDesc.asciiText ?? '', hdr.nCells);
    } else {
        const connRaw = await decompressBlock(ctx, info, connDesc, headerType, hasCompressor);
        const offsetsRaw = await decompressBlock(ctx, info, offsetsDesc, headerType, hasCompressor);
        rawConn = decodeConnectivity(connRaw, connDesc.type);
        rawOffsets = decodeConnectivity(offsetsRaw, offsetsDesc.type);
    }

    const { connectivity, triangleCellIndex } = buildTriangles(rawConn, rawOffsets, hdr.nCells);
    const numberOfTriangles = connectivity.length / 3;

    // PointData arrays (all numberOfComponents supported; multi-component values are interleaved)
    const pointData = new Map<string, VtpScalarArray>();
    for (let i = 0; i < hdr.pointDataArrays.length; i++) {
        const desc = hdr.pointDataArrays[i];
        const totalCount = hdr.nPoints * desc.numberOfComponents;
        ctx.update({ message: `Point data: ${desc.name}...`, current: 30 + Math.floor(i / Math.max(1, hdr.pointDataArrays.length) * 30), max: 100 });
        let values: TypedArray;
        if (desc.format === 'ascii') {
            values = parseAsciiFloat64(desc.asciiText ?? '', totalCount);
        } else {
            const raw = await decompressBlock(ctx, info, desc, headerType, hasCompressor);
            values = decodeTyped(raw, desc.type, totalCount);
        }
        pointData.set(desc.name, { desc, values });
    }

    // CellData arrays (all numberOfComponents supported; multi-component values are interleaved)
    const cellData = new Map<string, VtpScalarArray>();
    for (let i = 0; i < hdr.cellDataArrays.length; i++) {
        const desc = hdr.cellDataArrays[i];
        const totalCount = hdr.nCells * desc.numberOfComponents;
        ctx.update({ message: `Cell data: ${desc.name}...`, current: 60 + Math.floor(i / Math.max(1, hdr.cellDataArrays.length) * 35), max: 100 });
        let values: TypedArray;
        if (desc.format === 'ascii') {
            values = parseAsciiFloat64(desc.asciiText ?? '', totalCount);
        } else {
            const raw = await decompressBlock(ctx, info, desc, headerType, hasCompressor);
            values = decodeTyped(raw, desc.type, totalCount);
        }
        cellData.set(desc.name, { desc, values });
    }

    return Result.success({
        numberOfPoints: hdr.nPoints,
        numberOfCells: hdr.nCells,
        positions,
        connectivity,
        numberOfTriangles,
        triangleCellIndex,
        pointData,
        cellData,
    });
}

export function parseVtp(data: Uint8Array) {
    return Task.create<Result<VtpFile>>('Parse VTP', async ctx => {
        return await parseInternal(data, ctx);
    });
}
