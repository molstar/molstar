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
    // Decode the opening tag text to read the encoding attribute.
    const tagText = new TextDecoder('ascii').decode(data.subarray(tagPos, pos + 1));
    const encMatch = tagText.match(/encoding="([^"]+)"/i);
    const encoding = encMatch ? encMatch[1].toLowerCase() : 'raw';
    pos++; // skip '>'
    while (pos < data.length && (data[pos] === 0x20 || data[pos] === 0x09 || data[pos] === 0x0A || data[pos] === 0x0D)) pos++;
    if (pos >= data.length || data[pos] !== 0x5F) {
        throw new Error(`Expected '_' before appended data, got 0x${data[pos]?.toString(16) ?? 'EOF'}`);
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

const _ATTR_RE: Record<string, RegExp> = {
    Name:               /Name=["']([^"']*)["']/,
    type:               /type=["']([^"']*)["']/,
    NumberOfComponents: /NumberOfComponents=["']([^"']*)["']/,
    offset:             /offset=["']([^"']*)["']/,
    RangeMin:           /RangeMin=["']([^"']*)["']/,
    RangeMax:           /RangeMax=["']([^"']*)["']/,
    NumberOfPoints:     /NumberOfPoints=["']([^"']*)["']/,
    NumberOfPolys:      /NumberOfPolys=["']([^"']*)["']/,
    byte_order:         /byte_order=["']([^"']*)["']/,
    compressor:         /compressor=["']([^"']*)["']/,
    header_type:        /header_type=["']([^"']*)["']/,
};
function attrValue(attrStr: string, name: string): string | undefined {
    const re = _ATTR_RE[name] ?? new RegExp(name + '=["\']([^"\']*)["\']');
    const m = attrStr.match(re);
    return m ? m[1] : undefined;
}

function extractSectionArrays(header: string, sectionName: string): VtpDataArrayDescriptor[] {
    const re = new RegExp(`<${sectionName}(?:\\s[^>]*)?>([\\s\\S]*?)</${sectionName}>`);
    const sm = header.match(re);
    if (!sm) return [];
    const results: VtpDataArrayDescriptor[] = [];
    // Match self-closing (<DataArray .../>) or opening tag + text content.
    // For inline binary the text node (base64) always precedes any child elements (<InformationKey>, etc.),
    // so [A-Za-z0-9+/=\s]* captures the base64 and stops at the first '<' without needing </DataArray>.
    const daRe = /<DataArray\s+([^>]*?)(?:\/>|>([A-Za-z0-9+/=\s]*))/g;
    let m;
    while ((m = daRe.exec(sm[1])) !== null) {
        const attrStr = m[1];
        const rawContent = m[2];
        const offsetStr = attrValue(attrStr, 'offset');
        const b64 = rawContent ? rawContent.replace(/\s/g, '') : '';
        // Skip if no data source at all
        if (offsetStr === undefined && b64.length === 0) continue;
        const desc: VtpDataArrayDescriptor = {
            name: attrValue(attrStr, 'Name') ?? '',
            type: attrValue(attrStr, 'type') ?? 'Float32',
            numberOfComponents: parseInt(attrValue(attrStr, 'NumberOfComponents') ?? '1', 10),
            offset: offsetStr !== undefined ? parseInt(offsetStr, 10) : -1,
            rangeMin: parseFloat(attrValue(attrStr, 'RangeMin') ?? '0'),
            rangeMax: parseFloat(attrValue(attrStr, 'RangeMax') ?? '1'),
            ...(b64.length > 0 ? { inlineBase64: b64 } : {}),
        };
        results.push(desc);
    }
    return results;
}

interface ParsedHeader {
    nPoints: number;
    nCells: number;
    byteOrder: string;
    compressor: string;
    headerType: string;
    pointDataArrays: VtpDataArrayDescriptor[];
    cellDataArrays: VtpDataArrayDescriptor[];
    pointsArrays: VtpDataArrayDescriptor[];
    polysArrays: VtpDataArrayDescriptor[];
}

function parseHeader(header: string): ParsedHeader {
    const pieceMatch = header.match(/<Piece\s([^>]*)>/);
    if (!pieceMatch) throw new Error('Cannot find Piece element in VTP header');
    const pieceAttrs = pieceMatch[1];
    const nPointsStr = attrValue(pieceAttrs, 'NumberOfPoints');
    const nPolysStr = attrValue(pieceAttrs, 'NumberOfPolys');
    if (nPointsStr === undefined || nPolysStr === undefined) {
        throw new Error('Cannot parse NumberOfPoints/NumberOfPolys from VTP Piece element');
    }
    const vtkFileMatch = header.match(/<VTKFile\s+([^>]*)>/);
    const vtkAttrs = vtkFileMatch ? vtkFileMatch[1] : '';
    return {
        nPoints: parseInt(nPointsStr, 10),
        nCells: parseInt(nPolysStr, 10),
        byteOrder: attrValue(vtkAttrs, 'byte_order') ?? 'LittleEndian',
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
    charOffset: number,
    headerType: 'UInt32' | 'UInt64' = 'UInt32'
): Promise<Uint8Array> {
    const hdrItemBytes = headerType === 'UInt64' ? 8 : 4;

    // 1. Read first 3 header items (nblocks, blockSize, lastSize)
    const minHdrBytes = 3 * hdrItemBytes;
    const minHdr = decodeBase64Slice(rawB64, charOffset, minHdrBytes);
    const dvMin = new DataView(minHdr.buffer, minHdr.byteOffset, minHdrBytes);
    // For UInt64 we read the low 32 bits (little-endian; high bits are 0 for practical sizes)
    const nblocks = dvMin.getUint32(0, true);
    const blockSize = dvMin.getUint32(hdrItemBytes, true);
    const lastSize = dvMin.getUint32(2 * hdrItemBytes, true);

    // 2. Read full header (adds nblocks * hdrItemBytes bytes for compressed sizes)
    const totalHdrBytes = (3 + nblocks) * hdrItemBytes;
    const hdrRaw = decodeBase64Slice(rawB64, charOffset, totalHdrBytes);
    const hdrChars = Math.ceil(totalHdrBytes / 3) * 4;

    if (nblocks === 0) return new Uint8Array(0);

    const dvHdr = new DataView(hdrRaw.buffer, hdrRaw.byteOffset, hdrRaw.byteLength);
    const compSizes = new Uint32Array(nblocks);
    let totalCompressed = 0;
    for (let i = 0; i < nblocks; i++) {
        compSizes[i] = dvHdr.getUint32((3 + i) * hdrItemBytes, true);
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
        if (!hasCompressor) {
            // Uncompressed inline: [nbytes (hdrItemBytes)][raw data]
            const hdrItemBytes = headerType === 'UInt64' ? 8 : 4;
            const binary = atob(desc.inlineBase64);
            const bytes = new Uint8Array(binary.length);
            for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
            return bytes.slice(hdrItemBytes);
        }
        // Compressed inline: same block structure as AppendedData base64, starting at char 0
        return decompressVtkBlockB64(ctx, desc.inlineBase64, 0, headerType);
    }
    if (!info) throw new Error(`No data source for array "${desc.name}"`);
    if (info.encoding === 'base64') {
        return decompressVtkBlockB64(ctx, info.rawB64, desc.offset, headerType);
    } else {
        return decompressVtkBlockRaw(ctx, info.rawData, info.rawByteStart, desc.offset, headerType, hasCompressor);
    }
}

// --- Typed array decoding ---

function decodeFloat32(raw: Uint8Array, count: number): Float32Array {
    if (raw.byteOffset % 4 === 0) return new Float32Array(raw.buffer, raw.byteOffset, count);
    const aligned = raw.slice(0, count * 4);
    return new Float32Array(aligned.buffer, 0, count);
}

function decodePositions(raw: Uint8Array, type: string, count: number): Float32Array {
    if (type === 'Float32') return decodeFloat32(raw, count);
    // Decode any other type through Float64 then downcast — precision loss is acceptable for visualization
    return new Float32Array(decodeToFloat64(raw, type, count));
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
    throw new Error(`Unsupported VTP connectivity type: "${type}". Expected Int32, UInt32, Int64, or UInt64.`);
}

function decodeInt64AsInt32(raw: Uint8Array, count: number): Int32Array {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const out = new Int32Array(count);
    for (let i = 0; i < count; i++) out[i] = dv.getUint32(i * 8, true);
    return out;
}

function decodeToFloat64(raw: Uint8Array, type: string, count: number): Float64Array {
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const out = new Float64Array(count);
    switch (type) {
        case 'Float32': { const f32 = decodeFloat32(raw, count); for (let i = 0; i < count; i++) out[i] = f32[i]; break; }
        case 'Float64': for (let i = 0; i < count; i++) out[i] = dv.getFloat64(i * 8, true); break;
        case 'Int8':   for (let i = 0; i < count; i++) out[i] = dv.getInt8(i); break;
        case 'UInt8':  for (let i = 0; i < count; i++) out[i] = dv.getUint8(i); break;
        case 'Int16':  for (let i = 0; i < count; i++) out[i] = dv.getInt16(i * 2, true); break;
        case 'UInt16': for (let i = 0; i < count; i++) out[i] = dv.getUint16(i * 2, true); break;
        case 'Int32':  for (let i = 0; i < count; i++) out[i] = dv.getInt32(i * 4, true); break;
        case 'UInt32': for (let i = 0; i < count; i++) out[i] = dv.getUint32(i * 4, true); break;
        case 'Int64': case 'UInt64':
            for (let i = 0; i < count; i++) {
                const lo = dv.getUint32(i * 8, true);
                const hi = dv.getInt32(i * 8 + 4, true);
                out[i] = hi * 0x100000000 + lo;
            }
            break;
        default: throw new Error(`Unsupported VTP scalar type: "${type}".`);
    }
    return out;
}

// --- Fan-triangulate connectivity from VTP Polys offsets ---
// offsets[i] = cumulative vertex count through cell i; cell i uses
// rawConn[offsets[i-1]..offsets[i]-1] (with offsets[-1] = 0).

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
    if (hdr.headerType !== 'UInt32' && hdr.headerType !== 'UInt64') {
        throw new Error(`Unsupported VTP header_type: "${hdr.headerType}". Only UInt32 and UInt64 are supported.`);
    }
    const headerType = hdr.headerType as 'UInt32' | 'UInt64';
    const hasCompressor = !!hdr.compressor;

    ctx.update({ message: 'Locating binary section...', current: 5, max: 100 });
    const info = isInline ? null : extractAppendedInfo(data);

    // Points positions
    const pointsDesc = hdr.pointsArrays.find(d => d.name === 'Points') ?? hdr.pointsArrays[0];
    if (!pointsDesc) throw new Error('VTP file missing Points DataArray');

    ctx.update({ message: 'Decompressing positions...', current: 10, max: 100 });
    const posRaw = await decompressBlock(ctx, info, pointsDesc, headerType, hasCompressor);
    const positions = decodePositions(posRaw, pointsDesc.type, hdr.nPoints * 3);

    // Polys connectivity + offsets
    const connDesc = hdr.polysArrays.find(d => d.name === 'connectivity');
    const offsetsDesc = hdr.polysArrays.find(d => d.name === 'offsets');
    if (!connDesc || !offsetsDesc) throw new Error('VTP file missing Polys connectivity/offsets DataArrays');

    ctx.update({ message: 'Decompressing topology...', current: 20, max: 100 });
    const connRaw = await decompressBlock(ctx, info, connDesc, headerType, hasCompressor);
    const offsetsRaw = await decompressBlock(ctx, info, offsetsDesc, headerType, hasCompressor);

    const rawConn = decodeConnectivity(connRaw, connDesc.type);
    const rawOffsets = decodeConnectivity(offsetsRaw, offsetsDesc.type);

    const { connectivity, triangleCellIndex } = buildTriangles(rawConn, rawOffsets, hdr.nCells);
    const numberOfTriangles = connectivity.length / 3;

    // PointData scalar arrays
    const pointData = new Map<string, VtpScalarArray>();
    const skippedPD = hdr.pointDataArrays.filter(d => d.numberOfComponents !== 1);
    if (skippedPD.length > 0) {
        console.warn(`VTP: skipping ${skippedPD.length} PointData array(s) with numberOfComponents > 1: ${skippedPD.map(d => d.name).join(', ')}`);
    }
    const scalarPD = hdr.pointDataArrays.filter(d => d.numberOfComponents === 1);
    for (let i = 0; i < scalarPD.length; i++) {
        const desc = scalarPD[i];
        ctx.update({ message: `Point data: ${desc.name}...`, current: 30 + Math.floor(i / Math.max(1, scalarPD.length) * 30), max: 100 });
        const raw = await decompressBlock(ctx, info, desc, headerType, hasCompressor);
        const values = decodeToFloat64(raw, desc.type, hdr.nPoints);
        pointData.set(desc.name, { desc, values });
    }

    // CellData scalar arrays
    const cellData = new Map<string, VtpScalarArray>();
    const skippedCD = hdr.cellDataArrays.filter(d => d.numberOfComponents !== 1);
    if (skippedCD.length > 0) {
        console.warn(`VTP: skipping ${skippedCD.length} CellData array(s) with numberOfComponents > 1: ${skippedCD.map(d => d.name).join(', ')}`);
    }
    const scalarCD = hdr.cellDataArrays.filter(d => d.numberOfComponents === 1);
    for (let i = 0; i < scalarCD.length; i++) {
        const desc = scalarCD[i];
        ctx.update({ message: `Cell data: ${desc.name}...`, current: 60 + Math.floor(i / Math.max(1, scalarCD.length) * 35), max: 100 });
        const raw = await decompressBlock(ctx, info, desc, headerType, hasCompressor);
        const values = decodeToFloat64(raw, desc.type, hdr.nCells);
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
