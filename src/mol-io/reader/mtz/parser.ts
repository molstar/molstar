/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '../../../mol-task';
import { ReaderResult as Result } from '../result';
import { SimpleBuffer } from '../../common/simple-buffer';
import { MtzFile, MtzHeader, MtzColumn, MtzDataset } from './schema';
import { asciiSlice } from '../../common/ascii';

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/**
 * Trim trailing spaces/nulls from an ASCII string extracted from an MTZ record.
 * MTZ header records are fixed at 80 bytes, padded with spaces.
 */
function trimAscii(s: string): string {
    return s.replace(/\s+$/, '').replace(/\0/g, '').trim();
}

/**
 * Detect whether an MTZ buffer was written in little-endian byte order.
 * Strategy:
 *  1. Read the 4-byte machine stamp at bytes 8–11 (no endian conversion needed
 *     since we check individual bytes).
 *  2. If stamp is unambiguous (0x44/0x41 = LE or 0x11/0x11 = BE) use it.
 *  3. Otherwise heuristically check if the header-pointer word at bytes 4–7
 *     produces a plausible value under each interpretation.
 */
function detectLittleEndian(buf: Uint8Array, fileSize: number): boolean {
    const s0 = buf[8], s1 = buf[9];
    // CCP4 machine stamp: 0x44 (68) and 0x41 (65) → IEEE little-endian
    if (s0 === 68 && s1 === 65) return true;
    // CCP4 machine stamp: 0x11 (17) and 0x11 (17) → big-endian
    if (s0 === 17 && s1 === 17) return false;

    // Heuristic: try both interpretations of the header-pointer word and pick
    // the one that gives a value inside [21, fileSize/4].
    const view = new DataView(buf.buffer, buf.byteOffset);
    const headerWordLE = view.getInt32(4, true);
    const headerWordBE = view.getInt32(4, false);
    const maxWord = (fileSize / 4) | 0;
    const leOk = headerWordLE >= 21 && headerWordLE <= maxWord;
    const beOk = headerWordBE >= 21 && headerWordBE <= maxWord;
    if (leOk && !beOk) return true;
    if (beOk && !leOk) return false;
    // Default to little-endian (most modern systems)
    return true;
}

/**
 * Parse a single SYMINF record payload.
 * Format: `nsym nprim ltype spgno 'spgname' pgname`
 * Example: `4  2 P 4 'P 21 21 21' PG222`
 */
function parseSYMINF(payload: string): { nsym: number, spgno: number, spgname: string } {
    const parts = payload.trim().split(/\s+/);
    const nsym = parseInt(parts[0], 10) || 1;
    const spgno = parseInt(parts[3], 10) || 1;
    // Space-group name may be quoted
    const qStart = payload.indexOf("'");
    const qEnd = payload.lastIndexOf("'");
    const spgname = qStart >= 0 && qEnd > qStart
        ? payload.slice(qStart + 1, qEnd)
        : (parts[4] ?? 'P 1');
    return { nsym, spgno, spgname };
}

/**
 * Parse a COL header record payload.
 * Format: `label type min max datasetId`
 * Column labels may contain spaces so we parse from the right.
 */
function parseCOL(payload: string, columnIndex: number): MtzColumn {
    // The last tokens are: type(1char) min max datasetId
    const parts = payload.trim().split(/\s+/);
    const n = parts.length;
    // datasetId, max, min, type are the last 4 tokens; label is everything before
    const datasetId = parseInt(parts[n - 1], 10) || 0;
    const max = parseFloat(parts[n - 2]);
    const min = parseFloat(parts[n - 3]);
    const type = parts[n - 4] ?? 'R';
    const label = parts.slice(0, n - 4).join(' ');
    return { label: trimAscii(label), type, min, max, datasetId, index: columnIndex };
}

/**
 * Parse the 6 float cell parameters from a space-separated payload.
 * Returns [a, b, c, alpha, beta, gamma].
 */
function parseCellPayload(payload: string): [number, number, number, number, number, number] {
    const p = payload.trim().split(/\s+/).map(Number);
    return [p[0] ?? 0, p[1] ?? 0, p[2] ?? 0, p[3] ?? 90, p[4] ?? 90, p[5] ?? 90];
}

// ---------------------------------------------------------------------------
// Main parser
// ---------------------------------------------------------------------------

/**
 * Parse a binary MTZ reflection file.
 *
 * MTZ binary layout (all values IEEE float32 / int32, endianness determined
 * from the machine stamp):
 *
 *   Bytes 0–3   : magic "MTZ "
 *   Bytes 4–7   : int32  = 1-indexed word position of the main header
 *   Bytes 8–11  : machine stamp (endianness indicator)
 *   Bytes 80+   : float32 reflection data, row-major (nRefl × nCols)
 *   Bytes hdr*4 : 80-byte ASCII header records until "END "
 *
 * @param buffer  Raw file bytes (Uint8Array or SimpleBuffer).
 * @param name    File name / identifier stored on the result.
 */
export function parseMtz(buffer: Uint8Array, name: string): Task<Result<MtzFile>> {
    return Task.create<Result<MtzFile>>('Parse MTZ', async ctx => {
        try {
            if (buffer.length < 80) {
                return Result.error<MtzFile>('MTZ file too small to be valid.');
            }

            // Validate magic bytes "MTZ "
            const magic = asciiSlice(buffer, 0, 4);
            if (magic.slice(0, 3) !== 'MTZ') {
                return Result.error<MtzFile>('Not an MTZ file: magic bytes mismatch.');
            }

            const littleEndian = detectLittleEndian(buffer, buffer.length);
            const view = new DataView(buffer.buffer, buffer.byteOffset);

            const readInt32 = (byteOffset: number) => view.getInt32(byteOffset, littleEndian);
            const readFloat32 = (byteOffset: number) => view.getFloat32(byteOffset, littleEndian);

            // Word position of main header (1-indexed; multiply by 4 for byte offset,
            // then subtract 4 because header is 1-indexed from word 1 = byte 0).
            const headerWordPos = readInt32(4); // 1-indexed word number
            const headerBytePos = (headerWordPos - 1) * 4;

            if (headerBytePos <= 0 || headerBytePos >= buffer.length) {
                return Result.error<MtzFile>(`MTZ header position out of range: word ${headerWordPos}`);
            }

            // ----------------------------------------------------------------
            // Parse header records (each record is exactly 80 bytes of ASCII)
            // ----------------------------------------------------------------
            await ctx.update({ message: 'Reading MTZ header...' });

            let version = '';
            let nColumns = 0;
            let nReflections = 0;
            let nBatches = 0;
            let title = '';
            let globalCell: [number, number, number, number, number, number] = [1, 1, 1, 90, 90, 90];
            let spaceGroupNumber = 1;
            let spaceGroupName = 'P 1';
            let nSymOps = 1;
            let missingNumber = NaN;
            const symops: string[] = [];
            const columns: MtzColumn[] = [];
            const datasets: MtzDataset[] = [];
            // Accumulate datasets by ID for incremental fill from PROJECT/CRYSTAL/DATASET
            const datasetById = new Map<number, MtzDataset>();

            let pos = headerBytePos;
            while (pos + 80 <= buffer.length) {
                const record = asciiSlice(buffer, pos, pos + 80);
                pos += 80;

                // Keyword is the first whitespace-delimited token; payload is everything after.
                // This supports both the full keywords written by gemmi/newer CCP4 library
                // (COLUMN, SYMINF, PROJECT, CRYSTAL, DATASET, DCELL, DWAVEL, TITLE) and the
                // 4-char abbreviations used by some older writers (COL, SYMI, PROJ, …).
                const spaceIdx = record.indexOf(' ');
                const keyword = (spaceIdx >= 0 ? record.slice(0, spaceIdx) : record).trim();
                const payload = spaceIdx >= 0 ? record.slice(spaceIdx + 1) : '';

                if (keyword === 'END') break;

                switch (keyword) {
                    case 'VERS':
                        version = trimAscii(payload);
                        break;
                    case 'TITLE': // gemmi / CCP4 library
                    case 'TITL': // legacy 4-char form
                        title = trimAscii(payload);
                        break;
                    case 'NCOL': {
                        const p = payload.trim().split(/\s+/);
                        nColumns = parseInt(p[0], 10) || 0;
                        nReflections = parseInt(p[1], 10) || 0;
                        nBatches = parseInt(p[2], 10) || 0;
                        break;
                    }
                    case 'CELL':
                        globalCell = parseCellPayload(payload);
                        break;
                    case 'SYMINF': // gemmi / CCP4 library
                    case 'SYMI': // legacy 4-char form
                    {
                        const info = parseSYMINF(payload);
                        nSymOps = info.nsym;
                        spaceGroupNumber = info.spgno;
                        spaceGroupName = info.spgname;
                        break;
                    }
                    case 'SYMM':
                        symops.push(trimAscii(payload));
                        break;
                    case 'COLUMN': // gemmi / CCP4 library
                    case 'COL': // legacy 4-char form
                    {
                        const col = parseCOL(payload, columns.length);
                        columns.push(col);
                        break;
                    }
                    case 'VALM': {
                        const v = parseFloat(payload.trim());
                        missingNumber = isNaN(v) ? NaN : v;
                        break;
                    }
                    case 'PROJECT': // gemmi / CCP4 library
                    case 'PROJ': // legacy 4-char form
                    {
                        // payload: "datasetId projectName"
                        const p = payload.trim().split(/\s+/);
                        const id = parseInt(p[0], 10) || 0;
                        if (!datasetById.has(id)) {
                            datasetById.set(id, { id, project: '', crystal: '', name: '' });
                        }
                        datasetById.get(id)!.project = p.slice(1).join(' ');
                        break;
                    }
                    case 'CRYSTAL': // gemmi / CCP4 library
                    case 'CRYS': // legacy 4-char form
                    {
                        const p = payload.trim().split(/\s+/);
                        const id = parseInt(p[0], 10) || 0;
                        if (!datasetById.has(id)) {
                            datasetById.set(id, { id, project: '', crystal: '', name: '' });
                        }
                        datasetById.get(id)!.crystal = p.slice(1).join(' ');
                        break;
                    }
                    case 'DATASET': // gemmi / CCP4 library
                    case 'DATA': // legacy 4-char form
                    {
                        const p = payload.trim().split(/\s+/);
                        const id = parseInt(p[0], 10) || 0;
                        if (!datasetById.has(id)) {
                            datasetById.set(id, { id, project: '', crystal: '', name: '' });
                        }
                        datasetById.get(id)!.name = p.slice(1).join(' ');
                        break;
                    }
                    case 'DCELL': // gemmi / CCP4 library
                    case 'DCEL': // legacy 4-char form
                    {
                        const p = payload.trim().split(/\s+/);
                        const id = parseInt(p[0], 10) || 0;
                        if (!datasetById.has(id)) {
                            datasetById.set(id, { id, project: '', crystal: '', name: '' });
                        }
                        datasetById.get(id)!.cell = parseCellPayload(p.slice(1).join(' '));
                        break;
                    }
                    case 'DWAVEL': // gemmi / CCP4 library
                    case 'DWAV': // legacy 4-char form
                    {
                        const p = payload.trim().split(/\s+/);
                        const id = parseInt(p[0], 10) || 0;
                        if (!datasetById.has(id)) {
                            datasetById.set(id, { id, project: '', crystal: '', name: '' });
                        }
                        datasetById.get(id)!.wavelength = parseFloat(p[1]) || undefined;
                        break;
                    }
                    // SORT, RESO, NDIF, BATCH, COLSRC, COLGRP, history lines — intentionally ignored
                }
            }

            // Sort datasets by ID so datasets[0] is always HKL_base (id=0)
            for (const ds of datasetById.values()) datasets.push(ds);
            datasets.sort((a, b) => a.id - b.id);

            if (nColumns === 0 || nReflections === 0) {
                return Result.error<MtzFile>(`MTZ header did not contain valid NCOL record (nColumns=${nColumns}, nReflections=${nReflections}).`);
            }

            // ----------------------------------------------------------------
            // Read reflection data
            // Starts at byte 80 (word 21, 1-indexed), laid out as nRefl × nCols
            // float32 values.
            // ----------------------------------------------------------------
            await ctx.update({ message: 'Reading MTZ reflection data...' });

            const dataByteStart = 80; // (21 - 1) * 4
            const nValues = nReflections * nColumns;
            const dataByteEnd = dataByteStart + nValues * 4;

            if (dataByteEnd > buffer.length) {
                return Result.error<MtzFile>(
                    `MTZ file truncated: expected ${dataByteEnd} bytes for ${nReflections} × ${nColumns} reflections, got ${buffer.length}.`
                );
            }

            const data = new Float32Array(nValues);

            if (littleEndian === SimpleBuffer.IsNativeEndianLittle) {
                // Fast path: copy directly
                const src = new Float32Array(buffer.buffer, buffer.byteOffset + dataByteStart, nValues);
                data.set(src);
            } else {
                // Byte-swap each 4-byte float
                for (let i = 0; i < nValues; i++) {
                    data[i] = readFloat32(dataByteStart + i * 4);
                }
            }

            const header: MtzHeader = {
                version,
                title,
                nColumns,
                nReflections,
                nBatches,
                cell: globalCell,
                spaceGroupNumber,
                spaceGroupName,
                nSymOps,
                symops,
                columns,
                datasets,
                missingNumber,
            };

            return Result.success({ name, header, data });
        } catch (e: any) {
            return Result.error<MtzFile>(e?.message ?? String(e));
        }
    });
}
