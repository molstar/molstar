/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ReaderResult as Result } from '../result';
import { Task, RuntimeContext } from '../../../mol-task';
import { ChunkedArray } from '../../../mol-data/util';
import { ObjFile } from './schema';
import { StringLike } from '../../common/string-like';
import { Tokenizer } from '../common/text/tokenizer';
import { parseInt, parseFloat } from '../common/text/number-parser';

// OBJ file format specification: http://www.martinreddy.net/gfx/3d/OBJ.spec

interface State {
    tokenizer: Tokenizer
    positions: ChunkedArray<number, 3>
    normals: ChunkedArray<number, 3>
    vertexColors: ChunkedArray<number, 3>
    positionIndices: ChunkedArray<number, 3>
    normalIndices: ChunkedArray<number, 3>
    faceGroups: ChunkedArray<number, 1>
    materialNames: string[]
    materialMap: Map<string, number>
    currentMaterialIdx: number
    hasVertexColors: boolean
    warnings: string[]
}

function State(data: StringLike): State {
    return {
        tokenizer: Tokenizer(data),
        positions: ChunkedArray.create(Float32Array, 3, 1024),
        normals: ChunkedArray.create(Float32Array, 3, 1024),
        vertexColors: ChunkedArray.create(Float32Array, 3, 1024),
        positionIndices: ChunkedArray.create(Int32Array, 3, 1024),
        normalIndices: ChunkedArray.create(Int32Array, 3, 1024),
        faceGroups: ChunkedArray.create(Int32Array, 1, 1024),
        materialNames: [],
        materialMap: new Map(),
        currentMaterialIdx: 0,
        hasVertexColors: false,
        warnings: []
    };
}

// Character codes used for keyword identification without materializing strings
const CC_v = 118; // 'v'
const CC_n = 110; // 'n'
const CC_f = 102; // 'f'
const CC_u = 117; // 'u'
const CC_s = 115; // 's'
const CC_e = 101; // 'e'
const CC_m = 109; // 'm'
const CC_t = 116; // 't'
const CC_l = 108; // 'l'
const CC_HASH = 35; // '#'
const CC_NEWLINE = 10; // '\n'
const CC_CR = 13; // '\r'
const CC_SPACE = 32; // ' '
const CC_TAB = 9; // '\t'
const CC_SLASH = 47; // '/'

/** Skip to the end of the current line without returning it. */
function skipLine(tokenizer: Tokenizer): void {
    const { data } = tokenizer;
    while (tokenizer.position < tokenizer.length) {
        const c = data.charCodeAt(tokenizer.position);
        if (c === CC_NEWLINE) { ++tokenizer.position; ++tokenizer.lineNumber; return; }
        if (c === CC_CR) {
            ++tokenizer.position;
            ++tokenizer.lineNumber;
            if (tokenizer.position < tokenizer.length && data.charCodeAt(tokenizer.position) === CC_NEWLINE) ++tokenizer.position;
            return;
        }
        ++tokenizer.position;
    }
}

/** Skip inline whitespace (space/tab only — does not cross newlines). */
function skipInlineWS(tokenizer: Tokenizer): void {
    const { data } = tokenizer;
    while (tokenizer.position < tokenizer.length) {
        const c = data.charCodeAt(tokenizer.position);
        if (c !== CC_SPACE && c !== CC_TAB) return;
        ++tokenizer.position;
    }
}

/**
 * Read one whitespace-delimited token on the current line.
 * Returns false when end-of-line / end-of-file is reached before any character.
 * Leaves tokenizer.tokenStart/tokenEnd set to the token boundaries.
 */
function readInlineToken(tokenizer: Tokenizer): boolean {
    skipInlineWS(tokenizer);
    const { data } = tokenizer;
    if (tokenizer.position >= tokenizer.length) return false;
    const c = data.charCodeAt(tokenizer.position);
    if (c === CC_NEWLINE || c === CC_CR || c === CC_HASH) return false;
    tokenizer.tokenStart = tokenizer.position;
    while (tokenizer.position < tokenizer.length) {
        const ch = data.charCodeAt(tokenizer.position);
        if (ch === CC_SPACE || ch === CC_TAB || ch === CC_NEWLINE || ch === CC_CR || ch === CC_HASH) break;
        ++tokenizer.position;
    }
    tokenizer.tokenEnd = tokenizer.position;
    return tokenizer.tokenEnd > tokenizer.tokenStart;
}

/**
 * Read up to `maxCount` face-vertex tokens from the current line into `facePos` / `faceNorm`.
 * Returns the number of tokens read.
 * Face vertex format: posIdx[/[texIdx][/normIdx]]  (all 1-based, may be negative).
 */
function readFaceTokens(
    tokenizer: Tokenizer,
    facePos: Int32Array, faceNorm: Int32Array,
    maxCount: number,
    posCount: number, normCount: number
): number {
    const { data } = tokenizer;
    let count = 0;
    while (count < maxCount && readInlineToken(tokenizer)) {
        const start = tokenizer.tokenStart;
        const end = tokenizer.tokenEnd;

        // Find first slash within [start, end)
        let slash1 = -1;
        for (let i = start; i < end; ++i) {
            if (data.charCodeAt(i) === CC_SLASH) { slash1 = i; break; }
        }

        let posIdx: number;
        let normIdx = -1;

        if (slash1 === -1) {
            // "v"
            const p = parseInt(data, start, end);
            posIdx = p < 0 ? posCount + p : p - 1;
        } else {
            const p = parseInt(data, start, slash1);
            posIdx = p < 0 ? posCount + p : p - 1;

            // Find second slash
            let slash2 = -1;
            for (let i = slash1 + 1; i < end; ++i) {
                if (data.charCodeAt(i) === CC_SLASH) { slash2 = i; break; }
            }

            if (slash2 !== -1 && slash2 + 1 < end) {
                // "v/vt/vn" or "v//vn"
                const n = parseInt(data, slash2 + 1, end);
                normIdx = n < 0 ? normCount + n : n - 1;
            }
            // else "v/vt" — no normal
        }

        facePos[count] = posIdx;
        faceNorm[count] = normIdx;
        ++count;
    }
    return count;
}

// Reusable scratch buffers for face vertex data (polygons up to MAX_FACE_VERTICES vertices)
const MAX_FACE_VERTICES = 256;
const _facePos = new Int32Array(MAX_FACE_VERTICES);
const _faceNorm = new Int32Array(MAX_FACE_VERTICES);

function handleUseMtl(state: State): void {
    const { tokenizer } = state;
    if (readInlineToken(tokenizer)) {
        const name = tokenizer.data.substring(tokenizer.tokenStart, tokenizer.tokenEnd);
        if (!state.materialMap.has(name)) {
            const idx = state.materialNames.length;
            state.materialMap.set(name, idx);
            state.materialNames.push(name);
        }
        state.currentMaterialIdx = state.materialMap.get(name)!;
    }
    skipLine(tokenizer);
}

function handleVertex(state: State): void {
    const { tokenizer } = state;
    let x = 0, y = 0, z = 0;
    if (readInlineToken(tokenizer)) x = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
    if (readInlineToken(tokenizer)) y = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
    if (readInlineToken(tokenizer)) z = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
    ChunkedArray.add3(state.positions, x, y, z);

    // Non-standard vertex color extension: `v x y z r g b` with r,g,b in [0, 1].
    // Only mark hasVertexColors when all three components are successfully read.
    let r = 1, g = 1, b = 1;
    if (readInlineToken(tokenizer)) {
        const r0 = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
        if (readInlineToken(tokenizer)) {
            const g0 = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
            if (readInlineToken(tokenizer)) {
                r = r0; g = g0;
                b = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
                state.hasVertexColors = true;
            }
        }
    }
    ChunkedArray.add3(state.vertexColors, r, g, b);
    skipLine(tokenizer);
}

function handleNormal(state: State): void {
    const { tokenizer } = state;
    let x = 0, y = 0, z = 0;
    if (readInlineToken(tokenizer)) x = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
    if (readInlineToken(tokenizer)) y = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
    if (readInlineToken(tokenizer)) z = parseFloat(tokenizer.data, tokenizer.tokenStart, tokenizer.tokenEnd);
    ChunkedArray.add3(state.normals, x, y, z);
    skipLine(tokenizer);
}

function handleFace(state: State): void {
    const { tokenizer } = state;
    const posCount = state.positions.elementCount;
    const normCount = state.normals.elementCount;

    const n = readFaceTokens(tokenizer, _facePos, _faceNorm, MAX_FACE_VERTICES, posCount, normCount);
    if (n < 3) {
        state.warnings.push(`Line ${tokenizer.lineNumber}: degenerate face with ${n} vertices, skipped`);
        skipLine(tokenizer);
        return;
    }
    // Warn if the polygon exceeded the scratch buffer capacity and was truncated.
    if (n === MAX_FACE_VERTICES && readInlineToken(tokenizer)) {
        state.warnings.push(`Line ${tokenizer.lineNumber}: face with more than ${MAX_FACE_VERTICES} vertices truncated`);
    }

    // Fan-triangulate: (0,1,2), (0,2,3), ...
    const p0 = _facePos[0], n0 = _faceNorm[0];
    const group = state.currentMaterialIdx;
    for (let i = 1; i < n - 1; ++i) {
        ChunkedArray.add3(state.positionIndices, p0, _facePos[i], _facePos[i + 1]);
        ChunkedArray.add3(state.normalIndices, n0, _faceNorm[i], _faceNorm[i + 1]);
        ChunkedArray.add(state.faceGroups, group);
    }
    skipLine(tokenizer);
}

async function parseInternal(data: StringLike, ctx: RuntimeContext): Promise<Result<ObjFile>> {
    const state = State(data);
    const { tokenizer } = state;
    const updateChunk = 100000;

    while (tokenizer.position < tokenizer.length) {
        // Skip full-line whitespace and newlines between lines
        const c0 = tokenizer.data.charCodeAt(tokenizer.position);
        if (c0 === CC_NEWLINE) { ++tokenizer.position; ++tokenizer.lineNumber; continue; }
        if (c0 === CC_CR) {
            ++tokenizer.position; ++tokenizer.lineNumber;
            if (tokenizer.position < tokenizer.length && tokenizer.data.charCodeAt(tokenizer.position) === CC_NEWLINE) ++tokenizer.position;
            continue;
        }
        if (c0 === CC_SPACE || c0 === CC_TAB) { skipInlineWS(tokenizer); continue; }
        if (c0 === CC_HASH) { skipLine(tokenizer); continue; }

        // Identify keyword by inspecting character codes — no string allocation
        const c1 = tokenizer.position + 1 < tokenizer.length ? tokenizer.data.charCodeAt(tokenizer.position + 1) : -1;

        if (c0 === CC_f && (c1 === CC_SPACE || c1 === CC_TAB)) {
            // "f " — face
            tokenizer.position += 2;
            handleFace(state);
        } else if (c0 === CC_u) {
            // Check for "usemtl "
            const p = tokenizer.position;
            if (
                p + 6 < tokenizer.length &&
                tokenizer.data.charCodeAt(p + 1) === CC_s &&
                tokenizer.data.charCodeAt(p + 2) === CC_e &&
                tokenizer.data.charCodeAt(p + 3) === CC_m &&
                tokenizer.data.charCodeAt(p + 4) === CC_t &&
                tokenizer.data.charCodeAt(p + 5) === CC_l &&
                (tokenizer.data.charCodeAt(p + 6) === CC_SPACE || tokenizer.data.charCodeAt(p + 6) === CC_TAB)
            ) {
                tokenizer.position += 7;
                handleUseMtl(state);
            } else {
                skipLine(tokenizer);
            }
        } else if (c0 === CC_v) {
            if (c1 === CC_SPACE || c1 === CC_TAB) {
                // "v " — vertex position
                tokenizer.position += 2;
                handleVertex(state);
            } else if (c1 === CC_n) {
                // "vn" — vertex normal
                tokenizer.position += 2;
                handleNormal(state);
            } else {
                // "vt", "vp", etc — skip
                skipLine(tokenizer);
            }
        } else {
            // "g", "o", "s", "mtllib", etc. — skip entire line
            skipLine(tokenizer);
        }

        if (ctx.shouldUpdate && tokenizer.lineNumber % updateChunk === 0) {
            await ctx.update({ message: 'Parsing OBJ', current: tokenizer.position, max: tokenizer.length });
        }
    }

    const posArr = ChunkedArray.compact(state.positions) as Float32Array;
    const normArr = ChunkedArray.compact(state.normals) as Float32Array;
    const posIdxArr = ChunkedArray.compact(state.positionIndices) as Int32Array;
    const normIdxArr = ChunkedArray.compact(state.normalIndices) as Int32Array;
    const faceGroupsArr = ChunkedArray.compact(state.faceGroups) as Int32Array;
    const vertexColorsArr = state.hasVertexColors
        ? ChunkedArray.compact(state.vertexColors) as Float32Array
        : new Float32Array(0);

    const result: ObjFile = {
        positions: posArr,
        normals: normArr,
        positionIndices: posIdxArr,
        normalIndices: normIdxArr,
        positionCount: state.positions.elementCount,
        normalCount: state.normals.elementCount,
        triangleCount: posIdxArr.length / 3,
        vertexColors: vertexColorsArr,
        materialNames: state.materialNames,
        faceGroups: faceGroupsArr
    };

    return Result.success(result, state.warnings);
}

export function parseObj(data: StringLike) {
    return Task.create<Result<ObjFile>>('Parse OBJ', async ctx => {
        return await parseInternal(data, ctx);
    });
}
