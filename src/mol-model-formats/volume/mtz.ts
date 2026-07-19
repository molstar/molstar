/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Volume } from '../../mol-model/volume';
import { Task } from '../../mol-task';
import { SpacegroupCell, Box3D, Spacegroup } from '../../mol-math/geometry';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { ModelFormat } from '../format';
import { CustomProperties } from '../../mol-model/custom-property';
import { degToRad } from '../../mol-math/misc';
import { MtzFile, MtzHeader } from '../../mol-io/reader/mtz/schema';
import { computeElectronDensityFromReflections } from './shared/reflections';

// ---------------------------------------------------------------------------
// Column-pair detection
// ---------------------------------------------------------------------------

/**
 * Describes a detected amplitude + phase column pair that can be used to
 * compute an electron density map.
 */
export interface MtzColumnPair {
    /** Human-readable label for the map type, e.g. '2fo-fc' or 'fo-fc'. */
    label: string;
    /** Column label of the amplitude (type F) column. */
    ampLabel: string;
    /** Column label of the phase (type P) column. */
    phiLabel: string;
}

/**
 * Ordered list of well-known F+P column label pairs used in common refinement
 * programs. The first pair matched in the file is preferred.
 */
const KNOWN_PAIRS: Array<{ label: string, amp: string, phi: string }> = [
    { label: '2fo-fc', amp: 'FWT', phi: 'PHWT' }, // Refmac / CCP4
    { label: 'fo-fc', amp: 'DELFWT', phi: 'DELPHWT' }, // Refmac / CCP4
    { label: 'fo-fc', amp: 'DELFWT', phi: 'PHDELWT' }, // Gemmi?
    { label: '2fo-fc', amp: '2FOFCWT', phi: 'PH2FOFCWT' }, // Phenix
    { label: 'fo-fc', amp: 'FOFCWT', phi: 'PHFOFCWT' }, // Phenix
    { label: '2fo-fc', amp: 'FWT', phi: 'PHIC' }, // some older files
];

/**
 * Scan the MTZ header for recognised amplitude+phase column pairs.
 * Returns an array of detected pairs in priority order (first = most preferred).
 */
export function detectMtzColumnPairs(header: MtzHeader): MtzColumnPair[] {
    const labelSet = new Set(header.columns.map(c => c.label));
    const result: MtzColumnPair[] = [];
    const seen = new Set<string>();

    for (const kp of KNOWN_PAIRS) {
        if (labelSet.has(kp.amp) && labelSet.has(kp.phi)) {
            const key = `${kp.amp}|${kp.phi}`;
            if (!seen.has(key)) {
                seen.add(key);
                result.push({ label: kp.label, ampLabel: kp.amp, phiLabel: kp.phi });
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// Volume computation
// ---------------------------------------------------------------------------

export type MtzMapParams = Partial<{
    label: string;
    entryId: string;
    /** Label of the amplitude column (type F). Defaults to 'FWT'. */
    ampLabel: string;
    /** Label of the phase column (type P). Defaults to 'PHWT'. */
    phiLabel: string;
}>;

/**
 * Compute a real-space electron density volume from an MTZ reflection file.
 *
 * The cell is taken from the first non-base dataset's DCELL record; if absent
 * the global CELL header record is used instead.
 *
 * @param source  Parsed MTZ file.
 * @param params  Optional overrides for labels and identifiers.
 */
export function volumeFromMtz(source: MtzFile, params?: MtzMapParams): Task<Volume> {
    return Task.create<Volume>('Compute MTZ Map', async ctx => {
        const ampLabel = params?.ampLabel ?? 'FWT';
        const phiLabel = params?.phiLabel ?? 'PHWT';

        const { header, data } = source;

        // 1. Determine cell parameters
        // Prefer the first non-base (id > 0) dataset DCELL; fall back to global.
        let rawCell = header.cell;
        for (const ds of header.datasets) {
            if (ds.id > 0 && ds.cell) { rawCell = ds.cell; break; }
        }
        const [a, b, c, alpha, beta, gamma] = rawCell;

        const cellSize = Vec3.create(a, b, c);
        const cellAngles = Vec3.create(degToRad(alpha), degToRad(beta), degToRad(gamma));

        const sgId: number | string = header.spaceGroupNumber || header.spaceGroupName || 'P 1';
        const sgCell = SpacegroupCell.create(sgId, cellSize, cellAngles);
        const sg = Spacegroup.create(sgCell);

        // 2. Locate required columns
        await ctx.update({ message: 'Locating MTZ columns...' });

        const hCol = header.columns.find(c => c.type === 'H' && c.label === 'H');
        const kCol = header.columns.find(c => c.type === 'H' && c.label === 'K');
        const lCol = header.columns.find(c => c.type === 'H' && c.label === 'L');
        const aCol = header.columns.find(c => c.label === ampLabel);
        const pCol = header.columns.find(c => c.label === phiLabel);

        if (!hCol || !kCol || !lCol) {
            throw new Error('MTZ file does not contain H, K, L index columns.');
        }
        if (!aCol) throw new Error(`MTZ amplitude column '${ampLabel}' not found.`);
        if (!pCol) throw new Error(`MTZ phase column '${phiLabel}' not found.`);

        const hIdx = hCol.index, kIdx = kCol.index, lIdx = lCol.index;
        const aIdx = aCol.index, pIdx = pCol.index;
        const nCols = header.nColumns;
        const nRefl = header.nReflections;
        const missing = header.missingNumber;
        const isMissing = (v: number) => isNaN(missing) ? !isFinite(v) : v === missing || !isFinite(v);

        // 3. Extract reflection arrays (filter missing observations)
        await ctx.update({ message: 'Extracting MTZ reflections...' });

        // Pre-allocate at full size; trim after.
        const hArr = new Int16Array(nRefl);
        const kArr = new Int16Array(nRefl);
        const lArr = new Int16Array(nRefl);
        const ampArr = new Float32Array(nRefl);
        const phiArr = new Float32Array(nRefl);

        let validCount = 0;
        for (let r = 0; r < nRefl; r++) {
            const base = r * nCols;
            const aVal = data[base + aIdx];
            const pVal = data[base + pIdx];
            if (isMissing(aVal) || isMissing(pVal)) continue;
            hArr[validCount] = Math.round(data[base + hIdx]);
            kArr[validCount] = Math.round(data[base + kIdx]);
            lArr[validCount] = Math.round(data[base + lIdx]);
            ampArr[validCount] = aVal;
            phiArr[validCount] = pVal;
            validCount++;
        }

        if (validCount === 0) {
            throw new Error(`MTZ file contains no valid reflections for columns '${ampLabel}' and '${phiLabel}'.`);
        }

        // Trim to valid count (subarray avoids copy)
        const hView = hArr.subarray(0, validCount);
        const kView = kArr.subarray(0, validCount);
        const lView = lArr.subarray(0, validCount);
        const ampView = ampArr.subarray(0, validCount);
        const phiView = phiArr.subarray(0, validCount);

        // 4–6. Shared FFT computation
        const { density, N0, N1, N2 } = await computeElectronDensityFromReflections(
            hView, kView, lView, ampView, phiView, sg, ctx
        );

        // 7. Compute statistics
        await ctx.update({ message: 'Finalizing volume...' });
        const totalSize = N0 * N1 * N2;
        let min = Infinity, max = -Infinity, sum = 0, sum2 = 0;
        for (let i = 0; i < totalSize; i++) {
            const v = density[i];
            if (v < min) min = v;
            if (v > max) max = v;
            sum += v;
            sum2 += v * v;
        }
        const mean = sum / totalSize;
        const sigma = Math.sqrt(Math.max(0, sum2 / totalSize - mean * mean));

        // 8. Create Tensor and Volume
        const tensorSpace = Tensor.Space([N0, N1, N2], [0, 1, 2], Float32Array);
        const tensorData = Tensor.create(tensorSpace, Tensor.Data1(density));

        return {
            label: params?.label,
            entryId: params?.entryId,
            grid: {
                transform: {
                    kind: 'spacegroup',
                    cell: sgCell,
                    fractionalBox: Box3D.create(Vec3.zero(), Vec3.create(1, 1, 1)),
                },
                cells: tensorData,
                stats: { min, max, mean, sigma },
                periodicity: 'xyz',
            },
            instances: [{ transform: Mat4.identity() }],
            sourceData: MtzFormat.create(source),
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
            _localPropertyData: Object.create(null),
        };
    });
}

// ---------------------------------------------------------------------------
// MtzFormat (ModelFormat wrapper)
// ---------------------------------------------------------------------------

export { MtzFormat };

type MtzFormat = ModelFormat<MtzFile>

namespace MtzFormat {
    export function is(x?: ModelFormat): x is MtzFormat {
        return x?.kind === 'mtz';
    }

    export function create(mtz: MtzFile): MtzFormat {
        return { kind: 'mtz', name: mtz.name, data: mtz };
    }
}
