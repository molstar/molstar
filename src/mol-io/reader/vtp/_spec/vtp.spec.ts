/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import * as fs from 'fs';
import * as path from 'path';
import { parseVtp } from '../parser';

const VTP_EXAMPLES_DIR = '/Users/ludovicautin/Dev/surface_morphometrics/surface_morphometrics_example/morphometrics';

function skipIfMissing(dir: string) {
    if (!fs.existsSync(dir)) return true;
    return false;
}

describe('VTP parser', () => {
    it('parses surface.vtp (geometry only, no attributes)', async () => {
        const filePath = path.join(VTP_EXAMPLES_DIR, 'YTC042_2_lam10_ts_004_labels_IMM.surface.vtp');
        if (!fs.existsSync(filePath)) { console.warn('Skipping: file not found'); return; }

        const data = new Uint8Array(fs.readFileSync(filePath));
        const result = await parseVtp(data).run();

        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        expect(vtp.numberOfPoints).toBe(100751);
        expect(vtp.numberOfCells).toBe(194757);
        expect(vtp.positions.length).toBe(vtp.numberOfPoints * 3);
        expect(vtp.connectivity.length).toBe(vtp.numberOfTriangles * 3);

        // Verify positions decode to the correct numerical range
        const posArr = vtp.positions;
        let minVal = Infinity, maxVal = -Infinity;
        for (let i = 0; i < posArr.length; i++) {
            if (posArr[i] < minVal) minVal = posArr[i];
            if (posArr[i] > maxVal) maxVal = posArr[i];
        }
        expect(minVal).toBeCloseTo(-1.49, 1);
        expect(maxVal).toBeCloseTo(383.73, 1);

        expect(vtp.pointData.size).toBe(0);
        expect(vtp.cellData.size).toBe(0);
    });

    it('parses AVV_rh9.vtp (with CellData attributes)', async () => {
        const filePath = path.join(VTP_EXAMPLES_DIR, 'YTC042_2_lam10_ts_004_labels_IMM.AVV_rh9.vtp');
        if (!fs.existsSync(filePath)) { console.warn('Skipping: file not found'); return; }

        const data = new Uint8Array(fs.readFileSync(filePath));
        const result = await parseVtp(data).run();

        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        expect(vtp.numberOfPoints).toBeGreaterThan(0);
        expect(vtp.numberOfCells).toBeGreaterThan(0);

        // Should have scalar cell attributes
        expect(vtp.cellData.size).toBeGreaterThan(0);
        const cellKeys = [...vtp.cellData.keys()];
        console.log('CellData arrays:', cellKeys);
        expect(cellKeys).toContain('mean_curvature');

        const mc = vtp.cellData.get('mean_curvature')!;
        expect(mc.values.length).toBe(vtp.numberOfCells);
        expect(mc.desc.rangeMin).toBeCloseTo(-2.16, 1);
        expect(mc.desc.rangeMax).toBeCloseTo(3.88, 1);
    });

    it('parses scaled_cleaned.vtp (with CellData attributes)', async () => {
        const filePath = path.join(VTP_EXAMPLES_DIR, 'YTC042_2_lam10_ts_004_labels_IMM.scaled_cleaned.vtp');
        if (!fs.existsSync(filePath)) { console.warn('Skipping: file not found'); return; }

        const data = new Uint8Array(fs.readFileSync(filePath));
        const result = await parseVtp(data).run();

        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        console.log('scaled_cleaned numberOfPoints:', vtp.numberOfPoints, 'numberOfCells:', vtp.numberOfCells);
        console.log('CellData:', [...vtp.cellData.keys()]);
        console.log('PointData:', [...vtp.pointData.keys()]);

        expect(vtp.positions.length).toBe(vtp.numberOfPoints * 3);
        expect(vtp.connectivity.length).toBe(vtp.numberOfTriangles * 3);

        // Connectivity indices should be valid
        const maxIdx = Math.max(...Array.from(vtp.connectivity.slice(0, Math.min(300, vtp.connectivity.length))));
        expect(maxIdx).toBeLessThan(vtp.numberOfPoints);
    });
});
