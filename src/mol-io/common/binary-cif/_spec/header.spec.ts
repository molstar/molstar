/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { binaryCifHasCategory, binaryCifHasColumn, decodeBinaryCifHeader, getBinaryCifHeader } from '../../binary-cif';
import { decodeMsgPack } from '../../msgpack/decode';
import type { EncodedFile } from '../../binary-cif';
import { CifWriter } from '../../../writer/cif';

function createTestBcif() {
    const rows = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    const encoder = CifWriter.createEncoder({ binary: true, encoderName: 'test-encoder' });
    encoder.startDataBlock('TEST_BLOCK');
    encoder.writeCategory({
        name: 'atom_site',
        instance: () => CifWriter.categoryInstance([
            CifWriter.Field.int<number>('id', i => i),
            CifWriter.Field.float<number>('Cartn_x', i => Math.sqrt(i + 1)),
            CifWriter.Field.str<number>('type_symbol', i => i % 2 === 0 ? 'C' : 'N'),
        ], { data: rows, rowCount: rows.length })
    });
    encoder.writeCategory({
        name: 'refln',
        instance: () => CifWriter.categoryInstance([
            CifWriter.Field.float<number>('pdbx_FWT', i => i),
        ], { data: rows, rowCount: rows.length })
    });
    return encoder.getData() as Uint8Array;
}

describe('binary-cif header', () => {
    const data = createTestBcif();

    it('reads file, block, category and column metadata', () => {
        const header = decodeBinaryCifHeader(data);
        expect(header.encoder).toBe('test-encoder');
        expect(header.version.length).toBeGreaterThan(0);
        expect(header.dataBlocks.length).toBe(1);
        expect(header.dataBlocks[0].header).toBe('TEST_BLOCK');
        expect(header.dataBlocks[0].categories.map(c => c.name)).toEqual(['_atom_site', '_refln']);
        expect(header.dataBlocks[0].categories[0].rowCount).toBe(10);
        expect(header.dataBlocks[0].categories[0].columnNames).toEqual(['id', 'Cartn_x', 'type_symbol']);
    });

    it('matches a full msgpack decode', () => {
        const header = decodeBinaryCifHeader(data);
        const file = decodeMsgPack(data) as EncodedFile;
        expect(header.version).toBe(file.version);
        expect(header.encoder).toBe(file.encoder);
        expect(header.dataBlocks.map(b => b.header)).toEqual(file.dataBlocks.map(b => b.header));
        for (let i = 0; i < file.dataBlocks.length; i++) {
            expect(header.dataBlocks[i].categories.map(c => c.name)).toEqual(file.dataBlocks[i].categories.map(c => c.name));
            expect(header.dataBlocks[i].categories.map(c => c.columnNames)).toEqual(file.dataBlocks[i].categories.map(c => c.columns.map(f => f.name)));
        }
    });

    it('reads a Uint8Array view with a nonzero byte offset', () => {
        const padded = new Uint8Array(data.length + 17);
        padded.set(data, 17);
        expect(decodeBinaryCifHeader(padded.subarray(17))).toEqual(decodeBinaryCifHeader(data));
    });

    it('memoizes per buffer instance', () => {
        expect(getBinaryCifHeader(data)).toBe(getBinaryCifHeader(data));
        expect(getBinaryCifHeader(data)).not.toBe(getBinaryCifHeader(data.slice()));
    });

    it('supports category and column lookup', () => {
        const header = decodeBinaryCifHeader(data);
        expect(binaryCifHasCategory(header, '_atom_site')).toBe(true);
        expect(binaryCifHasCategory(header, '_volume_data_3d')).toBe(false);
        expect(binaryCifHasColumn(header, '_refln', 'pdbx_FWT')).toBe(true);
        expect(binaryCifHasColumn(header, '_refln', 'pdbx_DELFWT')).toBe(false);
        expect(binaryCifHasColumn(header, '_atom_site', 'pdbx_FWT')).toBe(false);
    });

    it('throws on malformed data', () => {
        expect(() => decodeBinaryCifHeader(new Uint8Array([0xc1]))).toThrow();
    });
});
