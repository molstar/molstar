/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { readFileSync } from 'fs';
import * as path from 'path';
import { parseHffData } from '../../../mol-io/reader/hff/parser';
import { _internals } from '../model';

function loadFixture(): ArrayBuffer {
    const buf = readFileSync(path.resolve(__dirname, '..', '..', '..', 'mol-io', 'reader', 'hff', '_spec', 'fixtures', 'sff_minimal.hff'));
    return buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength);
}

describe('hff mesh assembly', () => {
    it('builds a mesh whose vertex count matches all segments', async () => {
        const data = await parseHffData(loadFixture());
        const { mesh, segmentByGroup } = _internals.buildMesh(data);

        expect(mesh.vertexCount).toBe(4);
        expect(mesh.triangleCount).toBe(4);
        expect(segmentByGroup).toHaveLength(1);
        expect(segmentByGroup[0].id).toBe(1);
    });

    it('assigns each vertex the correct segment group', async () => {
        const data = await parseHffData(loadFixture());
        const { mesh } = _internals.buildMesh(data);

        const groups = mesh.groupBuffer.ref.value;
        // The fixture has one segment, so all vertices belong to group 0.
        for (let i = 0; i < mesh.vertexCount; i++) {
            expect(groups[i]).toBe(0);
        }
    });

    it('preserves vertex coordinates and triangle indices through the build', async () => {
        const data = await parseHffData(loadFixture());
        const { mesh } = _internals.buildMesh(data);

        const v = mesh.vertexBuffer.ref.value;
        expect(Array.from(v.subarray(0, 12))).toEqual([
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
        ]);

        const i = mesh.indexBuffer.ref.value;
        expect(Array.from(i.subarray(0, 12))).toEqual([
            0, 1, 2,
            0, 2, 3,
            0, 3, 1,
            1, 3, 2,
        ]);
    });
});
