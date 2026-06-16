/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { readFileSync } from 'fs';
import * as path from 'path';
import { parseHffData } from '../parser';

function loadFixture(name: string): ArrayBuffer {
    const buf = readFileSync(path.resolve(__dirname, 'fixtures', name));
    return buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength);
}

describe('hff parser', () => {
    it('reads top-level metadata, transforms, segments, and a mesh', async () => {
        const data = await parseHffData(loadFixture('sff_minimal.hff'));

        expect(data.name).toBe('minimal');
        expect(data.details).toBe('a tetrahedron');
        expect(data.primaryDescriptor).toBe('mesh_list');
        expect(data.version).toBe('0.8.0.dev1');

        expect(data.transforms).toHaveLength(1);
        expect(data.transforms[0].rows).toBe(3);
        expect(data.transforms[0].cols).toBe(4);
        // identity 3x4 → 1 on the diagonal
        expect(Array.from(data.transforms[0].data)).toEqual([
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
        ]);

        expect(data.segments).toHaveLength(1);
        const seg = data.segments[0];
        expect(seg.id).toBe(1);
        expect(seg.parentId).toBe(0);
        expect(seg.colour).toEqual([1, 0, 0, 1]);
        expect(seg.biologicalAnnotation?.name).toBe('tetra');
        expect(seg.biologicalAnnotation?.description).toBe('test segment');
        expect(seg.biologicalAnnotation?.numberOfInstances).toBe(1);

        expect(seg.meshes).toHaveLength(1);
        const m = seg.meshes[0];

        expect(m.vertices.count).toBe(4);
        expect(m.vertices.mode).toBe('float32');
        expect(m.vertices.data).toBeInstanceOf(Float32Array);
        expect(Array.from(m.vertices.data as Float32Array)).toEqual([
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
        ]);

        expect(m.triangles.count).toBe(4);
        expect(m.triangles.mode).toBe('uint32');
        expect(m.triangles.data).toBeInstanceOf(Uint32Array);
        expect(Array.from(m.triangles.data as Uint32Array)).toEqual([
            0, 1, 2,
            0, 2, 3,
            0, 3, 1,
            1, 3, 2,
        ]);

        expect(m.normals).toBeUndefined();
    });
});
