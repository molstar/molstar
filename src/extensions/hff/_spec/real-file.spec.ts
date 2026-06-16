/**
 * Smoke test against a real EMDB-SFF .hff file from the user's local disk.
 * Skipped when the file is not present.
 *
 * Set HFF_FIXTURE to override the default path.
 */
import { existsSync, readFileSync } from 'fs';
import { parseHffData } from '../../../mol-io/reader/hff/parser';
import { _internals } from '../model';

// Point HFF_FIXTURE at a real .hff file to enable this spec, e.g.
//   HFF_FIXTURE=/path/to/sample.hff npx jest src/extensions/hff/_spec/real-file.spec.ts
const PATH = process.env.HFF_FIXTURE;
const have = !!PATH && existsSync(PATH);
const maybeIt = have ? it : it.skip;

describe('hff real file (optional)', () => {
    maybeIt('parses and builds a non-empty mesh', async () => {
        const buf = readFileSync(PATH!);
        const ab = buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength);
        const t0 = Date.now();
        const sff = await parseHffData(ab);
        const t1 = Date.now();
        const built = _internals.buildMesh(sff);
        const t2 = Date.now();

        // eslint-disable-next-line no-console
        console.log(`HFF: ${sff.name}; ${sff.segments.length} segments, ${built.mesh.vertexCount} verts, ${built.mesh.triangleCount} tris; parse ${t1 - t0}ms, build ${t2 - t1}ms`);

        expect(sff.primaryDescriptor).toBe('mesh_list');
        expect(sff.segments.length).toBeGreaterThan(0);
        expect(built.mesh.vertexCount).toBeGreaterThan(0);
        expect(built.mesh.triangleCount).toBeGreaterThan(0);
        for (const s of sff.segments) {
            expect(s.colour.length).toBe(4);
            expect(s.meshes.length).toBeGreaterThan(0);
        }
    });
});
