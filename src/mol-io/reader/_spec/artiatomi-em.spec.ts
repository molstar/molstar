/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseArtiatomiEm, EmHeader_Size, EmDataType, EmMachineCoding } from '../artiatomi/em';
import { createParticleListFromArtiatomiEm, getArtiatomiMotivelistTomogramIds, ArtiatomiMotivelistRowCount } from '../../../mol-model-formats/particles/em';

/** Build a minimal EM motivelist buffer with PC little-endian encoding. */
function buildMotivelistBuffer(particles: number[][]): Uint8Array {
    const n = particles.length;
    const headerBytes = EmHeader_Size;
    const dataBytes = ArtiatomiMotivelistRowCount * n * 8; // double (8 bytes per value)
    const buf = new Uint8Array(headerBytes + dataBytes);
    const view = new DataView(buf.buffer);

    // Header
    view.setUint8(0, EmMachineCoding.PC); // machine coding: PC (little-endian)
    view.setUint8(1, 0);
    view.setUint8(2, 0);
    view.setUint8(3, EmDataType.Double); // data type: double
    view.setInt32(4, ArtiatomiMotivelistRowCount, true); // DimX = 20
    view.setInt32(8, n, true); // DimY = n (number of particles)
    view.setInt32(12, 1, true); // DimZ = 1

    // Data: particle i, field j → offset (i * 20 + j) * 8
    for (let p = 0; p < n; p++) {
        for (let f = 0; f < ArtiatomiMotivelistRowCount; f++) {
            const value = particles[p][f] ?? 0;
            view.setFloat64(headerBytes + (p * ArtiatomiMotivelistRowCount + f) * 8, value, true);
        }
    }

    return buf;
}

/** Construct a particle row array (20 elements, 0-based). */
function makeRow(overrides: Partial<Record<number, number>>): number[] {
    const row = new Array(ArtiatomiMotivelistRowCount).fill(0) as number[];
    for (const [k, v] of Object.entries(overrides)) {
        const i = Number(k);
        row[i] = v || 0;
    }
    return row;
}

test('parses Artiatomi EM header correctly', async () => {
    const buf = buildMotivelistBuffer([makeRow({ 7: 100, 8: 200, 9: 300 })]);
    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    const { header } = parsed.result;
    expect(header.machineCoding).toBe(EmMachineCoding.PC);
    expect(header.dataType).toBe(EmDataType.Double);
    expect(header.dimX).toBe(ArtiatomiMotivelistRowCount);
    expect(header.dimY).toBe(1);
    expect(header.dimZ).toBe(1);
    expect(header.isLittleEndian).toBe(true);
});

test('parses motivelist data values', async () => {
    const buf = buildMotivelistBuffer([
        makeRow({ 4: 5, 7: 100, 8: 200, 9: 300, 10: 10, 11: 20, 12: 30, 16: 90, 17: 0, 18: 45, 19: 1 }),
        makeRow({ 4: 7, 7: 400, 8: 500, 9: 600, 16: 0, 17: 0, 18: 0, 19: 2 }),
    ]);
    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    const { data, header } = parsed.result;
    expect(header.dimY).toBe(2);

    // Particle 0
    expect(data[0 * 20 + 4]).toBeCloseTo(5); // tomo
    expect(data[0 * 20 + 7]).toBeCloseTo(100); // x
    expect(data[0 * 20 + 9]).toBeCloseTo(300); // z
    expect(data[0 * 20 + 16]).toBeCloseTo(90); // phi
    expect(data[0 * 20 + 18]).toBeCloseTo(45); // theta

    // Particle 1
    expect(data[1 * 20 + 4]).toBeCloseTo(7); // tomo
    expect(data[1 * 20 + 7]).toBeCloseTo(400); // x
});

test('returns an error for a buffer that is too short', async () => {
    const tiny = new Uint8Array(100);
    const parsed = await parseArtiatomiEm(tiny).run();
    expect(parsed.isError).toBe(true);
});

test('returns an error for invalid dimensions', async () => {
    const buf = new Uint8Array(EmHeader_Size);
    const view = new DataView(buf.buffer);
    view.setUint8(0, EmMachineCoding.PC);
    view.setUint8(3, EmDataType.Double);
    view.setInt32(4, -1, true); // invalid DimX
    view.setInt32(8, 1, true);
    view.setInt32(12, 1, true);
    const parsed = await parseArtiatomiEm(buf).run();
    expect(parsed.isError).toBe(true);
});

test('getArtiatomiMotivelistTomogramIds returns sorted unique tomo numbers', async () => {
    const buf = buildMotivelistBuffer([
        makeRow({ 4: 3 }),
        makeRow({ 4: 1 }),
        makeRow({ 4: 3 }),
        makeRow({ 4: 7 }),
    ]);
    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    const ids = getArtiatomiMotivelistTomogramIds(parsed.result);
    expect(ids).toEqual([1, 3, 7]);
});

test('createParticleListFromArtiatomiEm builds correct coordinates with pixelSize', async () => {
    const buf = buildMotivelistBuffer([
        makeRow({ 4: 1, 7: 100, 8: 200, 9: 300, 10: 10, 11: 20, 12: 0, 16: 0, 17: 0, 18: 0 }),
    ]);
    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    const list = createParticleListFromArtiatomiEm(parsed.result, { pixelSize: 2 });
    expect(list.count).toBe(1);
    // coord = (x - dx) * pixelSize = (100 - 10) * 2 = 180
    expect(list.coordinates[0]).toBeCloseTo(180);
    // coord = (y - dy) * pixelSize = (200 - 20) * 2 = 360
    expect(list.coordinates[1]).toBeCloseTo(360);
    // coord = (z - dz) * pixelSize = (300 - 0) * 2 = 600
    expect(list.coordinates[2]).toBeCloseTo(600);
});

test('createParticleListFromArtiatomiEm filters by tomogram number', async () => {
    const buf = buildMotivelistBuffer([
        makeRow({ 4: 1, 7: 10, 8: 20, 9: 30 }),
        makeRow({ 4: 2, 7: 40, 8: 50, 9: 60 }),
        makeRow({ 4: 1, 7: 70, 8: 80, 9: 90 }),
    ]);
    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    const list = createParticleListFromArtiatomiEm(parsed.result, { pixelSize: 1, tomos: [1] });
    expect(list.count).toBe(2);
    expect(list.label).toBe('Particles (tomo 1)');
});

test('createParticleListFromArtiatomiEm assigns identity quaternion for zero angles', async () => {
    const buf = buildMotivelistBuffer([makeRow({ 7: 1, 8: 2, 9: 3, 16: 0, 17: 0, 18: 0 })]);
    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    const list = createParticleListFromArtiatomiEm(parsed.result, { pixelSize: 1 });
    expect(list.rotations).toBeDefined();
    // Identity quaternion: [0, 0, 0, 1] (x, y, z, w)
    expect(list.rotations![0]).toBeCloseTo(0); // x
    expect(list.rotations![1]).toBeCloseTo(0); // y
    expect(list.rotations![2]).toBeCloseTo(0); // z
    expect(list.rotations![3]).toBeCloseTo(1); // w
});

test('createParticleListFromArtiatomiEm throws when no particles match filter', async () => {
    const buf = buildMotivelistBuffer([makeRow({ 4: 1 })]);
    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    expect(() => createParticleListFromArtiatomiEm(parsed.result, { pixelSize: 1, tomos: [99] }))
        .toThrow(/No Artiatomi motivelist particles matched/);
});

test('createParticleListFromArtiatomiEm throws for non-motivelist EM file', async () => {
    // Build a 10 x 5 x 1 EM file (DimX ≠ 20)
    const n = 5;
    const headerBytes = EmHeader_Size;
    const buf = new Uint8Array(headerBytes + 10 * n * 8);
    const view = new DataView(buf.buffer);
    view.setUint8(0, EmMachineCoding.PC);
    view.setUint8(3, EmDataType.Double);
    view.setInt32(4, 10, true);
    view.setInt32(8, n, true);
    view.setInt32(12, 1, true);

    const parsed = await parseArtiatomiEm(buf).run();
    if (parsed.isError) throw new Error(parsed.message);

    expect(() => createParticleListFromArtiatomiEm(parsed.result, { pixelSize: 1 }))
        .toThrow(/does not appear to be a motivelist/);
});
