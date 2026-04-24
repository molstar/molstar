/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Augment Agent
 */

import { utf8ByteCount, utf8Write } from '../../common/utf8';
import { SyncRuntimeContext } from '../../../mol-task/execution/synchronous';
import { Structure } from '../../../mol-model/structure';
import { trajectoryFromZml } from '../../../mol-model-formats/structure/zml';
import { encodeZml } from '../../writer/zml';
import { zip } from '../../../mol-util/zip/zip';
import { parseZml } from '../zml/parser';
import { ZmlMolSys, ZmlV1FileNames, ZmlV3FileNames, ZmlV3Header } from '../zml/schema';

function toBytes(s: string): Uint8Array<ArrayBuffer> {
    const bytes = new Uint8Array(utf8ByteCount(s));
    utf8Write(bytes, 0, s);
    return bytes as Uint8Array<ArrayBuffer>;
}

function float64Bytes(values: number[]): Uint8Array<ArrayBuffer> {
    const buf = new ArrayBuffer(values.length * 8);
    new Float64Array(buf).set(values);
    return new Uint8Array(buf) as Uint8Array<ArrayBuffer>;
}

function waterMolSys(): ZmlMolSys {
    return {
        name: 'water',
        molecules: [{
            name: 'A', chain_id: 'A', comtype: 'Comtype.UNKNOWN',
            residues: [{
                id: 1, name: 'HOH', insertion_code: '',
                atoms: [
                    { id: 1, name: 'O', atomic_number: 8, formal_charge: 0, props: {} },
                    { id: 2, name: 'H1', atomic_number: 1, formal_charge: 0, props: {} },
                    { id: 3, name: 'H2', atomic_number: 1, formal_charge: 0, props: {} },
                ],
                props: {},
            }],
            props: {},
        }],
        bonds: [[0, 1, 1], [0, 2, 1]],
        box: null,
        props: {},
    };
}

async function buildArchive(molsys: ZmlMolSys, positions: number[]): Promise<Uint8Array> {
    const entries: { [k: string]: Uint8Array<ArrayBuffer> } = {
        [ZmlV1FileNames.version]: toBytes('1'),
        [ZmlV1FileNames.molsys]: toBytes(JSON.stringify(molsys)),
        [ZmlV1FileNames.positions]: float64Bytes(positions),
    };
    const buf = await zip(SyncRuntimeContext, entries);
    return new Uint8Array(buf);
}

function int16Bytes(values: number[]): Uint8Array<ArrayBuffer> {
    const buf = new ArrayBuffer(values.length * 2);
    new Int16Array(buf).set(values);
    return new Uint8Array(buf) as Uint8Array<ArrayBuffer>;
}
function int32Bytes(values: number[]): Uint8Array<ArrayBuffer> {
    const buf = new ArrayBuffer(values.length * 4);
    new Int32Array(buf).set(values);
    return new Uint8Array(buf) as Uint8Array<ArrayBuffer>;
}
function int8Bytes(values: number[]): Uint8Array<ArrayBuffer> {
    const buf = new Uint8Array(values.length);
    new Int8Array(buf.buffer).set(values);
    return buf as Uint8Array<ArrayBuffer>;
}

async function buildV3Archive(): Promise<{ data: Uint8Array, positions: number[] }> {
    const header: ZmlV3Header = {
        name: 'water',
        box: null,
        props: {},
        atom_names: ['O', 'H1', 'H2'],
        atom_props: [{}, {}, {}],
        res_names: ['HOH'],
        res_icodes: [''],
        res_props: [{}],
        mol_names: ['A'],
        mol_chain_ids: ['A'],
        mol_props: [{}],
    };
    const positions = [0, 0, 0, 0.96, 0, 0, -0.24, 0.93, 0];
    const entries: { [k: string]: Uint8Array<ArrayBuffer> } = {
        [ZmlV3FileNames.version]: toBytes('3'),
        [ZmlV3FileNames.header]: toBytes(JSON.stringify(header)),
        [ZmlV3FileNames.Z]: int16Bytes([8, 1, 1]),
        [ZmlV3FileNames.formalCharges]: int16Bytes([0, 0, 0]),
        [ZmlV3FileNames.atomIds]: int32Bytes([1, 2, 3]),
        [ZmlV3FileNames.resIds]: int32Bytes([1]),
        [ZmlV3FileNames.resAtomOffsets]: int32Bytes([0, 3]),
        [ZmlV3FileNames.molResOffsets]: int32Bytes([0, 1]),
        [ZmlV3FileNames.bonds]: int32Bytes([0, 1, 0, 2]),
        [ZmlV3FileNames.bondOrders]: int8Bytes([1, 1]),
        [ZmlV3FileNames.positions]: float64Bytes(positions),
    };
    const buf = await zip(SyncRuntimeContext, entries);
    return { data: new Uint8Array(buf), positions };
}

describe('zml reader', () => {
    it('parses a minimal V1 archive', async () => {
        const molsys = waterMolSys();
        const positions = [0, 0, 0, 0.96, 0, 0, -0.24, 0.93, 0];
        const data = await buildArchive(molsys, positions);
        const parsed = await parseZml(data, 'test').run();
        if (parsed.isError) throw new Error(parsed.message);
        const zml = parsed.result;

        expect(zml.atomCount).toBe(3);
        expect(zml.molsys.molecules.length).toBe(1);
        expect(zml.molsys.molecules[0].residues[0].atoms[0].name).toBe('O');
        expect(zml.molsys.bonds.length).toBe(2);
        expect(Array.from(zml.positions.slice(0, 3))).toEqual([0, 0, 0]);
        expect(zml.positions[3]).toBeCloseTo(0.96);
    });

    it('rejects unsupported versions', async () => {
        const entries: { [k: string]: Uint8Array<ArrayBuffer> } = {
            [ZmlV1FileNames.version]: toBytes('2'),
            'molsys.pkl': toBytes('irrelevant'),
        };
        const buf = await zip(SyncRuntimeContext, entries);
        const parsed = await parseZml(new Uint8Array(buf), 'bad').run();
        expect(parsed.isError).toBe(true);
        if (parsed.isError) expect(parsed.message).toMatch(/version/i);
    });

    it('parses a minimal V3 archive', async () => {
        const { data, positions } = await buildV3Archive();
        const parsed = await parseZml(data, 'v3').run();
        if (parsed.isError) throw new Error(parsed.message);
        const zml = parsed.result;

        expect(zml.atomCount).toBe(3);
        expect(zml.molsys.name).toBe('water');
        expect(zml.molsys.molecules.length).toBe(1);
        expect(zml.molsys.molecules[0].chain_id).toBe('A');
        const atoms = zml.molsys.molecules[0].residues[0].atoms;
        expect(atoms.map(a => a.name)).toEqual(['O', 'H1', 'H2']);
        expect(atoms.map(a => a.atomic_number)).toEqual([8, 1, 1]);
        expect(zml.molsys.bonds).toEqual([[0, 1, 1], [0, 2, 1]]);
        for (let i = 0; i < positions.length; i++) {
            expect(zml.positions[i]).toBeCloseTo(positions[i], 6);
        }
    });

    it('builds a Trajectory from a V3 archive', async () => {
        const { data } = await buildV3Archive();
        const parsed = await parseZml(data, 'v3').run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromZml(parsed.result).run();
        expect(trajectory.frameCount).toBe(1);
        const model = trajectory.representative;
        expect(model.atomicHierarchy.atoms._rowCount).toBe(3);
        expect(model.atomicHierarchy.atoms.type_symbol.value(0)).toBe('O');
        expect(model.atomicHierarchy.atoms.type_symbol.value(1)).toBe('H');
    });

    it('rejects positions with wrong length', async () => {
        const molsys = waterMolSys();
        const data = await buildArchive(molsys, [0, 0, 0]); // only 1 atom worth for 3 atoms
        const parsed = await parseZml(data, 'bad').run();
        expect(parsed.isError).toBe(true);
    });

    it('builds a Trajectory with correct elements', async () => {
        const molsys = waterMolSys();
        const data = await buildArchive(molsys, [0, 0, 0, 0.96, 0, 0, -0.24, 0.93, 0]);
        const parsed = await parseZml(data, 'traj').run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromZml(parsed.result).run();
        expect(trajectory.frameCount).toBe(1);
        const model = trajectory.representative;
        expect(model.atomicHierarchy.atoms._rowCount).toBe(3);
        expect(model.atomicHierarchy.atoms.type_symbol.value(0)).toBe('O');
        expect(model.atomicHierarchy.atoms.type_symbol.value(1)).toBe('H');
    });
});

describe('zml writer', () => {
    it('roundtrips through encodeZml', async () => {
        const originalMolSys = waterMolSys();
        const originalPositions = [0, 0, 0, 0.96, 0, 0, -0.24, 0.93, 0];
        const data = await buildArchive(originalMolSys, originalPositions);
        const parsed = await parseZml(data, 'rt').run();
        if (parsed.isError) throw new Error(parsed.message);
        const trajectory = await trajectoryFromZml(parsed.result).run();
        const structure = Structure.ofModel(trajectory.representative);

        const encoded = await encodeZml(SyncRuntimeContext, 'rt', structure);
        const reparsed = await parseZml(encoded, 'rt').run();
        if (reparsed.isError) throw new Error(reparsed.message);
        const rt = reparsed.result;

        expect(rt.atomCount).toBe(3);
        // Atom ordering should round-trip.
        const outAtoms = rt.molsys.molecules[0].residues[0].atoms.map(a => a.name);
        expect(outAtoms).toEqual(['O', 'H1', 'H2']);
        expect(rt.molsys.bonds.length).toBe(2);
        for (let i = 0; i < originalPositions.length; i++) {
            expect(rt.positions[i]).toBeCloseTo(originalPositions[i], 4);
        }
    });
});
