/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Augment Agent
 */

import { utf8Read } from '../../common/utf8';
import { RuntimeContext, Task } from '../../../mol-task';
import { unzip } from '../../../mol-util/zip/zip';
import { ReaderResult as Result } from '../result';
import {
    ZmlAtom, ZmlBond, ZmlFile, ZmlMol, ZmlMolSys, ZmlResidue,
    ZmlSupportedVersions, ZmlV1FileNames, ZmlV3FileNames, ZmlV3Header,
} from './schema';

type ZipEntries = { [k: string]: Uint8Array | { size: number, csize: number } };

function getEntry(entries: ZipEntries, fname: string): Uint8Array | undefined {
    const entry = entries[fname];
    return entry instanceof Uint8Array ? entry : undefined;
}

function countAtoms(molsys: ZmlMolSys): number {
    let n = 0;
    for (const mol of molsys.molecules) {
        for (const res of mol.residues) n += res.atoms.length;
    }
    return n;
}

function alignedCopy(raw: Uint8Array): ArrayBuffer {
    const buf = new ArrayBuffer(raw.byteLength);
    new Uint8Array(buf).set(raw);
    return buf;
}

function readFloat64Positions(raw: Uint8Array, atomCount: number, entryName: string): Float64Array {
    const expected = atomCount * 3 * 8;
    if (raw.byteLength !== expected) {
        throw new Error(
            `ZML ${entryName} has ${raw.byteLength} bytes but expected ${expected} ` +
            `(atomCount=${atomCount}, 3 coords, float64).`
        );
    }
    return new Float64Array(alignedCopy(raw));
}

function readTypedArray<T extends Int8Array | Int16Array | Int32Array>(
    raw: Uint8Array, Ctor: { new (buf: ArrayBuffer): T, BYTES_PER_ELEMENT: number }, entryName: string
): T {
    if (raw.byteLength % Ctor.BYTES_PER_ELEMENT !== 0) {
        throw new Error(
            `ZML ${entryName} byte length ${raw.byteLength} is not a multiple of ` +
            `${Ctor.BYTES_PER_ELEMENT}.`
        );
    }
    return new Ctor(alignedCopy(raw));
}

async function parseV1(entries: ZipEntries, name: string): Promise<Result<ZmlFile>> {
    const molsysBytes = getEntry(entries, ZmlV1FileNames.molsys);
    if (!molsysBytes) return Result.error(`Missing ${ZmlV1FileNames.molsys} in ZML archive.`);
    const posBytes = getEntry(entries, ZmlV1FileNames.positions);
    if (!posBytes) return Result.error(`Missing ${ZmlV1FileNames.positions} in ZML archive.`);

    let molsys: ZmlMolSys;
    try {
        molsys = JSON.parse(utf8Read(molsysBytes));
    } catch (e) {
        return Result.error(`Failed to parse molsys.json: ${e instanceof Error ? e.message : String(e)}`);
    }
    if (!molsys || !Array.isArray(molsys.molecules) || !Array.isArray(molsys.bonds)) {
        return Result.error('Malformed molsys.json: expected MolSys object with molecules and bonds.');
    }

    const atomCount = countAtoms(molsys);
    let positions: Float64Array;
    try {
        positions = readFloat64Positions(posBytes, atomCount, ZmlV1FileNames.positions);
    } catch (e) {
        return Result.error(e instanceof Error ? e.message : String(e));
    }

    return Result.success({ name, molsys, positions, atomCount });
}

function buildMolSysFromV3(header: ZmlV3Header, atomIds: Int32Array, Z: Int16Array, fc: Int16Array,
    resIds: Int32Array, resAtomOffsets: Int32Array, molResOffsets: Int32Array,
    bondsIJ: Int32Array, bondOrders: Int8Array): ZmlMolSys {

    const numAtoms = header.atom_names.length;
    if (Z.length !== numAtoms || fc.length !== numAtoms || atomIds.length !== numAtoms) {
        throw new Error(`ZML V3 atom columns length mismatch (expected ${numAtoms}).`);
    }
    const numRes = header.res_names.length;
    if (resIds.length !== numRes || resAtomOffsets.length !== numRes + 1) {
        throw new Error(`ZML V3 residue columns length mismatch (expected ${numRes}).`);
    }
    const numMols = header.mol_names.length;
    if (molResOffsets.length !== numMols + 1) {
        throw new Error(`ZML V3 molecule offsets length mismatch (expected ${numMols + 1}).`);
    }

    const atoms: ZmlAtom[] = new Array(numAtoms);
    for (let i = 0; i < numAtoms; i++) {
        atoms[i] = {
            id: atomIds[i],
            name: header.atom_names[i],
            atomic_number: Z[i],
            formal_charge: fc[i],
            props: header.atom_props[i] ?? {},
        };
    }

    const residues: (ZmlResidue & { atoms: ZmlAtom[] })[] = new Array(numRes);
    for (let r = 0; r < numRes; r++) {
        const a0 = resAtomOffsets[r];
        const a1 = resAtomOffsets[r + 1];
        residues[r] = {
            id: resIds[r],
            name: header.res_names[r],
            insertion_code: header.res_icodes[r] ?? '',
            atoms: atoms.slice(a0, a1),
            props: header.res_props[r] ?? {},
        };
    }

    const molecules: ZmlMol[] = new Array(numMols);
    for (let m = 0; m < numMols; m++) {
        const r0 = molResOffsets[m];
        const r1 = molResOffsets[m + 1];
        molecules[m] = {
            name: header.mol_names[m],
            chain_id: header.mol_chain_ids[m],
            comtype: 'Comtype.UNKNOWN',
            residues: residues.slice(r0, r1),
            props: header.mol_props[m] ?? {},
        };
    }

    const bondCount = bondOrders.length;
    if (bondsIJ.length !== bondCount * 2) {
        throw new Error(`ZML V3 bonds length ${bondsIJ.length} != 2 * ${bondCount}.`);
    }
    const bonds: ZmlBond[] = new Array(bondCount);
    for (let k = 0; k < bondCount; k++) {
        bonds[k] = [bondsIJ[2 * k], bondsIJ[2 * k + 1], bondOrders[k]];
    }

    return { name: header.name, molecules, bonds, box: header.box, props: header.props ?? {} };
}

async function parseV3(entries: ZipEntries, name: string): Promise<Result<ZmlFile>> {
    const headerBytes = getEntry(entries, ZmlV3FileNames.header);
    if (!headerBytes) return Result.error(`Missing ${ZmlV3FileNames.header} in ZML archive.`);

    let header: ZmlV3Header;
    try {
        header = JSON.parse(utf8Read(headerBytes));
    } catch (e) {
        return Result.error(`Failed to parse ${ZmlV3FileNames.header}: ${e instanceof Error ? e.message : String(e)}`);
    }
    if (!header || !Array.isArray(header.atom_names) || !Array.isArray(header.mol_names)) {
        return Result.error(`Malformed ${ZmlV3FileNames.header}: missing atom_names/mol_names.`);
    }

    const required: ReadonlyArray<keyof typeof ZmlV3FileNames> = [
        'Z', 'formalCharges', 'atomIds', 'resIds', 'resAtomOffsets', 'molResOffsets',
        'bonds', 'bondOrders', 'positions',
    ];
    const bufs: Partial<Record<keyof typeof ZmlV3FileNames, Uint8Array>> = {};
    for (const key of required) {
        const fname = ZmlV3FileNames[key];
        const raw = getEntry(entries, fname);
        if (!raw) return Result.error(`Missing ${fname} in ZML archive.`);
        bufs[key] = raw;
    }

    let molsys: ZmlMolSys;
    let positions: Float64Array;
    try {
        const atomIds = readTypedArray(bufs.atomIds!, Int32Array, ZmlV3FileNames.atomIds);
        const Z = readTypedArray(bufs.Z!, Int16Array, ZmlV3FileNames.Z);
        const fc = readTypedArray(bufs.formalCharges!, Int16Array, ZmlV3FileNames.formalCharges);
        const resIds = readTypedArray(bufs.resIds!, Int32Array, ZmlV3FileNames.resIds);
        const resAtomOffsets = readTypedArray(bufs.resAtomOffsets!, Int32Array, ZmlV3FileNames.resAtomOffsets);
        const molResOffsets = readTypedArray(bufs.molResOffsets!, Int32Array, ZmlV3FileNames.molResOffsets);
        const bondsIJ = readTypedArray(bufs.bonds!, Int32Array, ZmlV3FileNames.bonds);
        const bondOrders = readTypedArray(bufs.bondOrders!, Int8Array, ZmlV3FileNames.bondOrders);
        molsys = buildMolSysFromV3(header, atomIds, Z, fc, resIds, resAtomOffsets, molResOffsets, bondsIJ, bondOrders);
        positions = readFloat64Positions(bufs.positions!, header.atom_names.length, ZmlV3FileNames.positions);
    } catch (e) {
        return Result.error(e instanceof Error ? e.message : String(e));
    }

    return Result.success({ name, molsys, positions, atomCount: header.atom_names.length });
}

async function parseInternal(ctx: RuntimeContext, data: Uint8Array, name: string): Promise<Result<ZmlFile>> {
    const buf = data.buffer.slice(data.byteOffset, data.byteOffset + data.byteLength);
    let entries: ZipEntries;
    try {
        entries = await unzip(ctx, buf) as ZipEntries;
    } catch (e) {
        return Result.error(`Invalid ZML file: ${e instanceof Error ? e.message : String(e)}`);
    }

    const versionBytes = getEntry(entries, 'version');
    const version = versionBytes ? parseInt(utf8Read(versionBytes).trim(), 10) : 0;
    if (!ZmlSupportedVersions.includes(version)) {
        return Result.error(
            `Unsupported ZML version ${version}. Supported versions: ` +
            `${ZmlSupportedVersions.join(', ')}.`
        );
    }

    if (version === 1) return parseV1(entries, name);
    return parseV3(entries, name);
}

export function parseZml(data: Uint8Array, name: string = 'ZML') {
    return Task.create<Result<ZmlFile>>('Parse ZML', async ctx => {
        await ctx.update('Parsing ZML archive...');
        return parseInternal(ctx, data, name);
    });
}
