/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Yujie Wu <yujie@atommap.com>
 */

import { utf8ByteCount, utf8Write } from '../common/utf8';
import { ZmlAtom, ZmlBond, ZmlMol, ZmlMolSys, ZmlResidue, ZmlV1FileNames } from '../reader/zml/schema';
import { Structure, StructureElement, StructureProperties, Unit } from '../../mol-model/structure';
import { AtomNumber } from '../../mol-model/structure/model/properties/atomic/measures';
import { RuntimeContext } from '../../mol-task';
import { zip } from '../../mol-util/zip/zip';

function stringToBytes(s: string): Uint8Array<ArrayBuffer> {
    const bytes = new Uint8Array(utf8ByteCount(s));
    utf8Write(bytes, 0, s);
    return bytes as Uint8Array<ArrayBuffer>;
}

function structureToMolSys(name: string, structure: Structure): { molsys: ZmlMolSys, positions: Float64Array } {
    const loc = StructureElement.Location.create(structure);
    const P = StructureProperties;

    const molecules: ZmlMol[] = [];
    const positions: number[] = [];
    /** "unitId:elementIndex" -> flat ZML atom index, for bond lookup. */
    const indexMap = new Map<string, number>();
    let flatIndex = 0;

    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;
        loc.unit = unit;
        const elements = unit.elements;

        let currentMol: { mol: ZmlMol & { residues: ZmlResidue[] }, chainKey: number } | undefined;
        let currentRes: { res: ZmlResidue & { atoms: ZmlAtom[] }, resKey: number } | undefined;

        for (let i = 0; i < elements.length; i++) {
            loc.element = elements[i];
            const chainKey = P.chain.key(loc);
            const resKey = P.residue.key(loc);

            if (!currentMol || currentMol.chainKey !== chainKey) {
                const asymId = P.chain.auth_asym_id(loc);
                const mol: ZmlMol & { residues: ZmlResidue[] } = {
                    name: asymId, chain_id: asymId, comtype: 'Comtype.UNKNOWN',
                    residues: [], props: {},
                };
                molecules.push(mol);
                currentMol = { mol, chainKey };
                currentRes = undefined;
            }

            if (!currentRes || currentRes.resKey !== resKey) {
                const res: ZmlResidue & { atoms: ZmlAtom[] } = {
                    id: P.residue.auth_seq_id(loc),
                    name: P.residue.auth_comp_id(loc),
                    insertion_code: P.residue.pdbx_PDB_ins_code(loc) || '',
                    atoms: [], props: {},
                };
                currentMol.mol.residues.push(res);
                currentRes = { res, resKey };
            }

            const typeSymbol = P.atom.type_symbol(loc);
            const atom: ZmlAtom = {
                id: P.atom.id(loc),
                name: P.atom.auth_atom_id(loc),
                atomic_number: AtomNumber(typeSymbol),
                formal_charge: P.atom.pdbx_formal_charge(loc) || 0,
                props: {},
            };
            currentRes.res.atoms.push(atom);
            positions.push(P.atom.x(loc), P.atom.y(loc), P.atom.z(loc));
            indexMap.set(`${unit.id}:${elements[i]}`, flatIndex);
            flatIndex++;
        }
    }

    const bonds = collectBonds(structure, indexMap);
    const molsys: ZmlMolSys = { name, molecules, bonds, box: null, props: {} };
    return { molsys, positions: Float64Array.from(positions) };
}

function collectBonds(structure: Structure, indexMap: Map<string, number>): ZmlBond[] {
    const bonds: ZmlBond[] = [];
    const added = new Set<string>();
    const push = (a: number, b: number, order: number) => {
        if (a === b) return;
        const lo = a < b ? a : b, hi = a < b ? b : a;
        const key = `${lo}:${hi}`;
        if (added.has(key)) return;
        added.add(key);
        bonds.push([lo, hi, order] as const);
    };

    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;
        const { a, b, edgeProps } = unit.bonds;
        const elements = unit.elements;
        for (let i = 0; i < a.length; i++) {
            const ai = indexMap.get(`${unit.id}:${elements[a[i]]}`);
            const bi = indexMap.get(`${unit.id}:${elements[b[i]]}`);
            if (ai === undefined || bi === undefined) continue;
            push(ai, bi, edgeProps.order[i] ?? 1);
        }
    }
    for (const e of structure.interUnitBonds.edges) {
        const ua = structure.unitMap.get(e.unitA), ub = structure.unitMap.get(e.unitB);
        const ai = indexMap.get(`${ua.id}:${ua.elements[e.indexA]}`);
        const bi = indexMap.get(`${ub.id}:${ub.elements[e.indexB]}`);
        if (ai === undefined || bi === undefined) continue;
        push(ai, bi, e.props.order ?? 1);
    }

    bonds.sort((x, y) => x[0] - y[0] || x[1] - y[1]);
    return bonds;
}

/** Encode a `Structure` as a ZML V1 archive (ZIP) and return its byte buffer. */
export async function encodeZml(runtime: RuntimeContext, name: string, structure: Structure): Promise<Uint8Array<ArrayBuffer>> {
    const { molsys, positions } = structureToMolSys(name, structure);
    const molsysBytes = stringToBytes(JSON.stringify(molsys));
    const posBuffer = new ArrayBuffer(positions.byteLength);
    new Float64Array(posBuffer).set(positions);
    const posBytes: Uint8Array<ArrayBuffer> = new Uint8Array(posBuffer);
    const versionBytes = stringToBytes('1');

    const entries: { [k: string]: Uint8Array<ArrayBuffer> } = {};
    entries[ZmlV1FileNames.molsys] = molsysBytes;
    entries[ZmlV1FileNames.positions] = posBytes;
    entries[ZmlV1FileNames.version] = versionBytes;
    const buf = await zip(runtime, entries);
    return new Uint8Array(buf) as Uint8Array<ArrayBuffer>;
}
