/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from '../../mol-util';
import { JSONCifDataBlock } from './model';
import { arrayMapAdd } from '../../mol-util/map';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { Table } from '../../mol-data/db';
import { MolstarBondSiteTypeId, MolstarBondSiteValueOrder } from '../../mol-model/structure/export/categories/molstar_bond_site';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { VdwRadius } from '../../mol-model/structure/model/properties/atomic';
import { ElementSymbol } from '../../mol-model/structure/model/types';

type Atom = Partial<Table.Row<mmCIF_Schema['atom_site']>>

export interface JSONCifLigandGraphBondProps {
    value_order: MolstarBondSiteValueOrder | undefined;
    type_id: MolstarBondSiteTypeId | undefined;
}

export interface JSONCifLigandGraphAtom {
    key: string;
    final_id: number | undefined;
    row: Atom;
}

export interface JSONCifLigandGraphBond {
    atom_1: JSONCifLigandGraphAtom;
    atom_2: JSONCifLigandGraphAtom;
    props: JSONCifLigandGraphBondProps;
}

export interface JSONCifLigandGraphData {
    block: JSONCifDataBlock;
    atomIdRemapping: Map<number, number>;
    addedAtomIds: number[];
    removedAtomIds: number[];
}

const _State = {
    p1: Vec3(),
    p2: Vec3(),
};

export class JSONCifLigandGraph {
    readonly atoms: JSONCifLigandGraphAtom[] = [];
    readonly atomsById: Map<number, JSONCifLigandGraphAtom> = new Map();
    readonly bondByKey: Map<string, JSONCifLigandGraphBond[]> = new Map();
    readonly removedAtomIds: Set<number> = new Set();

    getAtomAtIndex(index: number) {
        return this.atoms[index];
    }

    getAtom(atomOrId: number | JSONCifLigandGraphAtom) {
        return typeof atomOrId === 'number' ? this.atomsById.get(atomOrId) : atomOrId;
    }

    getBonds(atomOrId: number | JSONCifLigandGraphAtom) {
        const atom = this.getAtom(atomOrId);
        if (!atom) return [];
        return this.bondByKey.get(atom.key) ?? [];
    }

    getAtomCoords(atomOrId: number | JSONCifLigandGraphAtom, out: Vec3 = Vec3()) {
        const atom = this.getAtom(atomOrId);
        if (!atom) return out;
        const { Cartn_x, Cartn_y, Cartn_z } = atom.row;
        return Vec3.set(out, Cartn_x!, Cartn_y!, Cartn_z!);
    }

    getBondDirection(bond: JSONCifLigandGraphBond, out: Vec3 = Vec3()) {
        const a1 = this.getAtomCoords(bond.atom_1, _State.p1);
        const a2 = this.getAtomCoords(bond.atom_2, _State.p2);
        const dir = Vec3.sub(out, a2, a1);
        return dir;
    }

    modifyAtom(id: number, data: Atom) {
        const atom = this.atomsById.get(id);
        if (!atom) return;
        atom.row = { ...atom.row, ...data };
    }

    addAtom(data: Omit<Atom, 'id'>) {
        const atom: JSONCifLigandGraphAtom = {
            key: UUID.create22(),
            final_id: undefined,
            row: { ...data, id: undefined },
        };
        this.atoms.push(atom);
        return atom;
    }

    transformCoords(xform: Mat4) {
        for (const a of this.atoms) {
            const p = this.getAtomCoords(a, _State.p1);
            Vec3.transformMat4(p, p, xform);
            a.row.Cartn_x = p[0];
            a.row.Cartn_y = p[1];
            a.row.Cartn_z = p[2];
        }
    }

    attachAtom(parent: number | JSONCifLigandGraphAtom, atom: Atom) {
        const p = this.getAtom(parent);
        if (!p) return;

        const c = this.getAtomCoords(p);
        const dir = this.getAddAtomDirection(p);
        const r = 2 / 5 * (VdwRadius(ElementSymbol(p.row.type_symbol ?? 'C')) + VdwRadius(ElementSymbol(atom.type_symbol ?? 'C')));
        const newAtom = this.addAtom({
            ...p.row,
            // NOTE: this is not correct for editing protein atoms
            // as they should have atom names from CCD, or at least the should be
            // unique. This should be fine for small ligand editing.
            auth_atom_id: atom.type_symbol,
            label_atom_id: atom.type_symbol,
            ...atom,
            Cartn_x: c[0] + dir[0] * r,
            Cartn_y: c[1] + dir[1] * r,
            Cartn_z: c[2] + dir[2] * r
        });
        this.addOrUpdateBond(p, newAtom, { value_order: 'sing', type_id: 'covale' });

        return newAtom;
    }

    removeAtom(atomOrId: number | JSONCifLigandGraphAtom) {
        const atom = this.getAtom(atomOrId);
        if (!atom) return;
        if (typeof atom.row.id === 'number') {
            this.removedAtomIds.add(atom.row.id);
            this.atomsById.delete(atom.row.id);
        }

        this.atoms.splice(this.atoms.indexOf(atom), 1);

        const bonds = this.bondByKey.get(atom.key);
        if (!bonds) return;
        this.bondByKey.delete(atom.key);

        for (const b of bonds) {
            const bBonds = this.bondByKey.get(b.atom_2.key);
            if (!bBonds) continue;
            this.bondByKey.set(b.atom_2.key, bBonds.filter(bb => bb.atom_2 !== atom));
        }
    }

    addOrUpdateBond(atom1: number | JSONCifLigandGraphAtom, atom2: number | JSONCifLigandGraphAtom, props: JSONCifLigandGraphBondProps) {
        const a1 = this.getAtom(atom1);
        const a2 = this.getAtom(atom2);
        if (!a1 || !a2) return;

        const ps = { ...props };
        this.removeBond(atom1, atom2);
        arrayMapAdd(this.bondByKey, a1.key, { atom_1: a1, atom_2: a2, props: ps });
        arrayMapAdd(this.bondByKey, a2.key, { atom_1: a2, atom_2: a1, props: ps });
    }

    // TODO: iterate on this?
    addBondPart(atom1: number | JSONCifLigandGraphAtom, atom2: number | JSONCifLigandGraphAtom, props: JSONCifLigandGraphBondProps) {
        const a1 = this.getAtom(atom1);
        const a2 = this.getAtom(atom2);
        if (!a1 || !a2) return;
        arrayMapAdd(this.bondByKey, a1.key, { atom_1: a1, atom_2: a2, props: { ...props } });
    }

    removeBond(atom1: number | JSONCifLigandGraphAtom, atom2: number | JSONCifLigandGraphAtom) {
        const a1 = this.getAtom(atom1);
        const a2 = this.getAtom(atom2);
        if (!a1 || !a2) return;
        const a1Bonds = this.bondByKey.get(a1.key);
        if (a1Bonds) {
            this.bondByKey.set(a1.key, a1Bonds.filter(b => b.atom_2 !== a2));
        }
        const a2Bonds = this.bondByKey.get(a2.key);
        if (a2Bonds) {
            this.bondByKey.set(a2.key, a2Bonds.filter(b => b.atom_2 !== a1));
        }
    }

    getData(): JSONCifLigandGraphData {
        const atomIdRemapping = new Map<number, number>();
        const addedAtomIds: number[] = [];
        const atoms: Atom[] = [];

        for (let i = 0; i < this.atoms.length; ++i) {
            const id = i + 1;

            if (this.atoms[i].row.id === undefined) {
                addedAtomIds.push(id);
            } else {
                atomIdRemapping.set(this.atoms[i].row.id!, id);
            }

            this.atoms[i].final_id = id;
            atoms.push({ ...this.atoms[i].row, id });
        }

        const block: JSONCifDataBlock = {
            ...this.data,
            categories: {
                ...this.data.categories,
                atom_site: {
                    ...this.data.categories['atom_site'],
                    rows: atoms,
                },
            }
        };

        const bonds: Record<string, any>[] = [];
        for (const a of this.atoms) {
            const xs = this.bondByKey.get(a.key);
            if (!xs) continue;
            for (const bb of xs) {
                if (a.final_id! >= bb.atom_2.final_id!) continue;

                bonds.push({
                    atom_id_1: a.final_id,
                    atom_id_2: bb.atom_2.final_id,
                    value_order: bb.props.value_order,
                    type_id: bb.props.type_id,
                });
            }
        }

        if (block.categories.molstar_bond_site) {
            block.categories['molstar_bond_site'] = {
                ...block.categories['molstar_bond_site'],
                rows: bonds
            };
        }

        return {
            block,
            atomIdRemapping,
            addedAtomIds,
            removedAtomIds: Array.from(this.removedAtomIds).sort((a, b) => a - b),
        };
    }

    private getAddAtomDirection(parent: JSONCifLigandGraphAtom) {
        // NOTE: this will not correctly handle all cases...

        let deltas: Vec3[] = [];
        const bonds = this.bondByKey.get(parent.key);
        if (!bonds?.length) return Vec3.create(1, 0, 0);

        const c = this.getAtomCoords(parent);
        for (const b of bonds) {
            const delta = Vec3.sub(Vec3(), this.getAtomCoords(b.atom_2), c);
            deltas.push(delta);
        }

        if (deltas.length === 1) {
            const ret = Vec3.negate(Vec3(), deltas[0]);
            Vec3.normalize(ret, ret);
            return ret;
        }

        if (deltas.length === 2) {
            const ret = Vec3.add(Vec3(), deltas[0], deltas[1]);
            Vec3.normalize(ret, ret);
            Vec3.negate(ret, ret);
            return ret;
        }

        // For now, just take the first three deltas and cross-product them
        deltas = deltas.slice(0, 3);
        const crossProducts: Vec3[] = [];
        for (let i = 0; i < deltas.length; ++i) {
            for (let j = i + 1; j < deltas.length; ++j) {
                const cross = Vec3.cross(Vec3(), deltas[i], deltas[j]);
                Vec3.normalize(cross, cross);
                crossProducts.push(cross);
            }
        }
        for (let i = 1; i < crossProducts.length; ++i) {
            Vec3.matchDirection(crossProducts[i], crossProducts[i], crossProducts[0]);
        }

        const avg = Vec3.create(0, 0, 0);
        for (const cp of crossProducts) {
            Vec3.add(avg, avg, cp);
        }
        Vec3.normalize(avg, avg);
        return avg;
    }

    constructor(private data: JSONCifDataBlock) {
        for (const row of data.categories['atom_site'].rows) {
            const atom: JSONCifLigandGraphAtom = {
                key: UUID.create22(),
                final_id: row.final_id,
                row: { ...row },
            };
            this.atoms.push(atom);
            this.atomsById.set(row.id, atom);
        }

        if (!data.categories.molstar_bond_site) return;

        for (const row of data.categories.molstar_bond_site.rows) {
            const atom_1 = this.atomsById.get(row.atom_id_1);
            const atom_2 = this.atomsById.get(row.atom_id_2);
            if (!atom_1 || !atom_2) continue;

            arrayMapAdd(this.bondByKey, atom_1.key, {
                atom_1: atom_1,
                atom_2: atom_2,
                props: {
                    value_order: row.value_order,
                    type_id: row.type_id,
                },
            });
            arrayMapAdd(this.bondByKey, atom_2.key, {
                atom_1: atom_2,
                atom_2: atom_1,
                props: {
                    value_order: row.value_order,
                    type_id: row.type_id,
                }
            });
        }
    }
}