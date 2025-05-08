/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Table } from '../../mol-data/db';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { MolstarBondSiteTypeId, MolstarBondSiteValueOrder } from '../../mol-model/structure/export/categories/molstar_bond_site';
import { UUID } from '../../mol-util';
import { arrayMapAdd } from '../../mol-util/map';
import { JSONCifDataBlock } from './model';

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
    readonly atomsByKey: Map<string, JSONCifLigandGraphAtom> = new Map();
    readonly atomsById: Map<number, JSONCifLigandGraphAtom> = new Map();
    /** Bond with the provided key is always atom_1 */
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

    modifyAtom(atomOrId: number | JSONCifLigandGraphAtom, data: Omit<Atom, 'id'>) {
        const atom = this.getAtom(atomOrId);
        if (!atom) return;
        atom.row = { ...atom.row, ...data, id: atom.row.id };
    }

    addAtom(data: Omit<Atom, 'id'>) {
        const atom: JSONCifLigandGraphAtom = {
            key: UUID.create22(),
            final_id: undefined,
            row: { ...data, id: undefined },
        };
        this.atomsByKey.set(atom.key, atom);
        this.atoms.push(atom);
        return atom;
    }

    removeAtom(atomOrId: number | JSONCifLigandGraphAtom) {
        const atom = this.getAtom(atomOrId);
        if (!atom) return;
        if (typeof atom.row.id === 'number') {
            this.removedAtomIds.add(atom.row.id);
            this.atomsById.delete(atom.row.id);
        }

        this.atoms.splice(this.atoms.indexOf(atom), 1);
        this.atomsByKey.delete(atom.key);

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

    private transformAtomCoords(xform: Mat4, atomOrId: number | JSONCifLigandGraphAtom) {
        const atom = this.getAtom(atomOrId);
        if (!atom) return;
        const p = this.getAtomCoords(atom, _State.p1);
        Vec3.transformMat4(p, p, xform);
        atom.row.Cartn_x = p[0];
        atom.row.Cartn_y = p[1];
        atom.row.Cartn_z = p[2];
    }

    transformCoords(xform: Mat4, atoms?: (number | JSONCifLigandGraphAtom)[]) {
        for (const a of atoms ?? this.atoms) {
            this.transformAtomCoords(xform, a);
        }
    }

    traverse<S>(
        atomOrId: number | JSONCifLigandGraphAtom,
        how: 'dfs' | 'bfs',
        state: S,
        visitAtom: (atom: JSONCifLigandGraphAtom, state: S, pred: JSONCifLigandGraphBond | undefined, graph: JSONCifLigandGraph) => void,
    ): S {
        const start = this.getAtom(atomOrId);
        if (!start) return state;

        const visited = new Set<string>();
        const pred = new Map<string, JSONCifLigandGraphBond>();
        const queue: string[] = [start.key];

        while (queue.length) {
            const key = how === 'bfs' ? queue.shift()! : queue.pop()!;
            if (visited.has(key)) continue;

            const a = this.atomsByKey.get(key)!;
            visited.add(a.key);
            visitAtom(a, state, pred.get(key), this);

            const bs = this.bondByKey.get(a.key);
            if (!bs?.length) continue;
            for (const b of bs) {
                if (visited.has(b.atom_2.key)) continue;
                queue.push(b.atom_2.key);
                pred.set(b.atom_2.key, b);
            }
        }
        return state;
    }

    getData(): JSONCifLigandGraphData {
        const atomIdRemapping = new Map<number, number>();
        const addedAtomIds: number[] = [];

        const sortedAtoms = this.atoms.map((a, i) => [a, i] as const);
        sortedAtoms.sort((a, b) => {
            const x = a[0].row.type_symbol;
            const y = b[0].row.type_symbol;
            if (x === 'H' && y !== 'H') return 1;
            if (x !== 'H' && y === 'H') return -1;
            return a[1] - b[1];
        });

        const atoms: Atom[] = [];

        for (let i = 0; i < sortedAtoms.length; ++i) {
            const a = sortedAtoms[i][0];
            const id = i + 1;

            if (a.row.id === undefined) {
                addedAtomIds.push(id);
            } else {
                atomIdRemapping.set(a.row.id!, id);
            }

            a.final_id = id;
            atoms.push({ ...a.row, id });
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
        for (const [a] of sortedAtoms) {
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
        bonds.sort((a, b) => {
            if (a.atom_id_1 !== b.atom_id_1) return a.atom_id_1 - b.atom_id_1;
            if (a.atom_id_2 !== b.atom_id_2) return a.atom_id_2 - b.atom_id_2;
            return 0;
        });

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

    constructor(private data: JSONCifDataBlock) {
        for (const row of data.categories['atom_site'].rows) {
            const atom: JSONCifLigandGraphAtom = {
                key: UUID.create22(),
                final_id: row.final_id,
                row: { ...row },
            };
            this.atoms.push(atom);
            this.atomsByKey.set(atom.key, atom);
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