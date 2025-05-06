/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from '@/mol-util';
import { JSONCifDataBlock } from './model';
import { arrayMapAdd } from '@/mol-util/map';
import { mmCIF_Schema } from '@/mol-io/reader/cif/schema/mmcif';
import { Table } from '@/mol-data/db';
import { MolstarBondSiteTypeId, MolstarBondSiteValueOrder } from '@/mol-model/structure/export/categories/molstar_bond_site';

type Atom = Partial<Table.Row<mmCIF_Schema['atom_site']>>

export interface LigandGraphBondProps {
    value_order: MolstarBondSiteValueOrder | undefined;
    type_id: MolstarBondSiteTypeId | undefined;
}

export interface AtomSiteRow {
    key: string;
    final_id: number | undefined;
    row: Atom;
}

export interface LigandGraphBond extends LigandGraphBondProps {
    atom_1: AtomSiteRow;
    atom_2: AtomSiteRow;
}

export interface JSONCifLigandGraphData {
    block: JSONCifDataBlock;
    atomIdRemapping: Map<number, number>;
    addedAtomIds: number[];
    removedAtomIds: number[];
}

export class JSONCifLigandGraph {
    readonly atoms: AtomSiteRow[] = [];
    readonly atomsById: Map<number, AtomSiteRow> = new Map();
    readonly bondByKey: Map<string, LigandGraphBond[]> = new Map();
    readonly removedAtomIds: Set<number> = new Set();

    getAtom(atomOrId: number | AtomSiteRow) {
        return typeof atomOrId === 'number' ? this.atomsById.get(atomOrId) : atomOrId;
    }

    modifyAtom(id: number, data: Atom) {
        const atom = this.atomsById.get(id);
        if (!atom) return;
        atom.row = { ...atom.row, ...data };
    }

    addAtom(data: Omit<Atom, 'id'>) {
        const atom: AtomSiteRow = {
            key: UUID.create22(),
            final_id: undefined,
            row: { ...data, id: undefined },
        };
        this.atoms.push(atom);
        return atom;
    }

    removeAtom(atomOrId: number | AtomSiteRow) {
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

    addOrUpdateBond(atom1: number | AtomSiteRow, atom2: number | AtomSiteRow, props: LigandGraphBondProps) {
        const a1 = this.getAtom(atom1);
        const a2 = this.getAtom(atom2);
        if (!a1 || !a2) return;

        this.removeBond(atom1, atom2);
        arrayMapAdd(this.bondByKey, a1.key, { atom_1: a1, atom_2: a2, ...props });
        arrayMapAdd(this.bondByKey, a2.key, { atom_1: a2, atom_2: a1, ...props });
    }

    removeBond(atom1: number | AtomSiteRow, atom2: number | AtomSiteRow) {
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

        for (let i = 0; i < this.atoms.length; ++i) {
            const id = i + 1;

            if (this.atoms[i].row.id === undefined) {
                addedAtomIds.push(id);
            } else {
                atomIdRemapping.set(this.atoms[i].row.id!, id);
            }

            this.atoms[i].final_id = id;
            this.atoms[i].row.id = id;
        }

        const block: JSONCifDataBlock = {
            ...this.data,
            categories: {
                ...this.data.categories,
                atom_site: {
                    ...this.data.categories['atom_site'],
                    rows: this.atoms.map(a => a.row),
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
                    value_order: bb.value_order,
                    type_id: bb.type_id,
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

    constructor(private data: JSONCifDataBlock) {
        for (const row of data.categories['atom_site'].rows) {
            const atom: AtomSiteRow = {
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
                value_order: row.value_order,
                type_id: row.type_id,
            });
            arrayMapAdd(this.bondByKey, atom_2.key, {
                atom_1: atom_2,
                atom_2: atom_1,
                value_order: row.value_order,
                type_id: row.type_id,
            });
        }
    }
}