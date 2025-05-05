/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { UUID } from '@/mol-util';
import { JSONCifDataBlock } from './model';
import { arrayMapAdd } from '@/mol-util/map';

export interface AtomSiteRow {
    key: string;
    final_id: number | undefined;
    row: Record<string, any>;
}

export interface Bond {
    atom_1: AtomSiteRow;
    atom_2: AtomSiteRow;
    value_order: string | undefined;
    type_id: string | undefined;
}

export class JSONCifLigandGraph {
    readonly atoms: AtomSiteRow[] = [];
    readonly atomsById: Map<number, AtomSiteRow> = new Map();
    readonly bondByKey: Map<string, Bond[]> = new Map();

    modifyAtom(id: number, data: Record<string, any>) {
        const atom = this.atomsById.get(id);
        if (!atom) return;
        atom.row = { ...atom.row, ...data };
    }

    removeAtom(id: number) {
        const atom = this.atomsById.get(id);
        if (!atom) return;
        this.atomsById.delete(id);
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

    getData(): JSONCifDataBlock {
        for (let i = 0; i < this.atoms.length; ++i) {
            this.atoms[i].final_id = i + 1;
            this.atoms[i].row['id'] = i + 1;
        }

        const ret: JSONCifDataBlock = {
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

        if (ret.categories.molstar_bond_site) {
            ret.categories['molstar_bond_site'] = {
                ...ret.categories['molstar_bond_site'],
                rows: bonds
            };
        }

        return ret;
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