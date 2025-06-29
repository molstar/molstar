/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getJSONCifCategory, JSONCifDataBlock } from '../../extensions/json-cif/model';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { MolstarBondSiteSchema, MolstarBondSiteTypeId, MolstarBondSiteValueOrder } from '../../mol-model/structure/export/categories/molstar_bond_site';

function padLeft(v: any, n = 3) {
    let s = `${v}`;
    while (s.length < n) s = ' ' + s;
    return s;
}

function padRight(v: any, n = 3) {
    let s = `${v}`;
    while (s.length < n) s = s + ' ';
    return s;
}

function mapMolChage(v: number) {
    switch (v) {
        case 3: return 1;
        case 2: return 2;
        case 1: return 3;
        case -1: return 5;
        case -2: return 6;
        case -3: return 7;
        default: return 0;
    }
}

function mapMolBondOrder(order: MolstarBondSiteValueOrder, type: MolstarBondSiteTypeId) {
    if (type !== 'covale') return 8;

    switch (order) {
        case 'sing': return 1;
        case 'doub': return 2;
        case 'trip': return 3;
        case 'arom': return 4;
        default: return 8;
    }
}

export function jsonCifToMolfile(data: JSONCifDataBlock, options?: { name?: string, comment?: string }) {
    // The method works in the sense that Mol* can re-open the file.
    // For production use, this will likely need more testing and tweaks (e.g., support for M CHG property).

    if (data.categories.atom_site === undefined || data.categories.molstar_bond_site === undefined) {
        throw new Error('The data block must contain atom_site and molstar_bond_site categories.');
    }

    const { atom_site: _atoms, molstar_bond_site: _bonds } = data.categories;

    const atoms = getJSONCifCategory<mmCIF_Schema['atom_site']>(data, 'atom_site')!;
    const bonds = getJSONCifCategory<MolstarBondSiteSchema['molstar_bond_site']>(data, 'molstar_bond_site')!;

    const lines = [
        `${options?.name ?? 'mol'}`,
        '  Molstar           3D',
        options?.comment ?? '',
        `${padLeft(atoms.rows.length)}${padLeft(bonds.rows.length)}  0  0  0  0  0  0  0  0 V2000`,
    ];

    const atomIdToIndex = new Map<number, number>();
    for (let i = 0; i < atoms.rows.length; ++i) {
        const a = atoms.rows[i];
        const { id, Cartn_x, Cartn_y, Cartn_z, type_symbol, pdbx_formal_charge } = a;
        atomIdToIndex.set(id, i + 1);
        const fields = [
            padLeft(Cartn_x.toFixed(4), 10),
            padLeft(Cartn_y.toFixed(4), 10),
            padLeft(Cartn_z.toFixed(4), 10),
            ' ',
            padRight(type_symbol, 2),
            '  0',
            padLeft(mapMolChage(pdbx_formal_charge), 3),
            '  0  0  0  0  0  0  0  0  0  0',
        ];
        lines.push(fields.join(''));
    }

    for (const b of bonds.rows) {
        const { atom_id_1, atom_id_2, value_order, type_id } = b;
        const fields = [
            padLeft(atomIdToIndex.get(atom_id_1)!, 3),
            padLeft(atomIdToIndex.get(atom_id_2)!, 3),
            padLeft(mapMolBondOrder(value_order, type_id), 3),
            '  0  0  0  0',
        ];
        lines.push(fields.join(''));
    }

    lines.push('M  END');
    return lines.join('\n');
}