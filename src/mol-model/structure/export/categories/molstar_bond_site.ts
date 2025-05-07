/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../../mol-data/db';
import { CifWriter } from '../../../../mol-io/writer/cif';
import { BondType } from '../../model/types';
import { StructureElement, StructureProperties, Unit } from '../../structure';
import { CifExportCategoryInfo, CifExportContext } from '../mmcif';
import { sortedCantorPairing } from '../../../../mol-data/util';

export function molstar_bond_site(ctx: CifExportContext): CifExportCategoryInfo | undefined {
    const entries = getEntries(ctx);
    if (entries.length === 0) return;
    return [Category, entries, { ignoreFilter: true }];
}

export type MolstarBondSiteValueOrder = 'sing' | 'doub' | 'trip' | 'quad' | 'arom';
export type MolstarBondSiteTypeId = 'covale' | 'disulf' | 'metalc' | 'hydrog';

export const MolstarBondSiteSchema = {
    molstar_bond_site: {
        atom_id_1: Column.Schema.int,
        atom_id_2: Column.Schema.int,
        value_order: Column.Schema.Aliased<MolstarBondSiteValueOrder>(Column.Schema.lstr),
        type_id: Column.Schema.Aliased<MolstarBondSiteTypeId>(Column.Schema.lstr),
    }
};

export type MolstarBondSiteSchema = typeof MolstarBondSiteSchema;

interface Entry {
    atom_id_1: number,
    atom_id_2: number,
    value_order?: MolstarBondSiteValueOrder,
    type_id?: MolstarBondSiteTypeId,
}

const Fields = CifWriter.fields<number, Entry[], keyof MolstarBondSiteSchema['molstar_bond_site']>()
    .int('atom_id_1', (i, xs) => xs[i].atom_id_1)
    .int('atom_id_2', (i, xs) => xs[i].atom_id_2)
    .str('value_order', (i, xs) => xs[i].value_order ?? '', {
        valueKind: (i, xs) => xs[i].value_order === undefined ? Column.ValueKinds.NotPresent : Column.ValueKinds.Present,
    })
    .str('type_id', (i, xs) => xs[i].type_id ?? '', {
        valueKind: (i, xs) => xs[i].type_id === undefined ? Column.ValueKinds.NotPresent : Column.ValueKinds.Present,
    })
    .getFields();


const Category: CifWriter.Category<Entry[]> = {
    name: 'molstar_bond_site',
    instance(entries: Entry[]) {
        return { fields: Fields, source: [{ data: entries, rowCount: entries.length }] };
    }
};

function assignValueOrder(order: number, flags: BondType.Flag, out: [MolstarBondSiteValueOrder | undefined, MolstarBondSiteTypeId | undefined]) {
    out[0] = undefined;
    out[1] = undefined;

    if (order === 1) out[0] = 'sing';
    else if (order === 2) out[0] = 'doub';
    else if (order === 3) out[0] = 'trip';
    else if (order === 4) out[0] = 'quad';
    if (BondType.is(flags, BondType.Flag.Aromatic)) out[0] = 'arom';

    if (BondType.is(flags, BondType.Flag.Disulfide)) out[1] = 'disulf';
    else if (BondType.is(flags, BondType.Flag.Covalent)) out[1] = 'covale';
    else if (BondType.is(flags, BondType.Flag.MetallicCoordination)) out[1] = 'metalc';
    else if (BondType.is(flags, BondType.Flag.HydrogenBond)) out[1] = 'hydrog';
}

function getEntries(ctx: CifExportContext) {
    const entries: Entry[] = [];
    const added = new Set<number>();

    const loc = StructureElement.Location.create();
    const { id: atom_id } = StructureProperties.atom;

    const info: [MolstarBondSiteValueOrder | undefined, MolstarBondSiteTypeId | undefined] = [undefined, undefined];

    const add = (a: number, b: number) => {
        const key = sortedCantorPairing(a, b);
        if (added.has(key)) return;
        added.add(key);
        if (a > b) {
            entries.push({ atom_id_1: b, atom_id_2: a, value_order: info[0], type_id: info[1] });
        } else {
            entries.push({ atom_id_1: a, atom_id_2: b, value_order: info[0], type_id: info[1] });
        }
    };

    for (const s of ctx.structures) {
        loc.structure = s;

        for (const u of s.units) {
            if (!Unit.isAtomic(u)) continue;

            const { elements } = u;
            const { a, b, edgeProps } = u.bonds;
            loc.unit = u;

            for (let i = 0; i < a.length; i++) {
                loc.element = elements[a[i]];
                const atom_id_1 = atom_id(loc);
                loc.element = elements[b[i]];
                const atom_id_2 = atom_id(loc);

                assignValueOrder(edgeProps.order[i], edgeProps.flags[i], info);
                add(atom_id_1, atom_id_2);
            }
        }

        s.interUnitBonds.edges.forEach((e) => {
            let u = s.unitMap.get(e.unitA);
            loc.unit = u;
            loc.element = u.elements[e.indexA];
            const atom_id_1 = atom_id(loc);

            u = s.unitMap.get(e.unitB);
            loc.unit = u;
            loc.element = u.elements[e.indexB];
            const atom_id_2 = atom_id(loc);

            assignValueOrder(e.props.order, e.props.flag, info);
            add(atom_id_1, atom_id_2);
        });
    }

    entries.sort((a, b) => {
        if (a.atom_id_1 !== b.atom_id_1) return a.atom_id_1 - b.atom_id_1;
        if (a.atom_id_2 !== b.atom_id_2) return a.atom_id_2 - b.atom_id_2;
        return 0;
    });

    return entries;
}