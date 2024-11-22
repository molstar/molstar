/**
 * Copyright (c) 2017-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomSite, BasicData } from './schema';
import { Column, Table } from '../../../mol-data/db';

export function getModelGroupName(model_id: number, data: BasicData) {
    const { ihm_model_group, ihm_model_group_link } = data;

    const link = Table.pickRow(ihm_model_group_link, i => ihm_model_group_link.model_id.value(i) === model_id);
    if (link) {
        const group = Table.pickRow(ihm_model_group, i => ihm_model_group.id.value(i) === link.group_id);
        if (group) return group.name;
    }
    return '';
}

//

function hasPresentValues(column: Column<any>) {
    for (let i = 0, il = column.rowCount; i < il; i++) {
        if (column.valueKind(i) === Column.ValueKinds.Present) return true;
    }
    return false;
}

function substUndefinedColumn<T extends Table<any>>(table: T, a: keyof T, b: keyof T) {
    if (!table[a].isDefined || !hasPresentValues(table[a])) table[a] = table[b];
    if (!table[b].isDefined || !hasPresentValues(table[b])) table[b] = table[a];
}

/** Fix possibly missing auth_/label_ columns */
export function getNormalizedAtomSite(atom_site: AtomSite) {
    const normalized = Table.ofColumns(atom_site._schema, atom_site);
    substUndefinedColumn(normalized, 'label_atom_id', 'auth_atom_id');
    substUndefinedColumn(normalized, 'label_comp_id', 'auth_comp_id');
    substUndefinedColumn(normalized, 'label_seq_id', 'auth_seq_id');
    substUndefinedColumn(normalized, 'label_asym_id', 'auth_asym_id');
    substUndefinedColumn(normalized, 'label_entity_id', 'label_asym_id');
    return normalized;
}