/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { AtomSite, BasicData } from './schema.js';
import { Column, Table } from '@molstar/core/data/db';

export function getModelGroupName(model_id: number, data: BasicData) {
  return getModelGroupNames(data)(model_id);
}

/**
 * The model group name of each model id: the name of the first `ihm_model_group` row with the group id of the model's
 * first `ihm_model_group_link` row, '' without one. Both tables are read once; `getModelGroupName` scanned them for
 * every model (perf-audit 04/P12).
 */
export function getModelGroupNames(data: BasicData): (model_id: number) => string {
  const { ihm_model_group, ihm_model_group_link } = data;
  const groupOfModel = new Map<number, number>();
  for (let i = 0, il = ihm_model_group_link._rowCount; i < il; ++i) {
    const modelId = ihm_model_group_link.model_id.value(i);
    if (!groupOfModel.has(modelId)) groupOfModel.set(modelId, ihm_model_group_link.group_id.value(i));
  }
  const nameOfGroup = new Map<number, string>();
  for (let i = 0, il = ihm_model_group._rowCount; i < il; ++i) {
    const id = ihm_model_group.id.value(i);
    if (!nameOfGroup.has(id)) nameOfGroup.set(id, ihm_model_group.name.value(i));
  }
  // NaN never matched the scans' `===`; a Map would match it
  return (model_id) => {
    const group = Number.isNaN(model_id) ? undefined : groupOfModel.get(model_id);
    return group === undefined || Number.isNaN(group) ? '' : (nameOfGroup.get(group) ?? '');
  };
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
