/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MVSAnnotationRow } from './schemas.js';
import { isDefined } from './utils.js';

/** Longest label made of the descriptions of individual rows */
const MaxListedLength = 60;
/** Most rows whose descriptions are listed one by one */
const MaxListedRows = 3;
/** Most chains listed by name when each row is just a chain */
const MaxListedChains = 6;

/** A short description of the criteria in one row, e.g. `Chain A, residues 10-50`. A row without criteria selects everything (`All`). */
export function describeRow(row: MVSAnnotationRow): string {
  const parts: string[] = [];
  if (isDefined(row.label_entity_id)) parts.push(`Entity ${row.label_entity_id}`);
  if (isDefined(row.label_asym_id)) parts.push(`Chain ${row.label_asym_id}`);
  if (isDefined(row.auth_asym_id)) parts.push(`Chain ${row.auth_asym_id} (auth)`);

  if (isDefined(row.label_seq_id)) parts.push(`residue ${row.label_seq_id}`);
  if (isDefined(row.auth_seq_id)) parts.push(`auth residue ${row.auth_seq_id}${row.pdbx_PDB_ins_code ?? ''}`);
  else if (isDefined(row.pdbx_PDB_ins_code)) parts.push(`insertion code ${row.pdbx_PDB_ins_code}`);
  if (isDefined(row.residue_index)) parts.push(`residue index ${row.residue_index}`);
  const labelRange = describeRange(row.beg_label_seq_id, row.end_label_seq_id);
  if (labelRange) parts.push(`residues ${labelRange}`);
  const authRange = describeRange(row.beg_auth_seq_id, row.end_auth_seq_id);
  if (authRange) parts.push(`auth residues ${authRange}`);

  if (isDefined(row.label_comp_id)) parts.push(row.label_comp_id);
  if (isDefined(row.auth_comp_id) && row.auth_comp_id !== row.label_comp_id) parts.push(`${row.auth_comp_id} (auth)`);

  if (isDefined(row.label_atom_id)) parts.push(`atom ${row.label_atom_id}`);
  if (isDefined(row.auth_atom_id)) parts.push(`auth atom ${row.auth_atom_id}`);
  if (isDefined(row.type_symbol)) parts.push(`${row.type_symbol.toUpperCase()} atoms`);
  if (isDefined(row.atom_id)) parts.push(`atom id ${row.atom_id}`);
  if (isDefined(row.atom_index)) parts.push(`atom index ${row.atom_index}`);

  if (isDefined(row.instance_id)) parts.push(`instance ${row.instance_id}`);

  if (parts.length === 0) return 'All';
  const text = parts.join(', ');
  return text[0].toUpperCase() + text.slice(1);
}

function describeRange(beg: number | undefined | null, end: number | undefined | null): string | undefined {
  if (isDefined(beg) && isDefined(end)) return beg === end ? `${beg}` : `${beg}-${end}`;
  if (isDefined(beg)) return `from ${beg}`;
  if (isDefined(end)) return `up to ${end}`;
  return undefined;
}

/** `true` if the only criterion of the row is the chain (`label_asym_id`) */
function isChainOnly(row: MVSAnnotationRow): boolean {
  if (!isDefined(row.label_asym_id)) return false;
  return Object.keys(row).every(
    (key) => key === 'label_asym_id' || key === 'group_id' || !isDefined((row as any)[key]),
  );
}

/** A short name for a selection made of `rows` (an element is selected if it satisfies any of them), for labels.
 * A few rows are described one by one (`Chain A`, `Chain A, residues 10-50`, `Chains A, B, C`),
 * more of them are only counted (`120 selections`, or `Chain A: 120 selections` if they are all in one chain). */
export function describeRows(rows: readonly MVSAnnotationRow[]): string {
  if (rows.length === 0) return 'Empty Selection';
  if (rows.length === 1) return describeRow(rows[0]);

  if (rows.every(isChainOnly)) {
    if (rows.length <= MaxListedChains) return `Chains ${rows.map((r) => r.label_asym_id).join(', ')}`;
    return `${rows.length} chains`;
  }

  if (rows.length <= MaxListedRows) {
    const listed = rows.map(describeRow).join(' + ');
    if (listed.length <= MaxListedLength) return listed;
  }

  const chain = rows[0].label_asym_id;
  if (isDefined(chain) && rows.every((r) => r.label_asym_id === chain))
    return `Chain ${chain}: ${rows.length} selections`;
  return `${rows.length} selections`;
}
