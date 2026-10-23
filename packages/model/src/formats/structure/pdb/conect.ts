/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortedCantorPairing } from '@molstar/core/data/util';
import { CifCategory, CifField } from '@molstar/io/reader/cif';
import type { mmCIF_Schema } from '@molstar/io/reader/cif/schema/mmcif';
import type { Tokens } from '@molstar/io/reader/common/text/tokenizer';
import { isMetal } from '@molstar/model/model/structure/model/properties/atomic/types';
import type { ElementSymbol } from '@molstar/model/model/structure/model/types';

export function parseConect(
  lines: Tokens,
  lineStart: number,
  lineEnd: number,
  sites: { [K in keyof mmCIF_Schema['atom_site']]?: CifField },
): CifCategory {
  const idMap: { [k: string]: number } = {};
  for (let i = 0, il = sites.id!.rowCount; i < il; ++i) {
    idMap[sites.id!.str(i)] = i;
  }

  const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);

  const id: string[] = [];
  const conn_type_id: string[] = [];
  const bondOrder: number[] = [];

  const ptnr1_label_asym_id: string[] = [];
  const ptnr1_label_seq_id: number[] = [];
  const ptnr1_auth_seq_id: number[] = [];
  const ptnr1_label_atom_id: string[] = [];
  const ptnr1_label_alt_id: string[] = [];
  const ptnr1_PDB_ins_code: string[] = [];

  const ptnr2_label_asym_id: string[] = [];
  const ptnr2_label_seq_id: number[] = [];
  const ptnr2_auth_seq_id: number[] = [];
  const ptnr2_label_atom_id: string[] = [];
  const ptnr2_label_alt_id: string[] = [];
  const ptnr2_PDB_ins_code: string[] = [];

  const pos = [11, 16, 21, 26];

  let k = 1;
  // One bond per atom pair, from the records of either atom: files list a bond on both atoms' records or on one, either
  // of them (perf-audit-2 04/B11: one listed on the higher serial's record only was dropped). A partner given several
  // times makes a double or triple bond (CONECT 1529 1528 1528 is a double bond): the order is the larger count of the
  // two atoms' records.
  const bondIndex = new Map<number, number>();
  /** per bond, the times the first and the second partner's records list it */
  const sideCounts: number[] = [];

  for (let i = lineStart; i < lineEnd; i++) {
    const line = getLine(i);
    const idxRecord = idMap[parseInt(line.substr(6, 5))];

    if (idxRecord === undefined) continue;

    for (let j = 0; j < 4; ++j) {
      const idB = parseInt(line.substr(pos[j], 5));
      if (Number.isNaN(idB)) continue;

      const idxPartner = idMap[idB];
      if (idxPartner === undefined) continue;

      const side = idxRecord <= idxPartner ? 0 : 1;
      const key = sortedCantorPairing(idxRecord, idxPartner);
      const b = bondIndex.get(key);
      if (b !== undefined) {
        sideCounts[2 * b + side] += 1;
        continue;
      }
      bondIndex.set(key, k - 1);
      sideCounts.push(1 - side, side);

      const idxA = Math.min(idxRecord, idxPartner);
      const idxB = Math.max(idxRecord, idxPartner);
      const atomIdA = sites.type_symbol!.str(idxA) as ElementSymbol;
      const atomIdB = sites.type_symbol!.str(idxB) as ElementSymbol;

      if (isMetal(atomIdA) || isMetal(atomIdB)) {
        id.push(`metalc${k}`);
        conn_type_id.push('metalc');
      } else {
        id.push(`covale${k}`);
        conn_type_id.push('covale');
      }

      ptnr1_label_asym_id.push(sites.label_asym_id!.str(idxA));
      ptnr1_label_seq_id.push(sites.label_seq_id!.int(idxA));
      ptnr1_auth_seq_id.push(sites.auth_seq_id!.int(idxA));
      ptnr1_label_atom_id.push(sites.label_atom_id!.str(idxA));
      ptnr1_label_alt_id.push(sites.label_alt_id!.str(idxA));
      ptnr1_PDB_ins_code.push(sites.pdbx_PDB_ins_code!.str(idxA));

      ptnr2_label_asym_id.push(sites.label_asym_id!.str(idxB));
      ptnr2_label_seq_id.push(sites.label_seq_id!.int(idxB));
      ptnr2_auth_seq_id.push(sites.auth_seq_id!.int(idxB));
      ptnr2_label_atom_id.push(sites.label_atom_id!.str(idxB));
      ptnr2_label_alt_id.push(sites.label_alt_id!.str(idxB));
      ptnr2_PDB_ins_code.push(sites.pdbx_PDB_ins_code!.str(idxB));

      k += 1;
    }
  }

  const struct_conn: Partial<CifCategory.Fields<mmCIF_Schema['struct_conn']>> = {
    id: CifField.ofStrings(id),
    conn_type_id: CifField.ofStrings(conn_type_id),

    ptnr1_label_asym_id: CifField.ofStrings(ptnr1_label_asym_id),
    ptnr1_label_seq_id: CifField.ofNumbers(ptnr1_label_seq_id),
    ptnr1_auth_seq_id: CifField.ofNumbers(ptnr1_auth_seq_id),
    ptnr1_label_atom_id: CifField.ofStrings(ptnr1_label_atom_id),
    pdbx_ptnr1_label_alt_id: CifField.ofStrings(ptnr1_label_alt_id),
    pdbx_ptnr1_PDB_ins_code: CifField.ofStrings(ptnr1_PDB_ins_code),

    ptnr2_label_asym_id: CifField.ofStrings(ptnr2_label_asym_id),
    ptnr2_label_seq_id: CifField.ofNumbers(ptnr2_label_seq_id),
    ptnr2_auth_seq_id: CifField.ofNumbers(ptnr2_auth_seq_id),
    ptnr2_label_atom_id: CifField.ofStrings(ptnr2_label_atom_id),
    pdbx_ptnr2_label_alt_id: CifField.ofStrings(ptnr2_label_alt_id),
    pdbx_ptnr2_PDB_ins_code: CifField.ofStrings(ptnr2_PDB_ins_code),
  };

  let hasMultipleBonds = false;
  for (let b = 0, il = k - 1; b < il; ++b) {
    const order = Math.max(sideCounts[2 * b], sideCounts[2 * b + 1]);
    bondOrder.push(order);
    if (order > 1) hasMultipleBonds = true;
  }

  if (hasMultipleBonds) {
    // the mmCIF value_order enumeration (perf-audit-2 04/B5: 'trpl' read back as order 1)
    const valueOrder = ['sing', 'doub', 'trip', 'quad'];
    struct_conn.pdbx_value_order = CifField.ofStrings(bondOrder.map((bo) => valueOrder[bo - 1] || 'sing'));
  }

  return CifCategory.ofFields('struct_conn', struct_conn);
}
