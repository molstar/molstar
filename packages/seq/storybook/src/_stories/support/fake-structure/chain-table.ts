/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ComponentExpressionT } from '@molstar/mvs-builder/tree/mvs/param-types';

// The fake structure viewer's own structure: a recorded chain table (a chain's ids, entity and operators, and per
// residue its label and author numbers, component and whether it was modeled), read from the mmCIF by the examples'
// recorder and held apart from anything of Seq*: the mapping crosses two independently held representations, as it
// does with a real viewer. Nothing here imports Seq*. A chain instance is a row; a cell is a residue on a row; the
// fake's words for places are spans of cells and MolViewSpec's selectors.

/** A residue of a chain instance. */
export interface ChainResidue {
  /** `label_seq_id` (1-based), null for a non-polymer's residue. */
  readonly label: number | null;
  readonly auth: string;
  /** `pdbx_PDB_ins_code`, empty for none. */
  readonly ins: string;
  readonly comp: string;
  readonly observed: boolean;
}

/** A chain as the mmCIF names it: its ids, entity, operators (a row each) and residues. */
export interface ChainEntry {
  readonly label_asym_id: string;
  readonly auth_asym_id: string;
  readonly entity: string;
  readonly operators: readonly string[];
  readonly residues: readonly ChainResidue[];
}

export interface ChainTable {
  readonly entry: string;
  readonly chains: readonly ChainEntry[];
}

/** A chain instance, drawn as a row: a chain under one operator. */
export interface Row {
  readonly index: number;
  readonly asym: string;
  readonly auth: string;
  readonly entity: string;
  readonly operator: string;
  /** Residues numbered by `label_seq_id` (a polymer); else one chip (a ligand). */
  readonly polymer: boolean;
  readonly residues: readonly ChainResidue[];
}

/** Cells `[start, end)` of a row. */
export interface Span {
  readonly row: number;
  readonly start: number;
  readonly end: number;
}

/** A cell: residue `index` of row `row`. */
export interface Cell {
  readonly row: number;
  readonly index: number;
}

/** The table's rows: its chains in order, each once for each of its operators (copies as rows). */
export function rowsOf(table: ChainTable): Row[] {
  const rows: Row[] = [];
  for (const c of table.chains)
    for (const operator of c.operators)
      rows.push({
        index: rows.length,
        asym: c.label_asym_id,
        auth: c.auth_asym_id,
        entity: c.entity,
        operator,
        polymer: c.residues.length > 0 && c.residues[0].label !== null,
        residues: c.residues,
      });
  return rows;
}

/** Spans sorted and merged where they touch or overlap. */
export function normalizeSpans(spans: readonly Span[]): Span[] {
  const sorted = spans.filter((s) => s.end > s.start).sort((a, b) => a.row - b.row || a.start - b.start);
  const out: Span[] = [];
  for (const s of sorted) {
    const last = out[out.length - 1];
    if (last && last.row === s.row && s.start <= last.end) {
      if (s.end > last.end) out[out.length - 1] = { row: s.row, start: last.start, end: s.end };
    } else out.push(s);
  }
  return out;
}

/** The observed cells of spans only (an unobserved residue has no atoms: it is never lit, selected or sent). */
export function observedSpans(rows: readonly Row[], spans: readonly Span[]): Span[] {
  const out: Span[] = [];
  for (const s of spans) {
    const res = rows[s.row].residues;
    let from = -1;
    for (let i = s.start; i <= s.end; i++) {
      const ok = i < s.end && res[i].observed;
      if (ok && from < 0) from = i;
      else if (!ok && from >= 0) {
        out.push({ row: s.row, start: from, end: i });
        from = -1;
      }
    }
  }
  return normalizeSpans(out);
}

export const spanOfCell = (c: Cell): Span => ({ row: c.row, start: c.index, end: c.index + 1 });

/** Whether a span set holds the cell. */
export function hasCell(spans: readonly Span[], c: Cell): boolean {
  return spans.some((s) => s.row === c.row && c.index >= s.start && c.index < s.end);
}

/** `a` with `b` toggled: the cells in exactly one of them. */
export function toggleSpans(rows: readonly Row[], a: readonly Span[], b: readonly Span[]): Span[] {
  const out: Span[] = [];
  for (const row of rows) {
    const n = row.residues.length;
    const bits = new Uint8Array(n);
    for (const s of a) if (s.row === row.index) bits.fill(1, s.start, s.end);
    for (const s of b) if (s.row === row.index) for (let i = s.start; i < s.end; i++) bits[i] ^= 1;
    let from = -1;
    for (let i = 0; i <= n; i++) {
      if (i < n && bits[i]) {
        if (from < 0) from = i;
      } else if (from >= 0) {
        out.push({ row: row.index, start: from, end: i });
        from = -1;
      }
    }
  }
  return out;
}

const inRange = (v: number, lo: number | undefined, hi: number | undefined) =>
  (lo === undefined || v >= lo) && (hi === undefined || v <= hi);

/**
 * The observed cells a MolViewSpec selector names: its chain (`label_asym_id`, `auth_asym_id`), `label_entity_id` and
 * `instance_id` pick the rows; a residue by `label_seq_id` (or `beg_` and `end_`), `auth_seq_id` or `label_comp_id`
 * picks the cells; atoms name their residue's cell. No residue field: the whole row.
 */
export function spansOfSelector(rows: readonly Row[], sel: ComponentExpressionT): Span[] {
  const out: Span[] = [];
  for (const row of rows) {
    if (sel.label_asym_id !== undefined && sel.label_asym_id !== row.asym) continue;
    if (sel.auth_asym_id !== undefined && sel.auth_asym_id !== row.auth) continue;
    if (sel.label_entity_id !== undefined && sel.label_entity_id !== row.entity) continue;
    if (sel.instance_id !== undefined && sel.instance_id !== row.operator) continue;
    const lo = sel.label_seq_id ?? sel.beg_label_seq_id;
    const hi = sel.label_seq_id ?? sel.end_label_seq_id;
    const alo = sel.auth_seq_id ?? sel.beg_auth_seq_id;
    const ahi = sel.auth_seq_id ?? sel.end_auth_seq_id;
    row.residues.forEach((r, i) => {
      if (lo !== undefined || hi !== undefined) {
        if (r.label === null || !inRange(r.label, lo, hi)) return;
      }
      if (alo !== undefined || ahi !== undefined) {
        if (!inRange(Number(r.auth), alo, ahi)) return;
        if (sel.pdbx_PDB_ins_code !== undefined && sel.pdbx_PDB_ins_code !== r.ins) return;
      }
      if (sel.label_comp_id !== undefined && sel.label_comp_id !== r.comp) return;
      out.push({ row: row.index, start: i, end: i + 1 });
    });
  }
  return observedSpans(rows, normalizeSpans(out));
}

/** The selectors of spans: a polymer's run by its label range, a ligand's chip by its chain. */
export function selectorsOfSpans(rows: readonly Row[], spans: readonly Span[]): ComponentExpressionT[] {
  const out: ComponentExpressionT[] = [];
  for (const s of normalizeSpans(spans)) {
    const row = rows[s.row];
    const base = { label_asym_id: row.asym, instance_id: row.operator };
    if (!row.polymer) out.push(base);
    else {
      const beg = row.residues[s.start].label!;
      const end = row.residues[s.end - 1].label!;
      out.push(
        beg === end ? { ...base, label_seq_id: beg } : { ...base, beg_label_seq_id: beg, end_label_seq_id: end },
      );
    }
  }
  return out;
}
