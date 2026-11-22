/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  CategoryColumn,
  defineIssues,
  type Issue,
  issueOf,
  type Node,
  printAddress,
  type Value,
} from '@molstar/core/repr/tree';
import { DefaultSpecBounds, type TooltipSchema } from '@molstar/seqx-core/lang';

// Tooltip tables as data (_spec/README.md "Resolve", Tooltip tables): a by-range tooltip's `ranges` table, inline
// `rows` or a `source`'s, is an annotation at the tooltip's address, read through the store as every data is (inline
// too: one entry over the whole system, indexed there). Tier 1 holds the tooltip resolved to references, never rows.
// An `address` table and MVS's schemas are read whole by the reader that shows them (`tooltipsAt`): no store entry.

/** A tooltip node with its address, as the resolve pass places data. */
export interface TooltipPlaced {
  readonly node: Node;
  readonly address: string;
}

/** Every tooltip of a spec, `views`' then each view's, in spec order, with its address. */
export function tooltipNodes(spec: Node): TooltipPlaced[] {
  const out: TooltipPlaced[] = [];
  const views: TooltipPlaced[] = [];
  for (const c of spec.children ?? []) {
    if (c.kind !== 'views' || c.id === undefined) continue;
    for (const n of c.children ?? []) {
      if (n.id === undefined) continue;
      if (n.kind === 'tooltip') out.push({ node: n, address: printAddress([c.id, n.id]) });
      else if (n.kind === 'view')
        for (const t of n.children ?? [])
          if (t.kind === 'tooltip' && t.id !== undefined)
            views.push({ node: t, address: printAddress([c.id, n.id, t.id]) });
    }
  }
  return [...out, ...views];
}

/** Whether a tooltip's table is a `ranges` table (by range, the schema absent or `ranges`): the store's. */
export function isRangesTable(n: Node): boolean {
  const p = n.params ?? {};
  return p.on === undefined && p.target !== undefined && (p.schema === undefined || p.schema === 'ranges');
}

/**
 * The tooltips' `ranges` tables as `ranges` data nodes at their tooltips' addresses: what the resolve pass makes into
 * annotations (an assembly's by its sequences' names, as data over an assembly is).
 */
export function tooltipTables(spec: Node): TooltipPlaced[] {
  const out: TooltipPlaced[] = [];
  for (const t of tooltipNodes(spec)) {
    if (!isRangesTable(t.node)) continue;
    const { target, rows, source } = t.node.params ?? {};
    const params = { target, ...(rows !== undefined ? { rows } : {}), ...(source !== undefined ? { source } : {}) };
    out.push({
      node: { kind: 'ranges', ...(t.node.id !== undefined ? { id: t.node.id } : {}), params },
      address: t.address,
    });
  }
  return out;
}

/**
 * The columns each MVS schema reads (MolViewSpec's annotation schemas): a row's other columns are ignored, as MVS
 * ignores them (a `beg_label_seq_id` under `chain` is whole chains).
 */
export const MvsSchemaColumns: Readonly<Record<Exclude<TooltipSchema, 'ranges' | 'address'>, readonly string[]>> = {
  whole_structure: [],
  entity: ['label_entity_id'],
  chain: ['label_entity_id', 'label_asym_id'],
  auth_chain: ['auth_asym_id'],
  residue: ['label_entity_id', 'label_asym_id', 'label_seq_id'],
  auth_residue: ['auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code'],
  residue_range: ['label_entity_id', 'label_asym_id', 'beg_label_seq_id', 'end_label_seq_id'],
  auth_residue_range: ['auth_asym_id', 'beg_auth_seq_id', 'end_auth_seq_id'],
  atom: ['label_entity_id', 'label_asym_id', 'label_seq_id', 'label_atom_id', 'type_symbol', 'atom_id', 'atom_index'],
  auth_atom: [
    'auth_asym_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code',
    'auth_atom_id',
    'type_symbol',
    'atom_id',
    'atom_index',
  ],
  all_atomic: [
    'label_entity_id',
    'label_asym_id',
    'label_seq_id',
    'label_atom_id',
    'auth_asym_id',
    'auth_seq_id',
    'pdbx_PDB_ins_code',
    'auth_atom_id',
    'type_symbol',
    'atom_id',
    'atom_index',
    'beg_label_seq_id',
    'end_label_seq_id',
    'beg_auth_seq_id',
    'end_auth_seq_id',
  ],
};

/**
 * The tooltips' `seq.resolve.tooltip.*` codes (_spec/README.md "Resolve", Tooltip tables): a catalog of their own,
 * loaded with the tooltips (a reader of a spec without one loads none of it).
 */
export const TooltipIssues = defineIssues('seq.resolve.tooltip', {
  'seq.resolve.tooltip.schema': {
    severity: 'error',
    description: 'A tooltip table’s columns are none of those its MolViewSpec schema locates by: it shows nothing',
    hint: 'Say the schema the table is written in, or give one of its columns',
  },
  'seq.resolve.tooltip.source': {
    severity: 'warning',
    description: 'A tooltip’s `address` or MolViewSpec table from a `source` is not read yet: it shows nothing',
    hint: 'Give the table inline (`rows`), or as a `ranges` table from the source',
  },
  'seq.resolve.tooltip.rows': {
    severity: 'warning',
    description:
      'A tooltip table read whole is past the bound (`SpecBounds.tooltipRows`): the rows past it are not read',
    hint: 'Give a `ranges` table (read by region), or fewer rows',
  },
});

/** An inline table's rows (its first column's length). */
export function tableLength(rows: Record<string, Value>): number {
  for (const c of Object.values(rows)) {
    if (c instanceof CategoryColumn) return c.length;
    if (Array.isArray(c) || ArrayBuffer.isView(c)) return (c as ArrayLike<unknown>).length;
  }
  return 0;
}

/**
 * What the check cannot say of the tooltips' tables (their schemas' columns load here): an MVS table off its schema,
 * an `address` or MVS table from a source (not read yet), a table read whole past `tooltipRows` (read to the bound).
 */
export function tooltipIssues(spec: Node): Issue[] {
  const out: Issue[] = [];
  const say = (code: keyof typeof TooltipIssues.codes, message: string, address: string, param: string) =>
    out.push(
      issueOf(TooltipIssues, code, message, {
        address,
        path: ['params', param],
        hint: TooltipIssues.codes[code].hint,
      }),
    );
  for (const t of tooltipNodes(spec)) {
    const p = t.node.params ?? {};
    const schema = (typeof p.schema === 'string' ? p.schema : 'ranges') as TooltipSchema;
    if (p.on !== undefined || p.target === undefined || schema === 'ranges') continue;
    if (p.source !== undefined && p.category === undefined && p.block === undefined) {
      say('seq.resolve.tooltip.source', `The ${schema} table from a source is not read yet`, t.address, 'source');
      continue;
    }
    const rows = p.rows as Record<string, Value> | undefined;
    if (!rows || typeof rows !== 'object') continue;
    const n = tableLength(rows);
    if (n > DefaultSpecBounds.tooltipRows)
      say('seq.resolve.tooltip.rows', `${n} rows, past ${DefaultSpecBounds.tooltipRows}`, t.address, 'rows');
    if (schema === 'address') continue;
    const locates = MvsSchemaColumns[schema];
    if (locates.length && !locates.some((c) => c in rows))
      say(
        'seq.resolve.tooltip.schema',
        `The table's columns (${Object.keys(rows).join(', ')}) do not fit the schema ${schema} (${locates.join(', ')})`,
        t.address,
        'rows',
      );
  }
  return out;
}
