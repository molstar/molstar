/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, printAddress } from '@molstar/core/repr/tree';

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
