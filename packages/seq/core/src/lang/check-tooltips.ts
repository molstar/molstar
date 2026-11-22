/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import type { Say } from './check-molecules.js';
import { MvsSchemaColumns, type TooltipSchema, TooltipDefaults, tooltipTemplate } from './tooltip-language.js';

// What the description cannot say of a tooltip (_spec/README.md, Language, Tooltips; `seq.lang.tooltip.*`): by data or
// by a table, never both; its template an f-string within the bound; a field it names in its table; a table fitting
// its schema; a CIF category said once (the structure runtime reads none yet).

/** The columns a `ranges` table's rows locate by, and the columns its items always have. */
const RangesColumns = ['start', 'end'];
const ItemFields = ['id', 'label', 'category', 'strand', 'start', 'end'];

/** A table's column names (an inline table's own; none when unknown). */
function columnsOf(rows: Value | undefined): string[] | undefined {
  if (!rows || typeof rows !== 'object' || Array.isArray(rows) || isColumn(rows)) return undefined;
  return Object.keys(rows as Record<string, Value>);
}

/** The fields a data node's items have, when they are known at check (its inline rows' columns); else undefined. */
function itemFields(data: Node | undefined): string[] | undefined {
  if (!data || data.kind !== 'ranges') return undefined;
  const cols = columnsOf(data.params?.rows);
  return cols && [...new Set([...ItemFields, ...cols])];
}

/**
 * The checks of one tooltip node; `target` names what a param names, `named` whether an address names a node of the
 * spec (by its leading ids), `maxText` the template's bound.
 */
export function checkTooltip(
  n: Node,
  target: (key: string) => Node | undefined,
  named: (address: string) => boolean,
  maxText: number,
  say: Say,
): void {
  const p = n.params ?? {};
  const on = p.on !== undefined;
  // an `on` naming nothing: that row inert, the rest of the spec opens (EC-250)
  const ons = Array.isArray(p.on) ? p.on : on ? [p.on] : [];
  ons.forEach((a, i) => {
    if (typeof a === 'string' && !named(a))
      say(
        'seq.lang.tooltip.on',
        `The tooltip's \`on\` names nothing in the spec: ${a}`,
        ...(Array.isArray(p.on) ? ['on', i] : ['on']),
      );
  });
  const table = p.rows !== undefined || p.source !== undefined;
  const ok = on
    ? !table && p.target === undefined
    : p.target !== undefined && (p.rows === undefined) !== (p.source === undefined);
  if (!ok)
    say(
      'seq.lang.tooltip.content',
      on
        ? 'A tooltip names `on` and a table: one of them'
        : 'A tooltip by range needs a `target` and one of `rows` and `source`',
      on ? 'on' : 'target',
    );
  // the template: an f-string within the bound (a 10,000-character template is refused here, EC-256)
  const field = typeof p.field === 'string' ? p.field : TooltipDefaults.field;
  const text = typeof p.text === 'string' ? p.text : TooltipDefaults.text;
  if (text.length > maxText) {
    say('seq.lang.tooltip.text', `The template is ${text.length} characters, past ${maxText}`, 'text');
    return;
  }
  const t = tooltipTemplate(text, field);
  if ('error' in t) {
    say('seq.lang.tooltip.text', `The template cannot be read: ${t.error}`, 'text');
    return;
  }
  // a CIF category: unavailable until the structure runtime reads one (said once, here)
  if (p.category !== undefined || p.block !== undefined)
    say('seq.lang.tooltip.category', 'A CIF category is not read yet: the tooltip shows nothing', 'category');
  // the table's columns, where known: the inline rows', or the data's items' (by data)
  const schema = (typeof p.schema === 'string' ? p.schema : 'ranges') as TooltipSchema;
  const columns = on ? itemFields(target('on')) : columnsOf(p.rows);
  if (!on && columns) {
    const locates = schema === 'ranges' ? RangesColumns : schema === 'address' ? ['address'] : MvsSchemaColumns[schema];
    const fits =
      schema === 'ranges'
        ? locates.every((c) => columns.includes(c))
        : !locates.length || locates.some((c) => columns.includes(c));
    if (!fits)
      say(
        'seq.lang.tooltip.schema',
        `The table's columns (${columns.join(', ')}) do not fit the schema ${schema} (${locates.join(', ')})`,
        'rows',
      );
  }
  if (columns)
    for (const f of t.fields)
      if (!columns.includes(f)) say('seq.lang.tooltip.field', `The template names ${f}, which the table lacks`, 'text');
}
