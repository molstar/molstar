/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { exactAddress, type Node, parseAddress, resolveAddress } from '@molstar/core/repr/tree';
import { FormatTemplate } from '@molstar/core/util/string-format';
import { followRef, type TooltipSchema, TooltipDefaults } from '@molstar/seqx-core/lang';
import { tooltipNodes } from '@molstar/seqx-core/resolve/tooltip-tables';

// A view's tooltips resolved (_spec/README.md "Resolved view", Tooltips): tier 1, references and the parsed template,
// never rows. `views`' tooltips before the view's own, each in spec order. By data: the data node or the layer `on`
// names (a ref standing for one followed); by range: the addresses `on` names, or a table at the tooltip's address
// (a `ranges` table the store's, an `address` or MVS table read whole by the reader that shows it).

export interface ResolvedTooltip {
  readonly address: string;
  /** By data: the data node whose items speak (`on` a data node, or a layer's data). */
  readonly data?: string;
  /** By data on a layer: that layer alone (its marks only). */
  readonly layer?: string;
  /** By range from addresses (`on` an address or a list): each a one-row table. */
  readonly addresses?: readonly string[];
  /** By range from a table: its system (`target`) and how its rows locate; a `ranges` table's rows are the store's. */
  readonly target?: string;
  readonly schema?: TooltipSchema;
  /** A table read from a source (its address). */
  readonly source?: string;
  /** The template (absent: it cannot be read, `seq.lang.tooltip.text`), the field `{}` reads. */
  readonly template?: FormatTemplate;
  readonly field: string;
  readonly label?: string;
  readonly format: 'markdown' | 'text';
  /** A CIF category: unavailable until the structure runtime reads one (`seq.lang.tooltip.category`). */
  readonly category?: true;
  /** The tooltip node (an inline `address` or MVS table is read from it, whole, by the reader). */
  readonly node: Node;
}

/** A tooltip's template parsed (MVS's `text_format`), undefined when it cannot be read (`seq.lang.tooltip.text`). */
export function tooltipTemplate(text: string): FormatTemplate | undefined {
  try {
    return FormatTemplate(text);
  } catch {
    return undefined;
  }
}

/** The data kinds whose items speak by data (each its own marks); another node `on` is an address (by range). */
const ItemKinds: ReadonlySet<string> = new Set(['ranges', 'pairs']);

/** The node an exact address names (its leading ids), a ref followed when it stands for one exact address. */
function nodeAt(spec: Node, address: string): { node: Node; address: string } | undefined {
  const f = followRef(spec, address);
  const to = 'addresses' in f && f.addresses.length === 1 ? exactAddress(f.addresses[0]) : exactAddress(address);
  if (to === undefined) return undefined;
  const ids: string[] = [];
  for (const s of parseAddress(to).pattern?.segments ?? []) {
    if (typeof s !== 'string') return undefined;
    ids.push(s);
  }
  const node = resolveAddress(spec, ids);
  return node ? { node, address: to } : undefined;
}

/** A view's tooltips: `views`' and its own (`view` its address), resolved over the spec. */
export function resolveTooltips(spec: Node, view: string): ResolvedTooltip[] {
  const out: ResolvedTooltip[] = [];
  for (const t of tooltipNodes(spec)) {
    // `views`' every view's (`/views/<id>`), a view's its own (`/views/<view>/<id>`)
    const own = t.address.slice(0, t.address.lastIndexOf('/'));
    if (own !== view && own.split('/').length > 2) continue;
    const p = t.node.params ?? {};
    const field = typeof p.field === 'string' ? p.field : TooltipDefaults.field;
    const template = tooltipTemplate(typeof p.text === 'string' ? p.text : TooltipDefaults.text);
    const base = {
      address: t.address,
      field,
      format: p.format === 'text' ? ('text' as const) : ('markdown' as const),
      node: t.node,
      ...(template ? { template } : {}),
      ...(typeof p.label === 'string' ? { label: p.label } : {}),
    };
    if (p.on !== undefined) {
      const ons = (Array.isArray(p.on) ? p.on : [p.on]).filter((a): a is string => typeof a === 'string');
      // one address naming a data node or a layer: by data; else by range, each address a row
      const one = ons.length === 1 ? nodeAt(spec, ons[0]) : undefined;
      if (one?.node.kind === 'layer') {
        const data = exactAddress(one.node.params?.data);
        out.push({ ...base, layer: one.address, ...(data !== undefined ? { data } : {}) });
      } else if (one && ItemKinds.has(one.node.kind) && exactAddress(ons[0]) !== undefined)
        out.push({ ...base, data: one.address });
      else out.push({ ...base, addresses: ons });
      continue;
    }
    const target = exactAddress(p.target);
    if (target === undefined) continue;
    const source = exactAddress(p.source);
    out.push({
      ...base,
      target,
      schema: (typeof p.schema === 'string' ? p.schema : 'ranges') as TooltipSchema,
      ...(source !== undefined ? { source } : {}),
      ...(p.category !== undefined || p.block !== undefined ? { category: true as const } : {}),
    });
  }
  return out;
}
