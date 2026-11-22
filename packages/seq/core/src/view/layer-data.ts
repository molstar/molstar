/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's data (_spec/README.md "Resolved view", Layer kinds): the node or query it names, the forms after its
// annotation split into its base and the user's filter, and over an assembly the plan on its section's sequence.

import { exactAddress, type Node, parseAddress, printAddress, type Value } from '@molstar/core/repr/tree';
import { followRef } from '@molstar/seqx-core/lang';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { LayerData } from './view.js';

/** What a layer's data says after its annotation: the track's base, the user's filter, the forms in the base. */
export interface DataForms {
  /** The leading forms as one (`item` alone, else `and` of them): rows they reject are not counted. */
  readonly base?: Node;
  readonly lead: readonly Node[];
  /** A trailing `(except …)`: counted. */
  readonly filter?: Node;
  /** Forms that are neither (a segment that is an id, a location form in the base): read as no filter, and said. */
  readonly unread?: boolean;
}

const BaseKinds = new Set(['item', 'intersect', 'merge']);

/**
 * The one rule for a layer's data after its annotation (_spec/README.md "Themes", Filters): a leading run of item forms
 * (`intersect` and `merge` of them too) is the base, a trailing `(except …)` the user's filter, a lone form of another kind the
 * filter (said unread by the layers); any other mix is `unread`.
 */
export function splitDataForms(rest: readonly (string | Node)[]): DataForms {
  if (!rest.length) return { lead: [] };
  if (rest.some((f) => typeof f === 'string')) return { lead: [], unread: true };
  const forms = rest as Node[];
  const last = forms[forms.length - 1];
  const filter = last.kind === 'except' ? last : undefined;
  const lead = filter ? forms.slice(0, -1) : forms;
  // a lone form that is none of these stays the layer's filter: the layers say they cannot read it (`seq.render.filter`)
  if (forms.length === 1 && !filter && !BaseKinds.has(last.kind)) return { lead: [], filter: last };
  if (!lead.every((f) => BaseKinds.has(f.kind))) return { lead: [], unread: true };
  const base =
    lead.length === 0 ? undefined : lead.length === 1 ? lead[0] : ({ kind: 'intersect', params: { of: lead } } as Node);
  return { lead, ...(base ? { base } : {}), ...(filter ? { filter } : {}) };
}

/**
 * A layer's data with a user-defined ref followed: a ref (a chain of them) naming one address stands for it
 * (_spec/README.md, Language, Refs); a set, or a ref that cannot be followed, is left as written (and names nothing).
 */
export function layerData(spec: Node, value: Value | undefined): Value | undefined {
  if (typeof value !== 'string') return value;
  const f = followRef(spec, value);
  return 'addresses' in f && f.addresses.length === 1 ? f.addresses[0] : value;
}

/** A system's length as a layer's data carries it (none when unknown). */
function lengthOf(model: StaticModel, system: string): { length?: number } {
  const length = model.systems.get(system)?.length;
  return length === undefined ? {} : { length };
}

/** A layer's data: the section's letters when it names none; an address's node; a query's first node. */
export function dataOf(
  model: StaticModel,
  value: Value | undefined,
  type: string,
  section: string,
): LayerData | undefined {
  if (value === undefined)
    return type === 'letters' || type === 'ruler'
      ? { kind: 'address', address: section, system: section, shape: 'letters' }
      : undefined;
  if (typeof value !== 'string') return undefined;
  const parsed = parseAddress(value).pattern;
  if (!parsed) return undefined;
  const exact = exactAddress(value);
  const ids: string[] = [];
  for (const s of parsed.segments) {
    if (typeof s !== 'string') break;
    ids.push(s);
  }
  if (!ids.length) return undefined;
  const prefix = printAddress(ids);
  const kind = exact === undefined ? 'query' : 'address';
  const address = exact ?? printAddress(parsed.segments);
  const a = model.annotations.get(prefix);
  // an alignment is its columns system too: a layer naming it reads its cells
  if (a?.shape === 'cells') return { kind, address, system: a.target, shape: 'cells' };
  if (model.systems.has(prefix))
    return { kind, address, system: prefix, ...(kind === 'address' ? { shape: 'letters' as const } : {}) };
  const rest = parsed.segments.slice(ids.length);
  // the forms after the annotation: a leading run of item forms is the layer's base, a trailing `(except …)` the user's
  // filter (`splitDataForms`, one parser for the viewer's filter too)
  const split = kind === 'query' ? splitDataForms(rest) : undefined;
  const base = split?.base;
  const filter = split?.filter;
  const unread = split?.unread;
  // over an assembly: its plan on the section's (or the block's) sequence, which reads nothing without one (Resolve,
  // data over an assembly)
  const plan = a?.bySystem?.get(section);
  // held pairs over an assembly (inline links): on each of their units' blocks, their ends on its sequence
  const onUnit = !plan && !!a?.units && model.assemblies.has(a.target) && a.units.includes(section);
  return a
    ? {
        kind,
        address,
        system: plan || onUnit ? section : a.target,
        ...(plan ? { access: plan, ...lengthOf(model, section) } : {}),
        shape: a.shape,
        ...(base ? { base } : {}),
        ...(filter ? { filter } : {}),
        ...(unread ? { unread: true } : {}),
      }
    : undefined;
}
