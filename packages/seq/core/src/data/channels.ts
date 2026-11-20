/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Coverage, FormatRuntime, OpenedSource } from './protocol.js';

// Values channels (_spec/README.md "Data", Store, Channels): a `values` node's `field` reads one column of a source
// opened once. The store holds a channel as a source of its own address, `<source>` and the field joined by a
// separator no address holds, so its tiles, its summaries and its entries are the field's alone, and every reader of a
// layer's source (render, describe, the query model) reads it by that address unchanged.

/** Between a source's address and a field in a channel's address: a control character no address or key holds. */
const Separator = '\u001f';

/** The address a source's values channel `field` is held under. */
export function channelAddress(source: string, field: string): string {
  return `${source}${Separator}${field}`;
}

/** The source and the field of a channel's address; undefined for a source's own. */
export function channelOf(address: string): { readonly source: string; readonly field: string } | undefined {
  const i = address.indexOf(Separator);
  return i < 0 ? undefined : { source: address.slice(0, i), field: address.slice(i + 1) };
}

/** The source an address reads: a channel's source, else itself. */
export function channelSource(address: string): string {
  const i = address.indexOf(Separator);
  return i < 0 ? address : address.slice(0, i);
}

/**
 * A source's channel: every read asks `field`; no nodes (the source's own attachment holds them, once); disposing it
 * leaves the source, which its own attachment disposes.
 */
export function fieldChannel(inner: OpenedSource, field: string): OpenedSource {
  const { nodes: _, ...coverage } = inner.coverage;
  return {
    coverage: coverage as Coverage,
    read: (request, signal) =>
      inner.read({ ...request, field }, signal).then((a) => (a.ok && a.nodes ? { ok: true, tables: a.tables } : a)),
    dispose() {},
  };
}

/**
 * The runtime each format is read with, and the formats two or more runtimes claim (_spec/README.md "Context",
 * Registry entries: one adapter a format). The same runtime twice is one; a conflicting format binds none, so nothing
 * depends on the order runtimes are listed in.
 */
export function runtimesByFormat(runtimes: readonly FormatRuntime[]): {
  readonly byFormat: ReadonlyMap<string, FormatRuntime>;
  readonly conflicts: ReadonlyMap<string, readonly string[]>;
} {
  const byFormat = new Map<string, FormatRuntime>();
  const conflicts = new Map<string, string[]>();
  for (const r of runtimes)
    for (const f of r.formats) {
      const had = byFormat.get(f);
      if (!had) {
        if (!conflicts.has(f)) byFormat.set(f, r);
        else if (!conflicts.get(f)!.includes(r.name)) conflicts.get(f)!.push(r.name);
      } else if (had !== r) {
        byFormat.delete(f);
        conflicts.set(f, [had.name, r.name]);
      }
    }
  return { byFormat, conflicts };
}

const ByFormat = new WeakMap<readonly FormatRuntime[], ReturnType<typeof runtimesByFormat>>();

/**
 * The one runtime of `runtimes` reading `format`: `undefined` when none does, `'conflict'` when several do (bound to
 * none, said). Kept per list.
 */
export function runtimeFor(runtimes: readonly FormatRuntime[], format: string): FormatRuntime | 'conflict' | undefined {
  let t = ByFormat.get(runtimes);
  if (!t) ByFormat.set(runtimes, (t = runtimesByFormat(runtimes)));
  return t.byFormat.get(format) ?? (t.conflicts.has(format) ? 'conflict' : undefined);
}
