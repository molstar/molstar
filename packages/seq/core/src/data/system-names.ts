/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// How a source's systems are matched to a requested name (_spec/README.md "Data", Systems by name): by the name or an
// alias the source lists, then by the requested name's own alias. An alignment's columns are named by the alignment,
// `<id>/columns`; `columns` is their alias, the name an alignment's own file gives them and older records sources
// used (kept for those sources: they are served as before).

/** The suffix of an alignment's columns' name, and the alias it keeps. */
const Columns = 'columns';
const ColumnsSuffix = `/${Columns}`;

/** The other names a requested system is served by: an alignment's columns (`<id>/columns`) also as `columns`. */
export function systemAliasOf(name: string): string | undefined {
  return name.endsWith(ColumnsSuffix) ? Columns : undefined;
}

/** The source's own name for a requested one, from its names map (its systems' names and aliases), the alias last. */
export function sourceNameOf(names: ReadonlyMap<string, string>, name: string): string | undefined {
  const own = names.get(name);
  if (own !== undefined) return own;
  const alias = systemAliasOf(name);
  return alias === undefined ? undefined : names.get(alias);
}

/** A coverage's system for a requested name: by its name or one of its aliases, then by the requested name's alias. */
export function coverageSystemOf<T extends { readonly name: string; readonly aliases?: readonly string[] }>(
  systems: readonly T[],
  name: string,
): T | undefined {
  const by = (n: string) => systems.find((s) => s.name === n || !!s.aliases?.includes(n));
  const own = by(name);
  if (own) return own;
  const alias = systemAliasOf(name);
  return alias === undefined ? undefined : by(alias);
}
