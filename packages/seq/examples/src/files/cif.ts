/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The record rules' view of an mmCIF file: a category's rows, each a record by item name, values as written (`.` and
// `?` kept). The file is parsed by `@molstar/io`'s text CIF parser in the record scripts (`scripts/cif.ts`, the one
// reader: the examples hold no parser of their own), so the rules stay synchronous and never run in a view.

/** A parsed file's categories. */
export interface Cif {
  /** The rows of a category (`_atom_site` or `atom_site`), each a record by item name; none when the file lacks it. */
  rows(category: string): readonly Record<string, string>[];
}

/** The rows of a category (`_atom_site`), each a record by item name; none when the file has no such category. */
export function cifRows(cif: Cif, category: string): readonly Record<string, string>[] {
  return cif.rows(category);
}

/** A category's rows as columns by item name (`atom_site`, no leading `_`): the same values, column-major. */
export function cifLoop(cif: Cif, category: string): Record<string, string[]> {
  const rows = cif.rows(category);
  const out: Record<string, string[]> = {};
  for (const name of Object.keys(rows[0] ?? {})) out[name] = rows.map((r) => r[name]);
  return out;
}
