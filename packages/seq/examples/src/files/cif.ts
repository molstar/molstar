/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The record scripts' mmCIF reader (the examples' one: `@molstar/io`'s parser would be a dependency of the package): a
// category's `loop_` or its single row, values whitespace-separated, quoted with ' or " (a quote inside a value when not
// followed by space), or a `;` text block. Enough for the PDB entries the datasets record; never run in a view.

/** The rows of a category (`_atom_site`), each a record by item name; none when the text has no such category. */
export function cifRows(text: string, category: string): Record<string, string>[] {
  const lines = text.split('\n');
  const prefix = `${category}.`;
  let i = lines.findIndex(
    (l, k) => l.startsWith(prefix) || (l.startsWith('loop_') && lines[k + 1]?.startsWith(prefix)),
  );
  if (i < 0) return [];
  const loop = lines[i].startsWith('loop_');
  if (loop) i++;
  const names: string[] = [];
  const tokens: string[] = [];
  const block = () => {
    let t = lines[i].slice(1);
    while (i + 1 < lines.length && !lines[++i].startsWith(';')) t += `\n${lines[i]}`;
    tokens.push(t);
  };
  if (loop) {
    for (; lines[i]?.startsWith(prefix); i++) names.push(lines[i].trim().slice(prefix.length));
    for (; i < lines.length; i++) {
      const l = lines[i];
      if (l.startsWith('#') || l.startsWith('loop_') || l.startsWith('_')) break;
      if (l.startsWith(';')) block();
      else tokens.push(...tokenize(l));
    }
  } else {
    for (; i < lines.length && (lines[i].startsWith(prefix) || lines[i].startsWith(';')); i++) {
      if (lines[i].startsWith(';')) {
        block();
        continue;
      }
      const m = /^\S+?\.(\S+)\s*(.*)$/.exec(lines[i])!;
      names.push(m[1]);
      if (m[2].trim()) tokens.push(...tokenize(m[2]));
    }
  }
  const rows: Record<string, string>[] = [];
  for (let r = 0; r + names.length <= tokens.length; r += names.length) {
    const row: Record<string, string> = {};
    names.forEach((n, k) => (row[n] = tokens[r + k]));
    rows.push(row);
    if (!loop) break;
  }
  return rows;
}

/** A category's rows as columns by item name (`atom_site`, no leading `_`): the same values, column-major. */
export function cifLoop(text: string, category: string): Record<string, string[]> {
  const rows = cifRows(text, `_${category}`);
  const out: Record<string, string[]> = {};
  for (const name of Object.keys(rows[0] ?? {})) out[name] = rows.map((r) => r[name]);
  return out;
}

/** A line's values: unquoted, a quote inside one kept when it is not followed by space (`O5'`). */
function tokenize(line: string): string[] {
  const out: string[] = [];
  const re = /'((?:[^']|'(?=\S))*)'(?=\s|$)|"((?:[^"]|"(?=\S))*)"(?=\s|$)|(\S+)/g;
  for (const m of line.matchAll(re)) out.push(m[1] ?? m[2] ?? m[3]);
  return out;
}
