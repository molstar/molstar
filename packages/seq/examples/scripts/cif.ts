/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The record scripts' mmCIF reader: `@molstar/io`'s text CIF parser (a dev dependency of the examples), its first data
// block as the rules' `Cif` (`src/files/cif.ts`): a category's rows by item name, `.` and `?` written back as they
// stand in the file, made once a category on first ask. The tests that read a recorded entry use it too.

import { parseCifText } from '@molstar/io/reader/cif/text/parser';
import type { Cif } from '../src/files/cif.js';

/** The written form of a value by the parser's kind: present (its text), not present (`.`), unknown (`?`). */
const Absent = ['', '.', '?'];

/** Parses a text mmCIF (its bytes or the text) into the rules' view of its first data block. */
export async function readCif(data: string | Uint8Array): Promise<Cif> {
  const parsed = await parseCifText(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const block = parsed.result.blocks[0];
  if (!block) throw new Error('no data block');
  const cache = new Map<string, Record<string, string>[]>();
  return {
    rows(category) {
      const name = category.startsWith('_') ? category.slice(1) : category;
      const cached = cache.get(name);
      if (cached) return cached;
      const rows: Record<string, string>[] = [];
      const c = block.categories[name];
      if (c) {
        const fields = c.fieldNames.map((f) => [f, c.getField(f)!] as const);
        for (let r = 0; r < c.rowCount; r++) {
          const row: Record<string, string> = {};
          for (const [f, field] of fields) {
            const kind = field.valueKind(r) as number;
            row[f] = kind === 0 ? field.str(r) : Absent[kind];
          }
          rows.push(row);
        }
      }
      cache.set(name, rows);
      return rows;
    },
  };
}
