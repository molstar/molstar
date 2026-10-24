/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CifField, CifFrame } from '@molstar/io/reader/cif';

// A ModelCIF's associated files (`ma_entry_associated_files`, and the files inside an archive,
// `ma_associated_archive_file_details`): Model Archive ships a model's PAE apart from it, as a CIF of its own, often in
// a zip with the MSA. Which file holds it is said by its content, "local pairwise QA scores".

/** Where a model's PAE is: a CIF at `url`, or the entry `path` of the zip at `url`. */
export interface AssociatedPae {
  readonly url: string;
  readonly path?: string;
}

const PairwiseContent = 'local pairwise qa scores';
const isCif = (format: string | undefined) => format === 'cif' || format === 'bcif';

/**
 * The associated file that holds the model's PAE, from the model's own file; undefined when it lists none. A relative
 * address is resolved against `base` (where the model came from), and dropped without one.
 */
export function associatedPaeOf(frame: CifFrame, base?: string): AssociatedPae | undefined {
  const files = frame.categories.ma_entry_associated_files;
  if (!files) return undefined;
  const field = (cat: typeof files, name: string) => cat.getField(name);
  const id = field(files, 'id');
  const url = field(files, 'file_url');
  const type = field(files, 'file_type');
  const format = field(files, 'file_format');
  const content = field(files, 'file_content');
  if (!url) return undefined;
  const lower = (f: CifField | undefined, i: number) => f?.str(i).toLowerCase();
  const resolved = (u: string) => {
    try {
      return new URL(u, base).href;
    } catch {
      return undefined;
    }
  };
  // A file of its own.
  for (let i = 0; i < files.rowCount; i++) {
    if (lower(type, i) === 'file' && isCif(lower(format, i)) && lower(content, i) === PairwiseContent) {
      const u = resolved(url.str(i));
      if (u) return { url: u };
    }
  }
  // An entry of an archive.
  const inner = frame.categories.ma_associated_archive_file_details;
  if (!inner || !id) return undefined;
  const archive = field(inner, 'archive_file_id');
  const path = field(inner, 'file_path');
  const innerFormat = field(inner, 'file_format');
  const innerContent = field(inner, 'file_content');
  if (!archive || !path) return undefined;
  for (let j = 0; j < inner.rowCount; j++) {
    if (!isCif(lower(innerFormat, j)) || lower(innerContent, j) !== PairwiseContent) continue;
    const of = archive.str(j);
    for (let i = 0; i < files.rowCount; i++) {
      if (id.str(i) !== of || lower(type, i) !== 'archive' || lower(format, i) !== 'zip') continue;
      const u = resolved(url.str(i));
      if (u) return { url: u, path: path.str(j) };
    }
  }
  return undefined;
}
