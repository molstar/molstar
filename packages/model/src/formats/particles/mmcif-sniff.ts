/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { binaryCifHasColumn, getBinaryCifHeader } from '@molstar/io/common/binary-cif';
import type { StringLike } from '@molstar/core/util/string-like';
import type { FileNameInfo } from '@molstar/core/util/file-info';

// Kept apart from mmcif.ts, which builds the particle lists: a host decides whether a file holds particles (and loads
// the particle code then) with this module alone.

const PETWORLD_CATEGORY = '_pdbx_struct_assembly_gen';
const PETWORLD_COLUMN = 'PDB_model_num';
const PETWORLD_FIELD = `${PETWORLD_CATEGORY}.${PETWORLD_COLUMN}`;

/**
 * Content-sniff a `.cif`/`.bcif` file (before full parsing) to detect whether it looks like
 * a CellPack or PetWorld mmCIF particle assembly, for use as a `DataFormatProvider.isApplicable`.
 * Returns `false` for ordinary/standard mmCIF structure files so those keep loading as
 * regular trajectories by default.
 */
export function looksLikeMmcifParticles(info: FileNameInfo, data: StringLike | Uint8Array): boolean {
  if (info.ext === 'bcif') {
    try {
      const bcifHeader = getBinaryCifHeader(data as Uint8Array);
      const header = bcifHeader.dataBlocks[0]?.header;
      if (header && header.toLowerCase().includes('cellpack')) return true;
      if (binaryCifHasColumn(bcifHeader, PETWORLD_CATEGORY, PETWORLD_COLUMN)) return true;
    } catch (e) {
      console.error(e);
    }
    return false;
  } else if (info.ext === 'cif') {
    const str = data as StringLike;
    const header = str.substring(0, 400);
    if (header.toLowerCase().includes('cellpack')) return true;
    if (str.includes(PETWORLD_FIELD)) return true;
  }
  return false;
}
