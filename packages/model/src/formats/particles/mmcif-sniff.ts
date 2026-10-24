/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { binaryCifHasColumn, getBinaryCifHeader } from '@molstar/io/common/binary-cif';
import type { StringLike } from '@molstar/core/util/string-like';
import type { FileNameInfo } from '@molstar/core/util/file-info';
import { utf8Read } from '@molstar/core/util/utf8';
import { Sniffable } from '@molstar/core/util/sniffable';

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
export function looksLikeMmcifParticles(info: FileNameInfo, data: StringLike | Uint8Array | Sniffable): boolean {
  const sniffable = Sniffable.of(data);
  if (info.ext === 'bcif') {
    // read whole: the header walks the MessagePack (a body holding text holds no BinaryCIF)
    const bytes = sniffable.heldBytes;
    if (!bytes) return false;
    try {
      const bcifHeader = getBinaryCifHeader(bytes);
      const header = bcifHeader.dataBlocks[0]?.header;
      if (header && header.toLowerCase().includes('cellpack')) return true;
      if (binaryCifHasColumn(bcifHeader, PETWORLD_CATEGORY, PETWORLD_COLUMN)) return true;
    } catch (e) {
      console.error(e);
    }
    return false;
  } else if (info.ext === 'cif') {
    // the first 400 bytes, then a search in the form the data is held in (`Sniffable`)
    const head = utf8Read(sniffable.header, 0, Math.min(400, sniffable.header.length));
    if (head.toLowerCase().includes('cellpack')) return true;
    if (sniffable.sniffText(PETWORLD_FIELD)) return true;
  }
  return false;
}
