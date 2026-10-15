/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import type { CifFile, CifFrame } from '@molstar/io/reader/cif/data-model';
import { toDatabase } from '@molstar/io/reader/cif/schema';
import { mmCIF_Schema, type mmCIF_Database } from '@molstar/io/reader/cif/schema/mmcif';
import type { ModelFormat } from '../format.js';

/**
 * The mmCIF model format: its data type and the kind guard. Code that only has to recognize models read from mmCIF
 * (core model properties, the mmCIF export, color themes) imports it from here, so it does not load the mmCIF parser and
 * the property providers it registers (`mmcif.ts`). Models of this format are created by `mmcif.ts` and `pdb.ts`, which
 * import that module.
 */
export { MmcifFormat };

type MmcifFormat = ModelFormat<MmcifFormat.Data>;

namespace MmcifFormat {
  export type Data = {
    db: mmCIF_Database;
    frame: CifFrame;
    file?: CifFile;
    /**
     * Original source format. Some formats, including PDB, are converted
     * to mmCIF before further processing.
     */
    source?: ModelFormat;
  };
  export function is(x?: ModelFormat): x is MmcifFormat {
    return x?.kind === 'mmCIF';
  }

  export function fromFrame(frame: CifFrame, db?: mmCIF_Database, source?: ModelFormat, file?: CifFile): MmcifFormat {
    if (!db) db = toDatabase<mmCIF_Schema, mmCIF_Database>(mmCIF_Schema, frame);
    return { kind: 'mmCIF', name: db._name, data: { db, file, frame, source } };
  }
}
