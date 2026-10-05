/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseCifText } from './cif/text/parser.js';
import { parseCifBinary } from './cif/binary/parser.js';
import type { CifFrame } from './cif/data-model.js';
import { toDatabaseCollection, toDatabase } from './cif/schema.js';
import { mmCIF_Schema, type mmCIF_Database } from './cif/schema/mmcif.js';
import { CCD_Schema, type CCD_Database } from './cif/schema/ccd.js';
import { BIRD_Schema, type BIRD_Database } from './cif/schema/bird.js';
import { dic_Schema, type dic_Database } from './cif/schema/dic.js';
import { DensityServer_Data_Schema, type DensityServer_Data_Database } from './cif/schema/density-server.js';
import { type CifCore_Database, CifCore_Schema, CifCore_Aliases } from './cif/schema/cif-core.js';
import { type Segmentation_Data_Database, Segmentation_Data_Schema } from './cif/schema/segmentation.js';
import { SF_Schema, type SF_Database } from './cif/schema/sf.js';
import { StringLike } from '@molstar/core/util/string-like';

export const CIF = {
  parse: (data: StringLike | Uint8Array) => (StringLike.is(data) ? parseCifText(data) : parseCifBinary(data)),
  parseText: parseCifText,
  parseBinary: parseCifBinary,
  toDatabaseCollection,
  toDatabase,
  schema: {
    mmCIF: (frame: CifFrame) => toDatabase<mmCIF_Schema, mmCIF_Database>(mmCIF_Schema, frame),
    CCD: (frame: CifFrame) => toDatabase<CCD_Schema, CCD_Database>(CCD_Schema, frame),
    BIRD: (frame: CifFrame) => toDatabase<BIRD_Schema, BIRD_Database>(BIRD_Schema, frame),
    dic: (frame: CifFrame) => toDatabase<dic_Schema, dic_Database>(dic_Schema, frame),
    cifCore: (frame: CifFrame) => toDatabase<CifCore_Schema, CifCore_Database>(CifCore_Schema, frame, CifCore_Aliases),
    densityServer: (frame: CifFrame) =>
      toDatabase<DensityServer_Data_Schema, DensityServer_Data_Database>(DensityServer_Data_Schema, frame),
    segmentation: (frame: CifFrame) =>
      toDatabase<Segmentation_Data_Schema, Segmentation_Data_Database>(Segmentation_Data_Schema, frame),
    SF: (frame: CifFrame) => toDatabase<SF_Schema, SF_Database>(SF_Schema, frame),
  },
};

export * from './cif/data-model.js';
