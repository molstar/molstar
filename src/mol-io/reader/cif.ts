/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseCifText } from './cif/text/parser';
import { parseCifBinary } from './cif/binary/parser';
import { CifFrame } from './cif/data-model';
import { toDatabaseCollection, toDatabase } from './cif/schema';
import { mmCIF_Schema, mmCIF_Database } from './cif/schema/mmcif';
import { CCD_Schema, CCD_Database } from './cif/schema/ccd';
import { BIRD_Schema, BIRD_Database } from './cif/schema/bird';
import { dic_Schema, dic_Database } from './cif/schema/dic';
import { DensityServer_Data_Schema, DensityServer_Data_Database } from './cif/schema/density-server';
import { CifCore_Database, CifCore_Schema, CifCore_Aliases } from './cif/schema/cif-core';
import { Segmentation_Data_Database, Segmentation_Data_Schema } from './cif/schema/segmentation';
import { StringLike } from '../common/string-like';


export const CIF = {
    parse: (data: StringLike | Uint8Array) => StringLike.is(data) ? parseCifText(data) : parseCifBinary(data),
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
        densityServer: (frame: CifFrame) => toDatabase<DensityServer_Data_Schema, DensityServer_Data_Database>(DensityServer_Data_Schema, frame),
        segmentation: (frame: CifFrame) => toDatabase<Segmentation_Data_Schema, Segmentation_Data_Database>(Segmentation_Data_Schema, frame),
    }
};

export * from './cif/data-model';