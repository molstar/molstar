/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from 'mol-model/structure/model/model';
import { Task } from 'mol-task';
import { ModelFormat } from './format';
import { _parse_mmCif } from './mmcif/parser';
import { CifFrame } from 'mol-io/reader/cif';

export function trajectoryFromMmCIF(frame: CifFrame): Task<Model.Trajectory> {
    return Task.create('Create mmCIF Model', ctx => _parse_mmCif(ModelFormat.mmCIF(frame), ctx));
}