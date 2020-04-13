/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PdbFile } from '../../mol-io/reader/pdb/schema';
import { pdbToMmCif } from './pdb/to-cif';
import { Model } from '../../mol-model/structure/model';
import { Task } from '../../mol-task';
import { MmcifFormat } from './mmcif';
import { createModels } from './basic/parser';

export function trajectoryFromPDB(pdb: PdbFile): Task<Model.Trajectory> {
    return Task.create('Parse PDB', async ctx => {
        await ctx.update('Converting to mmCIF');
        const cif = await pdbToMmCif(pdb);
        const format = MmcifFormat.fromFrame(cif);
        return createModels(format.data.db, format, ctx);
    });
}
