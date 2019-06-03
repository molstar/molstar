/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PdbFile } from '../../mol-io/reader/pdb/schema';
import { pdbToMmCif } from './pdb/to-cif';
import { Model } from '../../mol-model/structure/model';
import { Task } from '../../mol-task';
import { ModelFormat } from './format';
import { _parse_mmCif } from './mmcif/parser';

export function trajectoryFromPDB(pdb: PdbFile): Task<Model.Trajectory> {
    return Task.create('Parse PDB', async ctx => {
        await ctx.update('Converting to mmCIF');
        const cif = await pdbToMmCif(pdb);
        const format = ModelFormat.mmCIF(cif);
        return _parse_mmCif(format, ctx);
    })
}
