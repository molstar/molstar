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
import { Column } from '../../mol-data/db';
import { AtomPartialCharge } from './property/partial-charge';

export function trajectoryFromPDB(pdb: PdbFile): Task<Model.Trajectory> {
    return Task.create('Parse PDB', async ctx => {
        await ctx.update('Converting to mmCIF');
        const cif = await pdbToMmCif(pdb);
        const format = MmcifFormat.fromFrame(cif);
        const models = await createModels(format.data.db, format, ctx);
        const partial_charge = cif.categories['atom_site']?.getField('partial_charge');
        if (partial_charge) {
            // TODO works only for single, unsorted model, to work generally
            //      would need to do model splitting again
            if (models.length === 1) {
                const srcIndex = models[0].atomicHierarchy.atoms.sourceIndex;
                const isIdentity = Column.isIdentity(srcIndex);
                const srcIndexArray = isIdentity ? void 0 : srcIndex.toArray({ array: Int32Array });

                const q = partial_charge.toFloatArray();
                const partialCharge = srcIndexArray
                    ? Column.ofFloatArray(Column.mapToArray(srcIndex, i => q[i], Float32Array))
                    : Column.ofFloatArray(q);

                AtomPartialCharge.Provider.set(models[0], {
                    data: partialCharge,
                    type: 'GASTEIGER' // from PDBQT
                });
            }
        }
        return models;
    });
}
