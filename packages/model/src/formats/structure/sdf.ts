/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SdfFileCompound } from '@molstar/io/reader/sdf/parser';
import type { Trajectory } from '@molstar/model/model/structure';
import { Task } from '@molstar/core/task';
import type { ModelFormat } from '../format.js';
import { getMolModels } from './mol.js';

export { SdfFormat };

type SdfFormat = ModelFormat<SdfFileCompound>

namespace SdfFormat {
    export function is(x?: ModelFormat): x is SdfFormat {
        return x?.kind === 'sdf';
    }

    export function create(mol: SdfFileCompound): SdfFormat {
        return { kind: 'sdf', name: mol.molFile.title, data: mol };
    }
}

export function trajectoryFromSdf(mol: SdfFileCompound): Task<Trajectory> {
    return Task.create('Parse SDF', ctx => getMolModels(mol.molFile, SdfFormat.create(mol), ctx));
}
