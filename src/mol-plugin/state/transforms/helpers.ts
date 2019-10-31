/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../../mol-model/structure';
import { ComputedSecondaryStructure } from '../../../mol-model-props/computed/secondary-structure';

/**
 * Attaches ComputedSecondaryStructure property when unavailable in sourceData and
 * when not an archival file (i.e. no database_2.database_id field)
 */
export async function ensureSecondaryStructure(s: Structure) {
    if (s.models.length === 1 && s.model && s.model.sourceData.kind === 'mmCIF') {
        if (!s.model.sourceData.data.struct_conf.id.isDefined && !s.model.sourceData.data.struct_sheet_range.id.isDefined &&
            !s.model.sourceData.data.database_2.database_id.isDefined
        ) {
            await ComputedSecondaryStructure.attachFromCifOrCompute(s)
        }
    }
}