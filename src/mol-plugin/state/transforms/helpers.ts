/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../../mol-model/structure';
import { ComputedSecondaryStructure } from '../../../mol-model-props/computed/secondary-structure';
import { PluginStateObject } from '../objects';
import { DistanceData } from '../../../mol-repr/shape/loci/distance';
import { LabelData } from '../../../mol-repr/shape/loci/label';
import { OrientationData } from '../../../mol-repr/shape/loci/orientation';
import { AngleData } from '../../../mol-repr/shape/loci/angle';
import { DihedralData } from '../../../mol-repr/shape/loci/dihedral';

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

export function getDistanceDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): DistanceData {
    const lociA = s[0].loci
    const lociB = s[1].loci
    return { pairs: [{ lociA, lociB }] }
}

export function getAngleDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): AngleData {
    const lociA = s[0].loci
    const lociB = s[1].loci
    const lociC = s[2].loci
    return { triples: [{ lociA, lociB, lociC }] }
}

export function getDihedralDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): DihedralData {
    const lociA = s[0].loci
    const lociB = s[1].loci
    const lociC = s[2].loci
    const lociD = s[3].loci
    return { quads: [{ lociA, lociB, lociC, lociD }] }
}

export function getLabelDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): LabelData {
    const loci = s[0].loci
    return { infos: [{ loci }] }
}

export function getOrientationDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): OrientationData {
    const loci = s[0].loci
    return { loci }
}