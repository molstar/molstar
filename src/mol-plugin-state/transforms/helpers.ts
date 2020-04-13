/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../objects';
import { DistanceData } from '../../mol-repr/shape/loci/distance';
import { LabelData } from '../../mol-repr/shape/loci/label';
import { OrientationData } from '../../mol-repr/shape/loci/orientation';
import { AngleData } from '../../mol-repr/shape/loci/angle';
import { DihedralData } from '../../mol-repr/shape/loci/dihedral';

export function getDistanceDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): DistanceData {
    const lociA = s[0].loci;
    const lociB = s[1].loci;
    return { pairs: [ { loci: [lociA, lociB] as const }] };
}

export function getAngleDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): AngleData {
    const lociA = s[0].loci;
    const lociB = s[1].loci;
    const lociC = s[2].loci;
    return { triples: [{ loci: [lociA, lociB, lociC] as const }] };
}

export function getDihedralDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): DihedralData {
    const lociA = s[0].loci;
    const lociB = s[1].loci;
    const lociC = s[2].loci;
    const lociD = s[3].loci;
    return { quads: [{ loci: [lociA, lociB, lociC, lociD] as const }] };
}

export function getLabelDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): LabelData {
    const loci = s[0].loci;
    return { infos: [{ loci }] };
}

export function getOrientationDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): OrientationData {
    const loci = s[0].loci;
    return { locis: [loci] };
}