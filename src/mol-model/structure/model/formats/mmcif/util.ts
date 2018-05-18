/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Model from '../../model'

export function findEntityIdByAsymId(model: Model, asymId: string) {
    if (model.sourceData.kind !== 'mmCIF') return ''
    const { struct_asym } = model.sourceData.data
    for (let i = 0, n = struct_asym._rowCount; i < n; ++i) {
        if (struct_asym.id.value(i) === asymId) return struct_asym.entity_id.value(i)
    }
    return ''
}

export function findAtomIndexByLabelName(model: Model, residueIndex: number, atomName: string, altLoc: string | null) {
    const { segmentMap, segments } = model.atomicHierarchy.residueSegments
    const idx = segmentMap[residueIndex]
    const { label_atom_id, label_alt_id } = model.atomicHierarchy.atoms;
    for (let i = segments[idx], n = segments[idx + 1]; i <= n; ++i) {
        if (label_atom_id.value(i) === atomName && (!altLoc || label_alt_id.value(i) === altLoc)) return i;
    }
    return -1;
}