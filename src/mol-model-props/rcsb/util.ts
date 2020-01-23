/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { AssemblySymmetry } from './assembly-symmetry';

export function getAssemblyIds(units: ReadonlyArray<Unit>) {
    const ids = new Set<string>()
    units.forEach(u => {
        if (u.conformation.operator.assembly) ids.add(u.conformation.operator.assembly.id)
    })
    return Array.from(ids.values())
}

export function getSymmetrySelectParam(structure?: Structure) {
    const param = PD.Select<number>(-1, [[-1, 'No Symmetries']])
    if (structure && structure.models[0].customProperties.has(AssemblySymmetry.Descriptor)) {
        const assemblySymmetry = AssemblySymmetry.get(structure.models[0])!
        const assemblyIds = getAssemblyIds(structure.units)
        const s = assemblySymmetry.getSymmetries(assemblyIds)
        if (s._rowCount) {
            const options: [number, string][] = []
            for (let i = 0, il = s._rowCount; i < il; ++i) {
                options.push([ s.id.value(i), `${s.assembly_id.value(i)}: ${s.symbol.value(i)} ${s.kind.value(i)}` ])
            }
            if (options.length) {
                param.options = options
                param.defaultValue = options[0][0]
            }
        }
    }
    return param
}