/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF, { CifBlock } from 'mol-io/reader/cif'
import { readUrlAs } from 'mol-util/read';
import { Model, Format, StructureSymmetry, Structure } from 'mol-model/structure';
// import { parse as parseObj } from 'mol-io/reader/obj/parser'

// export async function getObjFromUrl(url: string) {
//     const data = await readUrlAs(url, false) as string
//     const comp = parseObj(data)
//     const parsed = await comp.run()
//     if (parsed.isError) throw parsed
//     return parsed.result
// }

export async function getCifFromUrl(url: string) {
    const data = await readUrlAs(url, false)
    const comp = CIF.parse(data)
    const parsed = await comp.run()
    if (parsed.isError) throw parsed
    return parsed.result.blocks[0]
}

export async function getModelsFromMmcif(cif: CifBlock) {
    return await Model.create(Format.mmCIF(cif)).run()
}

export async function getStructureFromModel(model: Model, assembly: string) {
    const assemblies = model.symmetry.assemblies
    if (assembly === '0') {
        return Structure.ofModel(model)
    } else if (assemblies.find(a => a.id === assembly)) {
        return await StructureSymmetry.buildAssembly(Structure.ofModel(model), assembly).run()
    }
}