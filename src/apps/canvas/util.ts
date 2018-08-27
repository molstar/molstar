/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF, { CifBlock } from 'mol-io/reader/cif'
import { readUrlAs } from "mol-util/read";
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

export async function getModelFromMmcif(cif: CifBlock) {
    const models = await Model.create(Format.mmCIF(cif)).run()
    return models[0]
}

export async function getStructureFromModel(model: Model, assembly = '0') {
    const assemblies = model.symmetry.assemblies
    if (assemblies.length) {
        return await StructureSymmetry.buildAssembly(Structure.ofModel(model), assembly).run()
    } else {
        return Structure.ofModel(model)
    }
}