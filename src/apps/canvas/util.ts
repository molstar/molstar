/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readUrl, readFile, readUrlAsBuffer, readFileAsBuffer } from 'mol-util/read';
import CIF, { CifBlock } from 'mol-io/reader/cif'
import { Model, Format, StructureSymmetry, Structure } from 'mol-model/structure';
import CCP4 from 'mol-io/reader/ccp4/parser'
import { FileHandle } from 'mol-io/common/file-handle';
import { Ccp4File } from 'mol-io/reader/ccp4/schema';
import { volumeFromCcp4 } from 'mol-model/volume/formats/ccp4';
// import { parse as parseObj } from 'mol-io/reader/obj/parser'

// export async function getObjFromUrl(url: string) {
//     const data = await readUrlAs(url, false) as string
//     const comp = parseObj(data)
//     const parsed = await comp.run()
//     if (parsed.isError) throw parsed
//     return parsed.result
// }

export async function getCifFromData(data: string | Uint8Array) {
    const comp = CIF.parse(data)
    const parsed = await comp.run()
    if (parsed.isError) throw parsed
    return parsed.result.blocks[0]
}

export async function getCifFromUrl(url: string, binary = false) {
    return getCifFromData(await readUrl(url, binary))
}

export async function getCifFromFile(file: File, binary = false) {
    return getCifFromData(await readFile(file, binary))
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

//

export async function getCcp4FromUrl(url: string) {
    return getCcp4FromData(await readUrlAsBuffer(url))
}

export async function getCcp4FromFile(file: File) {
    return getCcp4FromData(await readFileAsBuffer(file))
}

export async function getCcp4FromData(data: Uint8Array) {
    const file = FileHandle.fromBuffer(data)
    const parsed = await CCP4(file).run()
    if (parsed.isError) throw parsed
    return parsed.result
}

export async function getVolumeFromCcp4(ccp4: Ccp4File) {
    return await volumeFromCcp4(ccp4).run()
}