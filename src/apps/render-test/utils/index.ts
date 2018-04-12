/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF from 'mol-io/reader/cif'
import { Run, Progress } from 'mol-task'
import { Structure } from 'mol-model/structure'

export function log(progress: Progress) {
    const p = progress.root.progress
    console.log(`${p.message} ${(p.current/p.max*100).toFixed(2)}%`)
}

export async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data)
    const parsed = await Run(comp, log, 100);
    if (parsed.isError) throw parsed;
    return parsed
}

export async function getStructuresFromPdbId(pdbid: string) {
    const data = await fetch(`https://files.rcsb.org/download/${pdbid}.cif`)
    const parsed = await parseCif(await data.text())
    return Structure.ofData({ kind: 'mmCIF', data: CIF.schema.mmCIF(parsed.result.blocks[0]) })
}