/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import CIF from 'mol-io/reader/cif'
import { Run, Progress } from 'mol-task'
import { VolumeData, parseDensityServerData } from 'mol-model/volume'
import { DensityServer_Data_Database } from 'mol-io/reader/cif/schema/density-server';

export async function downloadCif(url: string, isBinary: boolean) {
    const data = await fetch(url);
    return parseCif(isBinary ? new Uint8Array(await data.arrayBuffer()) : await data.text());
}

export async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data)
    const parsed = await Run(comp, Progress.format);
    if (parsed.isError) throw parsed;
    return parsed.result
}

export type Volume = { source: DensityServer_Data_Database, volume: VolumeData }

export async function getVolumeFromEmdId(emdid: string): Promise<Volume> {
    const cif = await downloadCif(`https://webchem.ncbr.muni.cz/DensityServer/em/emd-${emdid}/cell?detail=4`, true)
    const data = CIF.schema.densityServer(cif.blocks[1])
    return { source: data, volume: await Run(parseDensityServerData(data)) }
}

export function resizeCanvas (canvas: HTMLCanvasElement, container: Element) {
    let w = window.innerWidth
    let h = window.innerHeight
    if (container !== document.body) {
        let bounds = container.getBoundingClientRect()
        w = bounds.right - bounds.left
        h = bounds.bottom - bounds.top
    }
    canvas.width = window.devicePixelRatio * w
    canvas.height = window.devicePixelRatio * h
    Object.assign(canvas.style, { width: `${w}px`, height: `${h}px` })
}