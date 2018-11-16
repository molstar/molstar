/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { BehaviorSubject } from 'rxjs';
import { App } from './app';
import { VolumeData } from 'mol-model/volume';
import { VolumeRepresentation } from 'mol-repr/volume/representation';
import { IsosurfaceRepresentation } from 'mol-repr/volume/isosurface-mesh';
import { DirectVolumeRepresentation } from 'mol-repr/volume/direct-volume';

export interface VolumeView {
    readonly app: App
    readonly viewer: Canvas3D

    readonly label: string
    readonly volume: VolumeData

    readonly active: { [k: string]: boolean }
    readonly volumeRepresentations: { [k: string]: VolumeRepresentation<any> }
    readonly updated: BehaviorSubject<null>

    setVolumeRepresentation(name: string, value: boolean): void
    destroy: () => void
}

interface VolumeViewProps {

}

export async function VolumeView(app: App, viewer: Canvas3D, volume: VolumeData, props: VolumeViewProps = {}): Promise<VolumeView> {
    const active: { [k: string]: boolean } = {
        isosurface: true,
        directVolume: false,
    }

    const volumeRepresentations: { [k: string]: VolumeRepresentation<any> } = {
        isosurface: IsosurfaceRepresentation(),
        directVolume: DirectVolumeRepresentation(),
    }

    const updated: BehaviorSubject<null> = new BehaviorSubject<null>(null)

    let label: string = 'Volume'

    async function setVolumeRepresentation(k: string, value: boolean) {
        active[k] = value
        await createVolumeRepr()
    }

    async function createVolumeRepr() {
        for (const k in volumeRepresentations) {
            if (active[k]) {
                await app.runTask(volumeRepresentations[k].createOrUpdate(app.reprCtx, {}, volume).run(
                    progress => app.log(progress)
                ), 'Create/update representation')
                viewer.add(volumeRepresentations[k])
            } else {
                viewer.remove(volumeRepresentations[k])
            }
        }

        // const center = Vec3.clone(volume.cell.size)
        // Vec3.scale(center, center, 0.5)
        // viewer.center(center)

        updated.next(null)
        viewer.requestDraw(true)
        console.log('stats', viewer.stats)
    }

    await createVolumeRepr()

    return {
        app,
        viewer,

        get label() { return label },
        volume,

        active,
        volumeRepresentations,
        setVolumeRepresentation,
        updated,

        destroy: () => {
            for (const k in volumeRepresentations) {
                viewer.remove(volumeRepresentations[k])
                volumeRepresentations[k].destroy()
            }
            viewer.requestDraw(true)
        }
    }
}