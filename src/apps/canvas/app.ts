/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Viewer from 'mol-view/viewer';
import { getCifFromUrl, getModelsFromMmcif, getCifFromFile, getCcp4FromUrl } from './util';
import { StructureView } from './structure-view';
import { BehaviorSubject } from 'rxjs';
import { CifBlock } from 'mol-io/reader/cif';
import { volumeFromCcp4 } from 'mol-model/volume/formats/ccp4';
import { VolumeRepresentation } from 'mol-geo/representation/volume';
import SurfaceVisual from 'mol-geo/representation/volume/surface';
import { VolumeIsoValue } from 'mol-model/volume';

export class App {
    viewer: Viewer
    container: HTMLDivElement | null = null;
    canvas: HTMLCanvasElement | null = null;
    structureView: StructureView | null = null;

    pdbIdLoaded: BehaviorSubject<StructureView | null> = new BehaviorSubject<StructureView | null>(null)

    initViewer(_canvas: HTMLCanvasElement, _container: HTMLDivElement) {
        this.canvas = _canvas
        this.container = _container

        try {
            this.viewer = Viewer.create(this.canvas, this.container)
            this.viewer.animate()
            return true
        } catch (e) {
            console.error(e)
            return false
        }
    }

    setStatus(msg: string) {

    }

    private taskCount = 0
    taskCountChanged = new BehaviorSubject({ count: 0, info: '' })

    private changeTaskCount(delta: number, info = '') {
        this.taskCount += delta
        this.taskCountChanged.next({ count: this.taskCount, info })
    }

    async runTask<T>(promise: Promise<T>, info: string) {
        this.changeTaskCount(1, info)
        let result: T
        try {
            result = await promise
        } finally {
            this.changeTaskCount(-1)
        }
        return result
    }

    //

    async loadMmcif(cif: CifBlock, assemblyId?: string) {
        const models = await this.runTask(getModelsFromMmcif(cif), 'Build models')
        this.structureView = await this.runTask(StructureView(this, this.viewer, models, { assemblyId }), 'Init structure view')
        this.pdbIdLoaded.next(this.structureView)
    }

    async loadPdbIdOrMmcifUrl(idOrUrl: string, options?: { assemblyId?: string, binary?: boolean }) {
        if (this.structureView) this.structureView.destroy();
        const url = idOrUrl.length <= 4 ? `https://files.rcsb.org/download/${idOrUrl}.cif` : idOrUrl;
        const cif = await this.runTask(getCifFromUrl(url, options ? !!options.binary : false), 'Load mmCIF from URL')
        this.loadMmcif(cif, options ? options.assemblyId : void 0)
    }

    async loadMmcifFile(file: File) {
        if (this.structureView) this.structureView.destroy();
        const binary = /\.bcif$/.test(file.name);
        const cif = await this.runTask(getCifFromFile(file, binary), 'Load mmCIF from file')
        this.loadMmcif(cif)
    }

    //

    async loadCcp4File() {
        const url = 'http://localhost:8091/ngl/data/betaGal.mrc'
        const ccp4 = await getCcp4FromUrl(url)
        console.log(ccp4)
        const volume = await volumeFromCcp4(ccp4).run()
        const volRepr = VolumeRepresentation(SurfaceVisual)
        await volRepr.createOrUpdate({
            isoValue: VolumeIsoValue.relative(volume.dataStats, 1)
        }, volume).run()
        this.viewer.add(volRepr)
        console.log('volRepr', volRepr)
        this.viewer.requestDraw(true)
    }
}