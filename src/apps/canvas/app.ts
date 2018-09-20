/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Viewer from 'mol-view/viewer';
import { getCifFromUrl, getModelsFromMmcif, getCifFromFile } from './util';
import { StructureView } from './structure-view';
import { BehaviorSubject } from 'rxjs';
import { CifBlock } from 'mol-io/reader/cif';

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

    async loadCif(cif: CifBlock, assemblyId?: string) {
        const models = await this.runTask(getModelsFromMmcif(cif), 'Build models')
        this.structureView = await this.runTask(StructureView(this, this.viewer, models, { assemblyId }), 'Init structure view')
        this.pdbIdLoaded.next(this.structureView)
    }

    async loadPdbIdOrUrl(idOrUrl: string, options?: { assemblyId?: string, binary?: boolean }) {
        if (this.structureView) this.structureView.destroy();
        const url = idOrUrl.length <= 4 ? `https://files.rcsb.org/download/${idOrUrl}.cif` : idOrUrl;
        const cif = await this.runTask(getCifFromUrl(url, options ? !!options.binary : false), 'Load mmCIF from URL')
        this.loadCif(cif, options ? options.assemblyId : void 0)
    }

    async loadCifFile(file: File) {
        if (this.structureView) this.structureView.destroy();
        const binary = /\.bcif$/.test(file.name);
        const cif = await this.runTask(getCifFromFile(file, binary), 'Load mmCIF from file')
        this.loadCif(cif)
    }
}