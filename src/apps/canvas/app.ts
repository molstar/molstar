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

    async loadCif(cif: CifBlock, assemblyId?: string) {
        const models = await getModelsFromMmcif(cif)
        this.structureView = await StructureView(this.viewer, models, { assemblyId })
        this.pdbIdLoaded.next(this.structureView)
    }

    async loadPdbIdOrUrl(idOrUrl: string, options?: { assemblyId?: string, binary?: boolean }) {
        if (this.structureView) this.structureView.destroy();
        const url = idOrUrl.length <= 4 ? `https://files.rcsb.org/download/${idOrUrl}.cif` : idOrUrl;
        const cif = await getCifFromUrl(url, options ? !!options.binary : false)
        this.loadCif(cif, options ? options.assemblyId : void 0)
    }

    async loadCifFile(file: File) {
        if (this.structureView) this.structureView.destroy();
        const binary = /\.bcif$/.test(file.name);
        const cif = await getCifFromFile(file, binary)
        this.loadCif(cif)
    }
}