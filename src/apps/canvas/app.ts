/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Viewer from 'mol-view/viewer';
import { getCifFromUrl, getModelsFromMmcif } from './util';
import { StructureView } from './structure-view';
import { BehaviorSubject } from 'rxjs';

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

    async loadPdbId(id: string, assemblyId?: string) {
        if (this.structureView) this.structureView.destroy()
        const cif = await getCifFromUrl(`https://files.rcsb.org/download/${id}.cif`)
        const models = await getModelsFromMmcif(cif)
        this.structureView = await StructureView(this.viewer, models, { assemblyId })
        this.pdbIdLoaded.next(this.structureView)
    }
}