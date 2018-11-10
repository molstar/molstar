/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Canvas3D from 'mol-canvas3d/canvas3d';
import { getCifFromUrl, getModelsFromMmcif, getCifFromFile, getCcp4FromUrl, getVolumeFromCcp4, getCcp4FromFile, getVolumeFromVolcif } from './util';
import { StructureView } from './structure-view';
import { BehaviorSubject } from 'rxjs';
import { CifBlock } from 'mol-io/reader/cif';
import { VolumeView } from './volume-view';
import { Ccp4File } from 'mol-io/reader/ccp4/schema';
import { Progress } from 'mol-task';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { StructureRepresentationRegistry } from 'mol-repr/structure/registry';
import { VolumeRepresentationRegistry } from 'mol-repr/volume/registry';

export class App {
    canvas3d: Canvas3D
    container: HTMLDivElement | null = null;
    canvas: HTMLCanvasElement | null = null;
    structureView: StructureView | null = null;
    volumeView: VolumeView | null = null;

    structureLoaded: BehaviorSubject<StructureView | null> = new BehaviorSubject<StructureView | null>(null)
    volumeLoaded: BehaviorSubject<VolumeView | null> = new BehaviorSubject<VolumeView | null>(null)

    colorThemeRegistry = new ColorTheme.Registry()
    sizeThemeRegistry = new SizeTheme.Registry()
    structureRepresentationRegistry = new StructureRepresentationRegistry()
    volumeRepresentationRegistry = new VolumeRepresentationRegistry()

    initViewer(_canvas: HTMLCanvasElement, _container: HTMLDivElement) {
        this.canvas = _canvas
        this.container = _container

        try {
            this.canvas3d = Canvas3D.create(this.canvas, this.container)
            this.canvas3d.animate()
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

    log(progress: Progress) {
        console.log(Progress.format(progress))
    }

    get reprCtx () {
        return {
            webgl: this.canvas3d.webgl,
            colorThemeRegistry: this.colorThemeRegistry,
            sizeThemeRegistry: this.sizeThemeRegistry
        }
    }

    //

    async loadMmcif(cif: CifBlock, assemblyId?: string) {
        const models = await this.runTask(getModelsFromMmcif(cif), 'Build models')
        this.structureView = await this.runTask(StructureView(this, this.canvas3d, models, { assemblyId }), 'Init structure view')
        this.structureLoaded.next(this.structureView)
    }

    async loadPdbIdOrMmcifUrl(idOrUrl: string, options?: { assemblyId?: string, binary?: boolean }) {
        if (this.structureView) this.structureView.destroy();
        const url = idOrUrl.length <= 4 ? `https://files.rcsb.org/download/${idOrUrl}.cif` : idOrUrl;
        const cif = await this.runTask(getCifFromUrl(url, options ? !!options.binary : false), 'Load mmCIF from URL')
        this.loadMmcif(cif.blocks[0], options ? options.assemblyId : void 0)
    }

    async loadMmcifFile(file: File) {
        if (this.structureView) this.structureView.destroy();
        const binary = /\.bcif$/.test(file.name);
        const cif = await this.runTask(getCifFromFile(file, binary), 'Load mmCIF from file')
        this.loadMmcif(cif.blocks[0])
    }

    //

    async loadCcp4(ccp4: Ccp4File) {
        const volume = await this.runTask(getVolumeFromCcp4(ccp4), 'Get Volume')
        this.volumeView = await this.runTask(VolumeView(this, this.canvas3d, volume), 'Init volume view')
        this.volumeLoaded.next(this.volumeView)
    }

    async loadCcp4File(file: File) {
        if (this.volumeView) this.volumeView.destroy();
        const ccp4 = await this.runTask(getCcp4FromFile(file), 'Load CCP4 from file')
        this.loadCcp4(ccp4)
    }

    async loadCcp4Url(url: string) {
        if (this.volumeView) this.volumeView.destroy();
        const ccp4 = await this.runTask(getCcp4FromUrl(url), 'Load CCP4 from URL')
        this.loadCcp4(ccp4)
    }

    //

    async loadVolcif(cif: CifBlock) {
        const volume = await this.runTask(getVolumeFromVolcif(cif), 'Get Volume')
        this.volumeView = await this.runTask(VolumeView(this, this.canvas3d, volume), 'Init volume view')
        this.volumeLoaded.next(this.volumeView)
    }

    async loadVolcifFile(file: File) {
        if (this.volumeView) this.volumeView.destroy();
        const binary = /\.bcif$/.test(file.name);
        const cif = await this.runTask(getCifFromFile(file, binary), 'Load volCif from file')
        this.loadVolcif(cif.blocks[1])
    }

    async loadVolcifUrl(url: string, binary?: boolean) {
        if (this.volumeView) this.volumeView.destroy();
        const cif = await this.runTask(getCifFromUrl(url, binary), 'Load volCif from URL')
        this.loadVolcif(cif.blocks[1])
    }
}