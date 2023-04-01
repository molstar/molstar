/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { Volume, Grid } from '../../../../mol-model/volume';
import { VolumeServerHeader, VolumeServerInfo } from './model';
import { Box3D } from '../../../../mol-math/geometry';
import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { Color } from '../../../../mol-util/color';
import { PluginBehavior } from '../../behavior';
import { LRUCache } from '../../../../mol-util/lru-cache';
import { urlCombine } from '../../../../mol-util/url';
import { CIF } from '../../../../mol-io/reader/cif';
import { volumeFromDensityServerData } from '../../../../mol-model-formats/volume/density-server';
import { PluginCommands } from '../../../commands';
import { StateSelection } from '../../../../mol-state';
import { StructureElement, Structure } from '../../../../mol-model/structure';
import { PluginContext } from '../../../context';
import { EmptyLoci, Loci, isEmptyLoci } from '../../../../mol-model/loci';
import { Asset } from '../../../../mol-util/assets';
import { GlobalModelTransformInfo } from '../../../../mol-model/structure/model/properties/global-transform';
import { distinctUntilChanged, filter, map, Observable, throttleTime } from 'rxjs';
import { Camera } from '../../../../mol-canvas3d/camera';
import { PluginCommand } from '../../../command';
import { SingleAsyncQueue } from '../../../../mol-util/single-async-queue';

export class VolumeStreaming extends PluginStateObject.CreateBehavior<VolumeStreaming.Behavior>({ name: 'Volume Streaming' }) { }

export namespace VolumeStreaming {
    export const RootTag = 'volume-streaming-info';

    export interface ChannelParams {
        isoValue: Volume.IsoValue,
        color: Color,
        wireframe: boolean,
        opacity: number
    }

    function channelParam(label: string, color: Color, defaultValue: Volume.IsoValue, stats: Grid['stats'], defaults: Partial<ChannelParams> = {}) {
        return PD.Group<ChannelParams>({
            isoValue: Volume.createIsoValueParam(defaults.isoValue ?? defaultValue, stats),
            color: PD.Color(defaults.color ?? color),
            wireframe: PD.Boolean(defaults.wireframe ?? false),
            opacity: PD.Numeric(defaults.opacity ?? 0.3, { min: 0, max: 1, step: 0.01 })
        }, { label, isExpanded: true });
    }

    const fakeSampling: VolumeServerHeader.Sampling = {
        byteOffset: 0,
        rate: 1,
        sampleCount: [1, 1, 1],
        valuesInfo: [{ mean: 0, min: -1, max: 1, sigma: 0.1 }, { mean: 0, min: -1, max: 1, sigma: 0.1 }]
    };

    export function createParams(options: { data?: VolumeServerInfo.Data, defaultView?: ViewTypes, channelParams?: DefaultChannelParams } = {}) {
        const { data, defaultView, channelParams } = options;
        const map = new Map<string, VolumeServerInfo.EntryData>();
        if (data) data.entries.forEach(d => map.set(d.dataId, d));
        const names = data ? data.entries.map(d => [d.dataId, d.dataId] as [string, string]) : [];
        const defaultKey = data ? data.entries[0].dataId : '';
        return {
            entry: PD.Mapped<EntryParams>(defaultKey, names, name => PD.Group(createEntryParams({ entryData: map.get(name)!, defaultView, structure: data && data.structure, channelParams }))),
        };
    }

    export type EntryParamDefinition = ReturnType<typeof createEntryParams>
    export type EntryParams = PD.Values<EntryParamDefinition>

    export function createEntryParams(options: { entryData?: VolumeServerInfo.EntryData, defaultView?: ViewTypes, structure?: Structure, channelParams?: DefaultChannelParams }) {
        const { entryData, defaultView, structure, channelParams = {} } = options;

        // fake the info
        const info = entryData || { kind: 'em', header: { sampling: [fakeSampling], availablePrecisions: [{ precision: 0, maxVoxels: 0 }] }, emDefaultContourLevel: Volume.IsoValue.relative(0) };
        const box = (structure && structure.boundary.box) || Box3D();

        return {
            view: PD.MappedStatic(defaultView || (info.kind === 'em' ? 'auto' : 'selection-box'), {
                'off': PD.Group<{}>({}),
                'box': PD.Group({
                    bottomLeft: PD.Vec3(box.min),
                    topRight: PD.Vec3(box.max),
                }, { description: 'Static box defined by cartesian coords.', isFlat: true }),
                'selection-box': PD.Group({
                    radius: PD.Numeric(5, { min: 0, max: 50, step: 0.5 }, { description: 'Radius in \u212B within which the volume is shown.' }),
                    bottomLeft: PD.Vec3(Vec3.create(0, 0, 0), {}, { isHidden: true }),
                    topRight: PD.Vec3(Vec3.create(0, 0, 0), {}, { isHidden: true }),
                }, { description: 'Box around focused element.', isFlat: true }),
                'camera-target': PD.Group({
                    radius: PD.Numeric(0.5, { min: 0, max: 1, step: 0.05 }, { description: 'Radius within which the volume is shown (relative to the field of view).' }),
                    // Minimal detail level for the inside of the zoomed region (real detail can be higher, depending on the region size)
                    dynamicDetailLevel: createDetailParams(info.header.availablePrecisions, 0, { label: 'Dynamic Detail' }),
                    bottomLeft: PD.Vec3(Vec3.create(0, 0, 0), {}, { isHidden: true }),
                    topRight: PD.Vec3(Vec3.create(0, 0, 0), {}, { isHidden: true }),
                }, { description: 'Box around camera target.', isFlat: true }),
                'cell': PD.Group<{}>({}),
                // Show selection-box if available and cell otherwise.
                'auto': PD.Group({
                    radius: PD.Numeric(5, { min: 0, max: 50, step: 0.5 }, { description: 'Radius in \u212B within which the volume is shown.' }),
                    selectionDetailLevel: createDetailParams(info.header.availablePrecisions, 6, { label: 'Selection Detail' }),
                    isSelection: PD.Boolean(false, { isHidden: true }),
                    bottomLeft: PD.Vec3(box.min, {}, { isHidden: true }),
                    topRight: PD.Vec3(box.max, {}, { isHidden: true }),
                }, { description: 'Box around focused element.', isFlat: true })
            }, { options: ViewTypeOptions, description: 'Controls what of the volume is displayed. "Off" hides the volume alltogether. "Bounded box" shows the volume inside the given box. "Around Interaction" shows the volume around the focused element/atom. "Whole Structure" shows the volume for the whole structure.' }),
            detailLevel: createDetailParams(info.header.availablePrecisions, 3),
            channels: info.kind === 'em'
                ? PD.Group({
                    'em': channelParam('EM', Color(0x638F8F), info.emDefaultContourLevel || Volume.IsoValue.relative(1), info.header.sampling[0].valuesInfo[0], channelParams['em'])
                }, { isFlat: true })
                : PD.Group({
                    '2fo-fc': channelParam('2Fo-Fc', Color(0x3362B2), Volume.IsoValue.relative(1.5), info.header.sampling[0].valuesInfo[0], channelParams['2fo-fc']),
                    'fo-fc(+ve)': channelParam('Fo-Fc(+ve)', Color(0x33BB33), Volume.IsoValue.relative(3), info.header.sampling[0].valuesInfo[1], channelParams['fo-fc(+ve)']),
                    'fo-fc(-ve)': channelParam('Fo-Fc(-ve)', Color(0xBB3333), Volume.IsoValue.relative(-3), info.header.sampling[0].valuesInfo[1], channelParams['fo-fc(-ve)']),
                }, { isFlat: true }),
        };
    }

    function createDetailParams(availablePrecisions: VolumeServerHeader.DetailLevel[], preferredPrecision: number, info?: PD.Info) {
        return PD.Select<number>(Math.min(preferredPrecision, availablePrecisions.length - 1),
            availablePrecisions.map((p, i) => [i, `${i + 1} [ ${Math.pow(p.maxVoxels, 1 / 3) | 0}^3 cells ]`] as [number, string]),
            {
                description: 'Determines the maximum number of voxels. Depending on the size of the volume options are in the range from 1 (0.52M voxels) to 7 (25.17M voxels).',
                ...info
            }
        );
    }

    export function copyParams(origParams: Params): Params {
        return {
            entry: {
                name: origParams.entry.name,
                params: {
                    detailLevel: origParams.entry.params.detailLevel,
                    channels: origParams.entry.params.channels,
                    view: {
                        name: origParams.entry.params.view.name,
                        params: { ...origParams.entry.params.view.params } as any,
                    }
                }
            }
        };
    }

    export const ViewTypeOptions = [['off', 'Off'], ['box', 'Bounded Box'], ['selection-box', 'Around Focus'], ['camera-target', 'Around Camera'], ['cell', 'Whole Structure'], ['auto', 'Auto']] as [ViewTypes, string][];

    export type ViewTypes = 'off' | 'box' | 'selection-box' | 'camera-target' | 'cell' | 'auto'

    export type ParamDefinition = ReturnType<typeof createParams>
    export type Params = PD.Values<ParamDefinition>


    type ChannelsInfo = { [name in ChannelType]?: { isoValue: Volume.IsoValue, color: Color, wireframe: boolean, opacity: number } }
    type ChannelsData = { [name in 'EM' | '2FO-FC' | 'FO-FC']?: Volume }

    export type ChannelType = 'em' | '2fo-fc' | 'fo-fc(+ve)' | 'fo-fc(-ve)'
    export const ChannelTypeOptions: [ChannelType, string][] = [['em', 'em'], ['2fo-fc', '2fo-fc'], ['fo-fc(+ve)', 'fo-fc(+ve)'], ['fo-fc(-ve)', 'fo-fc(-ve)']];
    export interface ChannelInfo {
        data: Volume,
        color: Color,
        wireframe: boolean,
        isoValue: Volume.IsoValue.Relative,
        opacity: number
    }
    export type Channels = { [name in ChannelType]?: ChannelInfo }

    export type DefaultChannelParams = { [name in ChannelType]?: Partial<ChannelParams> }

    export class Behavior extends PluginBehavior.WithSubscribers<Params> {
        private cache = LRUCache.create<{ data: ChannelsData, asset: Asset.Wrapper }>(25);
        public params: Params = {} as any;
        private lastLoci: StructureElement.Loci | EmptyLoci = EmptyLoci;
        private ref: string = '';
        public infoMap: Map<string, VolumeServerInfo.EntryData>;
        private updateQueue: SingleAsyncQueue;
        private cameraTargetObservable = this.plugin.canvas3d!.didDraw!.pipe(
            throttleTime(500, undefined, { 'leading': true, 'trailing': true }),
            map(() => this.plugin.canvas3d?.camera.getSnapshot()),
            distinctUntilChanged((a, b) => this.isCameraTargetSame(a, b)),
            filter(a => a !== undefined),
        ) as Observable<Camera.Snapshot>;
        private cameraTargetSubscription?: PluginCommand.Subscription = undefined;

        channels: Channels = {};

        public get info() {
            return this.infoMap.get(this.params.entry.name)!;
        }

        private async queryData(box?: Box3D) {
            let url = urlCombine(this.data.serverUrl, `${this.info.kind}/${this.info.dataId.toLowerCase()}`);

            if (box) {
                const { min: a, max: b } = box;
                url += `/box`
                    + `/${a.map(v => Math.round(1000 * v) / 1000).join(',')}`
                    + `/${b.map(v => Math.round(1000 * v) / 1000).join(',')}`;
            } else {
                url += `/cell`;
            }

            let detail = this.params.entry.params.detailLevel;
            if (this.params.entry.params.view.name === 'auto' && this.params.entry.params.view.params.isSelection) {
                detail = this.params.entry.params.view.params.selectionDetailLevel;
            }
            if (this.params.entry.params.view.name === 'camera-target' && box) {
                detail = this.decideDetail(box, this.params.entry.params.view.params.dynamicDetailLevel);
            }

            url += `?detail=${detail}`;

            const entry = LRUCache.get(this.cache, url);
            if (entry) return entry.data;

            const urlAsset = Asset.getUrlAsset(this.plugin.managers.asset, url);
            const asset = await this.plugin.runTask(this.plugin.managers.asset.resolve(urlAsset, 'binary'));
            const data = await this.parseCif(asset.data);
            if (!data) return;

            const removed = LRUCache.set(this.cache, url, { data, asset });
            if (removed) removed.asset.dispose();
            return data;
        }

        private async parseCif(data: Uint8Array): Promise<ChannelsData | undefined> {
            const parsed = await this.plugin.runTask(CIF.parseBinary(data));
            if (parsed.isError) {
                this.plugin.log.error('VolumeStreaming, parsing CIF: ' + parsed.toString());
                return;
            }
            if (parsed.result.blocks.length < 2) {
                this.plugin.log.error('VolumeStreaming: Invalid data.');
                return;
            }

            const ret: ChannelsData = {};
            for (let i = 1; i < parsed.result.blocks.length; i++) {
                const block = parsed.result.blocks[i];

                const densityServerCif = CIF.schema.densityServer(block);
                const volume = await this.plugin.runTask(volumeFromDensityServerData(densityServerCif));
                (ret as any)[block.header as any] = volume;
            }
            return ret;
        }

        private async updateParams(box: Box3D | undefined, autoIsSelection: boolean = false) {
            const newParams = copyParams(this.params);
            const viewType = newParams.entry.params.view.name;
            if (viewType !== 'off' && viewType !== 'cell') {
                newParams.entry.params.view.params.bottomLeft = box?.min || Vec3.zero();
                newParams.entry.params.view.params.topRight = box?.max || Vec3.zero();
            }
            if (viewType === 'auto') {
                newParams.entry.params.view.params.isSelection = autoIsSelection;
            }

            const state = this.plugin.state.data;
            const update = state.build().to(this.ref).update(newParams);

            await PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotUpdateCurrent: true } });
        }

        private getStructureRoot() {
            return this.plugin.state.data.select(StateSelection.Generators.byRef(this.ref).rootOfType(PluginStateObject.Molecule.Structure))[0];
        }

        register(ref: string): void {
            this.ref = ref;

            this.subscribeObservable(this.plugin.state.events.object.removed, o => {
                if (!PluginStateObject.Molecule.Structure.is(o.obj) || !StructureElement.Loci.is(this.lastLoci)) return;
                if (this.lastLoci.structure === o.obj.data) {
                    this.lastLoci = EmptyLoci;
                }
            });

            this.subscribeObservable(this.plugin.state.events.object.updated, o => {
                if (!PluginStateObject.Molecule.Structure.is(o.oldObj) || !StructureElement.Loci.is(this.lastLoci)) return;
                if (this.lastLoci.structure === o.oldObj.data) {
                    this.lastLoci = EmptyLoci;
                }
            });

            this.subscribeObservable(this.plugin.managers.structure.focus.behaviors.current, (entry) => {
                if (!this.plugin.state.data.tree.children.has(this.ref)) return;

                const loci = entry ? entry.loci : EmptyLoci;

                switch (this.params.entry.params.view.name) {
                    case 'auto':
                        this.updateAuto(loci);
                        break;
                    case 'selection-box':
                        this.updateSelectionBox(loci);
                        break;
                    default:
                        this.lastLoci = loci;
                        break;
                }
            });
        }

        unregister() {
            let entry = this.cache.entries.first;
            while (entry) {
                entry.value.data.asset.dispose();
                entry = entry.next;
            }
        }

        private isCameraTargetSame(a?: Camera.Snapshot, b?: Camera.Snapshot): boolean {
            if (!a || !b) return false;
            const targetSame = Vec3.equals(a.target, b.target);
            const sqDistA = Vec3.squaredDistance(a.target, a.position);
            const sqDistB = Vec3.squaredDistance(b.target, b.position);
            const distanceSame = Math.abs(sqDistA - sqDistB) / sqDistA < 1e-3;
            return targetSame && distanceSame;
        }
        private cameraTargetDistance(snapshot: Camera.Snapshot): number {
            return Vec3.distance(snapshot.target, snapshot.position);
        }

        private _invTransform: Mat4 = Mat4();
        private getBoxFromLoci(loci: StructureElement.Loci | EmptyLoci): Box3D {
            if (Loci.isEmpty(loci) || isEmptyLoci(loci)) {
                return Box3D();
            }

            const parent = this.plugin.helpers.substructureParent.get(loci.structure, true);
            if (!parent) return Box3D();
            const root = this.getStructureRoot();
            if (!root || root.obj?.data !== parent.obj?.data) return Box3D();

            const transform = GlobalModelTransformInfo.get(root.obj?.data.models[0]!);
            if (transform) Mat4.invert(this._invTransform, transform);

            const extendedLoci = StructureElement.Loci.extendToWholeResidues(loci);
            const box = StructureElement.Loci.getBoundary(extendedLoci, transform && !Number.isNaN(this._invTransform[0]) ? this._invTransform : void 0).box;

            if (StructureElement.Loci.size(extendedLoci) === 1) {
                Box3D.expand(box, box, Vec3.create(1, 1, 1));
            }

            return box;
        }

        private updateAuto(loci: StructureElement.Loci | EmptyLoci) {
            this.updateQueue.enqueue(async () => {
                this.lastLoci = loci;
                if (isEmptyLoci(loci)) {
                    await this.updateParams(this.info.kind === 'x-ray' ? this.data.structure.boundary.box : void 0, false);
                } else {
                    await this.updateParams(this.getBoxFromLoci(loci), true);
                }
            });
        }

        private updateSelectionBox(loci: StructureElement.Loci | EmptyLoci) {
            this.updateQueue.enqueue(async () => {
                if (Loci.areEqual(this.lastLoci, loci)) {
                    this.lastLoci = EmptyLoci;
                } else {
                    this.lastLoci = loci;
                }
                const box = this.getBoxFromLoci(this.lastLoci);
                await this.updateParams(box);
            });
        }

        private updateCameraTarget(snapshot: Camera.Snapshot) {
            this.updateQueue.enqueue(async () => {
                const origManualReset = this.plugin.canvas3d?.props.camera.manualReset;
                try {
                    if (!origManualReset) this.plugin.canvas3d?.setProps({ camera: { manualReset: true } });
                    const box = this.boxFromCameraTarget(snapshot, true);
                    await this.updateParams(box);
                } finally {
                    if (!origManualReset) this.plugin.canvas3d?.setProps({ camera: { manualReset: origManualReset } });
                }
            });
        }

        private boxFromCameraTarget(snapshot: Camera.Snapshot, boundByBoundarySize: boolean): Box3D {
            const target = snapshot.target;
            const distance = this.cameraTargetDistance(snapshot);
            const top = Math.tan(0.5 * snapshot.fov) * distance;
            let radius = top;
            const viewport = this.plugin.canvas3d?.camera.viewport;
            if (viewport && viewport.width > viewport.height) {
                radius *= viewport.width / viewport.height;
            }
            const relativeRadius = this.params.entry.params.view.name === 'camera-target' ? this.params.entry.params.view.params.radius : 0.5;
            radius *= relativeRadius;
            let radiusX, radiusY, radiusZ;
            if (boundByBoundarySize) {
                const bBoxSize = Vec3.zero();
                Box3D.size(bBoxSize, this.data.structure.boundary.box);
                radiusX = Math.min(radius, 0.5 * bBoxSize[0]);
                radiusY = Math.min(radius, 0.5 * bBoxSize[1]);
                radiusZ = Math.min(radius, 0.5 * bBoxSize[2]);
            } else {
                radiusX = radiusY = radiusZ = radius;
            }
            return Box3D.create(
                Vec3.create(target[0] - radiusX, target[1] - radiusY, target[2] - radiusZ),
                Vec3.create(target[0] + radiusX, target[1] + radiusY, target[2] + radiusZ)
            );
        }

        private decideDetail(box: Box3D, baseDetail: number): number {
            const cellVolume = this.info.kind === 'x-ray'
                ? Box3D.volume(this.data.structure.boundary.box)
                : this.info.header.spacegroup.size.reduce((a, b) => a * b, 1);
            const boxVolume = Box3D.volume(box);
            let ratio = boxVolume / cellVolume;
            const maxDetail = this.info.header.availablePrecisions.length - 1;
            let detail = baseDetail;
            while (ratio <= 0.5 && detail < maxDetail) {
                ratio *= 2;
                detail += 1;
            }
            // console.log(`Decided dynamic detail: ${detail}, (base detail: ${baseDetail}, box/cell volume ratio: ${boxVolume / cellVolume})`);
            return detail;
        }

        async update(params: Params) {
            const switchedToSelection = params.entry.params.view.name === 'selection-box' && this.params && this.params.entry && this.params.entry.params && this.params.entry.params.view && this.params.entry.params.view.name !== 'selection-box';

            this.params = params;
            let box: Box3D | undefined = void 0, emptyData = false;

            if (params.entry.params.view.name !== 'camera-target' && this.cameraTargetSubscription) {
                this.cameraTargetSubscription.unsubscribe();
                this.cameraTargetSubscription = undefined;
            }

            switch (params.entry.params.view.name) {
                case 'off':
                    emptyData = true;
                    break;
                case 'box':
                    box = Box3D.create(params.entry.params.view.params.bottomLeft, params.entry.params.view.params.topRight);
                    emptyData = Box3D.volume(box) < 0.0001;
                    break;
                case 'selection-box': {
                    if (switchedToSelection) {
                        box = this.getBoxFromLoci(this.lastLoci) || Box3D();
                    } else {
                        box = Box3D.create(Vec3.clone(params.entry.params.view.params.bottomLeft), Vec3.clone(params.entry.params.view.params.topRight));
                    }
                    const r = params.entry.params.view.params.radius;
                    emptyData = Box3D.volume(box) < 0.0001;
                    Box3D.expand(box, box, Vec3.create(r, r, r));
                    break;
                }
                case 'camera-target':
                    if (!this.cameraTargetSubscription) {
                        this.cameraTargetSubscription = this.subscribeObservable(this.cameraTargetObservable, (e) => this.updateCameraTarget(e));
                    }
                    box = this.boxFromCameraTarget(this.plugin.canvas3d!.camera.getSnapshot(), true);
                    break;
                case 'cell':
                    box = this.info.kind === 'x-ray'
                        ? this.data.structure.boundary.box
                        : void 0;
                    break;
                case 'auto':
                    box = params.entry.params.view.params.isSelection || this.info.kind === 'x-ray'
                        ? Box3D.create(Vec3.clone(params.entry.params.view.params.bottomLeft), Vec3.clone(params.entry.params.view.params.topRight))
                        : void 0;

                    if (box) {
                        emptyData = Box3D.volume(box) < 0.0001;
                        if (params.entry.params.view.params.isSelection) {
                            const r = params.entry.params.view.params.radius;
                            Box3D.expand(box, box, Vec3.create(r, r, r));
                        }
                    }

                    break;
            }

            const data = emptyData ? {} : await this.queryData(box);

            if (!data) return false;

            const info = params.entry.params.channels as ChannelsInfo;

            if (this.info.kind === 'x-ray') {
                this.channels['2fo-fc'] = this.createChannel(data['2FO-FC'] || Volume.One, info['2fo-fc'], this.info.header.sampling[0].valuesInfo[0]);
                this.channels['fo-fc(+ve)'] = this.createChannel(data['FO-FC'] || Volume.One, info['fo-fc(+ve)'], this.info.header.sampling[0].valuesInfo[1]);
                this.channels['fo-fc(-ve)'] = this.createChannel(data['FO-FC'] || Volume.One, info['fo-fc(-ve)'], this.info.header.sampling[0].valuesInfo[1]);
            } else {
                this.channels['em'] = this.createChannel(data['EM'] || Volume.One, info['em'], this.info.header.sampling[0].valuesInfo[0]);
            }

            return true;
        }

        private createChannel(data: Volume, info: ChannelsInfo['em'], stats: Grid['stats']): ChannelInfo {
            const i = info!;
            return {
                data,
                color: i.color,
                wireframe: i.wireframe,
                opacity: i.opacity,
                isoValue: i.isoValue.kind === 'relative' ? i.isoValue : Volume.IsoValue.toRelative(i.isoValue, stats)
            };
        }

        getDescription() {
            if (this.params.entry.params.view.name === 'selection-box') return 'Selection';
            if (this.params.entry.params.view.name === 'camera-target') return 'Camera';
            if (this.params.entry.params.view.name === 'box') return 'Static Box';
            if (this.params.entry.params.view.name === 'cell') return 'Cell';
            return '';
        }

        constructor(public plugin: PluginContext, public data: VolumeServerInfo.Data) {
            super(plugin, {} as any);

            this.infoMap = new Map<string, VolumeServerInfo.EntryData>();
            this.data.entries.forEach(info => this.infoMap.set(info.dataId, info));
            this.updateQueue = new SingleAsyncQueue();
        }
    }
}
