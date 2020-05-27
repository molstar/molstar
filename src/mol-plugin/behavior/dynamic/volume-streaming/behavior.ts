/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { Volume, Grid } from '../../../../mol-model/volume';
import { createIsoValueParam } from '../../../../mol-repr/volume/isosurface';
import { VolumeServerHeader, VolumeServerInfo } from './model';
import { Box3D } from '../../../../mol-math/geometry';
import { Vec3 } from '../../../../mol-math/linear-algebra';
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
            isoValue: createIsoValueParam(typeof defaults.isoValue !== 'undefined' ? defaults.isoValue : defaultValue, stats),
            color: PD.Color(typeof defaults.color !== 'undefined' ? defaults.color : color),
            wireframe: PD.Boolean(typeof defaults.wireframe !== 'undefined' ? defaults.wireframe : false),
            opacity: PD.Numeric(typeof defaults.opacity !== 'undefined' ? defaults.opacity : 0.3, { min: 0, max: 1, step: 0.01 })
        }, { label, isExpanded: true });
    }

    const fakeSampling: VolumeServerHeader.Sampling = {
        byteOffset: 0,
        rate: 1,
        sampleCount: [1, 1, 1],
        valuesInfo: [{ mean: 0, min: -1, max: 1, sigma: 0.1 }, { mean: 0, min: -1, max: 1, sigma: 0.1 }]
    };

    export function createParams(options: { data?: VolumeServerInfo.Data, defaultView?: ViewTypes, channelParams?: DefaultChannelParams } = { }) {
        const { data, defaultView, channelParams } = options;
        const map = new Map<string, VolumeServerInfo.EntryData>();
        if (data) data.entries.forEach(d => map.set(d.dataId, d));
        const names = data ? data.entries.map(d => [d.dataId, d.dataId] as [string, string]) : [];
        const defaultKey = data ? data.entries[0].dataId : '';
        return {
            entry: PD.Mapped<EntryParams>(defaultKey, names, name => PD.Group(createEntryParams({ entryData: map.get(name)!, defaultView, structure: data && data.structure, channelParams }))),
        };
    }

    export type EntryParamDefinition = typeof createEntryParams extends (...args: any[]) => (infer T) ? T : never
    export type EntryParams = EntryParamDefinition extends PD.Params ? PD.Values<EntryParamDefinition> : {}

    export function createEntryParams(options: { entryData?: VolumeServerInfo.EntryData, defaultView?: ViewTypes, structure?: Structure, channelParams?: DefaultChannelParams }) {
        const { entryData, defaultView, structure, channelParams = { } } = options;

        // fake the info
        const info = entryData || { kind: 'em', header: { sampling: [fakeSampling], availablePrecisions: [{ precision: 0, maxVoxels: 0 }] }, emDefaultContourLevel: Volume.IsoValue.relative(0) };
        const box = (structure && structure.boundary.box) || Box3D.empty();

        return {
            view: PD.MappedStatic(defaultView || (info.kind === 'em' ? 'cell' : 'selection-box'), {
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
                'cell': PD.Group({}),
                // 'auto': PD.Group({  }), // TODO based on camera distance/active selection/whatever, show whole structure or slice.
            }, { options: ViewTypeOptions, description: 'Controls what of the volume is displayed. "Off" hides the volume alltogether. "Bounded box" shows the volume inside the given box. "Around Interaction" shows the volume around the focused element/atom. "Whole Structure" shows the volume for the whole structure.' }),
            detailLevel: PD.Select<number>(Math.min(3, info.header.availablePrecisions.length - 1),
                info.header.availablePrecisions.map((p, i) => [i, `${i + 1} [ ${Math.pow(p.maxVoxels, 1 / 3) | 0}^3 cells ]`] as [number, string]), { description: 'Determines the maximum number of voxels. Depending on the size of the volume options are in the range from 0 (0.52M voxels) to 6 (25.17M voxels).' }),
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

    export const ViewTypeOptions = [['off', 'Off'], ['box', 'Bounded Box'], ['selection-box', 'Around Focus'], ['cell', 'Whole Structure']] as [ViewTypes, string][];

    export type ViewTypes = 'off' | 'box' | 'selection-box' | 'cell'

    export type ParamDefinition = typeof createParams extends (...args: any[]) => (infer T) ? T : never
    export type Params = ParamDefinition extends PD.Params ? PD.Values<ParamDefinition> : {}

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
        public infoMap: Map<string, VolumeServerInfo.EntryData>

        channels: Channels = {}

        public get info () {
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
            url += `?detail=${this.params.entry.params.detailLevel}`;

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

        private updateDynamicBox(box: Box3D) {
            if (this.params.entry.params.view.name !== 'selection-box') return;

            const state = this.plugin.state.data;
            const newParams: Params = {
                ...this.params,
                entry: {
                    name: this.params.entry.name,
                    params: {
                        ...this.params.entry.params,
                        view: {
                            name: 'selection-box' as 'selection-box',
                            params: {
                                radius: this.params.entry.params.view.params.radius,
                                bottomLeft: box.min,
                                topRight: box.max
                            }
                        }
                    }
                }
            };
            const update = state.build().to(this.ref).update(newParams);

            PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotUpdateCurrent: true } });
        }

        private getStructureRoot() {
            return this.plugin.state.data.select(StateSelection.Generators.byRef(this.ref).rootOfType([PluginStateObject.Molecule.Structure]))[0];
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

                if (this.params.entry.params.view.name !== 'selection-box') {
                    this.lastLoci = loci;
                } else {
                    this.updateInteraction(loci);
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

        private getBoxFromLoci(loci: StructureElement.Loci | EmptyLoci): Box3D {
            if (Loci.isEmpty(loci)) {
                return Box3D.empty();
            }

            const parent = this.plugin.helpers.substructureParent.get(loci.structure, true);
            if (!parent) return Box3D.empty();
            const root = this.getStructureRoot();
            if (!root || root.obj?.data !== parent.obj?.data) return Box3D.empty();

            const extendedLoci = StructureElement.Loci.extendToWholeResidues(loci);
            const box = StructureElement.Loci.getBoundary(extendedLoci).box;
            if (StructureElement.Loci.size(extendedLoci) === 1) {
                Box3D.expand(box, box, Vec3.create(1, 1, 1));
            }
            return box;
        }

        private updateInteraction(loci: StructureElement.Loci | EmptyLoci) {
            if (Loci.areEqual(this.lastLoci, loci)) {
                this.lastLoci = EmptyLoci;
                this.updateDynamicBox(Box3D.empty());
                return;
            }

            this.lastLoci = loci;

            if (isEmptyLoci(loci)) {
                this.updateDynamicBox(Box3D.empty());
                return;
            }

            const box = this.getBoxFromLoci(loci);
            this.updateDynamicBox(box);
        }

        async update(params: Params) {
            const switchedToSelection = params.entry.params.view.name === 'selection-box' && this.params && this.params.entry && this.params.entry.params && this.params.entry.params.view && this.params.entry.params.view.name !== 'selection-box';

            this.params = params;

            let box: Box3D | undefined = void 0, emptyData = false;

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
                        box = this.getBoxFromLoci(this.lastLoci) || Box3D.empty();
                    } else {
                        box = Box3D.create(Vec3.clone(params.entry.params.view.params.bottomLeft), Vec3.clone(params.entry.params.view.params.topRight));
                    }
                    const r = params.entry.params.view.params.radius;
                    emptyData = Box3D.volume(box) < 0.0001;
                    Box3D.expand(box, box, Vec3.create(r, r, r));
                    break;
                }
                case 'cell':
                    box = this.info.kind === 'x-ray'
                        ? this.data.structure.boundary.box
                        : void 0;
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
            if (this.params.entry.params.view.name === 'box') return 'Static Box';
            if (this.params.entry.params.view.name === 'cell') return 'Cell';
            return '';
        }

        constructor(public plugin: PluginContext, public data: VolumeServerInfo.Data) {
            super(plugin, {} as any);

            this.infoMap = new Map<string, VolumeServerInfo.EntryData>();
            this.data.entries.forEach(info => this.infoMap.set(info.dataId, info));
        }
    }
}