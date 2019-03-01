/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from 'mol-plugin/behavior';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { VolumeServerInfo, VolumeServerHeader } from './model';
import { createIsoValueParam } from 'mol-repr/volume/isosurface';
import { VolumeIsoValue, VolumeData } from 'mol-model/volume';
import { Color } from 'mol-util/color';
import { Vec3 } from 'mol-math/linear-algebra';
import { PluginContext } from 'mol-plugin/context';
import { LRUCache } from 'mol-util/lru-cache';
import CIF from 'mol-io/reader/cif';
import { Box3D } from 'mol-math/geometry';
import { urlCombine } from 'mol-util/url';
import { volumeFromDensityServerData } from 'mol-model-formats/volume/density-server';

export class VolumeStreaming extends PluginStateObject.CreateBehavior<VolumeStreaming.Behavior>({ name: 'Volume Streaming' }) { }

export namespace VolumeStreaming {
    function channelParam(label: string, color: Color, defaultValue: VolumeIsoValue, stats: VolumeData['dataStats']) {
        return PD.Group({
            isoValue: createIsoValueParam(defaultValue, stats),
            color: PD.Color(color),
            opacity: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 })
        }, { label, isExpanded: true });
    }

    const fakeSampling: VolumeServerHeader.Sampling = {
        byteOffset: 0,
        rate: 1,
        sampleCount: [1, 1, 1],
        valuesInfo: [{ mean: 0, min: -1, max: 1, sigma: 0.1 }, { mean: 0, min: -1, max: 1, sigma: 0.1 }]
    };

    export function createParams(data?: VolumeServerInfo.Data) {
        // fake the info
        const info = data || { kind: 'em', header: { sampling: [fakeSampling], availablePrecisions: [{ precision: 0, maxVoxels: 0 }] }, emDefaultContourLevel: VolumeIsoValue.relative(0) };

        return {
            view: PD.MappedStatic('box', {
                'box': PD.Group({
                    bottomLeft: PD.Vec3(Vec3.create(-22.4, -33.4, -21.6)),
                    topRight: PD.Vec3(Vec3.create(-7.1, -10, -0.9)),
                    autoUpdate: PD.Boolean(true, { description: 'Update the box when user clicks an element.' })
                }, { description: 'Static box defined by cartesian coords.', isFlat: true }),
                'cell': PD.Group({}),
                // 'auto': PD.Group({  }), // based on camera distance/active selection/whatever, show whole structure or slice.
            }, { options: [['box', 'Bounded Box'], ['cell', 'Whole Structure']] }),
            detailLevel: PD.Select<number>(Math.min(1, info.header.availablePrecisions.length - 1),
                info.header.availablePrecisions.map((p, i) => [i, `${i + 1} (${Math.pow(p.maxVoxels, 1 / 3) | 0}^3)`] as [number, string])),
            channels: info.kind === 'em'
                ? PD.Group({
                    'em': channelParam('EM', Color(0x638F8F), info.emDefaultContourLevel || VolumeIsoValue.relative(1), info.header.sampling[0].valuesInfo[0])
                }, { isFlat: true })
                : PD.Group({
                    '2fo-fc': channelParam('2Fo-Fc', Color(0x3362B2), VolumeIsoValue.relative(1.5), info.header.sampling[0].valuesInfo[0]),
                    'fo-fc(+ve)': channelParam('Fo-Fc(+ve)', Color(0x33BB33), VolumeIsoValue.relative(3), info.header.sampling[0].valuesInfo[1]),
                    'fo-fc(-ve)': channelParam('Fo-Fc(-ve)', Color(0xBB3333), VolumeIsoValue.relative(-3), info.header.sampling[0].valuesInfo[1]),
                }, { isFlat: true })
        };
    }

    type RT = typeof createParams extends (...args: any[]) => (infer T) ? T : never
    export type Params = RT extends PD.Params ? PD.Values<RT> : {}

    type ChannelsInfo = { [name in ChannelType]?: { isoValue: VolumeIsoValue, color: Color, opacity: number } }
    type ChannelsData = { [name in 'EM' | '2FO-FC' | 'FO-FC']?: VolumeData }

    export type ChannelType = 'em' | '2fo-fc' | 'fo-fc(+ve)' | 'fo-fc(-ve)'
    export const ChannelTypeOptions: [ChannelType, string][] = [['em', 'em'], ['2fo-fc', '2fo-fc'], ['fo-fc(+ve)', 'fo-fc(+ve)'], ['fo-fc(-ve)', 'fo-fc(-ve)']]
    interface ChannelInfo {
        data: VolumeData,
        color: Color,
        isoValue: VolumeIsoValue.Relative,
        opacity: number
    }
    export type Channels = { [name in ChannelType]?: ChannelInfo }

    export class Behavior implements PluginBehavior<{}> {
        private cache = LRUCache.create<ChannelsData>(25);
        private params: Params = {} as any;
        private ref: string = '';

        channels: Channels = {}

        private async queryData(box?: Box3D) {
            let url = urlCombine(this.info.serverUrl, `${this.info.kind}/${this.info.dataId}`);

            if (box) {
                const { min: a, max: b } = box;
                url += `/box`
                    + `/${a.map(v => Math.round(1000 * v) / 1000).join(',')}`
                    + `/${b.map(v => Math.round(1000 * v) / 1000).join(',')}`;
            } else {
                url += `/cell`;
            }
            url += `?detail=${this.params.detailLevel}`;

            let data = LRUCache.get(this.cache, url);
            if (data) {
                return data;
            }

            const cif = await this.plugin.runTask(this.plugin.fetch({ url, type: 'binary' }));
            data = await this.parseCif(cif as Uint8Array);
            if (!data) {
                return;
            }

            LRUCache.set(this.cache, url, data);
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
                const volume = await this.plugin.runTask(await volumeFromDensityServerData(densityServerCif));
                (ret as any)[block.header as any] = volume;
            }
            return ret;
        }

        register(ref: string): void {
            this.ref = ref;
        }

        async update(params: Params) {
            this.params = params;

            let box: Box3D | undefined = void 0;

            switch (params.view.name) {
                case 'box':
                    box = Box3D.create(params.view.params.bottomLeft, params.view.params.topRight);
                    break;
                case 'cell':
                    box = this.info.kind === 'x-ray'
                        ? this.info.structure.boundary.box
                        : void 0;
                    break;
            }

            const data = await this.queryData(box);

            if (!data) return false;

            const info = params.channels as ChannelsInfo;

            if (this.info.kind === 'x-ray') {
                this.channels['2fo-fc'] = this.createChannel(data['2FO-FC']!, info['2fo-fc'], this.info.header.sampling[0].valuesInfo[0]);
                this.channels['fo-fc(+ve)'] = this.createChannel(data['FO-FC']!, info['fo-fc(+ve)'], this.info.header.sampling[0].valuesInfo[1]);
                this.channels['fo-fc(-ve)'] = this.createChannel(data['FO-FC']!, info['fo-fc(-ve)'], this.info.header.sampling[0].valuesInfo[1]);
            } else {
                this.channels['em'] = this.createChannel(data['EM']!, info['em'], this.info.header.sampling[0].valuesInfo[0]);
            }

            return true;
        }

        private createChannel(data: VolumeData, info: ChannelsInfo['em'], stats: VolumeData['dataStats']): ChannelInfo {
            const i = info!;
            return {
                data,
                color: i.color,
                opacity: i.opacity,
                isoValue: i.isoValue.kind === 'relative' ? i.isoValue : VolumeIsoValue.toRelative(i.isoValue, stats)
            };
        }

        unregister(): void {
            // throw new Error('Method not implemented.');
        }

        constructor(public plugin: PluginContext, public info: VolumeServerInfo.Data) {

        }
    }
}