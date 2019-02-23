/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import CIF from 'mol-io/reader/cif';
import { Box3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import { volumeFromDensityServerData } from 'mol-model-formats/volume/density-server';
import { VolumeData } from 'mol-model/volume';
import { PluginContext } from 'mol-plugin/context';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { IsoValueParam } from 'mol-repr/volume/isosurface';
import { Color } from 'mol-util/color';
import { ColorNames } from 'mol-util/color/tables';
import { LRUCache } from 'mol-util/lru-cache';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginBehavior } from '../behavior';

export namespace VolumeStreamingBehavior {
    function channelParam(label: string) { return PD.Group({ color: PD.Color(Color(ColorNames.teal)), isoValue: IsoValueParam }, { label }); }

    export const Params = {
        id: PD.Text(''),
        levels: PD.MappedStatic('em', {
            'em': channelParam('EM'),
            'x-ray': PD.Group({
                '2fo-fc': channelParam('2Fo-Fc'),
                'fo-fc(+ve)': channelParam('Fo-Fc(+ve)'),
                'fo-fc(-ve)': channelParam('Fo-Fc(-ve)'),
            })
        }),
        box: PD.MappedStatic('static-box', {
            'static-box': PD.Group({
                bottomLeft: PD.Vec3(Vec3.create(0, 0, 0)),
                topRight: PD.Vec3(Vec3.create(1, 1, 1))
            }, { description: 'Static box defined by cartesian coords.' }),
            // 'around-selection': PD.Group({ radius: PD.Numeric(5, { min: 0, max: 10 }) }),
        }),
        detailLevel: PD.Numeric(3, { min: 0, max: 7 }),
        serverUrl: PD.Text('https://webchem.ncbr.muni.cz/DensityServer'),
    }
    export type Params = PD.Values<typeof Params>

    export type ChannelData = { [name in 'EM' | '2FO-FC' | 'FO-FC']?: VolumeData }
    export type LevelType = 'em' | '2fo-fc' | 'fo-fc(+ve)' | 'fo-fc(-ve)'

    export class Behavior implements PluginBehavior<Params> {
        private cache = LRUCache.create<ChannelData>(25);

        currentData: ChannelData = { }

        private async queryData(box?: Box3D) {
            let url = `${this.params.serverUrl}/${this.params.levels.name}/${this.params.id}`

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

            const cif = await this.ctx.runTask(this.ctx.fetch(url, 'binary'));
            data = await this.parseCif(cif as Uint8Array);
            if (!data) {
                return;
            }

            LRUCache.set(this.cache, url, data);
            return data;
        }

        private async parseCif(data: Uint8Array): Promise<ChannelData | undefined> {
            const parsed = await this.ctx.runTask(CIF.parseBinary(data));
            if (parsed.isError) {
                this.ctx.log.error('VolumeStreaming, parsing CIF: ' + parsed.toString());
                return;
            }
            if (parsed.result.blocks.length < 2) {
                this.ctx.log.error('VolumeStreaming: Invalid data.');
                return;
            }

            const ret: ChannelData = { };
            for (let i = 1; i < parsed.result.blocks.length; i++) {
                const block = parsed.result.blocks[i];

                const densityServerCif = CIF.schema.densityServer(block);
                const volume = this.ctx.runTask(await volumeFromDensityServerData(densityServerCif));
                (ret as any)[block.header as any] = volume;
            }
            return ret;
        }

        register(ref: string): void {
            this.update(this.params);
        }

        async update(params: Params): Promise<boolean> {
            this.params = params;

            const box: Box3D = Box3D.create(params.box.params.bottomLeft, params.box.params.topRight);
            const data = await this.queryData(box);
            this.currentData = data || { };

            return true;
        }

        unregister(): void {
        }

        constructor(public ctx: PluginContext, public params: Params) {
        }
    }

    export class Obj extends PluginStateObject.CreateBehavior<Behavior>({ name: 'Volume Streaming' }) { }
}