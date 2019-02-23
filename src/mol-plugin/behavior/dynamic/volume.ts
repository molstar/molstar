/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import CIF from 'mol-io/reader/cif';
import { Box3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import { volumeFromDensityServerData } from 'mol-model-formats/volume/density-server';
import { VolumeData, VolumeIsoValue } from 'mol-model/volume';
import { PluginContext } from 'mol-plugin/context';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { createIsoValueParam } from 'mol-repr/volume/isosurface';
import { Color } from 'mol-util/color';
import { LRUCache } from 'mol-util/lru-cache';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginBehavior } from '../behavior';

export namespace VolumeStreaming {
    function channelParam(label: string, color: Color, defaultValue: number) {
        return PD.Group({
            color: PD.Color(color),
            isoValue: createIsoValueParam(VolumeIsoValue.relative(defaultValue))
        }, { label });
    }

    export const Params = {
        id: PD.Text('1tqn'),
        levels: PD.MappedStatic('x-ray', {
            'em': channelParam('EM', Color(0x638F8F), 1.5),
            'x-ray': PD.Group({
                '2fo-fc': channelParam('2Fo-Fc', Color(0x3362B2), 1.5),
                'fo-fc(+ve)': channelParam('Fo-Fc(+ve)', Color(0x33BB33), 3),
                'fo-fc(-ve)': channelParam('Fo-Fc(-ve)', Color(0xBB3333), -3),
            })
        }),
        box: PD.MappedStatic('static-box', {
            'static-box': PD.Group({
                bottomLeft: PD.Vec3(Vec3.create(-22.4, -33.4, -21.6)),
                topRight: PD.Vec3(Vec3.create(-7.1, -10, -0.9))
            }, { description: 'Static box defined by cartesian coords.' }),
            // 'around-selection': PD.Group({ radius: PD.Numeric(5, { min: 0, max: 10 }) }),
            'cell': PD.Group({  }),
            // 'auto': PD.Group({  }), // based on camera distance/active selection/whatever, show whole structure or slice.
        }),
        detailLevel: PD.Numeric(3, { min: 0, max: 7 }),
        serverUrl: PD.Text('https://webchem.ncbr.muni.cz/DensityServer'),
    }
    export type Params = PD.Values<typeof Params>

    export type ChannelData = { [name in 'EM' | '2FO-FC' | 'FO-FC']?: VolumeData }
    export type LevelType = 'em' | '2fo-fc' | 'fo-fc(+ve)' | 'fo-fc(-ve)'

    export class Behavior implements PluginBehavior<Params> {
        // TODO: have special value for "cell"?
        private cache = LRUCache.create<ChannelData>(25);
        // private ref: string = '';

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
                const volume = await this.ctx.runTask(await volumeFromDensityServerData(densityServerCif));
                (ret as any)[block.header as any] = volume;
            }
            return ret;
        }

        register(ref: string): void {
            // TODO: register camera movement/loci so that "around selection box works"
            // alternatively, and maybe a better solution, write a global behavior that modifies this node from the outside
            // this.ref = ref;
            this.update(this.params);
        }

        async update(params: Params): Promise<boolean> {
            this.params = params;

            let box: Box3D | undefined = void 0;

            switch (params.box.name) {
                case 'static-box':
                    box = Box3D.create(params.box.params.bottomLeft, params.box.params.topRight);
                    break;
                case 'cell':
                    box = this.params.levels.name === 'x-ray'
                        ? void 0 // TODO get bounding box of the model (how to solve assemblies)
                        : void 0;
                    break;
            }

            const data = await this.queryData(box);
            this.currentData = data || { };

            return true;
        }

        unregister(): void {
            // TODO unsubscribe to events
        }

        constructor(public ctx: PluginContext, public params: Params) {
        }
    }

    export class Obj extends PluginStateObject.CreateBehavior<Behavior>({ name: 'Volume Streaming' }) { }
}