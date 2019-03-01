/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from 'mol-plugin/behavior';
import { PluginStateObject } from 'mol-plugin/state/objects';

export class VolumeStreaming extends PluginStateObject.CreateBehavior<VolumeStreaming.Behavior>({ name: 'Volume Streaming' }) { }

export namespace VolumeStreaming {

    export interface BehaviorParams {
        
    }

    export class Behavior implements PluginBehavior<{}> {
        register(ref: string): void {
            throw new Error('Method not implemented.');
        }

        unregister(): void {
            throw new Error('Method not implemented.');
        }
    }
}