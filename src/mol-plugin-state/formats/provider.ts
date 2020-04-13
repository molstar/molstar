/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import msgpackDecode from '../../mol-io/common/msgpack/decode';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectRef } from '../../mol-state';
import { FileInfo } from '../../mol-util/file-info';
import { PluginStateObject } from '../objects';

export interface DataFormatProvider<P = any, R = any, V = any> {
    label: string,
    description: string,
    category?: string,
    stringExtensions?: string[],
    binaryExtensions?: string[],
    isApplicable?(info: FileInfo, data: string | Uint8Array): boolean,
    parse(plugin: PluginContext, data: StateObjectRef<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, params?: P): Promise<R>,
    visuals?(plugin: PluginContext, data: R): Promise<V> | undefined
}

export function DataFormatProvider<P extends DataFormatProvider>(provider: P): P { return provider; }

type cifVariants = 'dscif' | 'coreCif' | -1
export function guessCifVariant(info: FileInfo, data: Uint8Array | string): cifVariants {
    if (info.ext === 'bcif') {
        try {
            // TODO: find a way to run msgpackDecode only once
            //      now it is run twice, here and during file parsing
            if (msgpackDecode(data as Uint8Array).encoder.startsWith('VolumeServer')) return 'dscif';
        } catch { }
    } else if (info.ext === 'cif') {
        const str = data as string;
        if (str.startsWith('data_SERVER\n#\n_density_server_result')) return 'dscif';
        if (str.includes('atom_site_fract_x') || str.includes('atom_site.fract_x')) return 'coreCif';
    }
    return -1;
}