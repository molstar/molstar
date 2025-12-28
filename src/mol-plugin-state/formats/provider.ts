/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { EncodedFile } from '../../mol-io/common/binary-cif';
import { decodeMsgPack } from '../../mol-io/common/msgpack/decode';
import { StringLike } from '../../mol-io/common/string-like';
import { PluginContext } from '../../mol-plugin/context';
import { StateObjectRef } from '../../mol-state';
import { FileNameInfo } from '../../mol-util/file-info';
import { PluginStateObject } from '../objects';


export interface DataFormatProvider<P = any, R = any, V = any> {
    label: string,
    description: string,
    category?: string,
    stringExtensions?: string[],
    binaryExtensions?: string[],
    isApplicable?(info: FileNameInfo, data: StringLike | Uint8Array): boolean,
    parse(plugin: PluginContext, data: StateObjectRef<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, params?: P): Promise<R>,
    visuals?(plugin: PluginContext, data: R): Promise<V> | undefined
}

export function DataFormatProvider<P extends DataFormatProvider>(provider: P): P { return provider; }

type CifVariants = 'dscif' | 'segcif' | 'coreCif' | -1
export function guessCifVariant(info: FileNameInfo, data: Uint8Array | StringLike): CifVariants {
    if (info.ext === 'bcif') {
        try {
            // TODO: find a way to run msgpackDecode only once
            //       now it is run twice, here and during file parsing
            const file = decodeMsgPack(data as Uint8Array) as EncodedFile;
            if (file.encoder.startsWith('VolumeServer')) return 'dscif';
            // Assumes volseg-volume-server only serves segments
            if (file.encoder.startsWith('volseg-volume-server')) return 'segcif';

            if (bcifHasCategory(file, 'volume_data_3d_info')) {
                if (bcifHasCategory(file, 'volume_data_3d')) return 'dscif';
                if (bcifHasCategory(file, 'segmentation_data_3d')) return 'segcif';
            }
        } catch (e) {
            console.error(e);
        }
    } else if (info.ext === 'cif') {
        const str = data as StringLike;
        if (str.startsWith('data_SERVER\n#\n_density_server_result')) return 'dscif';
        if (str.startsWith('data_SERVER\n#\ndata_SEGMENTATION_DATA')) return 'segcif';

        if (cifHasCategory(str, 'volume_data_3d_info')) {
            if (cifHasCategory(str, 'volume_data_3d')) return 'dscif';
            if (cifHasCategory(str, 'segmentation_data_3d')) return 'segcif';
        }

        if (str.includes('atom_site_fract_x') || str.includes('atom_site.fract_x')) return 'coreCif';
    }
    return -1;
}

function cifHasCategory(file: StringLike, categoryName: string): boolean {
    return file.includes(`_${categoryName}.`);
}

function bcifHasCategory(file: EncodedFile, categoryName: string): boolean {
    for (const block of file.dataBlocks) {
        for (const category of block.categories) {
            if (category.name === categoryName) return true;
        }
    }
    return false;
}