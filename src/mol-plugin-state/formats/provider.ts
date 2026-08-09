/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { binaryCifHasCategory, binaryCifHasColumn, getBinaryCifHeader } from '../../mol-io/common/binary-cif';
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
    /**
     * Controls the order in which `DataFormatRegistry.auto` tries providers that share the
     * same extension, independent of registration order. Higher values are tried first.
     * Defaults to 0. Use a positive value for providers with a restrictive `isApplicable`
     * that should take precedence over more generic fallback providers (e.g. a specialized
     * variant of a shared file extension).
     */
    priority?: number,
    isApplicable?(info: FileNameInfo, data: StringLike | Uint8Array): boolean,
    parse(plugin: PluginContext, data: StateObjectRef<PluginStateObject.Data.Binary | PluginStateObject.Data.String>, params?: P): Promise<R>,
    visuals?(plugin: PluginContext, data: R): Promise<V> | undefined
}

export function DataFormatProvider<P extends DataFormatProvider>(provider: P): P { return provider; }

type CifVariants = 'dscif' | 'segcif' | 'sfcif' | 'coreCif' | -1
export function guessCifVariant(info: FileNameInfo, data: Uint8Array | StringLike): CifVariants {
    if (info.ext === 'bcif') {
        try {
            const header = getBinaryCifHeader(data as Uint8Array);
            if (header.encoder.startsWith('VolumeServer')) return 'dscif';
            // Assumes volseg-volume-server only serves segments
            if (header.encoder.startsWith('volseg-volume-server')) return 'segcif';

            if (binaryCifHasCategory(header, '_volume_data_3d_info')) {
                if (binaryCifHasCategory(header, '_volume_data_3d')) return 'dscif';
                if (binaryCifHasCategory(header, '_segmentation_data_3d')) return 'segcif';
            }
            if (binaryCifHasCategory(header, '_refln')) {
                if (binaryCifHasColumn(header, '_refln', 'pdbx_FWT') || binaryCifHasColumn(header, '_refln', 'pdbx_DELFWT')) {
                    return 'sfcif';
                }
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

        // Structure factor CIF: has _refln category with map coefficients
        if (str.includes('_refln.pdbx_FWT') || str.includes('_refln.pdbx_DELFWT')) return 'sfcif';

        if (str.includes('atom_site_fract_x') || str.includes('atom_site.fract_x')) return 'coreCif';
    }
    return -1;
}

function cifHasCategory(file: StringLike, categoryName: string): boolean {
    return file.includes(`_${categoryName}.`);
}
