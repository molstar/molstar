/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer, StateTransform } from '../../mol-state';
import { PluginContext } from '../../mol-plugin/context';
import { Download, ReadFile, DownloadBlob, RawData } from '../transforms/data';
import { getFileInfo } from '../../mol-util/file-info';

export class DataBuilder {
    private get dataState() {
        return this.plugin.state.data;
    }

    rawData(params: StateTransformer.Params<RawData>, options?: Partial<StateTransform.Options>) {
        const data = this.dataState.build().toRoot().apply(RawData, params, options);
        return data.commit({ revertOnError: true });
    }

    download(params: StateTransformer.Params<Download>, options?: Partial<StateTransform.Options>) {
        const data = this.dataState.build().toRoot().apply(Download, params, options);
        return data.commit({ revertOnError: true });
    }

    downloadBlob(params: StateTransformer.Params<DownloadBlob>, options?: Partial<StateTransform.Options>) {
        const data = this.dataState.build().toRoot().apply(DownloadBlob, params, options);
        return data.commit({ revertOnError: true });
    }

    async readFile(params: StateTransformer.Params<ReadFile>, options?: Partial<StateTransform.Options>) {
        const data = await this.dataState.build().toRoot().apply(ReadFile, params, options).commit({ revertOnError: true });
        const fileInfo = getFileInfo(params.file?.file || '');
        return { data: data, fileInfo };
    }

    constructor(public plugin: PluginContext) {
    }
}