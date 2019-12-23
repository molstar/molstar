/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer, StateTransform, StateObjectSelector, StateObjectCell } from '../../../mol-state';
import { PluginContext } from '../../context';
import { Download, ReadFile } from '../transforms/data';
import { getFileInfo } from '../../../mol-util/file-info';
import { DataFormatProvider } from './data/provider';
import { BuiltInDataFormats } from './data/formats';
import { objectForEach } from '../../../mol-util/object';

export class DataManager {
    readonly formats: DataFormatProvider[] = [];

    addFormat(p: DataFormatProvider) {
        this.formats.push(p);
    }

    get dataState() {
        return this.plugin.state.dataState;
    }

    async download(params: StateTransformer.Params<Download>, options?: Partial<StateTransform.Options>) {
        const data = this.dataState.build().toRoot().apply(Download, params, options);
        await this.plugin.runTask(this.dataState.updateTree(data));
        return { data: data.selector };
    }

    async readFile(params: StateTransformer.Params<ReadFile>, options?: Partial<StateTransform.Options>) {
        const data = this.dataState.build().toRoot().apply(ReadFile, params, options);
        const fileInfo = getFileInfo(params.file);
        await this.plugin.runTask(this.dataState.updateTree(data));
        return { data: data.selector, fileInfo };
    }

    async parse<K extends keyof BuiltInDataFormats, P extends BuiltInDataFormats[K]>(provider: K, data: StateObjectSelector<DataFormatProvider.Data<P>> | StateObjectCell | string, params?: DataFormatProvider.Params<P>): Promise<DataFormatProvider.Ret<P>>
    async parse<P extends DataFormatProvider>(provider: P, data: StateObjectSelector<DataFormatProvider.Data<P>> | StateObjectCell | string, params?: DataFormatProvider.Params<P>): Promise<DataFormatProvider.Ret<P>>
    async parse<P extends DataFormatProvider>(providerOrBuildIn: P | string, data: StateObjectSelector<DataFormatProvider.Data<P>> | StateObjectCell | StateTransform.Ref, params?: DataFormatProvider.Params<P>) {
        const provider: P = typeof providerOrBuildIn === 'string' ? BuiltInDataFormats[providerOrBuildIn as keyof BuiltInDataFormats] as unknown as P : providerOrBuildIn as P;
        const cell = StateObjectCell.resolve(this.dataState, data);
        if (!cell) {
            throw new Error('Could not resolve data cell.');
        }
        return provider.apply({ state: this.dataState, plugin: this.plugin }, cell, params);
    }

    // async test() {
    //     const { data } = await this.download({ url: '' });
    //     const cif = await this.parse('mmcif', data);
    // }

    constructor(public plugin: PluginContext) {
        objectForEach(BuiltInDataFormats, f => this.formats.push(f));
    }
}