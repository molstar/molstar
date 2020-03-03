/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { State, StateObjectCell } from '../../../mol-state';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateObject } from '../../objects';
import { FileInfo } from '../../../mol-util/file-info';

export { DataFormatProvider }

type DataFormatProvider<Id extends string = string, D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String = PluginStateObject.Data.Binary | PluginStateObject.Data.String, P = any, S = {}> =
    DataFormatProvider.Base<Id, P> & DataFormatProvider.Definition<D, P, S>

function DataFormatProvider<Id extends string, P = {}>(provider: DataFormatProvider.Base<Id, P>) {
    return function<D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String, S>(p: DataFormatProvider.Definition<D, P, S>): DataFormatProvider<Id, D, P, S> {
        return { ...provider, ...p };
    };
}

namespace DataFormatProvider {
    export type Ret<P extends DataFormatProvider> = P extends DataFormatProvider<any, any, any, infer T> ? T : never
    export type Data<P extends DataFormatProvider> = P extends DataFormatProvider<any, infer T, any, any> ? T : never
    export type Params<P extends DataFormatProvider> = P extends DataFormatProvider<any, any, infer T, any> ? T : never

    export interface Base<Id extends string = string, P = any> {
        id: Id,
        display: { name: string, group?: string, description?: string },
        extensions: { text?: string[], binary?: string[] },
        params?(plugin: PluginContext, data: string | Uint8Array, info?: FileInfo): PD.Def<P>
        // TODO: default representation
        // defaultRepresentation?: RepresenatationProvider
    }
    export interface Definition<D extends PluginStateObject.Data.Binary | PluginStateObject.Data.String = PluginStateObject.Data.Binary | PluginStateObject.Data.String, P = any, S = any> {
        isApplicable?(plugin: PluginContext, data: string | Uint8Array, info?: FileInfo): boolean,
        apply(ctx: { state: State, plugin: PluginContext }, data: StateObjectCell<D>, params: P): Promise<S>
    }
}

// interface DataSourceProvider<Id extends string = string, Format extends DataFormatProvider = DataFormatProvider, P = any, S = {}> {
//     id: Id,
//     display: { name: string, group?: string, description?: string },
//     format: Format,
//     apply(ctx: { ctx: RuntimeContext, state: State, plugin: PluginContext }, params: P): Promise<S>,
//     params(plugin: PluginContext): PD.Def<P>,
// }
// function DataSourceProvider<Id extends string, Format extends DataFormatProvider, P, S>(provider: DataSourceProvider<Id, Format, P, S>) { return provider; }
