/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from '../../mol-task';
import { CustomPropertyDescriptor } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueBox } from '../../mol-util';
import { OrderedMap } from 'immutable';

type AjaxTask = import('../../mol-util/data-source').AjaxTask

export { CustomProperty }

namespace CustomProperty {
    export interface Context {
        runtime: RuntimeContext
        fetch: AjaxTask
    }

    export interface Container<P, V> {
        readonly props: P
        readonly data: ValueBox<V | undefined>
    }

    export interface Provider<Data, Params extends PD.Params, Value> {
        readonly label: string
        readonly descriptor: CustomPropertyDescriptor
        readonly getParams: (data: Data) => Params
        readonly isApplicable: (data: Data) => boolean
        readonly attach: (ctx: Context, data: Data, props?: Partial<PD.Values<Params>>) => Promise<void>
        readonly get: (data: Data) => ValueBox<Value | undefined>
        readonly set: (data: Data, props: PD.Values<Params>, value?: Value) => void
    }

    export class Registry<Data> {
        private providers = OrderedMap<string, Provider<Data, any, any>>().asMutable()
        private defaultAutoAttachValues = new Map<string, boolean>()

        /** Get params for all applicable property providers */
        getParams(data?: Data) {
            const params: PD.Params = {}
            if (data) {
            const values = this.providers.values();
                while (true) {
                    const v = values.next()
                    if (v.done) break
                    const provider = v.value
                    if (!provider.isApplicable(data)) continue
                    params[provider.descriptor.name] = PD.Group({
                        autoAttach: PD.Boolean(this.defaultAutoAttachValues.get(provider.descriptor.name)!),
                        ...provider.getParams(data),
                    }, { label: v.value.label })
                }
            }
            return params
        }

        setDefaultAutoAttach(name: string, value: boolean) {
            this.defaultAutoAttachValues.set(name, value)
        }

        get(name: string) {
            const prop = this.providers.get(name);
            if (!prop) throw new Error(`Custom property '${name}' is not registered.`)
            return this.providers.get(name)
        }

        register(provider: Provider<Data, any, any>, defaultAutoAttach: boolean) {
            this.providers.set(provider.descriptor.name, provider)
            this.defaultAutoAttachValues.set(provider.descriptor.name, defaultAutoAttach)
        }

        unregister(name: string) {
            this.providers.delete(name)
            this.defaultAutoAttachValues.delete(name)
        }
    }
}