/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from '../../mol-task';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueBox } from '../../mol-util';
import { OrderedMap } from 'immutable';
import { AssetManager, Asset } from '../../mol-util/assets';

export { CustomProperty };

namespace CustomProperty {
    export interface Context {
        runtime: RuntimeContext
        assetManager: AssetManager
    }

    export type Data<V> = { value: V, assets?: Asset.Wrapper[] }

    export interface Container<P, V> {
        readonly props: P
        readonly data: ValueBox<V | undefined>
    }

    export interface Provider<Data, Params extends PD.Params, Value> {
        readonly label: string
        readonly descriptor: CustomPropertyDescriptor
        /** hides property in ui and always attaches */
        readonly isHidden?: boolean
        readonly getParams: (data: Data) => Params
        readonly defaultParams: Params
        readonly isApplicable: (data: Data) => boolean
        readonly attach: (ctx: Context, data: Data, props?: Partial<PD.Values<Params>>, addRef?: boolean) => Promise<void>
        readonly ref: (data: Data, add: boolean) => void
        readonly get: (data: Data) => ValueBox<Value | undefined>
        readonly set: (data: Data, props: PD.Values<Params>, value?: Value) => void
        readonly props: (data: Data) => PD.Values<Params>
    }

    export class Registry<Data> {
        private providers = OrderedMap<string, Provider<Data, any, any>>().asMutable()
        private defaultAutoAttachValues = new Map<string, boolean>()

        /** Get params for all applicable property providers */
        getParams(data?: Data) {
            const propertiesParams: PD.Params = {};
            const autoAttachOptions: [string, string][] = [];
            const autoAttachDefault: string[] = [];
            if (data) {
                const values = this.providers.values();
                while (true) {
                    const v = values.next();
                    if (v.done) break;

                    const provider = v.value;
                    if (!provider.isApplicable(data)) continue;

                    if (!provider.isHidden) {
                        autoAttachOptions.push([provider.descriptor.name, provider.label]);
                        if (this.defaultAutoAttachValues.get(provider.descriptor.name)) {
                            autoAttachDefault.push(provider.descriptor.name);
                        }
                    }

                    propertiesParams[provider.descriptor.name] = PD.Group({
                        ...provider.getParams(data)
                    }, { label: provider.label, isHidden: provider.isHidden });
                }
            }
            return {
                autoAttach: PD.MultiSelect(autoAttachDefault, autoAttachOptions),
                properties: PD.Group(propertiesParams, { isFlat: true })
            };
        }

        setDefaultAutoAttach(name: string, value: boolean) {
            this.defaultAutoAttachValues.set(name, value);
        }

        get(name: string) {
            const prop = this.providers.get(name);
            if (!prop) {
                throw new Error(`Custom property '${name}' is not registered.`);
            }
            return this.providers.get(name);
        }

        register(provider: Provider<Data, any, any>, defaultAutoAttach: boolean) {
            this.providers.set(provider.descriptor.name, provider);
            this.defaultAutoAttachValues.set(provider.descriptor.name, defaultAutoAttach);
        }

        unregister(name: string) {
            this.providers.delete(name);
            this.defaultAutoAttachValues.delete(name);
        }
    }
}