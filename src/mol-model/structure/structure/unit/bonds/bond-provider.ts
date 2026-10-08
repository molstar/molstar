/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { CustomPropertyDescriptor } from '../../../../custom-property';
import type { Model } from '../../../model/model';
import { ParamDefinition as PD } from '../../../../../mol-util/param-definition';
import type { Structure } from '../../structure';
import type { Unit } from '../../unit';
import type { IntraUnitBonds } from './data';

/**
 * A Structure-scoped source of atomic intra-unit bonds, fixed when a Unit is created.
 *
 * The selected provider owns the full computation. It can delegate to the
 * exported default computation, replace it, or build a graph from scratch.
 */
export interface BondProvider {
    getBonds(unit: Unit.Atomic): IntraUnitBonds | undefined
}

export namespace BondProvider {
    export interface Context {
        structure?: Structure
    }

    export interface Provider<P extends PD.Params = any> {
        readonly name: string
        readonly label: string
        readonly getParams: (model: Model) => P
        readonly isApplicable: (model: Model) => boolean
        readonly factory: (model: Model, props: PD.Values<P>, context: Context) => BondProvider | undefined
    }
}

export class BondProviderRegistry {
    private readonly providers = new Map<string, BondProvider.Provider>();

    get list(): ReadonlyArray<BondProvider.Provider> {
        return Array.from(this.providers.values());
    }

    add(provider: BondProvider.Provider) {
        if (this.providers.has(provider.name)) {
            throw new Error(`Bond provider '${provider.name}' is already registered.`);
        }
        this.providers.set(provider.name, provider);
    }

    remove(provider: BondProvider.Provider) {
        if (this.providers.get(provider.name) !== provider) return false;
        return this.providers.delete(provider.name);
    }

    get<P extends PD.Params = any>(name: string): BondProvider.Provider<P> | undefined {
        return this.providers.get(name) as BondProvider.Provider<P> | undefined;
    }

    has(name: string): boolean {
        return this.providers.has(name);
    }

    getApplicable(model: Model): ReadonlyArray<BondProvider.Provider> {
        return this.list.filter(provider => provider.isApplicable(model));
    }
}

export namespace ModelBondProvider {
    export const Descriptor = CustomPropertyDescriptor({ name: 'molstar_model_bond_provider' });

    export type Props = PD.NamedParams

    interface Container {
        readonly data: {
            readonly value: Props | undefined
        }
    }

    function getContainer(model: Model): Container | undefined {
        if (!model.customProperties.hasReference(Descriptor)) return undefined;
        return model._dynamicPropertyData[Descriptor.name] as Container | undefined;
    }

    export function get(model: Model): Props | undefined {
        return getContainer(model)?.data.value;
    }
}
