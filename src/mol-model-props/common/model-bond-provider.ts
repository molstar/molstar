/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Model } from '../../mol-model/structure';
import { BondProviderRegistry, ModelBondProvider } from '../../mol-model/structure/structure/unit/bonds/bond-provider';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CustomModelProperty } from './custom-model-property';
import { CustomProperty } from './custom-property';

const DefaultModelBondProviderParams = {
    provider: PD.Mapped<any>('', [['', 'None']], () => PD.Value<any>({})),
};

export type ModelBondProviderParams = typeof DefaultModelBondProviderParams

function getModelBondProviderParams(registry: BondProviderRegistry, model: Model): ModelBondProviderParams {
    const providers = registry.getApplicable(model);
    if (providers.length === 0) return DefaultModelBondProviderParams;

    const options = providers.map(provider => [provider.name, provider.label] as [string, string]);
    return {
        provider: PD.Mapped<any>(
            options[0][0],
            options,
            name => PD.Group(registry.get(name)?.getParams(model) ?? {})
        ),
    };
}

export function createModelBondProviderProperty(registry: BondProviderRegistry): CustomModelProperty.Provider<ModelBondProviderParams, ModelBondProvider.Props | undefined> {
    return CustomModelProperty.createProvider({
        label: 'Bond Provider',
        descriptor: ModelBondProvider.Descriptor,
        type: 'dynamic',
        defaultParams: DefaultModelBondProviderParams,
        getParams: model => getModelBondProviderParams(registry, model),
        isApplicable: model => registry.getApplicable(model).length > 0,
        obtain: async (_ctx: CustomProperty.Context, model: Model, props: PD.Values<ModelBondProviderParams>) => {
            const provider = registry.get(props.provider.name);
            return {
                value: provider ? props.provider : undefined,
            };
        },
    });
}
