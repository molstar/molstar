/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementIndex, Model, CustomPropertyDescriptor } from '../../mol-model/structure';
import { StructureElement } from '../../mol-model/structure/structure';
import { Location } from '../../mol-model/location';
import { ThemeDataContext } from '../../mol-theme/theme';
import { ColorTheme, LocationColor } from '../../mol-theme/color';
import { Color } from '../../mol-util/color';
import { Loci } from '../../mol-model/loci';
import { OrderedSet } from '../../mol-data/int';
import { CustomModelProperty } from './custom-model-property';
import { CustomProperty } from './custom-property';
import { LociLabelProvider } from '../../mol-plugin/util/loci-label-manager';

export { CustomElementProperty };

interface CustomElementProperty<T> {
    propertyProvider: CustomModelProperty.Provider<{}, CustomElementProperty.Value<T>>
    colorThemeProvider?: ColorTheme.Provider<{}>
    labelProvider?: LociLabelProvider
}

namespace CustomElementProperty {
    export type Value<T> = Map<ElementIndex, T>

    export interface Builder<T> {
        label: string
        name: string
        getData(model: Model, ctx?: CustomProperty.Context): Value<T> | Promise<Value<T>>
        coloring?: {
            getColor: (p: T) => Color
            defaultColor: Color
        }
        getLabel?: (p: T) => string | undefined
        isApplicable?: (data: Model) => boolean,
        type?: 'dynamic' | 'static',
    }

    export function create<T>(builder: Builder<T>): CustomElementProperty<T> {
        const modelProperty = createModelProperty(builder)
        return {
            propertyProvider: modelProperty,
            colorThemeProvider: builder.coloring?.getColor && createColorThemeProvider(modelProperty, builder.coloring.getColor, builder.coloring.defaultColor),
            labelProvider: builder.getLabel && createLabelProvider(modelProperty, builder.getLabel)
        }
    }

    function createModelProperty<T>(builder: Builder<T>) {
        return CustomModelProperty.createProvider({
            label: builder.label,
            descriptor: CustomPropertyDescriptor({
                name: builder.name,
            }),
            type: builder.type || 'dynamic',
            defaultParams: {},
            getParams: (data: Model) => ({}),
            isApplicable: (data: Model) => !!builder.isApplicable?.(data),
            obtain: async (ctx: CustomProperty.Context, data: Model) => {
                return await builder.getData(data, ctx)
            }
        })
    }

    function createColorThemeProvider<T>(modelProperty: CustomModelProperty.Provider<{}, Value<T>>, getColor: (p: T) => Color, defaultColor: Color) {

        function Coloring(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
            let color: LocationColor;

            const property = ctx.structure && modelProperty.get(ctx.structure.models[0])
            const contextHash = property?.version

            if (property?.value && ctx.structure) {
                const data = property.value
                color = (location: Location) => {
                    if (StructureElement.Location.is(location)) {
                        const e = data.get(location.element);
                        if (typeof e !== 'undefined') return getColor(e);
                    }
                    return defaultColor;
                }
            } else {
                color = () => defaultColor;
            }

            return {
                factory: Coloring,
                granularity: 'group',
                color: color,
                props: props,
                contextHash,
                description: `Assign element colors based on '${modelProperty.label}' data.`
            };
        }

        return {
            label: modelProperty.label,
            factory: Coloring,
            getParams: () => ({}),
            defaultValues: {},
            isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && !!modelProperty.get(ctx.structure.models[0]).value
        }
    }

    function createLabelProvider<T>(modelProperty: CustomModelProperty.Provider<{}, Value<T>>, getLabel: (p: T) => string | undefined) {
        return function(loci: Loci): string | undefined {
            if (loci.kind === 'element-loci') {
                const e = loci.elements[0];
                if (!e) return
                const data = modelProperty.get(e.unit.model).value
                const element = e.unit.elements[OrderedSet.start(e.indices)]
                const value = data?.get(element)
                if (value === undefined) return
                return getLabel(value);
            }
            return
        }
    }
}