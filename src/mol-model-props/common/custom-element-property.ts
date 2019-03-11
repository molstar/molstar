/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ElementIndex, Model, ModelPropertyDescriptor } from 'mol-model/structure';
import { StructureElement } from 'mol-model/structure/structure';
import { Location } from 'mol-model/location';
import { CustomPropertyRegistry } from './custom-property-registry';
import { Task } from 'mol-task';
import { ThemeDataContext, ThemeProvider } from 'mol-theme/theme';
import { ColorTheme, LocationColor } from 'mol-theme/color';
import { Color } from 'mol-util/color';
import { TableLegend } from 'mol-util/color/tables';
import { Loci } from 'mol-model/loci';
import { OrderedSet } from 'mol-data/int';

export { CustomElementProperty };

namespace CustomElementProperty {
    export interface CreateParams<T> {
        isStatic: boolean,
        name: string,
        autoAttach?: boolean,
        display: string,
        attachableTo?: (model: Model) => boolean,
        getData(model: Model): Map<ElementIndex, T> | Promise<Map<ElementIndex, T>>,
        format?(e: T): string | undefined,
        coloring?: {
            getColor: (e: T) => Color,
            defaultColor: Color
        }
    }

    export function create<T>(params: CreateParams<T>) {
        const name = params.name;

        const Descriptor = ModelPropertyDescriptor({
            isStatic: params.isStatic,
            name: params.name,
        });

        function attach(model: Model) {
            return Task.create(`Attach ${params.display}`, async () => {
                try {
                    if (model.customProperties.has(Descriptor)) return true;

                    const data = await params.getData(model);

                    if (params.isStatic) {
                        model._staticPropertyData[name] = data;
                    } else {
                        model._dynamicPropertyData[name] = data;
                    }

                    model.customProperties.add(Descriptor);

                    return true;
                } catch (e) {
                    console.warn('Attach Property', e);
                    return false;
                }
            })
        }

        function getStatic(e: StructureElement) { return e.unit.model._staticPropertyData[name].get(e.element); }
        function getDynamic(e: StructureElement) { return e.unit.model._staticPropertyData[name].get(e.element); }

        const propertyProvider: CustomPropertyRegistry.Provider = {
            option: [name, params.display],
            descriptor: Descriptor,
            defaultSelected: !!params.autoAttach,
            attachableTo: params.attachableTo || (() => true),
            attach
        };

        const get = params.isStatic ? getStatic : getDynamic;

        function has(model: Model) { return model.customProperties.has(Descriptor); }

        function Coloring(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
            let color: LocationColor;
            const getColor = params.coloring!.getColor;
            const defaultColor = params.coloring!.defaultColor;

            if (ctx.structure && !ctx.structure.isEmpty && has(ctx.structure.models[0])) {
                color = (location: Location) => {
                    if (StructureElement.isLocation(location)) {
                        const e = get(location);
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
                description: 'Assign element colors based on the provided data.',
                legend: TableLegend([])
            };
        }

        const colorTheme: ThemeProvider<ColorTheme<{}>, {}> = {
            label: params.display,
            factory: Coloring,
            getParams: () => ({}),
            defaultValues: {},
            isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && !ctx.structure.isEmpty && has(ctx.structure.models[0])
        }

        function LabelProvider(loci: Loci): string | undefined {
            if (loci.kind === 'element-loci') {
                const e = loci.elements[0];
                if (!has(e.unit.model)) return void 0;
                return params.format!(get(StructureElement.create(e.unit, e.unit.elements[OrderedSet.getAt(e.indices, 0)])));
            }
            return void 0;
        }

        return {
            Descriptor,
            attach,
            get,
            propertyProvider,
            colorTheme: params.coloring ? colorTheme : void 0,
            labelProvider: params.format ? LabelProvider : ((loci: Loci) => void 0)
        };
    }
}