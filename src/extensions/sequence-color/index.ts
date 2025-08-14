/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { Loci } from '../../mol-model/loci';
import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ElementSet, Selector, SelectorParams } from '../mvs/components/selector';


/** Allows coloring residues in sequence panel */
export const SequenceColor = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'sequence-color',
    category: 'misc',
    display: {
        name: 'Sequence Color',
        description: 'Sequence Color extension',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private readonly provider = SequenceColorProvider;

        register(): void {
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);
        }
        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }
        unregister() {
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
    })
});


/** Parameter definition for custom structure property "SequenceColor" */
export type SequenceColorParams = typeof SequenceColorParams
export const SequenceColorParams = {
    colors: PD.ObjectList(
        {
            color: PD.Color(ColorNames.grey, { description: 'Color to apply' }),
            selector: SelectorParams,
        },
        obj => Color.toHexStyle(obj.color)
    ),
};

/** Parameter values of custom structure property "SequenceColor" */
export type SequenceColorProps = PD.Values<SequenceColorParams>

/** Values of custom structure property "SequenceColor" (and for its params at the same type) */
export type SequenceColorData = { selector: Selector, color: Color, elementSet?: ElementSet }[]


/** Provider for custom structure property "SequenceColor" */
export const SequenceColorProvider: CustomStructureProperty.Provider<SequenceColorParams, SequenceColorData> = CustomStructureProperty.createProvider({
    label: 'Sequence Color',
    descriptor: CustomPropertyDescriptor<any, any>({
        name: 'sequence-color',
    }),
    type: 'root',
    defaultParams: SequenceColorParams,
    getParams: (data: Structure) => SequenceColorParams,
    isApplicable: (data: Structure) => data.root === data,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<SequenceColorProps>) => {
        const fullProps = { ...PD.getDefaultValues(SequenceColorParams), ...props };
        const value = fullProps.colors.map(t => ({
            selector: t.selector,
            color: t.color,
        } satisfies SequenceColorData[number]));
        return { value: value } satisfies CustomProperty.Data<SequenceColorData>;
    },
});

export function sequenceColorForLocation(location: StructureElement.Location): Color | undefined {
    const colorData = SequenceColorProvider.get(location.structure).value;
    if (!colorData) return undefined;
    for (let i = colorData.length - 1; i >= 0; i--) { // last color matters
        const item = colorData[i];
        const elements = item.elementSet ??= ElementSet.fromSelector(location.structure, item.selector);
        if (ElementSet.has(elements, location)) {
            return item.color;
        }
    }
    return undefined;
}

export function sequenceColorForLoci(loci: Loci): Color | undefined {
    if (!StructureElement.Loci.is(loci)) return undefined;
    const location = StructureElement.Loci.getFirstLocation(loci);
    if (location === undefined) return undefined;
    return sequenceColorForLocation(location);
}
