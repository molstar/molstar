/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { CustomStructureProperty } from '../../mol-model-props/common/custom-structure-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { ElementIndex, Structure, StructureElement } from '../../mol-model/structure';
import { CustomStructureProperties } from '../../mol-plugin-state/transforms/model';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ElementSet, Selector, SelectorParams } from '../mvs/components/selector';
import { SequenceColoring } from '../../mol-plugin/util/sequence-coloring';


/** Allows coloring residues in sequence panel */
export const SequenceColor = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'sequence-color',
    category: 'misc',
    display: {
        name: 'Sequence Color',
        description: 'Sequence Color extension, allows assigning custom residue colors to be shown in the sequence panel',
    },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        register(): void {
            this.ctx.customStructureProperties.register(SequenceColorPropertyProvider, this.params.autoAttach);
            this.ctx.customSequenceColoringRegistry.register(SequenceColoringProvider);
        }
        update(p: { autoAttach: boolean }) {
            const updated = this.params.autoAttach !== p.autoAttach;
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(SequenceColorPropertyProvider.descriptor.name, this.params.autoAttach);
            return updated;
        }
        unregister() {
            this.ctx.customSequenceColoringRegistry.unregister(SequenceColoringProvider.name);
            this.ctx.customStructureProperties.unregister(SequenceColorPropertyProvider.descriptor.name);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(true),
    })
});


const SequenceColoringProvider: SequenceColoring.Provider = {
    name: 'custom-struct-prop-sequence-color',
    color: sequenceColorForLoci,
    subcribeForUpdates: (plugin, requestUpdate) => {
        const sub = plugin.state.events.cell.stateUpdated.subscribe(s => {
            if (s.cell.transform.transformer === CustomStructureProperties) {
                requestUpdate(s.cell.obj?.data);
            }
        });
        return sub;
    },
};


/** Parameter definition for custom structure property "SequenceColor" */
export type SequenceColorParams = typeof SequenceColorParams
export const SequenceColorParams = {
    colors: PD.ObjectList(
        {
            color: PD.Color(ColorNames.grey, { description: 'Color to apply to a substructure' }),
            selector: SelectorParams,
        },
        obj => Color.toHexStyle(obj.color),
        { description: 'List of substructure-color assignments' }
    ),
};

/** Parameter values of custom structure property "SequenceColor" */
export type SequenceColorProps = PD.Values<SequenceColorParams>

/** Values of custom structure property "SequenceColor" */
export interface SequenceColorData {
    items: {
        selector: Selector,
        color: Color,
        elementSet?: ElementSet,
    }[],
    colorCache: {
        [unitId: number]: {
            [elemIdx: ElementIndex]: Color | undefined,
        },
    },
}

/** Provider for custom structure property "SequenceColor" */
export const SequenceColorPropertyProvider: CustomStructureProperty.Provider<SequenceColorParams, SequenceColorData> = CustomStructureProperty.createProvider({
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
        const items = fullProps.colors.map(t => ({
            selector: t.selector,
            color: t.color,
        } satisfies SequenceColorData['items'][number]));
        return { value: { items, colorCache: {} } } satisfies CustomProperty.Data<SequenceColorData>;
    },
});

/** Get color assigned to a loci in custom structure property "SequenceColor" */
export function sequenceColorForLoci(loci: StructureElement.Loci): Color | undefined {
    const colorData = SequenceColorPropertyProvider.get(loci.structure).value;
    if (!colorData || colorData.items.length === 0) return undefined;
    const location = StructureElement.Loci.getFirstLocation(loci);
    if (location === undefined) return undefined;
    return sequenceColorForLocation(colorData, location);
}

function sequenceColorForLocation(colorData: SequenceColorData, location: StructureElement.Location): Color | undefined {
    const unitCache = colorData.colorCache[location.unit.id] ??= {};
    if (!(location.element in unitCache)) { // not using ??= pattern here, because cache may contain undefineds
        unitCache[location.element] = findSequenceColorForLocation(colorData, location);
    }
    return unitCache[location.element];
}

function findSequenceColorForLocation(colorData: SequenceColorData, location: StructureElement.Location): Color | undefined {
    for (let i = colorData.items.length - 1; i >= 0; i--) { // last color matters
        const item = colorData.items[i];
        const elements = item.elementSet ??= ElementSet.fromSelector(location.structure, item.selector);
        if (ElementSet.has(elements, location)) {
            return item.color;
        }
    }
    return undefined;
}
