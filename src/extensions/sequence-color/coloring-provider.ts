/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { StructureElement } from '../../mol-model/structure';
import { CustomStructureProperties } from '../../mol-plugin-state/transforms/model';
import { SequenceColoring } from '../../mol-plugin/util/sequence-coloring';
import { Color } from '../../mol-util/color';
import { ElementSet } from '../mvs/components/selector';
import { SequenceColorProperty } from './prop';


export const SequenceColoringProvider: SequenceColoring.Provider = {
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

/** Get color assigned to a loci in custom structure property "SequenceColor" */
function sequenceColorForLoci(loci: StructureElement.Loci): Color | undefined {
    const colorData = SequenceColorProperty.Provider.get(loci.structure).value;
    if (!colorData || colorData.items.length === 0) return undefined;
    const location = StructureElement.Loci.getFirstLocation(loci);
    if (location === undefined) return undefined;
    return sequenceColorForLocation(colorData, location);
}

function sequenceColorForLocation(colorData: SequenceColorProperty.Data, location: StructureElement.Location): Color | undefined {
    const unitCache = colorData.colorCache[location.unit.id] ??= {};
    if (!(location.element in unitCache)) { // not using ??= pattern here, because cache may contain undefineds
        unitCache[location.element] = findSequenceColorForLocation(colorData, location);
    }
    return unitCache[location.element];
}

function findSequenceColorForLocation(colorData: SequenceColorProperty.Data, location: StructureElement.Location): Color | undefined {
    for (let i = colorData.items.length - 1; i >= 0; i--) { // last color matters
        const item = colorData.items[i];
        const elements = item.elementSet ??= ElementSet.fromSelector(location.structure, item.selector);
        if (ElementSet.has(elements, location)) {
            return item.color;
        }
    }
    return undefined;
}
