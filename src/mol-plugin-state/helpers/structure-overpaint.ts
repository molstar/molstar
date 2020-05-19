/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObjectCell, StateSelection, StateTransform } from '../../mol-state';
import { Overpaint } from '../../mol-theme/overpaint';
import { Color } from '../../mol-util/color';
import { StructureComponentRef } from '../manager/structure/hierarchy-state';
import { EmptyLoci, Loci } from '../../mol-model/loci';

type OverpaintEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, overpaint?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle>>) => Promise<void>
const OverpaintManagerTag = 'overpaint-controls';

export async function setStructureOverpaint(plugin: PluginContext, components: StructureComponentRef[], color: Color | -1, lociGetter: (structure: Structure) => Promise<StructureElement.Loci | EmptyLoci>, types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, overpaintCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;

        const structure = repr.obj!.data.source.data;
        // always use the root structure to get the loci so the overpaint
        // stays applicable as long as the root structure does not change
        const loci = await lociGetter(structure.root);
        if (Loci.isEmpty(loci)) return;

        const layer = {
            bundle: StructureElement.Bundle.fromLoci(loci),
            color: color === -1 ? Color(0) : color,
            clear: color === -1
        };

        if (overpaintCell) {
            const bundleLayers = [...overpaintCell.params!.values.layers, layer];
            const filtered = getFilteredBundle(bundleLayers, structure);
            update.to(overpaintCell).update(Overpaint.toBundle(filtered));
        } else {
            const filtered = getFilteredBundle([layer], structure);
            update.to(repr.transform.ref)
                .apply(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, Overpaint.toBundle(filtered), { tags: OverpaintManagerTag });
        }
    });
}

export async function clearStructureOverpaint(plugin: PluginContext, components: StructureComponentRef[], types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, overpaintCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;
        if (overpaintCell) {
            update.delete(overpaintCell.transform.ref);
        }
    });
}

async function eachRepr(plugin: PluginContext, components: StructureComponentRef[], callback: OverpaintEachReprCallback) {
    const state = plugin.state.data;
    const update = state.build();
    for (const c of components) {
        for (const r of c.representations) {
            const overpaint = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, r.cell.transform.ref).withTag(OverpaintManagerTag));
            await callback(update, r.cell, overpaint[0]);
        }
    }

    return update.commit({ doNotUpdateCurrent: true });
}

/** filter overpaint layers for given structure */
function getFilteredBundle(layers: Overpaint.BundleLayer[], structure: Structure) {
    const overpaint = Overpaint.ofBundle(layers, structure.root);
    const merged = Overpaint.merge(overpaint);
    return Overpaint.filter(merged, structure);
}