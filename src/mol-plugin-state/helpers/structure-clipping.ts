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
import { StructureComponentRef } from '../manager/structure/hierarchy-state';
import { EmptyLoci, Loci } from '../../mol-model/loci';
import { Clipping } from '../../mol-theme/clipping';

type ClippingEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, clipping?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.ClippingStructureRepresentation3DFromBundle>>) => Promise<void>
const ClippingManagerTag = 'clipping-controls';

export async function setStructureClipping(plugin: PluginContext, components: StructureComponentRef[], groups: Clipping.Groups, lociGetter: (structure: Structure) => Promise<StructureElement.Loci | EmptyLoci>, types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, clippingCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;

        const structure = repr.obj!.data.source.data;
        // always use the root structure to get the loci so the clipping
        // stays applicable as long as the root structure does not change
        const loci = await lociGetter(structure.root);
        if (Loci.isEmpty(loci)) return;

        const layer = {
            bundle: StructureElement.Bundle.fromLoci(loci),
            groups
        };

        if (clippingCell) {
            const bundleLayers = [...clippingCell.params!.values.layers, layer];
            const filtered = getFilteredBundle(bundleLayers, structure);
            update.to(clippingCell).update(Clipping.toBundle(filtered));
        } else {
            const filtered = getFilteredBundle([layer], structure);
            update.to(repr.transform.ref)
                .apply(StateTransforms.Representation.ClippingStructureRepresentation3DFromBundle, Clipping.toBundle(filtered), { tags: ClippingManagerTag });
        }
    });
}

async function eachRepr(plugin: PluginContext, components: StructureComponentRef[], callback: ClippingEachReprCallback) {
    const state = plugin.state.data;
    const update = state.build();
    for (const c of components) {
        for (const r of c.representations) {
            const clipping = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.ClippingStructureRepresentation3DFromBundle, r.cell.transform.ref).withTag(ClippingManagerTag));
            await callback(update, r.cell, clipping[0]);
        }
    }

    return update.commit({ doNotUpdateCurrent: true });
}

/** filter clipping layers for given structure */
function getFilteredBundle(layers: Clipping.BundleLayer[], structure: Structure) {
    const clipping = Clipping.ofBundle(layers, structure.root);
    const merged = Clipping.merge(clipping);
    return Clipping.filter(merged, structure);
}