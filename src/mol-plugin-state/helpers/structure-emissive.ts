/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObjectCell, StateSelection, StateTransform } from '../../mol-state';
import { StructureComponentRef } from '../manager/structure/hierarchy-state';
import { EmptyLoci, isEmptyLoci, Loci } from '../../mol-model/loci';
import { Emissive } from '../../mol-theme/emissive';

type EmissiveEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, emissive?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.EmissiveStructureRepresentation3DFromBundle>>) => Promise<void>
const EmissiveManagerTag = 'emissive-controls';

export async function setStructureEmissive(plugin: PluginContext, components: StructureComponentRef[], value: number, lociGetter: (structure: Structure) => Promise<StructureElement.Loci | EmptyLoci>, types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, emissiveCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;

        const structure = repr.obj!.data.sourceData;
        // always use the root structure to get the loci so the emissive
        // stays applicable as long as the root structure does not change
        const loci = await lociGetter(structure.root);
        if (Loci.isEmpty(loci) || isEmptyLoci(loci)) return;

        const layer = {
            bundle: StructureElement.Bundle.fromLoci(loci),
            value,
        };

        if (emissiveCell) {
            const bundleLayers = [...emissiveCell.params!.values.layers, layer];
            const filtered = getFilteredBundle(bundleLayers, structure);
            update.to(emissiveCell).update(Emissive.toBundle(filtered));
        } else {
            const filtered = getFilteredBundle([layer], structure);
            update.to(repr.transform.ref)
                .apply(StateTransforms.Representation.EmissiveStructureRepresentation3DFromBundle, Emissive.toBundle(filtered), { tags: EmissiveManagerTag });
        }
    });
}

export async function clearStructureEmissive(plugin: PluginContext, components: StructureComponentRef[], types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, emissiveCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;
        if (emissiveCell) {
            update.delete(emissiveCell.transform.ref);
        }
    });
}

async function eachRepr(plugin: PluginContext, components: StructureComponentRef[], callback: EmissiveEachReprCallback) {
    const state = plugin.state.data;
    const update = state.build();
    for (const c of components) {
        for (const r of c.representations) {
            const emissive = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.EmissiveStructureRepresentation3DFromBundle, r.cell.transform.ref).withTag(EmissiveManagerTag));
            await callback(update, r.cell, emissive[0]);
        }
    }

    return update.commit({ doNotUpdateCurrent: true });
}

/** filter emissive layers for given structure */
function getFilteredBundle(layers: Emissive.BundleLayer[], structure: Structure) {
    const emissive = Emissive.ofBundle(layers, structure.root);
    const merged = Emissive.merge(emissive);
    return Emissive.filter(merged, structure) as Emissive<StructureElement.Loci>;
}