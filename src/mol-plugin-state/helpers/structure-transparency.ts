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
import { Transparency } from '../../mol-theme/transparency';

type TransparencyEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, transparency?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.TransparencyStructureRepresentation3DFromBundle>>) => Promise<void>
const TransparencyManagerTag = 'transparency-controls';

export async function setStructureTransparency(plugin: PluginContext, components: StructureComponentRef[], value: number, lociGetter: (structure: Structure) => Promise<StructureElement.Loci | EmptyLoci>, types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, transparencyCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;

        const structure = repr.obj!.data.source.data;
        // always use the root structure to get the loci so the transparency
        // stays applicable as long as the root structure does not change
        const loci = await lociGetter(structure.root);
        if (Loci.isEmpty(loci)) return;

        const layer = {
            bundle: StructureElement.Bundle.fromLoci(loci),
            value,
        };

        if (transparencyCell) {
            const bundleLayers = [...transparencyCell.params!.values.layers, layer];
            const filtered = getFilteredBundle(bundleLayers, structure);
            update.to(transparencyCell).update(Transparency.toBundle(filtered));
        } else {
            const filtered = getFilteredBundle([layer], structure);
            update.to(repr.transform.ref)
                .apply(StateTransforms.Representation.TransparencyStructureRepresentation3DFromBundle, Transparency.toBundle(filtered), { tags: TransparencyManagerTag });
        }
    });
}

export async function clearStructureTransparency(plugin: PluginContext, components: StructureComponentRef[], types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, transparencyCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;
        if (transparencyCell) {
            update.delete(transparencyCell.transform.ref);
        }
    });
}

async function eachRepr(plugin: PluginContext, components: StructureComponentRef[], callback: TransparencyEachReprCallback) {
    const state = plugin.state.data;
    const update = state.build();
    for (const c of components) {
        for (const r of c.representations) {
            const transparency = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.TransparencyStructureRepresentation3DFromBundle, r.cell.transform.ref).withTag(TransparencyManagerTag));
            await callback(update, r.cell, transparency[0]);
        }
    }

    return update.commit({ doNotUpdateCurrent: true });
}

/** filter transparency layers for given structure */
function getFilteredBundle(layers: Transparency.BundleLayer[], structure: Structure) {
    const transparency = Transparency.ofBundle(layers, structure.root);
    const merged = Transparency.merge(transparency);
    return Transparency.filter(merged, structure);
}