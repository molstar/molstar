/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginContext } from '../../mol-plugin/context';
import { StateBuilder, StateObjectCell, StateSelection, StateTransform } from '../../mol-state';
import { Substance } from '../../mol-theme/substance';
import { StructureComponentRef } from '../manager/structure/hierarchy-state';
import { EmptyLoci, isEmptyLoci, Loci } from '../../mol-model/loci';
import { Material } from '../../mol-util/material';

type SubstanceEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, substance?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.SubstanceStructureRepresentation3DFromBundle>>) => Promise<void>
const SubstanceManagerTag = 'substance-controls';

export async function setStructureSubstance(plugin: PluginContext, components: StructureComponentRef[], material: Material | undefined, lociGetter: (structure: Structure) => Promise<StructureElement.Loci | EmptyLoci>, types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, substanceCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;

        const structure = repr.obj!.data.sourceData;
        // always use the root structure to get the loci so the substance
        // stays applicable as long as the root structure does not change
        const loci = await lociGetter(structure.root);
        if (Loci.isEmpty(loci) || isEmptyLoci(loci)) return;

        const layer = {
            bundle: StructureElement.Bundle.fromLoci(loci),
            material: material ?? Material(),
            clear: !material
        };

        if (substanceCell) {
            const bundleLayers = [...substanceCell.params!.values.layers, layer];
            const filtered = getFilteredBundle(bundleLayers, structure);
            update.to(substanceCell).update(Substance.toBundle(filtered));
        } else {
            const filtered = getFilteredBundle([layer], structure);
            update.to(repr.transform.ref)
                .apply(StateTransforms.Representation.SubstanceStructureRepresentation3DFromBundle, Substance.toBundle(filtered), { tags: SubstanceManagerTag });
        }
    });
}

export async function clearStructureSubstance(plugin: PluginContext, components: StructureComponentRef[], types?: string[]) {
    await eachRepr(plugin, components, async (update, repr, substanceCell) => {
        if (types && types.length > 0 && !types.includes(repr.params!.values.type.name)) return;
        if (substanceCell) {
            update.delete(substanceCell.transform.ref);
        }
    });
}

async function eachRepr(plugin: PluginContext, components: StructureComponentRef[], callback: SubstanceEachReprCallback) {
    const state = plugin.state.data;
    const update = state.build();
    for (const c of components) {
        for (const r of c.representations) {
            const substance = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.SubstanceStructureRepresentation3DFromBundle, r.cell.transform.ref).withTag(SubstanceManagerTag));
            await callback(update, r.cell, substance[0]);
        }
    }

    return update.commit({ doNotUpdateCurrent: true });
}

/** filter substance layers for given structure */
function getFilteredBundle(layers: Substance.BundleLayer[], structure: Structure) {
    const substance = Substance.ofBundle(layers, structure.root);
    const merged = Substance.merge(substance);
    return Substance.filter(merged, structure) as Substance<StructureElement.Loci>;
}