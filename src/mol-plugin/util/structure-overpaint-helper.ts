/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../mol-plugin/state/objects';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StateSelection, StateObjectCell, StateTransform, StateBuilder } from '../../mol-state';
import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginContext } from '../context';
import { Color } from '../../mol-util/color';

type OverpaintEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, overpaint?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle>>) => void
const OverpaintManagerTag = 'overpaint-controls'

export class StructureOverpaintHelper {
    private async eachRepr(callback: OverpaintEachReprCallback) {
        const state = this.plugin.state.dataState;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3D));

        const update = state.build();
        for (const r of reprs) {
            const overpaint = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, r.transform.ref).withTag(OverpaintManagerTag))
            callback(update, r, overpaint[0])
        }

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }));
    }

    async set(color: Color | -1, lociGetter: (structure: Structure) => StructureElement.Loci, types?: string[]) {
        await this.eachRepr((update, repr, overpaint) => {
            if (types && !types.includes(repr.params!.values.type.name)) return

            // TODO merge overpaint layers, delete shadowed ones
            // TODO filter overpaint layers for given structure

            const structure = repr.obj!.data.source.data
            // always use the root structure to get the loci so the overpaint
            // stays applicable as long as the root structure does not change
            const loci = lociGetter(structure.root)
            if (StructureElement.Loci.isEmpty(loci)) return

            const layer = {
                bundle: StructureElement.Bundle.fromLoci(loci),
                color: color === -1 ? Color(0) : color,
                clear: color === -1
            }

            if (overpaint) {
                update.to(overpaint).update({ layers: [ ...overpaint.params!.values.layers, layer ], alpha: 1 })
            } else {
                update.to(repr.transform.ref)
                    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, { layers: [ layer ], alpha: 1 }, { tags: OverpaintManagerTag });
            }
        })
    }

    constructor(private plugin: PluginContext) {

    }
}