/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../mol-plugin/state/objects';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StateSelection, StateObjectCell, StateTransform, StateBuilder } from '../../mol-state';
import { Structure, StructureElement, StructureSelection, QueryContext } from '../../mol-model/structure';
import { PluginContext } from '../context';
import { Color } from '../../mol-util/color';
import { Overpaint } from '../../mol-theme/overpaint';
import Expression from '../../mol-script/language/expression';
import { compile } from '../../mol-script/runtime/query/compiler';

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

    async set(color: Color | -1, lociGetter: (structure: Structure) => StructureElement.Loci, types?: string[], alpha = 1) {
        await this.eachRepr((update, repr, overpaintCell) => {
            if (types && !types.includes(repr.params!.values.type.name)) return

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

            if (overpaintCell) {
                const bundleLayers = [ ...overpaintCell.params!.values.layers, layer ]
                const filtered = getFilteredBundle(bundleLayers, structure)
                update.to(overpaintCell).update(Overpaint.toBundle(filtered, alpha))
            } else {
                const filtered = getFilteredBundle([ layer ], structure)
                update.to(repr.transform.ref)
                    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, Overpaint.toBundle(filtered, alpha), { tags: OverpaintManagerTag });
            }
        })
    }

    async setFromExpression(color: Color | -1, expression: Expression, types?: string[], alpha = 1) {
        return this.set(color, (structure) => {
            const compiled = compile<StructureSelection>(expression)
            const result = compiled(new QueryContext(structure))
            return StructureSelection.toLociWithSourceUnits(result)
        }, types, alpha)
    }

    constructor(private plugin: PluginContext) {

    }
}

/** filter overpaint layers for given structure */
function getFilteredBundle(layers: Overpaint.BundleLayer[], structure: Structure) {
    const overpaint = Overpaint.ofBundle(layers, 1, structure.root)
    const merged = Overpaint.merge(overpaint)
    return Overpaint.filter(merged, structure)
}