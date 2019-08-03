/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../mol-plugin/state/objects';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StateSelection, StateObjectCell, StateTransform, StateBuilder } from '../../mol-state';
import { Structure, StructureElement } from '../../mol-model/structure';
import { isEmptyLoci, EmptyLoci } from '../../mol-model/loci';
import { PluginContext } from '../context';
import { Color } from '../../mol-util/color';
import { MolScriptBuilder } from '../../mol-script/language/builder';
import { formatMolScript } from '../../mol-script/language/expression-formatter';

type OverpaintEachReprCallback = (update: StateBuilder.Root, repr: StateObjectCell<PluginStateObject.Molecule.Structure.Representation3D, StateTransform<typeof StateTransforms.Representation.StructureRepresentation3D>>, rootStructure: Structure, overpaint?: StateObjectCell<any, StateTransform<typeof StateTransforms.Representation.OverpaintStructureRepresentation3D>>) => void
const OverpaintManagerTag = 'overpaint-controls'

export function getExpression(loci: StructureElement.Loci | EmptyLoci) {
    const scriptExpression = isEmptyLoci(loci)
        ? MolScriptBuilder.struct.generator.empty()
        : StructureElement.Loci.toScriptExpression(loci)
    return formatMolScript(scriptExpression)
}

export class StructureOverpaintHelper {
    private async eachRepr(callback: OverpaintEachReprCallback) {
        const state = this.plugin.state.dataState;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3D));

        const update = state.build();
        for (const r of reprs) {
            const overpaint = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.OverpaintStructureRepresentation3D, r.transform.ref).withTag(OverpaintManagerTag));

            const structure = r.obj!.data.source.data
            const rootStructure = structure.parent || structure

            callback(update, r, rootStructure, overpaint[0])
        }

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }));
    }

    async set(color: Color | -1, types?: string[]) {
        await this.eachRepr((update, repr, rootStructure, overpaint) => {
            if (types && !types.includes(repr.params!.values.type.name)) return

            const loci = this.plugin.helpers.structureSelectionManager.get(rootStructure)
            if (isEmptyLoci(loci) || loci.elements.length === 0) return
            const expression = getExpression(loci)

            const layer = {
                script: { language: 'mol-script', expression },
                color: color === -1 ? Color(0) : color,
                clear: color === -1
            }

            if (overpaint) {
                update.to(overpaint).update({ layers: [ ...overpaint.params!.values.layers, layer ], alpha: 1 })
            } else {
                update.to(repr.transform.ref)
                    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3D, { layers: [ layer ], alpha: 1 }, { tags: OverpaintManagerTag });
            }
        })
    }

    add(color: Color, types?: string[]) {
        this.set(color, types)
    }

    clear(types?: string[]) {
        this.set(-1, types)
    }

    clearAll() {
        this.eachRepr((update, repr, rootStructure, overpaint) => {
            if (overpaint) update.delete(overpaint.transform.ref)
        })
    }

    constructor(private plugin: PluginContext) {

    }
}