/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { StateSelection } from '../../mol-state';
import { PluginStateObject } from '../state/objects';
import { QueryContext, StructureSelection, QueryFn, Queries as _Queries } from '../../mol-model/structure';
import { compile } from '../../mol-script/runtime/query/compiler';
import { ButtonsType } from '../../mol-util/input/input-observer';
import { EmptyLoci } from '../../mol-model/loci';
import { PluginContext } from '../context';

export const StructureSelectionQueries = {
    all: () => compile<StructureSelection>(MS.struct.generator.all()),
    polymers: () => _Queries.internal.atomicSequence(),
    water: () => _Queries.internal.water(),
    ligands: () => _Queries.internal.atomicHet(),
    coarse: () => _Queries.internal.spheres(),
}

export class StructureSelectionHelper {
    select(query: QueryFn<StructureSelection>) {
        const state = this.plugin.state.dataState
        const structures = state.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure))

        for (const so of structures) {
            const s = so.obj!.data
            const result = query(new QueryContext(s))
            const loci = StructureSelection.toLoci2(result)

            // TODO use better API when available
            this.plugin.interactivity.lociSelections.apply({
                current: { loci },
                buttons: ButtonsType.Flag.Secondary,
                modifiers: { shift: false, alt: false, control: true, meta: false }
            })
        }
    }

    clearSelection() {
        // TODO use better API when available
        this.plugin.interactivity.lociSelections.apply({
            current: { loci: EmptyLoci },
            buttons: ButtonsType.Flag.Secondary,
            modifiers: { shift: false, alt: false, control: true, meta: false }
        })
    }

    constructor(private plugin: PluginContext) {

    }
}