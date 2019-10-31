/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject as PSO } from '../../mol-plugin/state/objects';
import { StateTransforms } from '../../mol-plugin/state/transforms';
import { StateTransformer, StateSelection, StateObjectCell, StateTransform, StateBuilder } from '../../mol-state';
import { StructureElement, Structure, StructureSelection, QueryContext } from '../../mol-model/structure';
import { PluginContext } from '../context';
import { StructureRepresentation3DHelpers } from '../state/transforms/representation';
import Expression from '../../mol-script/language/expression';
import { compile } from '../../mol-script/runtime/query/compiler';
import { StructureSelectionQueries as Q } from '../util/structure-selection-helper';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { VisualQuality } from '../../mol-geo/geometry/base';

type StructureTransform = StateObjectCell<PSO.Molecule.Structure, StateTransform<StateTransformer<any, PSO.Molecule.Structure, any>>>
type RepresentationTransform = StateObjectCell<PSO.Molecule.Structure.Representation3D, StateTransform<StateTransformer<any, PSO.Molecule.Structure.Representation3D, any>>>
const RepresentationManagerTag = 'representation-controls'

export function getRepresentationManagerTag(type: string) {
    return `${RepresentationManagerTag}-${type}`
}

function getCombinedLoci(mode: SelectionModifier, loci: StructureElement.Loci, currentLoci: StructureElement.Loci): StructureElement.Loci {
    switch (mode) {
        case 'add': return StructureElement.Loci.union(loci, currentLoci)
        case 'remove': return StructureElement.Loci.subtract(currentLoci, loci)
        case 'only': return loci
    }
}

type SelectionModifier = 'add' | 'remove' | 'only'

export class StructureRepresentationHelper {
    getRepresentationStructure(rootRef: string, type: string) {
        const state = this.plugin.state.dataState
        const selections = state.select(StateSelection.Generators.ofType(PSO.Molecule.Structure, rootRef).withTag(getRepresentationManagerTag(type)));
        return selections.length > 0 ? selections[0] : undefined
    }

    getRepresentation(rootRef: string, type: string) {
        const reprStructure = this.getRepresentationStructure(rootRef, type)
        if (!reprStructure) return
        const state = this.plugin.state.dataState
        const selections = state.select(StateSelection.Generators.ofType(PSO.Molecule.Structure.Representation3D, reprStructure.transform.ref))
        return selections.length > 0 ? selections[0] : undefined
    }

    private async _set(modifier: SelectionModifier, type: string, loci: StructureElement.Loci, structure: StructureTransform, props = {}) {
        const state = this.plugin.state.dataState
        const update = state.build()
        const s = structure.obj!.data

        const reprStructure = this.getRepresentationStructure(structure.transform.ref, type)

        if (reprStructure) {
            const currentLoci = StructureElement.Bundle.toLoci(reprStructure.params!.values.bundle, s)
            const combinedLoci = getCombinedLoci(modifier, loci, currentLoci)

            update.to(reprStructure).update({
                ...reprStructure.params!.values,
                bundle: StructureElement.Bundle.fromLoci(combinedLoci)
            })
        } else {
            const combinedLoci = getCombinedLoci(modifier, loci, StructureElement.Loci.none(s))
            const params = StructureRepresentation3DHelpers.getDefaultParams(this.plugin, type as any, s)

            const p = params.type.params
            Object.assign(p, props)
            if (p.ignoreHydrogens !== undefined) p.ignoreHydrogens = this._ignoreHydrogens
            if (p.quality !== undefined) p.quality = this._quality

            update.to(structure.transform.ref)
                .apply(
                    StateTransforms.Model.StructureSelectionFromBundle,
                    { bundle: StructureElement.Bundle.fromLoci(combinedLoci), label: type },
                    { tags: [ RepresentationManagerTag, getRepresentationManagerTag(type) ] }
                )
                .apply( StateTransforms.Representation.StructureRepresentation3D, params)
        }

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }))
    }

    async set(modifier: SelectionModifier, type: string, lociGetter: (structure: Structure) => StructureElement.Loci, props = {}) {
        const state = this.plugin.state.dataState;
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))

        for (const structure of structures) {
            const s = structure.obj!.data
            const loci = lociGetter(s)
            await this._set(modifier, type, loci, structure, props)
        }
    }

    async setFromExpression(modifier: SelectionModifier, type: string, expression: Expression, props = {}) {
        return this.set(modifier, type, (structure) => {
            const compiled = compile<StructureSelection>(expression)
            const result = compiled(new QueryContext(structure))
            return StructureSelection.toLociWithSourceUnits(result)
        }, props)
    }

    async clear() {
        const { registry } = this.plugin.structureRepresentation
        const state = this.plugin.state.dataState;
        const update = state.build()
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))
        const bundle = StructureElement.Bundle.Empty

        for (const structure of structures) {
            for (let i = 0, il = registry.types.length; i < il; ++i) {
                const type = registry.types[i][0]
                const reprStructure = this.getRepresentationStructure(structure.transform.ref, type)
                if (reprStructure) {
                    update.to(reprStructure).update({ ...reprStructure.params!.values, bundle })
                }
            }
        }
        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }))
    }

    async eachRepresentation(callback: (repr: RepresentationTransform, update: StateBuilder.Root) => void) {
        const { registry } = this.plugin.structureRepresentation
        const state = this.plugin.state.dataState;
        const update = state.build()
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))
        for (const structure of structures) {
            for (let i = 0, il = registry.types.length; i < il; ++i) {
                const repr = this.getRepresentation(structure.transform.ref, registry.types[i][0])
                if (repr) callback(repr, update)
            }
        }
        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }))
    }

    private _ignoreHydrogens = false
    get ignoreHydrogens () { return this._ignoreHydrogens }
    async setIgnoreHydrogens(ignoreHydrogens: boolean) {
        if (ignoreHydrogens === this._ignoreHydrogens) return
        await this.eachRepresentation((repr, update) => {
            if (repr.params && repr.params.values.type.params.ignoreHydrogens !== undefined) {
                const { name, params } = repr.params.values.type
                update.to(repr.transform.ref).update(
                    StateTransforms.Representation.StructureRepresentation3D,
                    props => ({ ...props, type: { name, params: { ...params, ignoreHydrogens }}})
                )
            }
        })
        this._ignoreHydrogens = ignoreHydrogens
    }

    private _quality = 'auto' as VisualQuality
    get quality () { return this._quality }
    async setQuality(quality: VisualQuality) {
        if (quality === this._quality) return
        await this.eachRepresentation((repr, update) => {
            if (repr.params && repr.params.values.type.params.quality !== undefined) {
                const { name, params } = repr.params.values.type
                update.to(repr.transform.ref).update(
                    StateTransforms.Representation.StructureRepresentation3D,
                    props => ({ ...props, type: { name, params: { ...params, quality }}})
                )
            }
        })
        this._quality = quality
    }

    async preset() {
        // TODO option to limit to specific structure
        const state = this.plugin.state.dataState;
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))

        if (structures.length === 0) return
        const s = structures[0].obj!.data

        if (s.elementCount < 50000) {
            await polymerAndLigand(this)
        } else if (s.elementCount < 200000) {
            await proteinAndNucleic(this)
        } else {
            if (s.unitSymmetryGroups[0].units.length > 10) {
                await capsid(this)
            } else {
                await coarseCapsid(this)
            }
        }
    }

    constructor(private plugin: PluginContext) {

    }
}

//

async function polymerAndLigand(r: StructureRepresentationHelper) {
    await r.clear()
    await r.setFromExpression('add', 'cartoon', Q.polymer.expression)
    await r.setFromExpression('add', 'carbohydrate', Q.branched.expression)
    await r.setFromExpression('add', 'ball-and-stick', MS.struct.modifier.union([
        MS.struct.combinator.merge([
            Q.ligandPlusConnected.expression,
            Q.branchedConnectedOnly.expression,
            Q.nonStandardPolymer.expression,
            Q.water.expression
        ])
    ]))
}

async function proteinAndNucleic(r: StructureRepresentationHelper) {
    await r.clear()
    await r.setFromExpression('add', 'cartoon', Q.protein.expression)
    await r.setFromExpression('add', 'gaussian-surface', Q.nucleic.expression)
}

async function capsid(r: StructureRepresentationHelper) {
    await r.clear()
    await r.setFromExpression('add', 'gaussian-surface', Q.polymer.expression, {
        smoothness: 0.5,
    })
}

async function coarseCapsid(r: StructureRepresentationHelper) {
    await r.clear()
    await r.setFromExpression('add', 'gaussian-surface', Q.trace.expression, {
        radiusOffset: 1,
        smoothness: 0.5,
        visuals: ['structure-gaussian-surface-mesh']
    })
}

export const StructureRepresentationPresets = {
    polymerAndLigand,
    proteinAndNucleic,
    capsid,
    coarseCapsid
}