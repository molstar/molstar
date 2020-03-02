/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

type ReprProps = {
    repr?: {},
    color?: string | [string, {}],
    size?: string | [string, {}],
}

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

    private getRepresentationParams(structure: Structure, type: string, repr: RepresentationTransform | undefined, props: ReprProps = {}) {
        const reprProps = {
            ...(repr?.params && repr.params.values.type.params),
            ignoreHydrogens: this._ignoreHydrogens,
            quality: this._quality,
            ...props.repr
        }
        const { themeCtx } =  this.plugin.structureRepresentation

        const p: StructureRepresentation3DHelpers.Props = {
            repr: [
                this.plugin.structureRepresentation.registry.get(type),
                () => reprProps
            ]
        }
        if (props.color) {
            const colorType = props.color instanceof Array ? props.color[0] : props.color
            const colorTheme = themeCtx.colorThemeRegistry.get(colorType)
            const colorProps = {
                ...(repr?.params && repr.params.values.colorTheme.name === colorType && repr.params.values.colorTheme.params),
                ...(props.color instanceof Array ? props.color[1] : {})
            }
            p.color = [colorTheme, () => colorProps]
        }
        if (props.size) {
            const sizeType = props.size instanceof Array ? props.size[0] : props.size
            const sizeTheme = themeCtx.sizeThemeRegistry.get(sizeType)
            const sizeProps = {
                ...(repr?.params && repr.params.values.sizeTheme.name === sizeType && repr.params.values.sizeTheme.params),
                ...(props.size instanceof Array ? props.size[1] : {})
            }
            p.size = [sizeTheme, () => sizeProps]
        }
        if (props.size) p.size = props.size

        return StructureRepresentation3DHelpers.createParams(this.plugin, structure, p)
    }

    private async _set(modifier: SelectionModifier, type: string, loci: StructureElement.Loci, structure: StructureTransform, props: ReprProps = {}) {
        const state = this.plugin.state.dataState
        const update = state.build()
        const s = structure.obj!.data

        const repr = this.getRepresentation(structure.transform.ref, type)
        const reprStructure = this.getRepresentationStructure(structure.transform.ref, type)
        const reprParams = this.getRepresentationParams(s.root, type, repr, props)

        if (reprStructure) {
            const currentLoci = StructureElement.Bundle.toLoci(reprStructure.params!.values.bundle, s)
            const combinedLoci = getCombinedLoci(modifier, loci, currentLoci)

            update.to(reprStructure).update({
                ...reprStructure.params!.values,
                bundle: StructureElement.Bundle.fromLoci(combinedLoci)
            })
            if (repr) update.to(repr).update(reprParams)
        } else {
            const combinedLoci = getCombinedLoci(modifier, loci, StructureElement.Loci.none(s))

            update.to(structure.transform.ref)
                .apply(
                    StateTransforms.Model.StructureSelectionFromBundle,
                    { bundle: StructureElement.Bundle.fromLoci(combinedLoci), label: type },
                    { tags: [ RepresentationManagerTag, getRepresentationManagerTag(type) ] }
                )
                .apply(StateTransforms.Representation.StructureRepresentation3D, reprParams)
        }

        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }))
    }

    async set(modifier: SelectionModifier, type: string, lociGetter: (structure: Structure) => StructureElement.Loci, props: ReprProps = {}) {
        const state = this.plugin.state.dataState;
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))

        for (const structure of structures) {
            const s = structure.obj!.data
            const loci = lociGetter(s)
            await this._set(modifier, type, loci, structure, props)
        }
    }

    async setFromExpression(modifier: SelectionModifier, type: string, expression: Expression, props: ReprProps = {}) {
        return this.set(modifier, type, (structure) => {
            const compiled = compile<StructureSelection>(expression)
            const result = compiled(new QueryContext(structure))
            return StructureSelection.toLociWithSourceUnits(result)
        }, props)
    }

    async eachStructure(callback: (structure: StructureTransform, type: string, update: StateBuilder.Root) => void) {
        const { registry } = this.plugin.structureRepresentation
        const state = this.plugin.state.dataState;
        const update = state.build()
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))

        for (const structure of structures) {
            for (let i = 0, il = registry.types.length; i < il; ++i) {
                const type = registry.types[i][0]
                const reprStructure = this.getRepresentationStructure(structure.transform.ref, type)
                if (reprStructure) callback(reprStructure, type, update)
            }
        }
        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }))
    }

    async clear() {
        const bundle = StructureElement.Bundle.Empty
        await this.eachStructure((structure, type, update) => {
            update.to(structure).update({ ...structure.params!.values, bundle })
        })
    }

    async clearExcept(exceptTypes: string[]) {
        const bundle = StructureElement.Bundle.Empty
        await this.eachStructure((structure, type, update) => {
            if (!exceptTypes.includes(type)) {
                update.to(structure).update({ ...structure.params!.values, bundle })
            }
        })
    }

    async eachRepresentation(callback: (repr: RepresentationTransform, type: string, update: StateBuilder.Root) => void) {
        const { registry } = this.plugin.structureRepresentation
        const state = this.plugin.state.dataState;
        const update = state.build()
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))
        for (const structure of structures) {
            for (let i = 0, il = registry.types.length; i < il; ++i) {
                const type = registry.types[i][0]
                const repr = this.getRepresentation(structure.transform.ref, type)
                if (repr) callback(repr, type, update)
            }
        }
        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }))
    }

    setRepresentationParams(repr: RepresentationTransform, type: string, update: StateBuilder.Root, props: ReprProps) {
        const state = this.plugin.state.dataState;
        const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))

        for (const structure of structures) {
            const s = structure.obj!.data
            const reprParams = this.getRepresentationParams(s.root, type, repr, props)
            update.to(repr).update(reprParams)
        }
    }

    async updateRepresentation(repr: RepresentationTransform, type: string, props: ReprProps) {
        const state = this.plugin.state.dataState;
        const update = state.build()
        this.setRepresentationParams(repr, type, update, props)
        await this.plugin.runTask(state.updateTree(update, { doNotUpdateCurrent: true }))
    }

    private _ignoreHydrogens = false
    get ignoreHydrogens () { return this._ignoreHydrogens }
    async setIgnoreHydrogens(ignoreHydrogens: boolean) {
        if (ignoreHydrogens === this._ignoreHydrogens) return
        await this.eachRepresentation((repr, type, update) => {
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
        await this.eachRepresentation((repr, type, update) => {
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

    constructor(private plugin: PluginContext) {

    }
}