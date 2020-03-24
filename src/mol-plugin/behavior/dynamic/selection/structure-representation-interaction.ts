/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InteractionsRepresentationProvider } from '../../../../mol-model-props/computed/representations/interactions';
import { InteractionTypeColorThemeProvider } from '../../../../mol-model-props/computed/themes/interaction-type';
import { EmptyLoci, isEmptyLoci, Loci } from '../../../../mol-model/loci';
import { Bond, Structure, StructureElement } from '../../../../mol-model/structure';
import { createStructureRepresentationParams } from '../../../../mol-plugin-state/helpers/structure-representation-params';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { StateTransforms } from '../../../../mol-plugin-state/transforms';
import { PluginBehavior } from '../../../../mol-plugin/behavior';
import { MolScriptBuilder as MS } from '../../../../mol-script/language/builder';
import { StateObjectCell, StateSelection, StateTransform } from '../../../../mol-state';
import { SizeTheme } from '../../../../mol-theme/size';
import { Binding } from '../../../../mol-util/binding';
import { ButtonsType, ModifiersKeys } from '../../../../mol-util/input/input-observer';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { PluginCommands } from '../../../commands';
import { PluginContext } from '../../../context';

const B = ButtonsType
const M = ModifiersKeys
const Trigger = Binding.Trigger

const DefaultStructureRepresentationInteractionBindings = {
    clickInteractionAroundOnly: Binding([Trigger(B.Flag.Secondary, M.create()), Trigger(B.Flag.Primary, M.create({ control: true }))], 'Structure element interaction', 'Click element using ${triggers}; on nothing/same element to hide'),
}

const StructureRepresentationInteractionParams = (plugin: PluginContext) => {
    const reprParams = StateTransforms.Representation.StructureRepresentation3D.definition.params!(void 0, plugin) as PD.Params;
    return {
        bindings: PD.Value(DefaultStructureRepresentationInteractionBindings, { isHidden: true }),
        // TODO: min = 0 to turn them off?
        expandRadius: PD.Numeric(5, { min: 1, max: 10, step: 1 }),
        focusParams: PD.Group(reprParams, {
            label: 'Focus',
            customDefault: createStructureRepresentationParams(plugin, void 0, { type: 'ball-and-stick', size: 'uniform' })
        }),
        surroundingsParams: PD.Group(reprParams, {
            label: 'Surroundings',
            customDefault: createStructureRepresentationParams(plugin, void 0, { type: 'ball-and-stick', color: 'element-symbol', size: 'uniform' })
        }),
        nciParams: PD.Group(reprParams, {
            label: 'Non-covalent Int.',
            customDefault: createStructureRepresentationParams(plugin, void 0, {
                type: InteractionsRepresentationProvider,
                color: InteractionTypeColorThemeProvider,
                size: SizeTheme.BuiltIn.uniform
            })
        })
    };
}

type StructureRepresentationInteractionProps = PD.ValuesFor<ReturnType<typeof StructureRepresentationInteractionParams>>

export enum StructureRepresentationInteractionTags {
    ResidueSel = 'structure-interaction-residue-sel',
    ResidueRepr = 'structure-interaction-residue-repr',
    SurrSel = 'structure-interaction-surr-sel',
    SurrRepr = 'structure-interaction-surr-repr',
    SurrNciRepr = 'structure-interaction-surr-nci-repr'
}

const TagSet: Set<StructureRepresentationInteractionTags> = new Set([StructureRepresentationInteractionTags.ResidueSel, StructureRepresentationInteractionTags.ResidueRepr, StructureRepresentationInteractionTags.SurrSel, StructureRepresentationInteractionTags.SurrRepr, StructureRepresentationInteractionTags.SurrNciRepr])

export class StructureRepresentationInteractionBehavior extends PluginBehavior.WithSubscribers<StructureRepresentationInteractionProps> {
    private get surrLabel() { return `[Focus +${this.params.expandRadius} Å Surrounding]`; }

    private ensureShape(cell: StateObjectCell<PluginStateObject.Molecule.Structure>) {
        const state = this.plugin.state.data, tree = state.tree;
        const builder = state.build();
        const refs = StateSelection.findUniqueTagsInSubtree(tree, cell.transform.ref, TagSet);

        // Selections
        if (!refs[StructureRepresentationInteractionTags.ResidueSel]) {
            refs[StructureRepresentationInteractionTags.ResidueSel] = builder
                .to(cell) // refs['structure-interaction-group'])
                .apply(StateTransforms.Model.StructureSelectionFromBundle,
                    { bundle: StructureElement.Bundle.Empty, label: '[Focus]' }, { tags: StructureRepresentationInteractionTags.ResidueSel }).ref;
        }

        if (!refs[StructureRepresentationInteractionTags.SurrSel]) {
            refs[StructureRepresentationInteractionTags.SurrSel] = builder
                .to(cell) // .to(refs['structure-interaction-group'])
                .apply(StateTransforms.Model.StructureSelectionFromExpression,
                    { expression: MS.struct.generator.empty(), label: this.surrLabel }, { tags: StructureRepresentationInteractionTags.SurrSel }).ref;
        }

        // Representations
        if (!refs[StructureRepresentationInteractionTags.ResidueRepr]) {
            refs[StructureRepresentationInteractionTags.ResidueRepr] = builder
                .to(refs['structure-interaction-residue-sel']!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.params.focusParams, { tags: StructureRepresentationInteractionTags.ResidueRepr }).ref;
        }

        if (!refs[StructureRepresentationInteractionTags.SurrRepr]) {
            refs[StructureRepresentationInteractionTags.SurrRepr] = builder
                .to(refs['structure-interaction-surr-sel']!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.params.surroundingsParams, { tags: StructureRepresentationInteractionTags.SurrRepr }).ref;
        }

        if (!refs[StructureRepresentationInteractionTags.SurrNciRepr]) {
            refs[StructureRepresentationInteractionTags.SurrNciRepr] = builder
                .to(refs['structure-interaction-surr-sel']!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.params.nciParams, { tags: StructureRepresentationInteractionTags.SurrNciRepr }).ref;
        }

        return { state, builder, refs };
    }

    private clear(root: StateTransform.Ref) {
        const state = this.plugin.state.data;

        const foci = state.select(StateSelection.Generators.byRef(root).subtree().withTag(StructureRepresentationInteractionTags.ResidueSel));
        const surrs = state.select(StateSelection.Generators.byRef(root).subtree().withTag(StructureRepresentationInteractionTags.SurrSel));
        if (foci.length === 0 && surrs.length === 0) return;

        const update = state.build();
        const bundle = StructureElement.Bundle.Empty;
        for (const f of foci) {
            update.to(f).update(StateTransforms.Model.StructureSelectionFromBundle, old => ({ ...old, bundle }));
        }

        const expression = MS.struct.generator.empty();
        for (const s of surrs) {
            update.to(s).update(StateTransforms.Model.StructureSelectionFromExpression, old => ({ ...old, expression }));
        }

        return PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });
    }

    register(ref: string): void {
        let lastLoci: Loci = EmptyLoci;

        this.subscribeObservable(this.plugin.events.state.object.removed, o => {
            if (!PluginStateObject.Molecule.Structure.is(o.obj) || !StructureElement.Loci.is(lastLoci)) return;
            if (lastLoci.structure === o.obj.data) {
                lastLoci = EmptyLoci;
            }
        });

        this.subscribeObservable(this.plugin.events.state.object.updated, o => {
            if (!PluginStateObject.Molecule.Structure.is(o.oldObj) || !StructureElement.Loci.is(lastLoci)) return;
            if (lastLoci.structure === o.oldObj.data) {
                lastLoci = EmptyLoci;
            }
        });

        this.subscribeObservable(this.plugin.behaviors.interaction.click, ({ current, button, modifiers }) => {
            const { clickInteractionAroundOnly } = this.params.bindings

            if (Binding.match(clickInteractionAroundOnly, button, modifiers)) {
                if (isEmptyLoci(current.loci)) {
                    this.clear(StateTransform.RootRef);
                    lastLoci = current.loci;
                    return;
                }

                let loci: StructureElement.Loci;
                if (StructureElement.Loci.is(current.loci)) {
                    loci = current.loci;
                } else if (Bond.isLoci(current.loci)) {
                    loci = Bond.toStructureElementLoci(current.loci);
                } else if (Structure.isLoci(current.loci)) {
                    loci = Structure.toStructureElementLoci(current.loci.structure);
                } else {
                    return;
                }

                if (StructureElement.Loci.isEmpty(loci)) return;

                const parent = this.plugin.helpers.substructureParent.get(loci.structure);
                if (!parent || !parent.obj) return;

                if (Loci.areEqual(lastLoci, loci)) {
                    lastLoci = EmptyLoci;
                    this.clear(parent.transform.ref);
                    return;
                }

                lastLoci = loci;

                const residueLoci = StructureElement.Loci.extendToWholeResidues(StructureElement.Loci.remap(loci, parent.obj!.data))
                const residueBundle = StructureElement.Bundle.fromLoci(residueLoci)

                const surroundings = MS.struct.modifier.includeSurroundings({
                    0: StructureElement.Bundle.toExpression(residueBundle),
                    radius: this.params.expandRadius,
                    'as-whole-residues': true
                });

                const { state, builder, refs } = this.ensureShape(parent);

                builder.to(refs[StructureRepresentationInteractionTags.ResidueSel]!).update(StateTransforms.Model.StructureSelectionFromBundle, old => ({ ...old, bundle: residueBundle }));
                builder.to(refs[StructureRepresentationInteractionTags.SurrSel]!).update(StateTransforms.Model.StructureSelectionFromExpression, old => ({ ...old, expression: surroundings, label: this.surrLabel }));

                PluginCommands.State.Update(this.plugin, { state, tree: builder, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });
            }
        });
    }

    async update(params: StructureRepresentationInteractionProps) {
        let oldRadius = this.params.expandRadius;
        this.params = params;

        const state = this.plugin.state.data;
        const builder = state.build();

        const all = StateSelection.Generators.root.subtree();
        for (const repr of state.select(all.withTag(StructureRepresentationInteractionTags.ResidueRepr))) {
            builder.to(repr).update(this.params.focusParams);
        }
        for (const repr of state.select(all.withTag(StructureRepresentationInteractionTags.SurrRepr))) {
            builder.to(repr).update(this.params.surroundingsParams);
        }
        for (const repr of state.select(all.withTag(StructureRepresentationInteractionTags.SurrNciRepr))) {
            builder.to(repr).update(this.params.nciParams);
        }

        await PluginCommands.State.Update(this.plugin, { state, tree: builder, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });

        // TODO: update properly
        if (params.expandRadius !== oldRadius) await this.clear(StateTransform.RootRef);

        return true;
    }
}

export const StructureRepresentationInteraction = PluginBehavior.create({
    name: 'create-structure-representation-interaction',
    display: { name: 'Structure Representation Interaction' },
    category: 'interaction',
    ctor: StructureRepresentationInteractionBehavior,
    params: (_, plugin) => StructureRepresentationInteractionParams(plugin)
});