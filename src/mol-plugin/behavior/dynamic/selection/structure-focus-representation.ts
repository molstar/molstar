/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InteractionsRepresentationProvider } from '../../../../mol-model-props/computed/representations/interactions';
import { InteractionTypeColorThemeProvider } from '../../../../mol-model-props/computed/themes/interaction-type';
import { StructureElement } from '../../../../mol-model/structure';
import { createStructureRepresentationParams } from '../../../../mol-plugin-state/helpers/structure-representation-params';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { StateTransforms } from '../../../../mol-plugin-state/transforms';
import { PluginBehavior } from '../../../behavior';
import { MolScriptBuilder as MS } from '../../../../mol-script/language/builder';
import { StateObjectCell, StateSelection, StateTransform } from '../../../../mol-state';
import { SizeTheme } from '../../../../mol-theme/size';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { PluginCommands } from '../../../commands';
import { PluginContext } from '../../../context';

const StructureFocusRepresentationParams = (plugin: PluginContext) => {
    const reprParams = StateTransforms.Representation.StructureRepresentation3D.definition.params!(void 0, plugin) as PD.Params;
    return {
        // TODO: min = 0 to turn them off?
        expandRadius: PD.Numeric(5, { min: 1, max: 10, step: 1 }),
        targetParams: PD.Group(reprParams, {
            label: 'Target',
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

type StructureFocusRepresentationProps = PD.ValuesFor<ReturnType<typeof StructureFocusRepresentationParams>>

export enum StructureFocusRepresentationTags {
    TargetSel = 'structure-focus-target-sel',
    TargetRepr = 'structure-focus-target-repr',
    SurrSel = 'structure-focus-surr-sel',
    SurrRepr = 'structure-focus-surr-repr',
    SurrNciRepr = 'structure-focus-surr-nci-repr'
}

const TagSet: Set<StructureFocusRepresentationTags> = new Set([StructureFocusRepresentationTags.TargetSel, StructureFocusRepresentationTags.TargetRepr, StructureFocusRepresentationTags.SurrSel, StructureFocusRepresentationTags.SurrRepr, StructureFocusRepresentationTags.SurrNciRepr])

export class StructureFocusRepresentationBehavior extends PluginBehavior.WithSubscribers<StructureFocusRepresentationProps> {
    private get surrLabel() { return `[Focus] Surroundings (${this.params.expandRadius} Å)`; }

    private ensureShape(cell: StateObjectCell<PluginStateObject.Molecule.Structure>) {
        const state = this.plugin.state.data, tree = state.tree;
        const builder = state.build();
        const refs = StateSelection.findUniqueTagsInSubtree(tree, cell.transform.ref, TagSet);

        // Selections
        if (!refs[StructureFocusRepresentationTags.TargetSel]) {
            refs[StructureFocusRepresentationTags.TargetSel] = builder
                .to(cell)
                .apply(StateTransforms.Model.StructureSelectionFromBundle,
                    { bundle: StructureElement.Bundle.Empty, label: '[Focus] Target' }, { tags: StructureFocusRepresentationTags.TargetSel }).ref;
        }

        if (!refs[StructureFocusRepresentationTags.SurrSel]) {
            refs[StructureFocusRepresentationTags.SurrSel] = builder
                .to(cell)
                .apply(StateTransforms.Model.StructureSelectionFromExpression,
                    { expression: MS.struct.generator.empty(), label: this.surrLabel }, { tags: StructureFocusRepresentationTags.SurrSel }).ref;
        }

        // Representations
        if (!refs[StructureFocusRepresentationTags.TargetRepr]) {
            refs[StructureFocusRepresentationTags.TargetRepr] = builder
                .to(refs[StructureFocusRepresentationTags.TargetSel]!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.params.targetParams, { tags: StructureFocusRepresentationTags.TargetRepr }).ref;
        }

        if (!refs[StructureFocusRepresentationTags.SurrRepr]) {
            refs[StructureFocusRepresentationTags.SurrRepr] = builder
                .to(refs[StructureFocusRepresentationTags.SurrSel]!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.params.surroundingsParams, { tags: StructureFocusRepresentationTags.SurrRepr }).ref;
        }

        if (!refs[StructureFocusRepresentationTags.SurrNciRepr]) {
            refs[StructureFocusRepresentationTags.SurrNciRepr] = builder
                .to(refs[StructureFocusRepresentationTags.SurrSel]!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.params.nciParams, { tags: StructureFocusRepresentationTags.SurrNciRepr }).ref;
        }

        return { state, builder, refs };
    }

    private clear(root: StateTransform.Ref) {
        const state = this.plugin.state.data;

        const foci = state.select(StateSelection.Generators.byRef(root).subtree().withTag(StructureFocusRepresentationTags.TargetSel));
        const surrs = state.select(StateSelection.Generators.byRef(root).subtree().withTag(StructureFocusRepresentationTags.SurrSel));
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

    private focus(loci: StructureElement.Loci) {
        const parent = this.plugin.helpers.substructureParent.get(loci.structure);
        if (!parent || !parent.obj) return;

        const residueLoci = StructureElement.Loci.extendToWholeResidues(StructureElement.Loci.remap(loci, parent.obj!.data))
        const residueBundle = StructureElement.Bundle.fromLoci(residueLoci)

        const surroundings = MS.struct.modifier.includeSurroundings({
            0: StructureElement.Bundle.toExpression(residueBundle),
            radius: this.params.expandRadius,
            'as-whole-residues': true
        });

        const { state, builder, refs } = this.ensureShape(parent);

        builder.to(refs[StructureFocusRepresentationTags.TargetSel]!).update(StateTransforms.Model.StructureSelectionFromBundle, old => ({ ...old, bundle: residueBundle }));
        builder.to(refs[StructureFocusRepresentationTags.SurrSel]!).update(StateTransforms.Model.StructureSelectionFromExpression, old => ({ ...old, expression: surroundings, label: this.surrLabel }));

        PluginCommands.State.Update(this.plugin, { state, tree: builder, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });
    }

    register(ref: string): void {
        this.subscribeObservable(this.plugin.managers.structure.focus.events.changed, () => {
            const entry = this.plugin.managers.structure.focus.current
            if (entry) this.focus(entry.loci)
            else this.clear(StateTransform.RootRef)
        });
    }

    async update(params: StructureFocusRepresentationProps) {
        let oldRadius = this.params.expandRadius;
        this.params = params;

        const state = this.plugin.state.data;
        const builder = state.build();

        const all = StateSelection.Generators.root.subtree();
        for (const repr of state.select(all.withTag(StructureFocusRepresentationTags.TargetRepr))) {
            builder.to(repr).update(this.params.targetParams);
        }
        for (const repr of state.select(all.withTag(StructureFocusRepresentationTags.SurrRepr))) {
            builder.to(repr).update(this.params.surroundingsParams);
        }
        for (const repr of state.select(all.withTag(StructureFocusRepresentationTags.SurrNciRepr))) {
            builder.to(repr).update(this.params.nciParams);
        }

        await PluginCommands.State.Update(this.plugin, { state, tree: builder, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });

        // TODO: update properly
        if (params.expandRadius !== oldRadius) await this.clear(StateTransform.RootRef);

        return true;
    }
}

export const StructureFocusRepresentation = PluginBehavior.create({
    name: 'create-structure-focus-representation',
    display: { name: 'Structure Focus Representation' },
    category: 'interaction',
    ctor: StructureFocusRepresentationBehavior,
    params: (_, plugin) => StructureFocusRepresentationParams(plugin)
});