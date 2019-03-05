/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from 'mol-model/structure';
import { PluginBehavior } from 'mol-plugin/behavior';
import { PluginCommands } from 'mol-plugin/command';
import { PluginContext } from 'mol-plugin/context';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { StateTransforms } from 'mol-plugin/state/transforms';
import { StructureRepresentation3DHelpers } from 'mol-plugin/state/transforms/representation';
import { BuiltInStructureRepresentations } from 'mol-repr/structure/registry';
import { MolScriptBuilder as MS } from 'mol-script/language/builder';
import { StateObjectCell, StateSelection } from 'mol-state';
import { BuiltInColorThemes } from 'mol-theme/color';
import { BuiltInSizeThemes } from 'mol-theme/size';
import { ColorNames } from 'mol-util/color/tables';
import { ButtonsType } from 'mol-util/input/input-observer';

type Params = { }

enum Tags {
    Group = 'structure-interaction-group',
    ResidueSel = 'structure-interaction-residue-sel',
    ResidueRepr = 'structure-interaction-residue-repr',
    SurrSel = 'structure-interaction-surr-sel',
    SurrRepr = 'structure-interaction-surr-repr'
}

const TagSet: Set<Tags> = new Set([Tags.Group, Tags.ResidueSel, Tags.ResidueRepr, Tags.SurrSel, Tags.SurrRepr])

export class StructureRepresentationInteractionBehavior extends PluginBehavior.WithSubscribers<Params> {

    private createResVisualParams(s: Structure) {
        return StructureRepresentation3DHelpers.createParams(this.plugin, s, {
            repr: BuiltInStructureRepresentations['ball-and-stick'],
            size: [BuiltInSizeThemes.uniform, () => ({ value: 0.85 } )]
        });
    }

    private createSurVisualParams(s: Structure) {
        return StructureRepresentation3DHelpers.createParams(this.plugin, s, {
            repr: BuiltInStructureRepresentations['ball-and-stick'],
            color: [BuiltInColorThemes.uniform, () => ({ value: ColorNames.gray })],
            size: [BuiltInSizeThemes.uniform, () => ({ value: 0.33 } )]
        });
    }

    private ensureShape(cell: StateObjectCell<PluginStateObject.Molecule.Structure>) {
        const state = this.plugin.state.dataState, tree = state.tree;
        const builder = state.build();
        const refs = StateSelection.findUniqueTagsInSubtree(tree, cell.transform.ref, TagSet);

        if (!refs['structure-interaction-group']) {
            refs['structure-interaction-group'] = builder.to(cell.transform.ref).group(StateTransforms.Misc.CreateGroup,
                { label: 'Current Interaction' }, { props: { tag: Tags.Group } }).ref;
        }

        // Selections
        if (!refs[Tags.ResidueSel]) {
            refs[Tags.ResidueSel] = builder.to(refs['structure-interaction-group']).apply(StateTransforms.Model.StructureSelection,
                { query: { } as any, label: 'Residue' }, { props: { tag: Tags.ResidueSel } }).ref;
        }

        if (!refs[Tags.SurrSel]) {
            refs[Tags.SurrSel] = builder.to(refs['structure-interaction-group']).apply(StateTransforms.Model.StructureSelection,
                { query: { } as any, label: 'Surroundings' }, { props: { tag: Tags.SurrSel } }).ref;
        }

        // Representations
        // TODO: ability to customize how it looks in the behavior params
        if (!refs[Tags.ResidueRepr]) {
            refs[Tags.ResidueRepr] = builder.to(refs['structure-interaction-residue-sel']!).apply(StateTransforms.Representation.StructureRepresentation3D,
                this.createResVisualParams(cell.obj!.data), { props: { tag: Tags.ResidueRepr } }).ref;
        }

        if (!refs[Tags.SurrRepr]) {
            refs[Tags.SurrRepr] = builder.to(refs['structure-interaction-surr-sel']!).apply(StateTransforms.Representation.StructureRepresentation3D,
                this.createSurVisualParams(cell.obj!.data), { props: { tag: Tags.SurrRepr } }).ref;
        }

        return { state, builder, refs };
    }

    private clear() {
        const state = this.plugin.state.dataState;
        const groups = state.select(StateSelection.Generators.root.subtree().filter(o => o.transform.props.tag === Tags.Group));
        if (groups.length === 0) return;

        const update = state.build();
        const query = MS.struct.generator.empty();
        for (const g of groups) {
            // TODO: update props of the group node to ghost

            const res = StateSelection.findTagInSubtree(state.tree, g.transform.ref, Tags.ResidueSel);
            const surr = StateSelection.findTagInSubtree(state.tree, g.transform.ref, Tags.SurrSel);
            if (res) update.to(res).update(StateTransforms.Model.StructureSelection, old => ({ ...old, query }));
            if (surr) update.to(surr).update(StateTransforms.Model.StructureSelection, old => ({ ...old, query }));

            // update.delete(g.transform.ref);
        }

        PluginCommands.State.Update.dispatch(this.plugin, { state, tree: update, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });
    }

    register(ref: string): void {
        // this.ref = ref;

        this.subscribeObservable(this.plugin.behaviors.canvas3d.click, ({ current, buttons, modifiers }) => {
            if (buttons !== ButtonsType.Flag.Secondary) return;

            if (current.loci.kind === 'empty-loci') {
                if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                    this.clear();
                    return;
                }
            }

            // TODO: support link loci as well?
            if (!StructureElement.isLoci(current.loci)) return;

            const parent = this.plugin.helpers.substructureParent.get(current.loci.structure);
            if (!parent || !parent.obj) return;

            const core = MS.struct.modifier.wholeResidues([
                StructureElement.Loci.toScriptExpression(current.loci)
            ]);

            const surroundings = MS.struct.modifier.exceptBy({
                0: MS.struct.modifier.includeSurroundings({
                    0: core,
                    radius: 5,
                    'as-whole-residues': true
                }),
                by: core
            });

            const { state, builder, refs } = this.ensureShape(parent);

            builder.to(refs[Tags.ResidueSel]!).update(StateTransforms.Model.StructureSelection, old => ({ ...old, query: core }));
            builder.to(refs[Tags.SurrSel]!).update(StateTransforms.Model.StructureSelection, old => ({ ...old, query: surroundings }));

            PluginCommands.State.Update.dispatch(this.plugin, { state, tree: builder, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });
        });
    }

    async update(params: Params) {
        return false;
    }

    constructor(public plugin: PluginContext) {
        super(plugin);
    }
}

export const StructureRepresentationInteraction = PluginBehavior.create({
    name: 'create-structure-representation-interaction',
    display: { name: 'Structure Representation Interaction' },
    category: 'interaction',
    ctor: StructureRepresentationInteractionBehavior
});