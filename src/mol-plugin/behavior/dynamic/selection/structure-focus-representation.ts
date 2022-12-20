/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { Material } from '../../../../mol-util/material';
import { Clip } from '../../../../mol-util/clip';

const StructureFocusRepresentationParams = (plugin: PluginContext) => {
    const reprParams = StateTransforms.Representation.StructureRepresentation3D.definition.params!(void 0, plugin) as PD.Params;
    return {
        expandRadius: PD.Numeric(5, { min: 1, max: 10, step: 1 }),
        targetParams: PD.Group(reprParams, {
            label: 'Target',
            customDefault: createStructureRepresentationParams(plugin, void 0, {
                type: 'ball-and-stick',
                size: 'physical',
                typeParams: { sizeFactor: 0.22, sizeAspectRatio: 0.73, adjustCylinderLength: true, xrayShaded: true, aromaticBonds: false, multipleBonds: 'off', excludeTypes: ['hydrogen-bond', 'metal-coordination'] },
            })
        }),
        surroundingsParams: PD.Group(reprParams, {
            label: 'Surroundings',
            customDefault: createStructureRepresentationParams(plugin, void 0, {
                type: 'ball-and-stick',
                size: 'physical',
                typeParams: { sizeFactor: 0.16, excludeTypes: ['hydrogen-bond', 'metal-coordination'] }
            })
        }),
        nciParams: PD.Group(reprParams, {
            label: 'Non-covalent Int.',
            customDefault: createStructureRepresentationParams(plugin, void 0, {
                type: InteractionsRepresentationProvider,
                color: InteractionTypeColorThemeProvider,
                size: SizeTheme.BuiltIn.uniform
            })
        }),
        components: PD.MultiSelect(FocusComponents, PD.arrayToOptions(FocusComponents)),
        excludeTargetFromSurroundings: PD.Boolean(false, { label: 'Exclude Target', description: 'Exclude the focus "target" from the surroudings component.' }),
        ignoreHydrogens: PD.Boolean(false),
        ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
        ignoreLight: PD.Boolean(false),
        material: Material.getParam(),
        clip: PD.Group(Clip.Params),
    };
};

const FocusComponents = ['target' as const, 'surroundings' as const, 'interactions' as const];

type StructureFocusRepresentationProps = PD.ValuesFor<ReturnType<typeof StructureFocusRepresentationParams>>

export enum StructureFocusRepresentationTags {
    TargetSel = 'structure-focus-target-sel',
    TargetRepr = 'structure-focus-target-repr',
    SurrSel = 'structure-focus-surr-sel',
    SurrRepr = 'structure-focus-surr-repr',
    SurrNciRepr = 'structure-focus-surr-nci-repr'
}

const TagSet: Set<StructureFocusRepresentationTags> = new Set([StructureFocusRepresentationTags.TargetSel, StructureFocusRepresentationTags.TargetRepr, StructureFocusRepresentationTags.SurrSel, StructureFocusRepresentationTags.SurrRepr, StructureFocusRepresentationTags.SurrNciRepr]);

class StructureFocusRepresentationBehavior extends PluginBehavior.WithSubscribers<StructureFocusRepresentationProps> {
    private get surrLabel() { return `[Focus] Surroundings (${this.params.expandRadius} Å)`; }

    private getReprParams(reprParams: PD.Values<PD.Params>) {
        return {
            ...reprParams,
            type: {
                name: reprParams.type.name,
                params: { ...reprParams.type.params, ignoreHydrogens: this.params.ignoreHydrogens, ignoreHydrogensVariant: this.params.ignoreHydrogensVariant, ignoreLight: this.params.ignoreLight, material: this.params.material, clip: this.params.clip }
            }
        };
    }

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

        const components = this.params.components;

        // Representations
        if (components.indexOf('target') >= 0 && !refs[StructureFocusRepresentationTags.TargetRepr]) {
            refs[StructureFocusRepresentationTags.TargetRepr] = builder
                .to(refs[StructureFocusRepresentationTags.TargetSel]!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.getReprParams(this.params.targetParams), { tags: StructureFocusRepresentationTags.TargetRepr }).ref;
        }

        if (components.indexOf('surroundings') >= 0 && !refs[StructureFocusRepresentationTags.SurrRepr]) {
            refs[StructureFocusRepresentationTags.SurrRepr] = builder
                .to(refs[StructureFocusRepresentationTags.SurrSel]!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.getReprParams(this.params.surroundingsParams), { tags: StructureFocusRepresentationTags.SurrRepr }).ref;
        }

        if (components.indexOf('interactions') >= 0 && !refs[StructureFocusRepresentationTags.SurrNciRepr] && cell.obj && InteractionsRepresentationProvider.isApplicable(cell.obj?.data)) {
            refs[StructureFocusRepresentationTags.SurrNciRepr] = builder
                .to(refs[StructureFocusRepresentationTags.SurrSel]!)
                .apply(StateTransforms.Representation.StructureRepresentation3D, this.getReprParams(this.params.nciParams), { tags: StructureFocusRepresentationTags.SurrNciRepr }).ref;
        }

        return { state, builder, refs };
    }

    private clear(root: StateTransform.Ref) {
        const state = this.plugin.state.data;
        this.currentSource = void 0;

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

    private currentSource: StructureElement.Loci | undefined = void 0;
    private async focus(sourceLoci: StructureElement.Loci) {
        const parent = this.plugin.helpers.substructureParent.get(sourceLoci.structure);
        if (!parent || !parent.obj) return;

        this.currentSource = sourceLoci;
        const loci = StructureElement.Loci.remap(sourceLoci, parent.obj!.data);

        const residueLoci = StructureElement.Loci.extendToWholeResidues(loci);
        const residueBundle = StructureElement.Bundle.fromLoci(residueLoci);

        const target = StructureElement.Bundle.toExpression(residueBundle);
        let surroundings = MS.struct.modifier.includeSurroundings({
            0: target,
            radius: this.params.expandRadius,
            'as-whole-residues': true
        });

        if (this.params.excludeTargetFromSurroundings) {
            surroundings = MS.struct.modifier.exceptBy({
                0: surroundings,
                by: target
            });
        }

        const { state, builder, refs } = this.ensureShape(parent);

        builder.to(refs[StructureFocusRepresentationTags.TargetSel]!).update(StateTransforms.Model.StructureSelectionFromBundle, old => ({ ...old, bundle: residueBundle }));
        builder.to(refs[StructureFocusRepresentationTags.SurrSel]!).update(StateTransforms.Model.StructureSelectionFromExpression, old => ({ ...old, expression: surroundings, label: this.surrLabel }));

        await PluginCommands.State.Update(this.plugin, { state, tree: builder, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });
    }

    register(ref: string): void {
        this.subscribeObservable(this.plugin.managers.structure.focus.behaviors.current, (entry) => {
            if (entry) this.focus(entry.loci);
            else this.clear(StateTransform.RootRef);
        });
    }

    async update(params: StructureFocusRepresentationProps) {
        const old = this.params;
        this.params = params;

        if (old.excludeTargetFromSurroundings !== params.excludeTargetFromSurroundings) {
            if (this.currentSource) {
                this.focus(this.currentSource);
            }
            return true;
        }

        const state = this.plugin.state.data;
        const builder = state.build();

        const all = StateSelection.Generators.root.subtree();

        const components = this.params.components;

        // TODO: create component if previously didnt exist
        let hasComponent = components.indexOf('target') >= 0;
        for (const repr of state.select(all.withTag(StructureFocusRepresentationTags.TargetRepr))) {
            if (!hasComponent) builder.delete(repr.transform.ref);
            else builder.to(repr).update(this.getReprParams(this.params.targetParams));
        }

        hasComponent = components.indexOf('surroundings') >= 0;
        for (const repr of state.select(all.withTag(StructureFocusRepresentationTags.SurrRepr))) {
            if (!hasComponent) builder.delete(repr.transform.ref);
            else builder.to(repr).update(this.getReprParams(this.params.surroundingsParams));
        }

        hasComponent = components.indexOf('interactions') >= 0;
        for (const repr of state.select(all.withTag(StructureFocusRepresentationTags.SurrNciRepr))) {
            if (!hasComponent) builder.delete(repr.transform.ref);
            else builder.to(repr).update(this.getReprParams(this.params.nciParams));
        }

        await PluginCommands.State.Update(this.plugin, { state, tree: builder, options: { doNotLogTiming: true, doNotUpdateCurrent: true } });

        if (params.expandRadius !== old.expandRadius) {
            if (this.currentSource) {
                this.focus(this.currentSource);
            }
            return true;
        }

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