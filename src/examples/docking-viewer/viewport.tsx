/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../../mol-plugin-ui/base';
import { Viewport, ViewportControls } from '../../mol-plugin-ui/viewport';
import { BackgroundTaskProgress } from '../../mol-plugin-ui/task';
import { LociLabels } from '../../mol-plugin-ui/controls';
import { Toasts } from '../../mol-plugin-ui/toast';
import { Button } from '../../mol-plugin-ui/controls/common';
import { StructureRepresentationPresetProvider, presetStaticComponent } from '../../mol-plugin-state/builder/structure/representation-preset';
import { StateObjectRef } from '../../mol-state';
import { StructureSelectionQueries, StructureSelectionQuery } from '../../mol-plugin-state/helpers/structure-selection-query';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { InteractionsRepresentationProvider } from '../../mol-model-props/computed/representations/interactions';
import { InteractionTypeColorThemeProvider } from '../../mol-model-props/computed/themes/interaction-type';
import { compile } from '../../mol-script/runtime/query/compiler';
import { StructureSelection, QueryContext, Structure } from '../../mol-model/structure';
import { PluginCommands } from '../../mol-plugin/commands';
import { PluginContext } from '../../mol-plugin/context';

function shinyStyle(plugin: PluginContext) {
    return PluginCommands.Canvas3D.SetSettings(plugin, { settings: {
        renderer: {
            ...plugin.canvas3d!.props.renderer,
            style: { name: 'plastic', params: {} },
        },
        postprocessing: {
            ...plugin.canvas3d!.props.postprocessing,
            occlusion: { name: 'off', params: {} },
            outline: { name: 'off', params: {} }
        }
    } });
}

function occlusionStyle(plugin: PluginContext) {
    return PluginCommands.Canvas3D.SetSettings(plugin, { settings: {
        renderer: {
            ...plugin.canvas3d!.props.renderer,
            style: { name: 'flat', params: {} }
        },
        postprocessing: {
            ...plugin.canvas3d!.props.postprocessing,
            occlusion: { name: 'on', params: {
                kernelSize: 8,
                bias: 0.8,
                radius: 64
            } },
            outline: { name: 'on', params: {
                scale: 1.0,
                threshold: 0.8
            } }
        }
    } });
}

const ligandPlusSurroundings = StructureSelectionQuery('Surrounding Residues (5 \u212B) of Ligand plus Ligand itself', MS.struct.modifier.union([
    MS.struct.modifier.includeSurroundings({
        0: StructureSelectionQueries.ligand.expression,
        radius: 5,
        'as-whole-residues': true
    })
]));

const ligandSurroundings = StructureSelectionQuery('Surrounding Residues (5 \u212B) of Ligand', MS.struct.modifier.union([
    MS.struct.modifier.exceptBy({
        0: ligandPlusSurroundings.expression,
        by: StructureSelectionQueries.ligand.expression
    })
]));

const PresetParams = {
    ...StructureRepresentationPresetProvider.CommonParams,
};

export const StructurePreset = StructureRepresentationPresetProvider({
    id: 'preset-structure',
    display: { name: 'Structure' },
    params: () => PresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const components = {
            ligand: await presetStaticComponent(plugin, structureCell, 'ligand'),
            polymer: await presetStaticComponent(plugin, structureCell, 'polymer'),
        };

        const { update, builder, typeParams, color } = StructureRepresentationPresetProvider.reprBuilder(plugin, params);
        const representations = {
            ligand: builder.buildRepresentation(update, components.ligand, { type: 'ball-and-stick', typeParams: { ...typeParams, sizeFactor: 0.26 }, color }, { tag: 'ligand' }),
            polymer: builder.buildRepresentation(update, components.polymer, { type: 'cartoon', typeParams: { ...typeParams }, color }, { tag: 'polymer' }),
        };

        await update.commit({ revertOnError: true });
        await shinyStyle(plugin);
        plugin.managers.interactivity.setProps({ granularity: 'residue' });

        return { components, representations };
    }
});

export const IllustrativePreset = StructureRepresentationPresetProvider({
    id: 'preset-illustrative',
    display: { name: 'Illustrative' },
    params: () => PresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const components = {
            all: await presetStaticComponent(plugin, structureCell, 'all')
        };

        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, params);
        const representations = {
            all: builder.buildRepresentation(update, components.all, { type: 'spacefill', typeParams: { ...typeParams }, color: 'illustrative' }, { tag: 'all' }),
        };

        await update.commit({ revertOnError: true });
        await occlusionStyle(plugin);
        plugin.managers.interactivity.setProps({ granularity: 'residue' });

        return { components, representations };
    }
});

const PocketPreset = StructureRepresentationPresetProvider({
    id: 'preset-pocket',
    display: { name: 'Pocket' },
    params: () => PresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const structure = structureCell?.obj?.data;
        if (!structureCell || !structure) return {};

        const components = {
            ligand: await presetStaticComponent(plugin, structureCell, 'ligand'),
            surroundings: await plugin.builders.structure.tryCreateComponentFromSelection(structureCell, ligandSurroundings, `surroundings`),
        };

        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, params);
        const representations = {
            ligand: builder.buildRepresentation(update, components.ligand, { type: 'ball-and-stick', typeParams: { ...typeParams, sizeFactor: 0.26 }, color: 'partial-charge' }, { tag: 'ligand' }),
            surroundings: builder.buildRepresentation(update, components.surroundings, { type: 'molecular-surface', typeParams: { ...typeParams, includeParent: true, quality: 'custom', resolution: 0.2, doubleSided: true }, color: 'partial-charge' }, { tag: 'surroundings' }),
        };

        await update.commit({ revertOnError: true });
        await shinyStyle(plugin);
        plugin.managers.interactivity.setProps({ granularity: 'element' });

        const compiled = compile<StructureSelection>(StructureSelectionQueries.ligand.expression);
        const result = compiled(new QueryContext(structure));
        const selection = StructureSelection.unionStructure(result);
        plugin.managers.camera.focusLoci(Structure.toStructureElementLoci(selection));

        return { components, representations };
    }
});

const InteractionsPreset = StructureRepresentationPresetProvider({
    id: 'preset-interactions',
    display: { name: 'Interactions' },
    params: () => PresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        const structure = structureCell?.obj?.data;
        if (!structureCell || !structure) return {};

        const components = {
            ligand: await presetStaticComponent(plugin, structureCell, 'ligand'),
            selection: await plugin.builders.structure.tryCreateComponentFromSelection(structureCell, ligandPlusSurroundings, `selection`)
        };

        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, params);
        const representations = {
            ligand: builder.buildRepresentation(update, components.ligand, { type: 'ball-and-stick', typeParams: { ...typeParams, sizeFactor: 0.26 }, color: 'partial-charge' }, { tag: 'ligand' }),
            ballAndStick: builder.buildRepresentation(update, components.selection, { type: 'ball-and-stick', typeParams: { ...typeParams, sizeFactor: 0.1, sizeAspectRatio: 1 }, color: 'partial-charge' }, { tag: 'ball-and-stick' }),
            interactions: builder.buildRepresentation(update, components.selection, { type: InteractionsRepresentationProvider, typeParams: { ...typeParams }, color: InteractionTypeColorThemeProvider }, { tag: 'interactions' }),
        };

        await update.commit({ revertOnError: true });
        await shinyStyle(plugin);
        plugin.managers.interactivity.setProps({ granularity: 'element' });

        const compiled = compile<StructureSelection>(StructureSelectionQueries.ligand.expression);
        const result = compiled(new QueryContext(structure));
        const selection = StructureSelection.unionStructure(result);
        plugin.managers.camera.focusLoci(Structure.toStructureElementLoci(selection));

        return { components, representations };
    }
});

export class ViewportComponent extends PluginUIComponent {
    structurePreset = () => {
        this.plugin.managers.structure.component.applyPreset(
            this.plugin.managers.structure.hierarchy.selection.structures,
            StructurePreset
        );
    }

    illustrativePreset = () => {
        this.plugin.managers.structure.component.applyPreset(
            this.plugin.managers.structure.hierarchy.selection.structures,
            IllustrativePreset
        );
    }

    pocketPreset = () => {
        this.plugin.managers.structure.component.applyPreset(
            this.plugin.managers.structure.hierarchy.selection.structures,
            PocketPreset
        );
    }

    interactionsPreset = () => {
        this.plugin.managers.structure.component.applyPreset(
            this.plugin.managers.structure.hierarchy.selection.structures,
            InteractionsPreset
        );
    }

    render() {
        const VPControls = this.plugin.spec.components?.viewport?.controls || ViewportControls;

        return <>
            <Viewport />
            <div className='msp-viewport-top-left-controls'>
                <div style={{ marginBottom: '4px' }}>
                    <Button onClick={this.structurePreset} >Structure</Button>
                </div>
                <div style={{ marginBottom: '4px' }}>
                    <Button onClick={this.illustrativePreset}>Illustrative</Button>
                </div>
                <div style={{ marginBottom: '4px' }}>
                    <Button onClick={this.pocketPreset}>Pocket</Button>
                </div>
                <div style={{ marginBottom: '4px' }}>
                    <Button onClick={this.interactionsPreset}>Interactions</Button>
                </div>
            </div>
            <VPControls />
            <BackgroundTaskProgress />
            <div className='msp-highlight-toast-wrapper'>
                <LociLabels />
                <Toasts />
            </div>
        </>;
    }
}