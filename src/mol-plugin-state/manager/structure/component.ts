/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement, StructureSelection } from '../../../mol-model/structure';
import { structureAreIntersecting, structureSubtract, structureUnion } from '../../../mol-model/structure/query/utils/structure-set';
import { PluginContext } from '../../../mol-plugin/context';
import { StateBuilder } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { UUID } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentationProvider } from '../../builder/structure/provider';
import { StructureComponentParams } from '../../helpers/structure-component';
import { setStructureOverpaint } from '../../helpers/structure-overpaint';
import { StructureSelectionQuery, StructureSelectionQueryOptions } from '../../helpers/structure-selection-query';
import { HierarchyRef, StructureComponentRef, StructureRef, StructureRepresentationRef } from './hierarchy-state';

export { StructureComponentManager };

class StructureComponentManager {
    applyPreset<P = any, S = {}>(structures: StructureRef[], provider: StructureRepresentationProvider<P, S>, params?: P): Promise<any>  {
        return this.plugin.runTask(this.dataState.transaction(async () => {
            await this.clearComponents(structures);
            for (const s of structures) {
                await this.plugin.builders.structure.representation.structurePreset(s.cell, provider, params);
            }
        }));
    }

    clear(structures: StructureRef[]) {
        return this.clearComponents(structures);
    }

    removeRepresentations(components: StructureComponentRef[], pivot: StructureRepresentationRef) {
        if (components.length === 0) return;
        const index = components[0].representations.indexOf(pivot);
        if (index < 0) return;

        const toRemove: HierarchyRef[] = [];
        for (const c of components) {
            if (index >= c.representations.length) continue;
            toRemove.push(c.representations[index]);
        }
        return this.plugin.managers.structure.hierarchy.remove(toRemove);
    }

    modify(action: StructureComponentManager.ModifyAction, structures?: ReadonlyArray<StructureRef>) {        
        return this.plugin.runTask(this.dataState.transaction(async () => {
            if (!structures) structures = this.plugin.managers.structure.hierarchy.state.currentStructures;
            if (structures.length === 0) return;

            switch (action.kind) {
                case 'add': await this.modifyAdd(action, structures); break;
                case 'merge': await this.modifyMerge(action, structures); break;
                case 'subtract': await this.modifySubtract(action, structures); break;
                case 'color': await this.modifyColor(action, structures); break;
            }
        }))
    }

    private async modifyAdd(params: StructureComponentManager.ModifyActionAdd, structures: ReadonlyArray<StructureRef>) {
        const componentKey = UUID.create22();
        for (const s of structures) {
            const component = await this.plugin.builders.structure.tryCreateQueryComponent({ 
                structure: s.cell,
                query: params.selection,
                key: componentKey,
                label: params.label,
            });
            if (params.representation === 'none' || !component) continue;
            await this.plugin.builders.structure.representation.addRepresentation(component, {
                repr: this.plugin.structureRepresentation.registry.get(params.representation)
            });
        }
    }

    private updateComponent(builder: StateBuilder.Root, component: StructureComponentRef, by: Structure, action: 'union' | 'subtract') {
        const structure = component.cell.obj?.data;
        if (!structure) return;
        if (!structureAreIntersecting(structure, by)) return;

        const parent = component.structure.cell.obj?.data!;
        const modified = action === 'union' ? structureUnion(parent, [structure, by]) : structureSubtract(structure, by);

        if (modified.elementCount === 0) {
            builder.delete(component.cell.transform.ref);
        } else {
            const bundle = StructureElement.Bundle.fromLoci(StructureSelection.toLociWithSourceUnits(StructureSelection.Singletons(parent, modified)));
            const params: StructureComponentParams = {
                type: { name: 'bundle', params: bundle },
                nullIfEmpty: true,
                label: component.cell.obj?.label!
            };
            builder.to(component.cell).update(params)
        }
    }

    private async modifyMerge(params: StructureComponentManager.ModifyActionMerge, structures: ReadonlyArray<StructureRef>) {
        return this.plugin.runTask(Task.create('Merge', async taskCtx => {
            const b = this.dataState.build();        
            for (const s of structures) {
                const by = await StructureSelectionQuery.getStructure(this.plugin, taskCtx, params.selection, s.cell.obj?.data!);
                for (const c of s.components) {
                    if (params.componentKey !== 'intersecting' && params.componentKey !== c.key) continue;
                    this.updateComponent(b, c, by, 'union');
                }
            }
            await this.dataState.updateTree(b).runInContext(taskCtx);
        }));
    }

    private async modifySubtract(params: StructureComponentManager.ModifyActionSubtract, structures: ReadonlyArray<StructureRef>) {
        return this.plugin.runTask(Task.create('Subtract', async taskCtx => {
            const b = this.dataState.build();        
            for (const s of structures) {
                const by = await StructureSelectionQuery.getStructure(this.plugin, taskCtx, params.selection, s.cell.obj?.data!);
                for (const c of s.components) {
                    if (params.componentKey !== 'intersecting' && params.componentKey !== c.key) continue;
                    this.updateComponent(b, c, by, 'subtract');
                }
            }
            await this.dataState.updateTree(b).runInContext(taskCtx);
        }));
    }

    private async modifyColor(params: StructureComponentManager.ModifyActionColor, structures: ReadonlyArray<StructureRef>) {
        const getLoci = (s: Structure) => this.plugin.managers.structure.selection.getLoci(s);
        for (const s of structures) {
            await setStructureOverpaint(this.plugin, s.components, params.action.name === 'color' ? params.action.params : -1, getLoci);
        }
    }

    private get dataState() {
        return this.plugin.state.dataState;
    }

    private clearComponents(structures: StructureRef[]) {
        const deletes = this.dataState.build();
        for (const s of structures) {
            for (const c of s.components) {
                deletes.delete(c.cell.transform.ref);
            }
            if (s.currentFocus) {
                if (s.currentFocus.focus) deletes.delete(s.currentFocus.focus.cell.transform.ref);
                if (s.currentFocus.surroundings) deletes.delete(s.currentFocus.surroundings.cell.transform.ref);
            }
        }
        return this.plugin.runTask(this.dataState.updateTree(deletes));
    }

    constructor(public plugin: PluginContext) {

    }
}

namespace StructureComponentManager {
    export type ActionType = 'add' | 'merge' | 'subtract' | 'color'

    const SelectionParam = PD.Select(StructureSelectionQueryOptions[1][0], StructureSelectionQueryOptions)

    function getComponentsOptions(plugin: PluginContext, custom: [string, string][], label?: string) {
        const types = [
            ...custom,
            ...plugin.managers.structure.hierarchy.componentGroups.map(g => [g[0].key!, g[0].cell.obj?.label])
        ] as [string, string][];
        return PD.Select(types[0][0], types, { label });
    }

    function getRepresentationTypes(plugin: PluginContext, pivot: StructureRef | undefined, custom: [string, string][], label?: string) {
        const types = [
            ...custom,
            ...(pivot?.cell.obj?.data
                ? plugin.structureRepresentation.registry.getApplicableTypes(pivot.cell.obj?.data!)
                : plugin.structureRepresentation.registry.types)
        ] as [string, string][];
        return PD.Select(types[0][0], types, { label });
    }

    export function getActionParams(plugin: PluginContext, action: ActionType) {
        switch (action) {
            case 'add': 
                return {
                    kind: PD.Value<ActionType>(action, { isHidden: true }),
                    selection: SelectionParam,
                    label: PD.Text(''),
                    representation: getRepresentationTypes(plugin, plugin.managers.structure.hierarchy.state.currentStructures[0], [['none', '< None >']])
                };
            case 'merge':
            case 'subtract':
                return {
                    kind: PD.Value<ActionType>(action, { isHidden: true }),
                    selection: SelectionParam,
                    componentKey: getComponentsOptions(plugin, [['intersecting', '< Intersecting >']], 'Target')
                };
            case 'color':
                // TODO: ability to reset
                return {
                    kind: PD.Value<ActionType>(action, { isHidden: true }),
                    action: PD.MappedStatic('color', {
                        color: PD.Color(ColorNames.black),
                        reset: PD.EmptyGroup()
                    }),
                    // TODO: filter by representation type
                    // representation: getRepresentationTypes(plugin, void 0, [['all', '< All >']])
                };
        }
    }

    export type ModifyActionAdd = { kind: 'add', selection: StructureSelectionQuery, label: string, representation: string }
    export type ModifyActionMerge = { kind: 'merge', selection: StructureSelectionQuery, componentKey: 'intersecting' | string }
    export type ModifyActionSubtract = { kind: 'subtract', selection: StructureSelectionQuery, componentKey: 'intersecting' | string }
    export type ModifyActionColor = { kind: 'color', action: { name: 'color', params: Color } | { name: 'reset', params: any } } //, representationType?: string }
    
    export type ModifyAction = ModifyActionAdd | ModifyActionMerge | ModifyActionSubtract | ModifyActionColor
}
