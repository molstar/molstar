/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { VisualQualityOptions } from '../../../mol-geo/geometry/base';
import { InteractionsProvider } from '../../../mol-model-props/computed/interactions';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { structureAreIntersecting, structureSubtract, structureUnion, structureIntersect } from '../../../mol-model/structure/query/utils/structure-set';
import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginContext } from '../../../mol-plugin/context';
import { StateBuilder, StateTransformer } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { UUID } from '../../../mol-util';
import { arraySetAdd } from '../../../mol-util/array';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentationProvider } from '../../builder/structure/provider';
import { PluginComponent } from '../../component';
import { StructureComponentParams } from '../../helpers/structure-component';
import { clearStructureOverpaint, setStructureOverpaint } from '../../helpers/structure-overpaint';
import { StructureSelectionQueries, StructureSelectionQuery, StructureSelectionQueryOptions } from '../../helpers/structure-selection-query';
import { CustomStructureProperties } from '../../transforms/model';
import { StructureRepresentation3D } from '../../transforms/representation';
import { HierarchyRef, StructureComponentRef, StructureRef, StructureRepresentationRef } from './hierarchy-state';

export { StructureComponentManager };

interface StructureComponentManagerState {
    options: StructureComponentManager.Options
}

class StructureComponentManager extends PluginComponent<StructureComponentManagerState> {
    readonly events = {
        optionsUpdated: this.ev<undefined>()
    }

    get currentStructures() {
        return this.plugin.managers.structure.hierarchy.state.current.structures;
    }

    get pivotStructure(): StructureRef | undefined {
        return this.currentStructures[0];
    }

    async setOptions(options: StructureComponentManager.Options) {
        const interactionChanged = options.interactions !== this.state.options.interactions;
        this.updateState({ options });
        this.events.optionsUpdated.next();

        const update = this.dataState.build();

        for (const s of this.currentStructures) {
            for (const c of s.components) {
                this.updateReprParams(update, c);
            }
            if (s.currentFocus?.focus) this.updateReprParams(update, s.currentFocus.focus);
            if (s.currentFocus?.surroundings) this.updateReprParams(update, s.currentFocus.surroundings);
        }

        return this.plugin.dataTransaction(async () => {
            await this.plugin.runTask(this.dataState.updateTree(update));
            if (interactionChanged) await this.updateInterationProps();
        });
    }

    private updateReprParams(update: StateBuilder.Root, component: StructureComponentRef) {
        const { showHydrogens, visualQuality: quality } = this.state.options;
        const ignoreHydrogens = !showHydrogens;
        for (const r of component.representations) {
            if (r.cell.transform.transformer !== StructureRepresentation3D) continue;

            const params = r.cell.transform.params as StateTransformer.Params<StructureRepresentation3D>;
            if (!!params.type.params.ignoreHydrogens !== ignoreHydrogens || params.type.params.quality !== quality) {
                update.to(r.cell).update(old => {
                    old.type.params.ignoreHydrogens = ignoreHydrogens;
                    old.type.params.quality = quality;
                });
            }
        }
    }

    private async updateInterationProps() {
        for (const s of this.currentStructures) {
            const interactionParams = InteractionsProvider.getParams(s.cell.obj?.data!);
            
            if (s.properties) {
                const params = s.properties.cell.transform.params;
                if (PD.areEqual(interactionParams, params, this.state.options.interactions)) continue;

                const b = this.dataState.build();
                b.to(s.properties.cell).update((old: StateTransformer.Params<CustomStructureProperties>) => {
                    arraySetAdd(old.autoAttach, InteractionsProvider.descriptor.name);
                    old.properties[InteractionsProvider.descriptor.name] = this.state.options.interactions;
                });
                await this.plugin.runTask(this.dataState.updateTree(b));
            } else {
                const pd = this.plugin.customStructureProperties.getParams(s.cell.obj?.data);
                const params = PD.getDefaultValues(pd);
                if (PD.areEqual(interactionParams, params.properties[InteractionsProvider.descriptor.name], this.state.options.interactions)) continue;

                arraySetAdd(params.autoAttach, InteractionsProvider.descriptor.name);
                params.properties[InteractionsProvider.descriptor.name] = this.state.options.interactions;
                await this.plugin.builders.structure.insertStructureProperties(s.cell, params);
            }
        }
    }

    applyPreset<P = any, S = {}>(structures: ReadonlyArray<StructureRef>, provider: StructureRepresentationProvider<P, S>, params?: P): Promise<any>  {
        return this.plugin.dataTransaction(async () => {
            await this.clearComponents(structures);
            for (const s of structures) {
                await this.plugin.builders.structure.representation.applyPreset(s.cell, provider, params);
            }
        }, { canUndo: true });
    }

    clear(structures: ReadonlyArray<StructureRef>) {
        return this.clearComponents(structures);
    }

    selectThis(components: ReadonlyArray<StructureComponentRef>) {
        const mng = this.plugin.managers.structure.selection;
        mng.clear();
        for (const c of components) {
            const loci =  Structure.toSubStructureElementLoci(c.structure.cell.obj!.data, c.cell.obj?.data!)
            mng.fromLoci('set', loci);
        }
    }

    canBeModified(ref: HierarchyRef) {
        return this.plugin.builders.structure.isComponentTransform(ref.cell);
    }

    modifyByCurrentSelection(components: ReadonlyArray<StructureComponentRef>, action: StructureComponentManager.ModifyAction) {
        return this.plugin.runTask(Task.create('Modify Component', async taskCtx => {
            const b = this.dataState.build();        
            for (const c of components) {
                if (!this.canBeModified(c)) continue;

                const selection = this.plugin.managers.structure.selection.getStructure(c.structure.cell.obj!.data);
                if (!selection || selection.elementCount === 0) continue;                
                this.modifyComponent(b, c, selection, action);
            }
            await this.dataState.updateTree(b, { canUndo: true }).runInContext(taskCtx);
        }));
    }

    toggleVisibility(components: ReadonlyArray<StructureComponentRef>) {
        if (components.length === 0) return;
        const isHidden = !components[0].cell.state.isHidden;
        for (const c of components) {
            setSubtreeVisibility(this.dataState, c.cell.transform.ref, isHidden);
        }
    }

    removeRepresentations(components: ReadonlyArray<StructureComponentRef>, pivot?: StructureRepresentationRef) {
        if (components.length === 0) return;

        const toRemove: HierarchyRef[] = [];
        if (pivot) {
            const index = components[0].representations.indexOf(pivot);
            if (index < 0) return;

            for (const c of components) {
                if (c.representations[index]) toRemove.push(c.representations[index]);
            }
        } else {
            for (const c of components) {
                for (const r of c.representations) {
                    toRemove.push(r);
                }
            }
        }

        return this.plugin.managers.structure.hierarchy.remove(toRemove, true);
    }

    updateRepresentations(components: ReadonlyArray<StructureComponentRef>, pivot: StructureRepresentationRef, params: StateTransformer.Params<StructureRepresentation3D>) {        
        if (components.length === 0) return Promise.resolve();

        const index = components[0].representations.indexOf(pivot);
        if (index < 0) return Promise.resolve();

        const update = this.dataState.build();

        for (const c of components) {
            const repr = c.representations[index];
            if (!repr) continue;
            update.to(repr.cell).update(params);
        }

        return this.plugin.runTask(this.dataState.updateTree(update, { canUndo: true }));
    }

    addRepresentation(components: ReadonlyArray<StructureComponentRef>, type: string) {
        if (components.length === 0) return;

        const { showHydrogens, visualQuality: quality } = this.state.options;
        const ignoreHydrogens = !showHydrogens;
        const typeParams = { ignoreHydrogens, quality };

        return this.plugin.dataTransaction(async () => {
            for (const component of components) {
                await this.plugin.builders.structure.representation.addRepresentation(component.cell, {
                    type: this.plugin.structureRepresentation.registry.get(type),
                    typeParams
                });
            }
        }, { canUndo: true });
    }

    async add(params: StructureComponentManager.AddParams, structures?: ReadonlyArray<StructureRef>) {
        return this.plugin.dataTransaction(async () => {
            const xs = structures || this.currentStructures;
            if (xs.length === 0) return;

            const componentKey = UUID.create22();
            for (const s of xs) {
                const component = await this.plugin.builders.structure.tryCreateQueryComponent({ 
                    structure: s.cell,
                    query: params.selection,
                    key: componentKey,
                    label: params.label || (params.selection === StructureSelectionQueries.current ? 'Custom Selection' : ''),
                });
                if (params.representation === 'none' || !component) continue;
                await this.plugin.builders.structure.representation.addRepresentation(component, {
                    type: this.plugin.structureRepresentation.registry.get(params.representation)
                });
            }
        }, { canUndo: true });
    }

    async applyColor(params: StructureComponentManager.ColorParams, structures?: ReadonlyArray<StructureRef>) {
        return this.plugin.dataTransaction(async () => {
            const xs = structures || this.currentStructures;
            if (xs.length === 0) return;
            const getLoci = (s: Structure) => this.plugin.managers.structure.selection.getLoci(s);

            for (const s of xs) {
                if (params.action.name === 'reset') {
                    await clearStructureOverpaint(this.plugin, s.components, params.representations);
                } else {
                    const p = params.action.params;
                    await setStructureOverpaint(this.plugin, s.components, p.color, getLoci, params.representations, p.opacity);
                }
            }
        }, { canUndo: true });
    }

    private modifyComponent(builder: StateBuilder.Root, component: StructureComponentRef, by: Structure, action: StructureComponentManager.ModifyAction) {
        const structure = component.cell.obj?.data;
        if (!structure) return;
        if ((action === 'subtract' || action === 'intersect') && !structureAreIntersecting(structure, by)) return;

        const parent = component.structure.cell.obj?.data!;
        const modified = action === 'union' 
            ? structureUnion(parent, [structure, by]) 
            : action === 'intersect'
            ? structureIntersect(structure, by)
            : structureSubtract(structure, by);

        if (modified.elementCount === 0) {
            builder.delete(component.cell.transform.ref);
        } else {
            const bundle = StructureElement.Bundle.fromSubStructure(parent, modified);
            const params: StructureComponentParams = {
                type: { name: 'bundle', params: bundle },
                nullIfEmpty: true,
                label: component.cell.obj?.label!
            };
            builder.to(component.cell).update(params)
        }
    }

    private get dataState() {
        return this.plugin.state.dataState;
    }

    private clearComponents(structures: ReadonlyArray<StructureRef>) {
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
        super({ options: PD.getDefaultValues(StructureComponentManager.OptionsParams) })
    }
}

namespace StructureComponentManager {
    export const OptionsParams = {
        showHydrogens: PD.Boolean(true, { description: 'Toggle display of hydrogen atoms in representations' }),
        visualQuality: PD.Select('auto', VisualQualityOptions, { description: 'Control the visual/rendering quality of representations' }),
        interactions: PD.Group(InteractionsProvider.defaultParams, { label: 'Non-covalent Interactions' }),
    }
    export type Options = PD.Values<typeof OptionsParams>
    
    const SelectionParam = PD.Select(StructureSelectionQueryOptions[1][0], StructureSelectionQueryOptions)
    
    export function getAddParams(plugin: PluginContext) {
        return {
            selection: SelectionParam,
            representation: getRepresentationTypesSelect(plugin, plugin.managers.structure.component.pivotStructure, [['none', '< None >']]),
            label: PD.Text('')
        };
    }
    export type AddParams = { selection: StructureSelectionQuery, label: string, representation: string }

    export function getColorParams(plugin: PluginContext, pivot: StructureRef | StructureComponentRef | undefined) {
        return {
            action: PD.MappedStatic('color', {
                color: PD.Group({
                    color: PD.Color(ColorNames.blue, { isExpanded: true }),
                    opacity: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
                }, { isFlat: true }),
                reset: PD.EmptyGroup()
            }),
            representations: PD.MultiSelect([], getRepresentationTypes(plugin, pivot), { emptyValue: 'All' })
        };
    }
    export type ColorParams = PD.Values<ReturnType<typeof getColorParams>>

    export function getRepresentationTypes(plugin: PluginContext, pivot: StructureRef | StructureComponentRef | undefined) {
        return pivot?.cell.obj?.data
            ? plugin.structureRepresentation.registry.getApplicableTypes(pivot.cell.obj?.data!)
            : plugin.structureRepresentation.registry.types;
    }

    function getRepresentationTypesSelect(plugin: PluginContext, pivot: StructureRef | undefined, custom: [string, string][], label?: string) {
        const types = [
            ...custom,
            ...getRepresentationTypes(plugin, pivot)
        ] as [string, string][];
        return PD.Select(types[0][0], types, { label });
    }

    export type ModifyAction = 'union' | 'subtract' | 'intersect'
}
