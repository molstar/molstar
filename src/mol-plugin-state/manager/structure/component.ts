/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { VisualQualityOptions } from '../../../mol-geo/geometry/base';
import { InteractionsProvider } from '../../../mol-model-props/computed/interactions';
import { Structure, StructureElement, StructureSelection } from '../../../mol-model/structure';
import { structureAreEqual, structureAreIntersecting, structureIntersect, structureSubtract, structureUnion } from '../../../mol-model/structure/query/utils/structure-set';
import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginContext } from '../../../mol-plugin/context';
import { StateBuilder, StateObjectRef, StateTransformer } from '../../../mol-state';
import { Task } from '../../../mol-task';
import { ColorTheme } from '../../../mol-theme/color';
import { SizeTheme } from '../../../mol-theme/size';
import { UUID } from '../../../mol-util';
import { ColorNames } from '../../../mol-util/color/names';
import { objectForEach } from '../../../mol-util/object';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentationPresetProvider } from '../../builder/structure/representation-preset';
import { StatefulPluginComponent } from '../../component';
import { StructureComponentParams } from '../../helpers/structure-component';
import { setStructureOverpaint } from '../../helpers/structure-overpaint';
import { createStructureColorThemeParams, createStructureSizeThemeParams } from '../../helpers/structure-representation-params';
import { StructureSelectionQueries, StructureSelectionQuery } from '../../helpers/structure-selection-query';
import { StructureRepresentation3D } from '../../transforms/representation';
import { StructureHierarchyRef, StructureComponentRef, StructureRef, StructureRepresentationRef } from './hierarchy-state';
import { Clipping } from '../../../mol-theme/clipping';
import { setStructureClipping } from '../../helpers/structure-clipping';
import { setStructureTransparency } from '../../helpers/structure-transparency';

export { StructureComponentManager };

interface StructureComponentManagerState {
    options: StructureComponentManager.Options
}

class StructureComponentManager extends StatefulPluginComponent<StructureComponentManagerState> {
    readonly events = {
        optionsUpdated: this.ev<undefined>()
    }

    get currentStructures() {
        return this.plugin.managers.structure.hierarchy.selection.structures;
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
        }

        return this.plugin.dataTransaction(async () => {
            await update.commit();
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
                const oldParams = s.properties.cell.transform.params?.properties[InteractionsProvider.descriptor.name];
                if (PD.areEqual(interactionParams, oldParams, this.state.options.interactions)) continue;

                await this.dataState.build().to(s.properties.cell)
                    .update(old => {
                        old.properties[InteractionsProvider.descriptor.name] = this.state.options.interactions;
                    })
                    .commit();
            } else {
                const pd = this.plugin.customStructureProperties.getParams(s.cell.obj?.data);
                const params = PD.getDefaultValues(pd);
                if (PD.areEqual(interactionParams, params.properties[InteractionsProvider.descriptor.name], this.state.options.interactions)) continue;
                params.properties[InteractionsProvider.descriptor.name] = this.state.options.interactions;
                await this.plugin.builders.structure.insertStructureProperties(s.cell, params);
            }
        }
    }

    applyPreset<P extends StructureRepresentationPresetProvider>(structures: ReadonlyArray<StructureRef>, provider: P, params?: StructureRepresentationPresetProvider.Params<P>): Promise<any>  {
        return this.plugin.dataTransaction(async () => {
            for (const s of structures) {
                const preset = await this.plugin.builders.structure.representation.applyPreset(s.cell, provider, params);
                await this.syncPreset(s, preset);
            }
        }, { canUndo: 'Preset' });
    }

    private syncPreset(root: StructureRef, preset?: StructureRepresentationPresetProvider.Result) {
        if (!preset || !preset.components) return this.clearComponents([root]);

        const keptRefs = new Set<string>();
        objectForEach(preset.components, c => {
            if (c) keptRefs.add(c.ref);
        });

        if (preset.representations) {
            objectForEach(preset.representations, r => {
                if (r) keptRefs.add(r.ref);
            });
        }

        if (keptRefs.size === 0) return this.clearComponents([root]);

        let changed = false;
        const update = this.dataState.build();

        const sync = (r: StructureHierarchyRef) => {
            if (!keptRefs.has(r.cell.transform.ref)) {
                changed = true;
                update.delete(r.cell);
            }
        };

        for (const c of root.components) {
            sync(c);
            for (const r of c.representations) sync(r);
            if (c.genericRepresentations) {
                for (const r of c.genericRepresentations) sync(r);
            }
        }

        if (root.genericRepresentations) {
            for (const r of root.genericRepresentations) {
                sync(r);
            }
        }

        if (changed) return update.commit();
    }

    clear(structures: ReadonlyArray<StructureRef>) {
        return this.clearComponents(structures);
    }

    selectThis(components: ReadonlyArray<StructureComponentRef>) {
        const mng = this.plugin.managers.structure.selection;
        mng.clear();
        for (const c of components) {
            const loci =  Structure.toSubStructureElementLoci(c.structure.cell.obj!.data, c.cell.obj?.data!);
            mng.fromLoci('set', loci);
        }
    }

    canBeModified(ref: StructureHierarchyRef) {
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
            await this.dataState.updateTree(b, { canUndo: 'Modify Selection' }).runInContext(taskCtx);
        }));
    }

    toggleVisibility(components: ReadonlyArray<StructureComponentRef>, reprPivot?: StructureRepresentationRef) {
        if (components.length === 0) return;

        if (!reprPivot) {
            const isHidden = !components[0].cell.state.isHidden;
            for (const c of components) {
                setSubtreeVisibility(this.dataState, c.cell.transform.ref, isHidden);
            }
        } else {
            const index = components[0].representations.indexOf(reprPivot);
            const isHidden = !reprPivot.cell.state.isHidden;

            for (const c of components) {
                // TODO: is it ok to use just the index here? Could possible lead to ugly edge cases, but perhaps not worth the trouble to "fix".
                const repr = c.representations[index];
                if (!repr) continue;
                setSubtreeVisibility(this.dataState, repr.cell.transform.ref, isHidden);
            }
        }
    }

    removeRepresentations(components: ReadonlyArray<StructureComponentRef>, pivot?: StructureRepresentationRef) {
        if (components.length === 0) return;

        const toRemove: StructureHierarchyRef[] = [];
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
            // TODO: is it ok to use just the index here? Could possible lead to ugly edge cases, but perhaps not worth the trouble to "fix".
            const repr = c.representations[index];
            if (!repr) continue;
            if (repr.cell.transform.transformer !== pivot.cell.transform.transformer) continue;

            update.to(repr.cell).update(params);
        }

        return update.commit({ canUndo: 'Update Representation' });
    }

    /**
     * To update theme for all selected structures, use
     *   plugin.dataTransaction(async () => {
     *      for (const s of structure.hierarchy.selection.structures) await updateRepresentationsTheme(s.componets, ...);
     *   }, { canUndo: 'Update Theme' });
     */
    updateRepresentationsTheme<C extends ColorTheme.BuiltIn, S extends SizeTheme.BuiltIn>(components: ReadonlyArray<StructureComponentRef>, params: StructureComponentManager.UpdateThemeParams<C, S>): Promise<any> | undefined
    updateRepresentationsTheme<C extends ColorTheme.BuiltIn, S extends SizeTheme.BuiltIn>(components: ReadonlyArray<StructureComponentRef>, params: (c: StructureComponentRef, r: StructureRepresentationRef) => StructureComponentManager.UpdateThemeParams<C, S>): Promise<any> | undefined
    updateRepresentationsTheme(components: ReadonlyArray<StructureComponentRef>, paramsOrProvider: StructureComponentManager.UpdateThemeParams<any, any> | ((c: StructureComponentRef, r: StructureRepresentationRef) => StructureComponentManager.UpdateThemeParams<any, any>)) {
        if (components.length === 0) return;

        const update = this.dataState.build();

        for (const c of components) {
            for (const repr of c.representations) {
                const old = repr.cell.transform.params;
                const params: StructureComponentManager.UpdateThemeParams<any, any> = typeof paramsOrProvider === 'function' ? paramsOrProvider(c, repr) : paramsOrProvider;

                const colorTheme = params.color === 'default'
                    ? createStructureColorThemeParams(this.plugin, c.structure.cell.obj?.data, old?.type.name)
                    : params.color
                        ? createStructureColorThemeParams(this.plugin, c.structure.cell.obj?.data, old?.type.name, params.color, params.colorParams)
                        : void 0;
                const sizeTheme = params.size === 'default'
                    ? createStructureSizeThemeParams(this.plugin, c.structure.cell.obj?.data, old?.type.name)
                    : params.color
                        ? createStructureSizeThemeParams(this.plugin, c.structure.cell.obj?.data, old?.type.name, params.size, params.sizeParams)
                        : void 0;

                if (colorTheme || sizeTheme) {
                    update.to(repr.cell).update(prev => {
                        if (colorTheme) prev.colorTheme = colorTheme;
                        if (sizeTheme) prev.sizeTheme = sizeTheme;
                    });
                }
            }
        }

        return update.commit({ canUndo: 'Update Theme' });
    }

    addRepresentation(components: ReadonlyArray<StructureComponentRef>, type: string) {
        if (components.length === 0) return;

        const { showHydrogens, visualQuality: quality } = this.state.options;
        const ignoreHydrogens = !showHydrogens;
        const typeParams = { ignoreHydrogens, quality };

        return this.plugin.dataTransaction(async () => {
            for (const component of components) {
                await this.plugin.builders.structure.representation.addRepresentation(component.cell, {
                    type: this.plugin.representation.structure.registry.get(type),
                    typeParams
                });
            }
        }, { canUndo: 'Add Representation' });
    }

    private tryFindComponent(structure: StructureRef, selection: StructureSelectionQuery) {
        if (structure.components.length === 0) return;

        return this.plugin.runTask(Task.create('Find Component', async taskCtx => {

            const data = structure.cell.obj?.data;
            if (!data) return;
            const sel = StructureSelection.unionStructure(await selection.getSelection(this.plugin, taskCtx, data));

            for (const c of structure.components) {
                const comp = c.cell.obj?.data;
                if (!comp || !c.cell.parent) continue;

                if (structureAreEqual(sel, comp)) return c.cell;
            }
        }));
    }

    async add(params: StructureComponentManager.AddParams, structures?: ReadonlyArray<StructureRef>) {
        return this.plugin.dataTransaction(async () => {
            const xs = structures || this.currentStructures;
            if (xs.length === 0) return;

            const { showHydrogens, visualQuality: quality } = this.state.options;
            const ignoreHydrogens = !showHydrogens;
            const typeParams = { ignoreHydrogens, quality };

            const componentKey = UUID.create22();
            for (const s of xs) {
                let component: StateObjectRef | undefined = void 0;

                if (params.options.checkExisting) {
                    component = await this.tryFindComponent(s, params.selection);
                }

                if (!component) {
                    component = await this.plugin.builders.structure.tryCreateComponentFromSelection(s.cell, params.selection, componentKey, {
                        label: params.options.label || (params.selection === StructureSelectionQueries.current ? 'Custom Selection' : ''),
                    });
                }

                if (params.representation === 'none' || !component) continue;
                await this.plugin.builders.structure.representation.addRepresentation(component, {
                    type: this.plugin.representation.structure.registry.get(params.representation),
                    typeParams
                });
            }
        }, { canUndo: 'Add Selection' });
    }

    async applyTheme(params: StructureComponentManager.ThemeParams, structures?: ReadonlyArray<StructureRef>) {
        return this.plugin.dataTransaction(async ctx => {
            const xs = structures || this.currentStructures;
            if (xs.length === 0) return;

            const getLoci = async (s: Structure) => StructureSelection.toLociWithSourceUnits(await params.selection.getSelection(this.plugin, ctx, s));
            for (const s of xs) {
                if (params.action.name === 'reset') {
                    await setStructureOverpaint(this.plugin, s.components, -1, getLoci, params.representations);
                } else if (params.action.name === 'color') {
                    const p = params.action.params;
                    await setStructureOverpaint(this.plugin, s.components, p.color, getLoci, params.representations);
                } else if (params.action.name === 'transparency') {
                    const p = params.action.params;
                    await setStructureTransparency(this.plugin, s.components, p.value, getLoci, params.representations);
                } else if (params.action.name === 'clipping') {
                    const p = params.action.params;
                    await setStructureClipping(this.plugin, s.components, Clipping.Groups.fromNames(p.excludeGroups), getLoci, params.representations);
                }
            }
        }, { canUndo: 'Apply Theme' });
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
            builder.to(component.cell).update(params);
        }
    }

    updateLabel(component: StructureComponentRef, label: string) {
        const params: StructureComponentParams = {
            type: component.cell.params?.values.type,
            nullIfEmpty: component.cell.params?.values.nullIfEmpty,
            label
        };
        this.dataState.build().to(component.cell).update(params).commit();
    }

    private get dataState() {
        return this.plugin.state.data;
    }

    private clearComponents(structures: ReadonlyArray<StructureRef>) {
        const deletes = this.dataState.build();
        for (const s of structures) {
            for (const c of s.components) {
                deletes.delete(c.cell.transform.ref);
            }
        }
        return deletes.commit({ canUndo: 'Clear Selections' });
    }

    constructor(public plugin: PluginContext) {
        super({ options: PD.getDefaultValues(StructureComponentManager.OptionsParams) });
    }
}

namespace StructureComponentManager {
    export const OptionsParams = {
        showHydrogens: PD.Boolean(true, { description: 'Toggle display of hydrogen atoms in representations' }),
        visualQuality: PD.Select('auto', VisualQualityOptions, { description: 'Control the visual/rendering quality of representations' }),
        interactions: PD.Group(InteractionsProvider.defaultParams, { label: 'Non-covalent Interactions' }),
    };
    export type Options = PD.Values<typeof OptionsParams>

    export function getAddParams(plugin: PluginContext, params?: { pivot?: StructureRef, allowNone: boolean, hideSelection?: boolean, checkExisting?: boolean }) {
        const { options } = plugin.query.structure.registry;
        params = {
            pivot: plugin.managers.structure.component.pivotStructure,
            allowNone: true,
            hideSelection: false,
            checkExisting: false,
            ...params
        };
        return {
            selection: PD.Select(options[1][0], options, { isHidden: params?.hideSelection }),
            representation: getRepresentationTypesSelect(plugin, params?.pivot, params?.allowNone ? [['none', '< Create Later >']] : []),
            options: PD.Group({
                label: PD.Text(''),
                checkExisting: PD.Boolean(!!params?.checkExisting, { help: () => ({ description: 'Checks if a selection with the specifield elements already exists to avoid creating duplicate components.' }) }),
            })
        };
    }
    export type AddParams = { selection: StructureSelectionQuery, options: { checkExisting: boolean, label: string }, representation: string }

    export function getThemeParams(plugin: PluginContext, pivot: StructureRef | StructureComponentRef | undefined) {
        const { options } = plugin.query.structure.registry;
        return {
            selection: PD.Select(options[1][0], options, { isHidden: false }),
            action: PD.MappedStatic('color', {
                color: PD.Group({
                    color: PD.Color(ColorNames.blue, { isExpanded: true }),
                }, { isFlat: true }),
                reset: PD.EmptyGroup({ label: 'Reset Color' }),
                transparency: PD.Group({
                    value: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
                }, { isFlat: true }),
                clipping: PD.Group({
                    excludeGroups: PD.MultiSelect([] as Clipping.Groups.Names[], PD.objectToOptions(Clipping.Groups.Names)),
                }, { isFlat: true }),
            }),
            representations: PD.MultiSelect([], getRepresentationTypes(plugin, pivot), { emptyValue: 'All' })
        };
    }
    export type ThemeParams = PD.Values<ReturnType<typeof getThemeParams>>

    export function getRepresentationTypes(plugin: PluginContext, pivot: StructureRef | StructureComponentRef | undefined) {
        return pivot?.cell.obj?.data
            ? plugin.representation.structure.registry.getApplicableTypes(pivot.cell.obj?.data!)
            : plugin.representation.structure.registry.types;
    }

    function getRepresentationTypesSelect(plugin: PluginContext, pivot: StructureRef | undefined, custom: [string, string][], label?: string) {
        const types = [
            ...custom,
            ...getRepresentationTypes(plugin, pivot)
        ] as [string, string][];
        return PD.Select(types[0][0], types, { label });
    }

    export type ModifyAction = 'union' | 'subtract' | 'intersect'

    export interface UpdateThemeParams<C extends ColorTheme.BuiltIn, S extends SizeTheme.BuiltIn> {
        /**
         * this works for any theme name (use 'name as any'), but code completion will break
         */
        color?: C | 'default',
        colorParams?: ColorTheme.BuiltInParams<C>,
        size?: S | 'default',
        sizeParams?: SizeTheme.BuiltInParams<S>
    }
}
