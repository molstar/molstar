/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { EmptyLoci } from 'mol-model/loci';
import { QueryContext, StructureElement, StructureSelection } from 'mol-model/structure';
import { PluginContext } from 'mol-plugin/context';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { Representation } from 'mol-repr/representation';
import Expression from 'mol-script/language/expression';
import { parseMolScript } from 'mol-script/language/parser';
import { compile } from 'mol-script/runtime/query/compiler';
import { transpileMolScript } from 'mol-script/script/mol-script/symbols';
import { StateObjectTracker, StateSelection } from 'mol-state';
import { labelFirst } from 'mol-theme/label';
import { Overpaint } from 'mol-theme/overpaint';
import { Color } from 'mol-util/color';
import { ColorNames } from 'mol-util/color/tables';
import { ButtonsType } from 'mol-util/input/input-observer';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginBehavior } from '../behavior';

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            let prev: Representation.Loci = { loci: EmptyLoci, repr: void 0 };
            const sel = this.ctx.helpers.structureSelection;

            this.subscribeObservable(this.ctx.behaviors.canvas3d.highlight, ({ current, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                if (StructureElement.isLoci(current.loci)) {
                    let loci: StructureElement.Loci = current.loci;
                    if (modifiers && modifiers.shift) {
                        loci = sel.tryGetRange(loci) || loci;
                    }

                    this.ctx.canvas3d.mark(prev, MarkerAction.RemoveHighlight);
                    const toHighlight = { loci, repr: current.repr };
                    this.ctx.canvas3d.mark(toHighlight, MarkerAction.Highlight);
                    prev = toHighlight;
                } else {
                    if (!Representation.Loci.areEqual(prev, current)) {
                        this.ctx.canvas3d.mark(prev, MarkerAction.RemoveHighlight);
                        this.ctx.canvas3d.mark(current, MarkerAction.Highlight);
                        prev = current;
                    }
                }

            });
        }
    },
    display: { name: 'Highlight Loci on Canvas' }
});

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            const sel = this.ctx.helpers.structureSelection;

            const toggleSel = (current: Representation.Loci<StructureElement.Loci>) => {
                if (sel.has(current.loci)) {
                    sel.remove(current.loci);
                    this.ctx.canvas3d.mark(current, MarkerAction.Deselect);
                } else {
                    sel.add(current.loci);
                    this.ctx.canvas3d.mark(current, MarkerAction.Select);
                }
            }

            this.subscribeObservable(this.ctx.behaviors.canvas3d.click, ({ current, buttons, modifiers }) => {
                if (!this.ctx.canvas3d) return;

                if (current.loci.kind === 'empty-loci') {
                    if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                        // clear the selection on Ctrl + Right-Click on empty
                        const sels = sel.clear();
                        for (const s of sels) this.ctx.canvas3d.mark({ loci: s }, MarkerAction.Deselect);
                    }
                } else if (StructureElement.isLoci(current.loci)) {
                    if (modifiers.control && buttons === ButtonsType.Flag.Secondary) {
                        // select only the current element on Ctrl + Right-Click
                        const old = sel.get(current.loci.structure);
                        this.ctx.canvas3d.mark({ loci: old }, MarkerAction.Deselect);
                        sel.set(current.loci);
                        this.ctx.canvas3d.mark(current, MarkerAction.Select);
                    } else if (modifiers.control && buttons === ButtonsType.Flag.Primary) {
                        // toggle current element on Ctrl + Left-Click
                        toggleSel(current as Representation.Loci<StructureElement.Loci>);
                    } else if (modifiers.shift && buttons === ButtonsType.Flag.Primary) {
                        // try to extend sequence on Shift + Left-Click
                        let loci: StructureElement.Loci = current.loci;
                        if (modifiers && modifiers.shift) {
                            loci = sel.tryGetRange(loci) || loci;
                        }
                        toggleSel({ loci, repr: current.repr });
                    }
                } else {
                    if (!ButtonsType.has(buttons, ButtonsType.Flag.Secondary)) return;
                    this.ctx.canvas3d.mark(current, MarkerAction.Toggle);
                }
            });
        }
    },
    display: { name: 'Select Loci on Canvas' }
});

export const DefaultLociLabelProvider = PluginBehavior.create({
    name: 'default-loci-label-provider',
    category: 'interaction',
    ctor: class implements PluginBehavior<undefined> {
        private f = labelFirst;
        register(): void { this.ctx.lociLabels.addProvider(this.f); }
        unregister() { this.ctx.lociLabels.removeProvider(this.f); }
        constructor(protected ctx: PluginContext) { }
    },
    display: { name: 'Provide Default Loci Label' }
});

type ColorMappings = { query: Expression, color: Color }[]
namespace ColorMappings {
    export function areEqual(colorMappingsA: ColorMappings, colorMappingsB: ColorMappings) {
        return false
    }
}

export namespace ColorRepresentation3D {
    export const Params = {
        layers: PD.ObjectList({
            query: PD.ScriptExpression({ language: 'mol-script', expression: '(sel.atom.atom-groups :residue-test (= atom.resname LYS))' }),
            color: PD.Color(ColorNames.blueviolet)
        }, e => `${Color.toRgbString(e.color)}`, {
            defaultValue: [
                {
                    query: {
                        language: 'mol-script',
                        expression: '(sel.atom.atom-groups :residue-test (= atom.resname LYS))'
                    },
                    color: ColorNames.blueviolet
                },
                {
                    query: {
                        language: 'mol-script',
                        expression: '(sel.atom.atom-groups :residue-test (= atom.resname ALA))'
                    },
                    color: ColorNames.chartreuse
                }
            ]
        }),
        alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity' }),
    }
    export type Params = PD.Values<typeof Params>

    export class Behavior implements PluginBehavior<Params> {
        private currentColorMappings: ColorMappings = [];
        private repr: StateObjectTracker<PluginStateObject.Molecule.Structure.Representation3D>;
        private structure: StateObjectTracker<PluginStateObject.Molecule.Structure>;

        private updateData() {
            const reprUpdated = this.repr.update();
            const strucUpdated = this.structure.update();
            return reprUpdated || strucUpdated;
        }

        register(ref: string): void {
            this.repr.setQuery(StateSelection.Generators.byRef(ref).ancestorOfType([PluginStateObject.Molecule.Structure.Representation3D]));
            this.structure.setQuery(StateSelection.Generators.byRef(ref).ancestorOfType([PluginStateObject.Molecule.Structure]));
            this.update(this.params);
        }

        update(params: Params): boolean {
            const { layers, alpha } = params
            const colorMappings: ColorMappings = []
            for (let i = 0, il = params.layers.length; i < il; ++i) {
                const { query, color } = layers[i]
                const parsed = parseMolScript(query.expression);
                if (parsed.length === 0) throw new Error('No query');
                colorMappings.push({ query: transpileMolScript(parsed[0]), color })
            }
            return this.applyMappings(colorMappings, alpha)
        }

        private applyMappings(colorMappings: ColorMappings, alpha: number): boolean {
            if (!this.updateData() && ColorMappings.areEqual(colorMappings, this.currentColorMappings)) return false;
            this.currentColorMappings = colorMappings;
            if (!this.repr.data || !this.structure.data) return true;

            const layers: Overpaint.Layer[] = []
            for (let i = 0, il = this.currentColorMappings.length; i < il; ++i) {
                const { query, color } = this.currentColorMappings[i]
                const compiled = compile<StructureSelection>(query);
                const result = compiled(new QueryContext(this.structure.data));
                const loci = StructureSelection.toLoci2(result)
                layers.push({ loci, color })
            }
            return this.applyLayers({ alpha, layers })
        }

        private applyLayers(overpaint: Overpaint): boolean {
            if (!this.repr.data) return true;
            this.repr.data.repr.setState({ overpaint })
            this.ctx.canvas3d.add(this.repr.data.repr);
            this.ctx.canvas3d.requestDraw(true);
            return true;
        }

        unregister(): void {
            this.applyLayers(Overpaint.Empty) // clear
            this.repr.cell = void 0;
            this.structure.cell = void 0;
        }

        constructor(private ctx: PluginContext, private params: Params) {
            this.repr = new StateObjectTracker(ctx.state.dataState);
            this.structure = new StateObjectTracker(ctx.state.dataState);
        }
    }

    export class Obj extends PluginStateObject.CreateBehavior<Behavior>({ name: 'Color Representation3D Behavior' }) { }
}