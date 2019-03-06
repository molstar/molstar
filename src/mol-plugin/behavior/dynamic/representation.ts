/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { EmptyLoci, EveryLoci } from 'mol-model/loci';
import { StructureUnitTransforms } from 'mol-model/structure/structure/util/unit-transforms';
import { PluginContext } from 'mol-plugin/context';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { StateObjectTracker, StateSelection } from 'mol-state';
import { labelFirst } from 'mol-theme/label';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginBehavior } from '../behavior';
import { Representation } from 'mol-repr/representation';
import { ButtonsType } from 'mol-util/input/input-observer';
import { StructureElement, StructureSelection, QueryContext } from 'mol-model/structure';
import { ColorNames } from 'mol-util/color/tables';
// import { MolScriptBuilder as MS } from 'mol-script/language/builder';
import Expression from 'mol-script/language/expression';
import { Color } from 'mol-util/color';
import { compile } from 'mol-script/runtime/query/compiler';
import { Overpaint } from 'mol-theme/overpaint';
import { parseMolScript } from 'mol-script/language/parser';
import { transpileMolScript } from 'mol-script/script/mol-script/symbols';

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

export namespace ExplodeRepresentation3D {
    export const Params = {
        t: PD.Numeric(0, { min: 0, max: 1, step: 0.01 })
    }
    export type Params = PD.Values<typeof Params>

    export class Behavior implements PluginBehavior<Params> {
        private currentT = 0;
        private repr: StateObjectTracker<PluginStateObject.Molecule.Representation3D>;
        private structure: StateObjectTracker<PluginStateObject.Molecule.Structure>;
        private transforms: StructureUnitTransforms;

        private updateData() {
            const reprUpdated = this.repr.update();
            const strUpdated = this.structure.update();
            if (strUpdated && this.structure.data) {
                this.transforms = new StructureUnitTransforms(this.structure.data);
            }
            return reprUpdated || strUpdated;
        }

        register(ref: string): void {
            this.repr.setQuery(StateSelection.Generators.byRef(ref).ancestorOfType([PluginStateObject.Molecule.Representation3D]));
            this.structure.setQuery(StateSelection.Generators.byRef(ref).rootOfType([PluginStateObject.Molecule.Structure]));
            this.update(this.params);
        }

        private centerVec = Vec3.zero();
        private transVec = Vec3.zero();
        private transMat = Mat4.zero();

        update(params: Params): boolean | Promise<boolean> {
            if (!this.updateData() && params.t === this.currentT) return false;
            this.currentT = params.t;
            if (!this.structure.data || !this.repr.data) return true;

            const structure = this.structure.data;
            const boundary = structure.boundary.sphere;
            const d = boundary.radius * params.t;

            for (let i = 0, _i = structure.units.length; i < _i; i++) {
                const u = structure.units[i];

                Vec3.transformMat4(this.centerVec, u.lookup3d.boundary.sphere.center, u.conformation.operator.matrix);
                Vec3.sub(this.transVec, this.centerVec, boundary.center);
                Vec3.setMagnitude(this.transVec, this.transVec, d);
                Mat4.fromTranslation(this.transMat, this.transVec)

                this.transforms.setTransform(this.transMat, u);
            }

            // TODO: should be be "auto updated"?
            // perhaps have Representation3D.setState(state, autoSync = false)?

            // TODO: where to handle unitTransforms composition?
            // Manually or inside the representation? "inside" would better compose with future additions.
            this.repr.data.setState({ unitTransforms: this.transforms });
            this.ctx.canvas3d.add(this.repr.data);
            this.ctx.canvas3d.requestDraw(true);

            return true;
        }

        unregister(): void {
            this.update({ t: 0 })
            this.repr.cell = void 0;
            this.structure.cell = void 0;
        }

        constructor(private ctx: PluginContext, private params: Params) {
            this.repr = new StateObjectTracker(ctx.state.dataState);
            this.structure = new StateObjectTracker(ctx.state.dataState);
        }
    }

    export class Obj extends PluginStateObject.CreateBehavior<Behavior>({ name: 'Explode Representation3D Behavior' }) { }
}

type ColorMappings = { query: Expression, color: Color }[]
namespace ColorMappings {
    export function areEqual(colorMappingsA: ColorMappings, colorMappingsB: ColorMappings) {
        return false
    }
}

export namespace ColorRepresentation3D {
    export const Params = {
        query: PD.ScriptExpression({ language: 'mol-script', expression: '(sel.atom.atom-groups :residue-test (= atom.resname LYS))' }),
        color: PD.Color(ColorNames.blueviolet)
        // colorMappings: PD.Value<ColorMappings>([{ query: MS.struct.generator.atomGroups({
        //     'residue-test': MS.core.rel.eq([MS.ammp('auth_comp_id'), 'ALA'])
        // }), color: ColorNames.greenyellow }], { isHidden: true }),
    }
    export type Params = PD.Values<typeof Params>

    export class Behavior implements PluginBehavior<Params> {
        private currentColorMappings: ColorMappings = [];
        private repr: StateObjectTracker<PluginStateObject.Molecule.Representation3D>;
        private structure: StateObjectTracker<PluginStateObject.Molecule.Structure>;

        private updateData() {
            const reprUpdated = this.repr.update();
            const strucUpdated = this.structure.update();
            return reprUpdated || strucUpdated;
        }

        register(ref: string): void {
            this.repr.setQuery(StateSelection.Generators.byRef(ref).ancestorOfType([PluginStateObject.Molecule.Representation3D]));
            this.structure.setQuery(StateSelection.Generators.byRef(ref).ancestorOfType([PluginStateObject.Molecule.Structure]));
            this.update(this.params);
        }

        update(params: Params): boolean {
            const parsed = parseMolScript(params.query.expression);
            if (parsed.length === 0) throw new Error('No query');
            const query = transpileMolScript(parsed[0]);

            return this.applyMappings([{ query, color: params.color }])
        }

        private applyMappings(colorMappings: ColorMappings): boolean {
            if (!this.updateData() && ColorMappings.areEqual(colorMappings, this.currentColorMappings)) return false;
            this.currentColorMappings = colorMappings;
            if (!this.repr.data || !this.structure.data) return true;

            const layers: Overpaint.Layers = [
                // unsets the overpaint
                // TODO do smarter by looking at the current color mappings
                { loci: EveryLoci, color: ColorNames.black, alpha: 0 }
            ]
            // console.log('currentColorMappings', this.currentColorMappings)
            for (let i = 0, il = this.currentColorMappings.length; i < il; ++i) {
                const { query, color } = this.currentColorMappings[i]
                const compiled = compile<StructureSelection>(query);
                const result = compiled(new QueryContext(this.structure.data));
                const loci = StructureSelection.toLoci2(result)
                layers.push({ loci, color, alpha: 1 })
            }
            return this.applyLayers(layers)
        }

        private applyLayers(layers: Overpaint.Layers): boolean {
            if (!this.repr.data) return true;
            this.repr.data.setOverpaint(layers)
            this.ctx.canvas3d.add(this.repr.data);
            this.ctx.canvas3d.requestDraw(true);
            return true;
        }

        unregister(): void {
            this.applyLayers([{ loci: EveryLoci, color: ColorNames.black, alpha: 0 }]) // clear
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