/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { EmptyLoci } from 'mol-model/loci';
import { StructureUnitTransforms } from 'mol-model/structure/structure/util/unit-transforms';
import { PluginContext } from 'mol-plugin/context';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { StateObjectTracker, StateSelection } from 'mol-state';
import { labelFirst } from 'mol-theme/label';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { PluginBehavior } from '../behavior';
import { Representation } from 'mol-repr/representation';

export const HighlightLoci = PluginBehavior.create({
    name: 'representation-highlight-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            let prev: Representation.Loci = { loci: EmptyLoci, repr: void 0 };

            this.subscribeObservable(this.ctx.events.canvas3d.highlight, ({ loci }) => {
                if (!this.ctx.canvas3d) return;

                if (!Representation.Loci.areEqual(prev, loci)) {
                    this.ctx.canvas3d.mark(prev, MarkerAction.RemoveHighlight);
                    this.ctx.canvas3d.mark(loci, MarkerAction.Highlight);
                    prev = loci;
                }
            });
        }
    },
    display: { name: 'Highlight Loci on Canvas', group: 'Representation' }
});

export const SelectLoci = PluginBehavior.create({
    name: 'representation-select-loci',
    category: 'interaction',
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            let prev = Representation.Loci.Empty;
            this.subscribeObservable(this.ctx.events.canvas3d.click, ({ loci: current }) => {
                if (!this.ctx.canvas3d) return;
                if (!Representation.Loci.areEqual(prev, current)) {
                    this.ctx.canvas3d.mark(prev, MarkerAction.Deselect);
                    this.ctx.canvas3d.mark(current, MarkerAction.Select);
                    prev = current;
                } else {
                    this.ctx.canvas3d.mark(current, MarkerAction.Toggle);
                }
                // this.ctx.canvas3d.mark(loci, MarkerAction.Toggle);
            });
        }
    },
    display: { name: 'Select Loci on Canvas', group: 'Representation' }
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
    display: { name: 'Provide Default Loci Label', group: 'Representation' }
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