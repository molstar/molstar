/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior';
import { ParamDefinition as PD } from '../../../../../mol-util/param-definition';
import { ValenceModelProvider } from '../../../../../mol-model-props/computed/valence-model';
import { Loci } from '../../../../../mol-model/loci';
import { PluginStateObject } from '../../../../../mol-plugin-state/objects';
import { StateSelection } from '../../../../../mol-state';
import { Structure, StructureElement } from '../../../../../mol-model/structure';
import { OrderedSet } from '../../../../../mol-data/int';
import { geometryLabel } from '../../../../../mol-model-props/computed/chemistry/geometry';
import { arraySetAdd } from '../../../../../mol-util/array';

export const ValenceModel = PluginBehavior.create<{ autoAttach: boolean, showTooltip: boolean }>({
    name: 'computed-valence-model-prop',
    category: 'custom-props',
    display: { name: 'Valence Model' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean, showTooltip: boolean }> {
        private provider = ValenceModelProvider

        private getStructures(structure: Structure) {
            const structures: Structure[] = [];
            const root = this.ctx.helpers.substructureParent.get(structure);
            if (root) {
                const state = this.ctx.state.data;
                const selections = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure, root.transform.ref));
                for (const s of selections) {
                    if (s.obj) arraySetAdd(structures, s.obj.data);
                }
            }
            return structures;
        }

        private labelProvider = {
            label: (loci: Loci): string | undefined => {
                if (!this.params.showTooltip) return void 0;

                switch (loci.kind) {
                    case 'element-loci':
                        if (loci.elements.length === 0) return void 0;

                        const labels: string[] = [];
                        const structures = this.getStructures(loci.structure);

                        for (const s of structures) {
                            const valenceModel = this.provider.get(s).value;
                            if (!valenceModel) continue;

                            const l = StructureElement.Loci.remap(loci, s);
                            if (l.elements.length !== 1) continue;

                            const e = l.elements[0];
                            if (OrderedSet.size(e.indices) !== 1) continue;

                            const vm = valenceModel.get(e.unit.id);
                            if (!vm) continue;

                            const idx = OrderedSet.start(e.indices);
                            const charge = vm.charge[idx];
                            const idealGeometry = vm.idealGeometry[idx];
                            const implicitH = vm.implicitH[idx];
                            const totalH = vm.totalH[idx];

                            labels.push(`Valence Model: <small>Charge</small> ${charge} | <small>Ideal Geometry</small> ${geometryLabel(idealGeometry)} | <small>Implicit H</small> ${implicitH} | <small>Total H</small> ${totalH}`);
                        }

                        return labels.length ? labels.join('<br/>') : undefined;

                    default: return void 0;
                }
            }
        }

        update(p: { autoAttach: boolean, showTooltip: boolean }) {
            let updated = (
                this.params.autoAttach !== p.autoAttach ||
                this.params.showTooltip !== p.showTooltip
            );
            this.params.autoAttach = p.autoAttach;
            this.params.showTooltip = p.showTooltip;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        register(): void {
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);
            this.ctx.managers.lociLabels.addProvider(this.labelProvider);
        }

        unregister() {
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
            this.ctx.managers.lociLabels.removeProvider(this.labelProvider);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false),
        showTooltip: PD.Boolean(true)
    })
});