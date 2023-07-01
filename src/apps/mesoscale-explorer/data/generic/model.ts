/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';
import { Structure } from '../../../../mol-model/structure';
import { PluginStateObject as SO, PluginStateTransform } from '../../../../mol-plugin-state/objects';
import { Task } from '../../../../mol-task';
import { StateObject } from '../../../../mol-state';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { SymmetryOperator } from '../../../../mol-math/geometry';
import { mergeUnits, partitionUnits } from '../util';

export { StructureFromGeneric };
type StructureFromGeneric = typeof StructureFromGeneric
const StructureFromGeneric = PluginStateTransform.BuiltIn({
    name: 'structure-from-generic',
    display: { name: 'Structure from Generic', description: 'Create a molecular structure from Generic models.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
    params: {
        transforms: PD.Value<Mat4[]>([]),
        label: PD.Optional(PD.Text('')),
        cellSize: PD.Numeric(500, { min: 0, max: 10000, step: 100 }),
    }
})({
    apply({ a, params }) {
        return Task.create('Build Structure', async ctx => {
            if (params.transforms.length === 0) return StateObject.Null;

            const model = a.data;
            const label = params.label || model.label;

            const base = Structure.ofModel(a.data);

            let structure: Structure;
            if (params.transforms.length === 1 && Mat4.isIdentity(params.transforms[0])) {
                const mergedUnits = partitionUnits(base.units, params.cellSize);
                structure = Structure.create(mergedUnits, { label });
            } else {
                const assembler = Structure.Builder({ label });
                const unit = mergeUnits(base.units, 0);
                for (let i = 0, il = params.transforms.length; i < il; ++i) {
                    const t = params.transforms[i];
                    const op = SymmetryOperator.create(`op-${i}`, t);
                    assembler.addWithOperator(unit, op);
                }
                structure = assembler.getStructure();
            }

            const props = { label, description: Structure.elementDescription(structure) };
            return new SO.Molecule.Structure(structure, props);
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});
