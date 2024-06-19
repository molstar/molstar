/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';
import { ElementIndex, Model, Structure, Unit } from '../../../../mol-model/structure';
import { PluginStateObject as SO, PluginStateTransform } from '../../../../mol-plugin-state/objects';
import { Task } from '../../../../mol-task';
import { StateObject, StateTransformer } from '../../../../mol-state';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { SymmetryOperator } from '../../../../mol-math/geometry';
import { mergeUnits, partitionUnits } from '../util';
import { Assembly, Symmetry } from '../../../../mol-model/structure/model/properties/symmetry';
import { ModelSymmetry } from '../../../../mol-model-formats/structure/property/symmetry';
import { SortedArray } from '../../../../mol-data/int';
import { GenericInstances, getTransforms } from './preset';
import { Asset } from '../../../../mol-util/assets';
import { PluginContext } from '../../../../mol-plugin/context';
import { deepEqual } from '../../../../mol-util';

function createModelChainMap(model: Model) {
    const builder = new Structure.StructureBuilder();
    const units = new Map<string, Unit>();

    const { label_asym_id, _rowCount } = model.atomicHierarchy.chains;
    const { offsets } = model.atomicHierarchy.chainAtomSegments;

    for (let i = 0; i < _rowCount; i++) {
        const elements = SortedArray.ofBounds(offsets[i] as ElementIndex, offsets[i + 1] as ElementIndex);
        const unit = builder.addUnit(Unit.Kind.Atomic, model, SymmetryOperator.Default, elements, Unit.Trait.FastBoundary);
        units.set(label_asym_id.value(i), unit);
    }

    return units;
}

function buildAssembly(model: Model, assembly: Assembly) {
    const coordinateSystem = SymmetryOperator.create(assembly.id, Mat4.identity(), { assembly: { id: assembly.id, operId: 0, operList: [] } });
    const assembler = Structure.Builder({
        coordinateSystem,
        label: model.label,
    });

    const units = createModelChainMap(model);

    for (const g of assembly.operatorGroups) {
        for (const oper of g.operators) {
            for (const id of g.asymIds!) {
                const u = units.get(id);
                if (u) {
                    assembler.addWithOperator(u, oper);
                } else {
                    console.log(`missing asymId '${id}'`);
                }
            }
        }
    }

    return assembler.getStructure();
}

const EmptyInstances: GenericInstances<Asset> = {
    positions: { data: [] },
    rotations: { variant: 'euler', data: [] }
};

export { StructureFromGeneric };
type StructureFromGeneric = typeof StructureFromGeneric
const StructureFromGeneric = PluginStateTransform.BuiltIn({
    name: 'structure-from-generic',
    display: { name: 'Structure from Generic', description: 'Create a molecular structure from Generic models.' },
    from: SO.Molecule.Model,
    to: SO.Molecule.Structure,
    params: {
        instances: PD.Value<GenericInstances<Asset>>(EmptyInstances),
        label: PD.Optional(PD.Text('')),
        description: PD.Optional(PD.Text('')),
        cellSize: PD.Numeric(500, { min: 0, max: 10000, step: 100 }),
    }
})({
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Build Structure', async ctx => {
            const transforms = await getTransforms(plugin, params.instances);
            if (transforms.length === 0) return StateObject.Null;

            const model = a.data;
            const label = params.label || model.label;

            const base = Structure.ofModel(a.data);

            let structure: Structure;
            if (transforms.length === 1 && Mat4.isIdentity(transforms[0])) {
                const symmetry = ModelSymmetry.Provider.get(model);
                const id = symmetry?.assemblies[0]?.id;
                const asm = Symmetry.findAssembly(model, id || '');
                if (asm) {
                    structure = buildAssembly(model, asm);
                } else {
                    const mergedUnits = partitionUnits(base.units, params.cellSize);
                    structure = Structure.create(mergedUnits, { label });
                }
            } else {
                const assembler = Structure.Builder({ label });
                const unit = mergeUnits(base.units, 0);
                for (let i = 0, il = transforms.length; i < il; ++i) {
                    const t = transforms[i];
                    const op = SymmetryOperator.create(`op-${i}`, t);
                    assembler.addWithOperator(unit, op);
                }
                structure = assembler.getStructure();
            }

            const props = { label, description: params.description || Structure.elementDescription(structure) };
            return new SO.Molecule.Structure(structure, props);
        });
    },
    update({ newParams, oldParams }, plugin: PluginContext) {
        if (deepEqual(newParams, oldParams)) {
            return StateTransformer.UpdateResult.Unchanged;
        }

        if (oldParams.instances) releaseInstances(plugin, oldParams.instances);
        return StateTransformer.UpdateResult.Recreate;
    },
    dispose({ b, params }, plugin: PluginContext) {
        b?.data.customPropertyDescriptors.dispose();
        if (params?.instances) releaseInstances(plugin, params.instances);
    }
});

function releaseInstances(plugin: PluginContext, instances: GenericInstances<Asset>) {
    if (!Array.isArray(instances.positions.data)) {
        plugin.managers.asset.release(instances.positions.data.file);
    }
    if (!Array.isArray(instances.rotations.data)) {
        plugin.managers.asset.release(instances.rotations.data.file);
    }
}
