/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from '../../../../mol-data/int';
import { SymmetryOperator } from '../../../../mol-math/geometry';
import { Mat4 } from '../../../../mol-math/linear-algebra';
import { ModelSymmetry } from '../../../../mol-model-formats/structure/property/symmetry';
import { CustomStructureProperty } from '../../../../mol-model-props/common/custom-structure-property';
import { ElementIndex, EntityIndex, Model, Structure, Unit } from '../../../../mol-model/structure';
import { Assembly, Symmetry } from '../../../../mol-model/structure/model/properties/symmetry';
import { PluginStateObject as PSO, PluginStateTransform } from '../../../../mol-plugin-state/objects';
import { PluginContext } from '../../../../mol-plugin/context';
import { StateTransformer } from '../../../../mol-state/transformer';
import { Task } from '../../../../mol-task';
import { deepEqual } from '../../../../mol-util';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { partitionUnits } from '../util';

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

export { MmcifAssembly };
type MmcifAssembly = typeof MmcifAssembly
const MmcifAssembly = PluginStateTransform.BuiltIn({
    name: 'mmcif-assembly',
    display: { name: 'Mmcif Assembly' },
    from: PSO.Molecule.Model,
    to: PSO.Molecule.Structure,
    params: {
        id: PD.Text('', { label: 'Asm Id', description: 'Assembly Id (use empty for the 1st assembly)' }),
    }
})({
    canAutoUpdate({ newParams }) {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('Build Structure', async ctx => {
            const model = a.data;

            let id = params.id;
            let asm: Assembly | undefined = void 0;
            const symmetry = ModelSymmetry.Provider.get(model);

            // if no id is specified, use the 1st assembly.
            if (!id && symmetry && symmetry.assemblies.length !== 0) {
                id = symmetry.assemblies[0].id;
            }

            if (!symmetry || symmetry.assemblies.length === 0) {
                plugin.log.warn(`Model '${model.entryId}' has no assembly, returning model structure.`);
            } else {
                asm = Symmetry.findAssembly(model, id || '');
                if (!asm) {
                    plugin.log.warn(`Model '${model.entryId}' has no assembly called '${id}', returning model structure.`);
                }
            }

            const base = Structure.ofModel(model);
            if (!asm) {
                const label = { label: 'Model', description: Structure.elementDescription(base) };
                return new PSO.Molecule.Structure(base, label);
            }

            const s = buildAssembly(model, asm);

            const objProps = { label: `Assembly ${id}`, description: Structure.elementDescription(s) };
            return new PSO.Molecule.Structure(s, objProps);
        });
    },
    update({ newParams, oldParams }) {
        return deepEqual(newParams, oldParams)
            ? StateTransformer.UpdateResult.Unchanged
            : StateTransformer.UpdateResult.Recreate;
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});

type UnitsByEntity = Map<EntityIndex, Unit[]>;
const UnitsByEntity = CustomStructureProperty.createSimple<UnitsByEntity>('units_by_entity', 'root');

function getUnitsByEntity(structure: Structure): UnitsByEntity {
    if (UnitsByEntity.get(structure).value) {
        return UnitsByEntity.get(structure).value!;
    }

    const atomicIndex = structure.model.atomicHierarchy.index;
    const spheresIndex = structure.model.coarseHierarchy.spheres;
    const map: UnitsByEntity = new Map();
    for (const ug of structure.unitSymmetryGroups) {
        const u = ug.units[0];
        let e: EntityIndex;
        if (Unit.isAtomic(u)) {
            e = atomicIndex.getEntityFromChain(u.chainIndex[u.elements[0]]);
        } else if (Unit.isSpheres(u)) {
            e = spheresIndex.getEntityFromChain(u.coarseElements.chainElementSegments.index[u.elements[0]]);
        } else {
            continue;
        }

        if (!map.has(e)) map.set(e, []);
        const entityUnits = map.get(e)!;

        for (let i = 0, il = ug.units.length; i < il; ++i) {
            entityUnits.push(ug.units[i]);
        }
    }

    UnitsByEntity.set(structure, { value: map }, map);
    return map;
}

export { MmcifStructure };
type MmcifStructure = typeof MmcifStructure
const MmcifStructure = PluginStateTransform.BuiltIn({
    name: 'mmcif-structure',
    display: { name: 'Mmcif Structure' },
    from: PSO.Root,
    to: PSO.Molecule.Structure,
    params: {
        structureRef: PD.Text(''),
        entityId: PD.Text(''),
        cellSize: PD.Numeric(500, { min: 0, max: 10000, step: 100 }),
    }
})({
    canAutoUpdate({ newParams }) {
        return true;
    },
    apply({ a, params, dependencies }) {
        return Task.create('Build Structure', async ctx => {
            const parent = dependencies![params.structureRef].data as Structure;
            const { entities } = parent.model;
            const idx = entities.getEntityIndex(params.entityId);

            const unitsByEntity = getUnitsByEntity(parent);
            const units = unitsByEntity.get(idx) || [];
            const unitCount = units.length;

            let structure: Structure;
            if (unitCount > 1 && units.every(u => u.conformation.operator.isIdentity)) {
                const mergedUnits = partitionUnits(units, params.cellSize);
                structure = Structure.create(mergedUnits);
            } else {
                structure = Structure.create(units);
            }
            // could also use _struct_ref.pdbx_db_accession to point to uniprot with _struct_ref.db_name == UNP
            const label = entities.data.pdbx_description.value(idx).join(', ') || 'model';
            const description = `*Entity id* ${entities.data.id.value(idx)} *src_method* ${entities.data.src_method.value(idx)} *type* ${entities.data.type.value(idx)}`;
            return new PSO.Molecule.Structure(structure, { label, description: description });
        });
    },
    update({ newParams, oldParams }) {
        return deepEqual(newParams, oldParams)
            ? StateTransformer.UpdateResult.Unchanged
            : StateTransformer.UpdateResult.Recreate;
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});
