/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
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
import { Task } from '../../../../mol-task';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';

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

function buildCellpackAssembly(model: Model, assembly: Assembly) {
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

export { CellpackAssembly };
type CellpackAssembly = typeof CellpackAssembly
const CellpackAssembly = PluginStateTransform.BuiltIn({
    name: 'cellpack-assembly',
    display: { name: 'Cellpack Assembly' },
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

            const s = buildCellpackAssembly(model, asm);

            const objProps = { label: `Assembly ${id}`, description: Structure.elementDescription(s) };
            return new PSO.Molecule.Structure(s, objProps);
        });
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
    const map: UnitsByEntity = new Map();
    for (const ug of structure.unitSymmetryGroups) {
        const u = ug.units[0] as Unit.Atomic;
        const e = atomicIndex.getEntityFromChain(u.chainIndex[u.elements[0]]);

        if (!map.has(e)) map.set(e, []);
        const entityUnits = map.get(e)!;

        for (let i = 0, il = ug.units.length; i < il; ++i) {
            entityUnits.push(ug.units[i]);
        }
    }

    UnitsByEntity.set(structure, { value: map }, map);
    return map;
}

export { CellpackStructure };
type CellpackStructure = typeof CellpackStructure
const CellpackStructure = PluginStateTransform.BuiltIn({
    name: 'cellpack-structure',
    display: { name: 'Cellpack Structure' },
    from: PSO.Root,
    to: PSO.Molecule.Structure,
    params: {
        structureRef: PD.Text(''),
        entityId: PD.Text('')
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

            const structure = Structure.create(units);
            const description_label = entities.data.pdbx_description.value(idx)[0] || 'model';
            const label = description_label.split('.').at(-1) || a.label;
            const description = entities.data.pdbx_parent_entity_id.value(idx) || label;
            return new PSO.Molecule.Structure(structure, { label, description: description }); // `${a.description}`
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});
