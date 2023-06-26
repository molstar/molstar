/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet, SortedArray } from '../../../../mol-data/int';
import { Box3D, GridLookup3D, PositionData, Sphere3D, SymmetryOperator } from '../../../../mol-math/geometry';
import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
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
                    console.log(`missing asymid '${id}'`);
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

function mergeUnits(units: Unit[], id: number): Unit {
    const u = units[0];

    let start = -1 as ElementIndex, end = -1 as ElementIndex;
    let elements = SortedArray.Empty as SortedArray<ElementIndex>;

    for (let i = 0, il = units.length; i < il; ++i) {
        const e = units[i].elements;
        if (SortedArray.isRange(e)) {
            if (end !== -1 && e[0] === end + 1) {
                // extend range
                end = e[e.length - 1];
            } else {
                if (end !== -1) {
                    // pending range
                    elements = SortedArray.union(elements, SortedArray.ofRange(start, end));
                }
                // new range
                start = e[0];
                end = e[e.length - 1];
            }
        } else {
            if (end !== -1) {
                // pending range
                elements = SortedArray.union(elements, SortedArray.ofRange(start, end));
                start = -1 as ElementIndex, end = -1 as ElementIndex;
            }
            elements = SortedArray.union(elements, e);
        }
    }

    if (end !== -1) {
        // pending range
        elements = SortedArray.union(elements, SortedArray.ofRange(start, end));
    }

    return Unit.create(id, id, 0, u.traits | Unit.Trait.MultiChain, u.kind, u.model, u.conformation.operator, elements);
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
                const mergedUnits: Unit[] = [];
                const cellSize = params.cellSize;

                const box = Box3D.setEmpty(Box3D());
                const x = new Float32Array(unitCount);
                const y = new Float32Array(unitCount);
                const z = new Float32Array(unitCount);
                const indices = OrderedSet.ofBounds(0, unitCount);

                for (let i = 0, il = unitCount; i < il; ++i) {
                    const v = units[i].boundary.sphere.center;
                    x[i] = v[0];
                    y[i] = v[1];
                    z[i] = v[2];
                    Box3D.add(box, v);
                }
                Box3D.expand(box, box, Vec3.create(1, 1, 1));

                const positionData: PositionData = { x, y, z, indices };
                const boundary = { box, sphere: Sphere3D.fromBox3D(Sphere3D(), box) };
                const lookup = GridLookup3D(positionData, boundary, Vec3.create(cellSize, cellSize, cellSize));

                const { array, offset, count } = lookup.buckets;

                for (let i = 0, il = offset.length; i < il; ++i) {
                    const start = offset[i];
                    const size = count[i];
                    const cellUnits: Unit[] = [];
                    for (let j = start, jl = start + size; j < jl; ++j) {
                        cellUnits.push(units[array[j]]);
                    }
                    mergedUnits.push(mergeUnits(cellUnits, i));
                }

                structure = Structure.create(mergedUnits);
            } else {
                structure = Structure.create(units);
            }

            const label = entities.data.pdbx_description.value(idx).join(', ') || 'model';
            return new PSO.Molecule.Structure(structure, { label, description: `${a.description}` });
        });
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});
