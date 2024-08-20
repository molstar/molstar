/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../mol-model/structure';
import { DSSPComputationParams, DSSPComputationProps, DefaultDSSPComputationProps, computeUnitDSSP } from './secondary-structure/dssp';
import { SecondaryStructure } from '../../mol-model/structure/model/properties/secondary-structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Unit } from '../../mol-model/structure/structure';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { CustomProperty } from '../common/custom-property';
import { ModelSecondaryStructure } from '../../mol-model-formats/structure/property/secondary-structure';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { Model } from '../../mol-model/structure/model';
import { computeUnitZhangSkolnik } from './secondary-structure/zhang-skolnik';

function getSecondaryStructureParams(_data?: Structure) {
    return {
        type: PD.MappedStatic('auto', {
            'auto': PD.EmptyGroup({ label: 'Automatic' }),
            'model': PD.EmptyGroup({ label: 'Model' }),
            'dssp': PD.Group(DSSPComputationParams, { label: 'DSSP', isFlat: true }),
            'zhang-skolnick': PD.EmptyGroup({ label: 'Zhang-Skolnick' }),
        }, { options: [['auto', 'Automatic'], ['model', 'Model'], ['dssp', 'DSSP'], ['zhang-skolnick', 'Zhang-Skolnick']] })
    };
}

export const SecondaryStructureParams = getSecondaryStructureParams();
export type SecondaryStructureParams = typeof SecondaryStructureParams
export type SecondaryStructureProps = PD.Values<SecondaryStructureParams>

/** Maps `unit.id` to `SecondaryStructure` */
export type SecondaryStructureValue = Map<number, SecondaryStructure>

export const SecondaryStructureProvider: CustomStructureProperty.Provider<SecondaryStructureParams, SecondaryStructureValue> = CustomStructureProperty.createProvider({
    label: 'Secondary Structure',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_computed_secondary_structure',
        // TODO `cifExport` and `symbol`
    }),
    type: 'root',
    defaultParams: SecondaryStructureParams,
    getParams: getSecondaryStructureParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<SecondaryStructureProps>) => {
        const p = { ...PD.getDefaultValues(SecondaryStructureParams), ...props };
        switch (p.type.name) {
            case 'auto': return { value: await computeAuto(data) };
            case 'dssp': return { value: await computeDssp(data, p.type.params) };
            case 'model': return { value: await computeModel(data) };
            case 'zhang-skolnick': return { value: await computeZhangSkolnik(data) };
        }
    }
});

async function computeAuto(structure: Structure): Promise<SecondaryStructureValue> {
    const map = new Map<number, SecondaryStructure>();
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const u = structure.unitSymmetryGroups[i].units[0];
        const m = u.model;
        if ((Model.isFromPdbArchive(m) && Model.isExperimental(m) && !Model.isCoarseGrained(m)) || Model.hasSecondaryStructure(m)) {
            const secondaryStructure = ModelSecondaryStructure.Provider.get(m);
            if (secondaryStructure) map.set(u.invariantId, secondaryStructure);
        } else if (Unit.isAtomic(u) && !Model.isCoarseGrained(m)) {
            const secondaryStructure = await computeUnitDSSP(u, DefaultDSSPComputationProps);
            map.set(u.invariantId, secondaryStructure);
        } else if (Unit.isAtomic(u)) {
            const secondaryStructure = await computeUnitZhangSkolnik(u);
            map.set(u.invariantId, secondaryStructure);
        }
    }
    return map;
}

async function computeDssp(structure: Structure, props: DSSPComputationProps): Promise<SecondaryStructureValue> {
    // TODO take inter-unit hbonds into account for bridge, ladder, sheet assignment
    const map = new Map<number, SecondaryStructure>();
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const u = structure.unitSymmetryGroups[i].units[0];
        if (Unit.isAtomic(u) && !Model.isCoarseGrained(u.model)) {
            const secondaryStructure = await computeUnitDSSP(u, props);
            map.set(u.invariantId, secondaryStructure);
        }
    }
    return map;
}

async function computeZhangSkolnik(structure: Structure): Promise<SecondaryStructureValue> {
    const map = new Map<number, SecondaryStructure>();
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const u = structure.unitSymmetryGroups[i].units[0];
        if (Unit.isAtomic(u)) {
            const secondaryStructure = await computeUnitZhangSkolnik(u);
            map.set(u.invariantId, secondaryStructure);
        }
    }
    return map;
}

async function computeModel(structure: Structure): Promise<SecondaryStructureValue> {
    const map = new Map<number, SecondaryStructure>();
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const u = structure.unitSymmetryGroups[i].units[0];
        if (Unit.isAtomic(u)) {
            const secondaryStructure = ModelSecondaryStructure.Provider.get(u.model);
            if (secondaryStructure) map.set(u.invariantId, secondaryStructure);
        }
    }
    return map;
}