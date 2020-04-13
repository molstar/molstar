/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Unit } from '../../../mol-model/structure';
import { StructureElement } from '../../../mol-model/structure/structure';
import { Elements, isHalogen } from '../../../mol-model/structure/model/properties/atomic/types';
import { ElementSymbol, BondType } from '../../../mol-model/structure/model/types';
import { eachBondedAtom, bondCount, typeSymbol, bondToElementCount } from './util';

function isAromatic(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    // TODO also extend unit.rings with geometry/composition-based aromaticity detection and use it here in addition
    const { offset, edgeProps } = unit.bonds;
    for (let i = offset[index], il = offset[index + 1]; i < il; ++i) {
        if (BondType.is(BondType.Flag.Aromatic, edgeProps.flags[i])) return true;
    }
    return false;
}

function bondToCarbonylCount(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let carbonylCount = 0;
    eachBondedAtom(structure, unit, index, (unit: Unit.Atomic, index: StructureElement.UnitIndex) => {
        if (isCarbonyl(structure, unit, index)) carbonylCount += 1;
    });
    return carbonylCount;
}

//

/**
 * Nitrogen in a quaternary amine
 */
export function isQuaternaryAmine(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return (
        typeSymbol(unit, index) === Elements.N &&
        bondCount(structure, unit, index) === 4 &&
        bondToElementCount(structure, unit, index, Elements.H) === 0
    );
}

/**
 * Nitrogen in a tertiary amine
 */
export function isTertiaryAmine(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex, idealValence: number) {
    return (
        typeSymbol(unit, index) === Elements.N &&
        bondCount(structure, unit, index) === 4 &&
        idealValence === 3
    );
}

/**
 * Nitrogen in an imide
 */
export function isImide(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let flag = false;
    if (typeSymbol(unit, index) === Elements.N &&
        (bondCount(structure, unit, index) - bondToElementCount(structure, unit, index, Elements.H)) === 2
    ) {
        flag = bondToCarbonylCount(structure, unit, index) === 2;
    }
    return flag;
}

/**
 * Nitrogen in an amide
 */
export function isAmide(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let flag = false;
    if (typeSymbol(unit, index) === Elements.N &&
        (bondCount(structure, unit, index) - bondToElementCount(structure, unit, index, Elements.H)) === 2
    ) {
        flag = bondToCarbonylCount(structure, unit, index) === 1;
    }
    return flag;
}

/**
 * Sulfur in a sulfonium group
 */
export function isSulfonium(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return (
        typeSymbol(unit, index) === Elements.S &&
        bondCount(structure, unit, index) === 3 &&
        bondToElementCount(structure, unit, index, Elements.H) === 0
    );
}

/**
 * Sulfur in a sulfonic acid or sulfonate group
 */
export function isSulfonicAcid(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return (
        typeSymbol(unit, index) === Elements.S &&
        bondToElementCount(structure, unit, index, Elements.O) === 3
    );
}

/**
 * Sulfur in a sulfate group
 */
export function isSulfate(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return (
        typeSymbol(unit, index) === Elements.S &&
        bondToElementCount(structure, unit, index, Elements.O) === 4
    );
}

/**
 * Phosphor in a phosphate group
 */
export function isPhosphate (structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return (
        typeSymbol(unit, index) === Elements.P &&
        bondToElementCount(structure, unit, index, Elements.O) === bondCount(structure, unit, index)
    );
}

/**
 * Halogen with one bond to a carbon
 */
export function isHalocarbon (structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return (
        isHalogen(typeSymbol(unit, index)) &&
        bondCount(structure, unit, index) === 1 &&
        bondToElementCount(structure, unit, index, Elements.C) === 1
    );
}

/**
 * Carbon in a carbonyl/acyl group
 *
 * TODO currently only checks intra bonds for group detection
 */
export function isCarbonyl(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let flag = false;
    if (typeSymbol(unit, index) === Elements.C) {
        const { offset, edgeProps, b } = unit.bonds;
        for (let i = offset[index], il = offset[index + 1]; i < il; ++i) {
            if (edgeProps.order[i] === 2 && typeSymbol(unit, b[i] as StructureElement.UnitIndex) === Elements.O) {
                flag = true;
                break;
            }
        }
    }
    return flag;
}

/**
 * Carbon in a carboxylate group
 */
export function isCarboxylate (structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let terminalOxygenCount = 0;
    if (
        typeSymbol(unit, index) === Elements.C &&
        bondToElementCount(structure, unit, index, Elements.O) === 2 &&
        bondToElementCount(structure, unit, index, Elements.C) === 1
    ) {
        eachBondedAtom(structure, unit, index, (unit: Unit.Atomic, index: StructureElement.UnitIndex) => {
            if (
                typeSymbol(unit, index) === Elements.O &&
                bondCount(structure, unit, index) - bondToElementCount(structure, unit, index, Elements.H) === 1
            ) {
                terminalOxygenCount += 1;
            }
        });
    }
    return terminalOxygenCount === 2;
}

/**
 * Carbon in a guanidine group
 */
export function isGuanidine (structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let terminalNitrogenCount = 0;
    if (
        typeSymbol(unit, index) === Elements.C &&
        bondCount(structure, unit, index) === 3 &&
        bondToElementCount(structure, unit, index, Elements.N) === 3
    ) {
        eachBondedAtom(structure, unit, index, (unit: Unit.Atomic, index: StructureElement.UnitIndex) => {
            if (
                bondCount(structure, unit, index) - bondToElementCount(structure, unit, index, Elements.H) === 1
            ) {
                terminalNitrogenCount += 1;
            }
        });
    }
    return terminalNitrogenCount === 2;
}

/**
 * Carbon in a acetamidine group
 */
export function isAcetamidine (structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let terminalNitrogenCount = 0;
    if (
        typeSymbol(unit, index) === Elements.C &&
        bondCount(structure, unit, index) === 3 &&
        bondToElementCount(structure, unit, index, Elements.N) === 2 &&
        bondToElementCount(structure, unit, index, Elements.C) === 1
    ) {
        eachBondedAtom(structure, unit, index, (unit: Unit.Atomic, index: StructureElement.UnitIndex) => {
            if (
                bondCount(structure, unit, index) - bondToElementCount(structure, unit, index, Elements.H) === 1
            ) {
                terminalNitrogenCount += 1;
            }
        });
    }
    return terminalNitrogenCount === 2;
}

const PolarElements = new Set<ElementSymbol>([ 'N', 'O', 'S', 'F', 'CL', 'BR', 'I' ] as ElementSymbol[]);
export function isPolar(element: ElementSymbol) { return PolarElements.has(element); }

export function hasPolarNeighbour (structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let flag = false;
    eachBondedAtom(structure, unit, index, (unit: Unit.Atomic, index: StructureElement.UnitIndex) => {
        if (isPolar(typeSymbol(unit, index))) flag = true;
    });
    return flag;
}

export function hasAromaticNeighbour (structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    let flag = false;
    eachBondedAtom(structure, unit, index, (unit: Unit.Atomic, index: StructureElement.UnitIndex) => {
        if (isAromatic(unit, index)) flag = true;
    });
    return flag;
}
