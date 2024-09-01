/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SecondaryStructure } from '../../../mol-model/structure/model/properties/secondary-structure';
import { SecondaryStructureType } from '../../../mol-model/structure/model/types';
import { Unit } from '../../../mol-model/structure';
import { ElementIndex, ResidueIndex } from '../../../mol-model/structure/model';
import { SortedArray } from '../../../mol-data/int';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';

const HelixDistances = [5.45, 5.18, 6.37];
const HelixDelta = 2.1;

const SheetDistances = [6.1, 10.4, 13.0];
const SheetDelta = 1.42;

const posA = Vec3();
const posB = Vec3();

function zhangSkolnickAtomicSS(unit: Unit.Atomic, residueIndices: SortedArray<ResidueIndex>, i: number, distances: number[], delta: number) {
    const c = unit.conformation;
    const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;

    for (let j = Math.max(0, i - 2); j <= i; ++j) {
        for (let k = 2; k < 5; ++k) {
            if (j + k >= residueIndices.length) return false;

            const rA = residueIndices[j];
            const rB = residueIndices[j + k];

            const aA = traceElementIndex[rA];
            const aB = traceElementIndex[rB];
            if (aA === -1 || aB === -1) return false;

            c.invariantPosition(aA as ElementIndex, posA);
            c.invariantPosition(aB as ElementIndex, posB);
            const d = Vec3.distance(posA, posB);

            if (Math.abs(d - distances[k - 2]) >= delta) return false;
        }
    }

    return true;
}

/**
 * Secondary-structure assignment based on Zhang and Skolnick's TM-align paper.
 * TM-align: a protein structure alignment algorithm based on the Tm-score (2005) NAR, 33(7) 2302-2309.
 *
 * While not as accurate as DSSP, it is faster and works for coarse-grained/backbone-only models.
 */
export async function computeUnitZhangSkolnik(unit: Unit.Atomic): Promise<SecondaryStructure> {
    const count = unit.proteinElements.length;
    const type = new Uint32Array(count) as unknown as SecondaryStructureType[];
    const keys: number[] = [];
    const elements: SecondaryStructure.Element[] = [];

    const { proteinElements, residueIndex } = unit;
    const residueCount = proteinElements.length;
    const unitProteinResidues = new Uint32Array(residueCount);
    for (let i = 0; i < residueCount; ++i) {
        const rI = residueIndex[proteinElements[i]];
        unitProteinResidues[i] = rI;
    }
    const residueIndices = SortedArray.ofSortedArray<ResidueIndex>(unitProteinResidues);
    const getIndex = (rI: ResidueIndex) => SortedArray.indexOf(residueIndices, rI);

    for (let i = 0, il = residueIndices.length; i < il; ++i) {
        let flag = SecondaryStructureType.Flag.None;
        if (zhangSkolnickAtomicSS(unit, residueIndices, i, HelixDistances, HelixDelta)) {
            flag = SecondaryStructureType.Flag.Helix;
        } else if (zhangSkolnickAtomicSS(unit, residueIndices, i, SheetDistances, SheetDelta)) {
            flag = SecondaryStructureType.Flag.Beta;
        }
        type[i] = flag;
        if (elements.length === 0 || flag !== getFlag(elements[elements.length - 1])) {
            elements[elements.length] = createElement(mapToKind(flag), flag);
        }
        keys[i] = elements.length - 1;
    }

    return SecondaryStructure(type, keys, elements, getIndex);
}

function createElement(kind: string, flag: SecondaryStructureType.Flag): SecondaryStructure.Element {
    if (kind === 'helix') {
        return { kind: 'helix', flags: flag } as SecondaryStructure.Helix;
    } else if (kind === 'sheet') {
        return { kind: 'sheet', flags: flag } as SecondaryStructure.Sheet;
    } else {
        return { kind: 'none' };
    }
}

function mapToKind(flag: SecondaryStructureType.Flag) {
    if (flag === SecondaryStructureType.Flag.Helix) {
        return 'helix';
    } else if (flag === SecondaryStructureType.Flag.Beta) {
        return 'sheet';
    } else {
        return 'none';
    }
}

function getFlag(element: SecondaryStructure.Element) {
    if (element.kind === 'helix') {
        return element.flags;
    } else if (element.kind === 'sheet') {
        return element.flags;
    } else {
        return SecondaryStructureType.Flag.None;
    }
}
