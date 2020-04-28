/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ModelCrossLinkRestraint } from './format';
import { Unit, StructureElement, Structure, Bond} from '../../../mol-model/structure';
import { PairRestraints, PairRestraint } from '../pair-restraints';
import { CustomStructureProperty } from '../../common/custom-structure-property';
import { CustomProperty } from '../../common/custom-property';
import { DataLocation } from '../../../mol-model/location';
import { DataLoci } from '../../../mol-model/loci';
import { Sphere3D } from '../../../mol-math/geometry';
import { CentroidHelper } from '../../../mol-math/geometry/centroid-helper';
import { bondLabel } from '../../../mol-theme/label';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';

export type CrossLinkRestraintValue = PairRestraints<CrossLinkRestraint>

export const CrossLinkRestraintProvider: CustomStructureProperty.Provider<{}, CrossLinkRestraintValue> = CustomStructureProperty.createProvider({
    label: 'Cross Link Restraint',
    descriptor: CustomPropertyDescriptor({
        name: 'integrative-cross-link-restraint',
        // TODO `cifExport` and `symbol`
    }),
    type: 'local',
    defaultParams: {},
    getParams: (data: Structure) => ({}),
    isApplicable: (data: Structure) => data.models.some(m => !!ModelCrossLinkRestraint.Provider.get(m)),
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<{}>) => {
        return { value: extractCrossLinkRestraints(data) };
    }
});

export { CrossLinkRestraint };

interface CrossLinkRestraint extends PairRestraint {
    readonly restraintType: 'harmonic' | 'upper bound' | 'lower bound'
    readonly distanceThreshold: number
    readonly psi: number
    readonly sigma1: number
    readonly sigma2: number
}

namespace CrossLinkRestraint {
    export enum Tag {
        CrossLinkRestraint = 'cross-link-restraint'
    }

    export function isApplicable(structure: Structure) {
        return structure.models.some(m => !!ModelCrossLinkRestraint.Provider.get(m));
    }

    const distVecA = Vec3(), distVecB = Vec3();
    export function distance(pair: CrossLinkRestraint) {
        pair.unitA.conformation.position(pair.unitA.elements[pair.indexA], distVecA);
        pair.unitB.conformation.position(pair.unitB.elements[pair.indexB], distVecB);
        return Vec3.distance(distVecA, distVecB);
    }

    type StructureCrossLinkRestraints = { readonly structure: Structure, readonly crossLinkRestraints: CrossLinkRestraintValue }

    export type Element = number
    export interface Location extends DataLocation<StructureCrossLinkRestraints, Element> {}

    export function Location(crossLinkRestraints: CrossLinkRestraintValue, structure: Structure, index?: number): Location {
        return DataLocation('cross-link-restraints', { structure, crossLinkRestraints }, index as any);
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'data-location' && x.tag === 'cross-link-restraints';
    }

    export function areLocationsEqual(locA: Location, locB: Location) {
        return (
            locA.data.structure === locB.data.structure &&
            locA.data.crossLinkRestraints === locB.data.crossLinkRestraints &&
            locA.element === locB.element
        );
    }

    function _label(crossLinkRestraints: CrossLinkRestraintValue, element: Element): string {
        const p = crossLinkRestraints.pairs[element];
        return `Cross Link Restraint | Type: ${p.restraintType} | Threshold: ${p.distanceThreshold} \u212B | Psi: ${p.psi} | Sigma 1: ${p.sigma1} | Sigma 2: ${p.sigma2} | Distance: ${distance(p).toFixed(2)} \u212B`;
    }

    export function locationLabel(location: Location): string {
        return _label(location.data.crossLinkRestraints, location.element);
    }

    export interface Loci extends DataLoci<StructureCrossLinkRestraints, Element> { }

    export function Loci(structure: Structure, crossLinkRestraints: CrossLinkRestraintValue, elements: ReadonlyArray<Element>): Loci {
        return DataLoci('cross-link-restraints', { structure, crossLinkRestraints }, elements,
            (boundingSphere) => getBoundingSphere(crossLinkRestraints, elements, boundingSphere),
            () => getLabel(structure, crossLinkRestraints, elements));
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'data-loci' && x.tag === 'interactions';
    }

    export function getBoundingSphere(crossLinkRestraints: CrossLinkRestraintValue, elements: ReadonlyArray<Element>, boundingSphere: Sphere3D) {
        return CentroidHelper.fromPairProvider(elements.length, (i, pA, pB) => {
            const p = crossLinkRestraints.pairs[elements[i]];
            p.unitA.conformation.position(p.unitA.elements[p.indexA], pA);
            p.unitB.conformation.position(p.unitB.elements[p.indexB], pB);
        }, boundingSphere);
    }

    export function getLabel(structure: Structure, crossLinkRestraints: CrossLinkRestraintValue, elements: ReadonlyArray<Element>) {
        const element = elements[0];
        if (element === undefined) return '';
        const p = crossLinkRestraints.pairs[element];
        return [
            _label(crossLinkRestraints, element),
            bondLabel(Bond.Location(structure, p.unitA, p.indexA, structure, p.unitB, p.indexB))
        ].join('</br>');
    }
}

//

function _addRestraints(map: Map<number, number>, unit: Unit, restraints: ModelCrossLinkRestraint) {
    const { elements } = unit;
    const elementCount = elements.length;
    const kind = unit.kind;

    for (let i = 0; i < elementCount; i++) {
        const e = elements[i];
        restraints.getIndicesByElement(e, kind).forEach(ri => map.set(ri, i));
    }
}

function extractInter(pairs: CrossLinkRestraint[], unitA: Unit, unitB: Unit) {
    if (unitA.model !== unitB.model) return;
    if (unitA.model.sourceData.kind !== 'mmCIF') return;

    const restraints = ModelCrossLinkRestraint.Provider.get(unitA.model);
    if (!restraints) return;

    const rA = new Map<number, StructureElement.UnitIndex>();
    const rB = new Map<number, StructureElement.UnitIndex>();
    _addRestraints(rA, unitA, restraints);
    _addRestraints(rB, unitB, restraints);

    rA.forEach((indexA, ri) => {
        const indexB = rB.get(ri);
        if (indexB !== undefined) {
            pairs.push(
                createCrossLinkRestraint(unitA, indexA, unitB, indexB, restraints, ri),
                createCrossLinkRestraint(unitB, indexB, unitA, indexA, restraints, ri)
            );
        }
    });
}

function extractIntra(pairs: CrossLinkRestraint[], unit: Unit) {
    if (unit.model.sourceData.kind !== 'mmCIF') return;

    const restraints = ModelCrossLinkRestraint.Provider.get(unit.model);
    if (!restraints) return;

    const { elements } = unit;
    const elementCount = elements.length;
    const kind = unit.kind;

    const r = new Map<number, StructureElement.UnitIndex[]>();

    for (let i = 0; i < elementCount; i++) {
        const e = elements[i];
        restraints.getIndicesByElement(e, kind).forEach(ri => {
            const il = r.get(ri);
            if (il) il.push(i as StructureElement.UnitIndex);
            else r.set(ri, [i as StructureElement.UnitIndex]);
        });
    }

    r.forEach((il, ri) => {
        if (il.length < 2) return;
        const [ indexA, indexB ] = il;
        pairs.push(
            createCrossLinkRestraint(unit, indexA, unit, indexB, restraints, ri),
            createCrossLinkRestraint(unit, indexB, unit, indexA, restraints, ri)
        );
    });
}

function createCrossLinkRestraint(unitA: Unit, indexA: StructureElement.UnitIndex, unitB: Unit, indexB: StructureElement.UnitIndex, restraints: ModelCrossLinkRestraint, row: number): CrossLinkRestraint {
    return {
        unitA, indexA, unitB, indexB,

        restraintType: restraints.data.restraint_type.value(row),
        distanceThreshold: restraints.data.distance_threshold.value(row),
        psi: restraints.data.psi.value(row),
        sigma1: restraints.data.sigma_1.value(row),
        sigma2: restraints.data.sigma_2.value(row),
    };
}

function extractCrossLinkRestraints(structure: Structure): PairRestraints<CrossLinkRestraint> {
    const pairs: CrossLinkRestraint[] = [];
    if (!structure.models.some(m => ModelCrossLinkRestraint.Provider.get(m))) {
        return new PairRestraints(pairs);
    }

    const n = structure.units.length;
    for (let i = 0; i < n; ++i) {
        const unitA = structure.units[i];
        extractIntra(pairs, unitA);
        for (let j = i + 1; j < n; ++j) {
            const unitB = structure.units[j];
            if (unitA.model === unitB.model) {
                extractInter(pairs, unitA, unitB);
            }
        }
    }

    return new PairRestraints(pairs);
}