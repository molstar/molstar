import Unit from '../../unit'
import { Vec3 } from 'mol-math/linear-algebra';
import { AccessibleSurfaceAreaComputationParameters, AccessibleSurfaceArea } from './data';
import { isHydrogen, getElementIdx } from '../links/common'; // TODO move these functions somewhere more common
import { MoleculeType, ElementSymbol, MaxAsa, DefaultMaxAsa } from 'mol-model/structure/model/types';
import { VdwRadius } from 'mol-model/structure/model/properties/atomic/measures';

const trigonalCarbonVdw = 1.76;
const tetrahedralCarbonVdw = 1.87;
const trigonalNitrogenVdw = 1.65;
const tetrahedralNitrogenVdw = 1.50;
/** deviating radii from the definition in types.ts */
const oxygenVdw = 1.4;
const sulfurVdw = 1.85;
const missingAccessibleSurfaceAreaValue = -1.0;

function _computeAccessibleSurfaceArea(unit: Unit.Atomic, params: AccessibleSurfaceAreaComputationParameters): AccessibleSurfaceArea {
    const ctx = initialize(unit, params);
    assignRadiusForHeavyAtoms(ctx);
    computePerResidue(ctx);
    normalizeAccessibleSurfaceArea(ctx);

    return {
        atomRadius: ctx.atomRadius,
        accessibleSurfaceArea: ctx.accessibleSurfaceArea,
        relativeAccessibleSurfaceArea: ctx.relativeAccessibleSurfaceArea,
        buried: void 0 // TODO impl - rasa < 0.16 - find Rost reference
    };
}

function normalizeAccessibleSurfaceArea(ctx: AccessibleSurfaceAreaContext) {
    const { residues, derived } = ctx.unit.model.atomicHierarchy;
    const { accessibleSurfaceArea, relativeAccessibleSurfaceArea } = ctx;

    for (let i = 0; i < residues.label_comp_id.rowCount; ++i) {
        // skip entities not part of a peptide chain
        if (derived.residue.moleculeType[i] !== MoleculeType.protein) continue;

        const maxAsa = (MaxAsa as any)[residues.label_comp_id.value(i)];
        relativeAccessibleSurfaceArea[i] = accessibleSurfaceArea[i] / (maxAsa === undefined ? DefaultMaxAsa : maxAsa);
    }
}

/**
 * notes on performance - scenario: compute for first 10 units of 3j3q @ 960 sphere points
 * lookup3d + refinement: ~5000ms
 * naive approach: ~5600ms - higher variance
 */
function computePerResidue(ctx: AccessibleSurfaceAreaContext) { // runs at roughly 5000 ms
    const { atomRadius, accessibleSurfaceArea, maxLookupRadius, spherePoints, cons } = ctx;
    const { probeSize } = ctx.params;
    const { elements: atoms, residueIndex } = ctx.unit;
    const { x, y, z } = ctx.unit.model.atomicConformation;
    const atomCount = ctx.unit.elements.length;
    const { lookup3d } = ctx.unit;

    const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i]);
    const a1Pos = Vec3.zero();
    const a2Pos = Vec3.zero();

    for (let _aI = 0; _aI < atomCount; ++_aI) {
        // map the atom index of this unit to the actual 'global' atom index
        const aI = atoms[_aI];
        const radii1 = atomRadius[aI];
        if (radii1 === missingAccessibleSurfaceAreaValue) continue;

        // find suitable neighbor candidates by lookup
        const { indices, count } = lookup3d.find(x[aI], y[aI], z[aI], maxLookupRadius);
        position(aI, a1Pos);

        // refine list by actual criterion
        const cutoff = probeSize + probeSize + radii1;
        const filteredIndicies = [];
        for (let ni = 0; ni < count; ni++) {
            const _bI = indices[ni];
            const bI = atoms[_bI];
            const radii2 = atomRadius[bI];
            if (bI === aI || radii2 === missingAccessibleSurfaceAreaValue) continue;

            const cutoff2 = (cutoff + radii2) * (cutoff + radii2);
            // accurately check for neighborhood
            position(bI, a2Pos);
            if (Vec3.squaredDistance(a1Pos, a2Pos) < cutoff2) {
                filteredIndicies[filteredIndicies.length] = bI;
            }
        }

        // test sphere points
        const r = probeSize + radii1;
        let accessiblePointCount = 0;
        for (let si = 0; si < spherePoints.length; ++si) {
            const spherePoint = spherePoints[si];
            const testPoint = [spherePoint[0] * r + a1Pos[0], spherePoint[1] * r + a1Pos[1], spherePoint[2] * r + a1Pos[2]] as Vec3;
            let accessible = true;

            for (let ni = 0; ni < filteredIndicies.length; ++ni) {
                const naI = filteredIndicies[ni];
                position(naI, a2Pos);
                const cutoff3 = (atomRadius[naI] + probeSize) * (atomRadius[naI] + probeSize);
                if (Vec3.squaredDistance(testPoint, a2Pos) < cutoff3) {
                    accessible = false;
                    break;
                }
            }

            if (accessible) ++accessiblePointCount;
        }

        const value = cons * accessiblePointCount * r * r;
        accessibleSurfaceArea[residueIndex[aI]] += value;
        // +30% computation by normalizing partial solutions
        // relativeAccessibleSurfaceArea[residueIndex[aI]] += value * (NormalizationFactors as any)[residueIndex[aI]];
    }
}

function assignRadiusForHeavyAtoms(ctx: AccessibleSurfaceAreaContext) {
    const atomCount = ctx.unit.elements.length;
    const { elements: atoms, residueIndex } = ctx.unit;
    const { moleculeType } = ctx.unit.model.atomicHierarchy.derived.residue;
    const { type_symbol, label_atom_id } = ctx.unit.model.atomicHierarchy.atoms;
    const { label_comp_id } = ctx.unit.model.atomicHierarchy.residues;

    for (let _aI = 0; _aI < atomCount; ++_aI) {
        const aI =  atoms[_aI];
        const raI = residueIndex[aI];
        const aeI = getElementIdx(type_symbol.value(aI)!);

        // skip hydrogen atoms
        if (isHydrogen(aeI)) {
            ctx.atomRadius[ctx.atomRadius.length] = missingAccessibleSurfaceAreaValue;
            continue;
        }

        // skip non-peptide groups
        if (moleculeType[raI] !== MoleculeType.protein) {
            ctx.atomRadius[ctx.atomRadius.length] = missingAccessibleSurfaceAreaValue;
            continue;
        }

        const atomId = label_atom_id.value(aI);
        const element = type_symbol.value(aI);
        const resn = label_comp_id.value(raI)!;

        ctx.atomRadius[aI] = determineRadius(atomId, element, resn);

        // while having atom->parent mapping at hand, initialize residue-level of information
        ctx.accessibleSurfaceArea[raI] = 0.0;
        ctx.relativeAccessibleSurfaceArea[raI] = 0.0;
    }
}

/**
 * Gets the van der Waals radius of the given atom following the values defined by Chothia (1976)
 * J.Mol.Biol.105,1-14. NOTE: the vdw values defined by the paper assume no Hydrogens and thus "inflates" slightly
 * the heavy atoms to account for Hydrogens. Thus this method cannot be used in a structure that contains Hydrogens!
 */
function determineRadius(atomId: string, element: ElementSymbol, compId: string): number {
    switch (element) {
        case 'O':
        return oxygenVdw;
        case 'S':
        return sulfurVdw;
        case 'N':
        return atomId === 'NZ' ? tetrahedralNitrogenVdw : trigonalNitrogenVdw;
        case 'C':
        switch (atomId) {
            case 'C': case 'CE1': case'CE2': case 'CE3': case 'CH2': case 'CZ': case 'CZ2': case 'CZ3':
            return trigonalCarbonVdw;
            case 'CA': case 'CB': case 'CE': case 'CG1': case 'CG2':
            return tetrahedralCarbonVdw;
            default:
            switch (compId) {
                case 'PHE': case 'TRP': case 'TYR': case 'HIS': case 'ASP': case 'ASN':
                return trigonalCarbonVdw;
                case 'PRO': case 'LYS': case 'ARG': case 'MET': case 'ILE': case 'LEU':
                return tetrahedralCarbonVdw;
                case 'GLU': case 'GLN':
                return atomId === 'CD' ? trigonalCarbonVdw : tetrahedralCarbonVdw;
            }
        }
    }
    return VdwRadius(element);
}

interface AccessibleSurfaceAreaContext {
    unit: Unit.Atomic,
    params: AccessibleSurfaceAreaComputationParameters,
    spherePoints: Vec3[],
    cons: number,
    atomRadius: number[],
    accessibleSurfaceArea: number[],
    relativeAccessibleSurfaceArea: number[],
    maxLookupRadius: number
}

function initialize(unit: Unit.Atomic, params: AccessibleSurfaceAreaComputationParameters): AccessibleSurfaceAreaContext {
    return {
        unit: unit,
        params: params,
        spherePoints: generateSpherePoints(params.numberOfSpherePoints),
        cons: 4.0 * Math.PI / params.numberOfSpherePoints,
        atomRadius: [],
        accessibleSurfaceArea: [],
        relativeAccessibleSurfaceArea: [],
        maxLookupRadius: 1.4 + 1.4 + 1.87 + 1.87
    }
}

// class Count {
//     static count = 10;
//     get count(): number {
//         return Count.count;
//     }
//     set count(v: number) {
//         Count.count = v;
//     }
// }

function computeAccessibleSurfaceArea(unit: Unit.Atomic, params?: Partial<AccessibleSurfaceAreaComputationParameters>): AccessibleSurfaceArea {
    // const count = new Count();
    // count.count = count.count - 1;
    // if (count.count > 0) {
    console.log(`computing accessible surface area for unit #${ unit.id + 1 }`);
    return _computeAccessibleSurfaceArea(unit, {
        numberOfSpherePoints: (params && params.numberOfSpherePoints) || 92 /* original value Shrake, A. & Rupley, J. A. (1973) J. Mol. Biol. 79: 92, BioJava: 960, EPPIC: 3000 */ ,
        /** 92: 1600ms, 960: 5000ms, 3000: 13000ms */
        probeSize: (params && params.probeSize) || 1.4
    });
    // } else {
    //     return {
    //         atomRadius: [],
    //         accessibleSurfaceArea: [],
    //         relativeAccessibleSurfaceArea: [],
    //         buried: void 0
    //     }
    // }
}

/** Creates a collection of points on a sphere by the Golden Section Spiral algorithm. */
function generateSpherePoints(numberOfSpherePoints: number): Vec3[] {
    const points: Vec3[] = [];
    const inc = Math.PI * (3.0 - Math.sqrt(5.0));
    const offset = 2.0 / numberOfSpherePoints;
    for (let k = 0; k < numberOfSpherePoints; ++k) {
        const y = k * offset - 1.0 + (offset / 2.0);
        const r = Math.sqrt(1.0 - y * y);
        const phi = k * inc;
        points[points.length] = [Math.cos(phi) * r, y, Math.sin(phi) * r] as Vec3;
    }
    return points;
}

export { computeAccessibleSurfaceArea, missingAccessibleSurfaceAreaValue }