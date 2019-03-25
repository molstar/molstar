/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import Unit from '../../unit'
import { Vec3 } from 'mol-math/linear-algebra';
import { AccessibleSurfaceAreaComputationParams, AccessibleSurfaceArea, SolventAccessibility } from './data';
import { isHydrogen, getElementIdx } from '../links/common'; // TODO these functions are relevant for many tasks: move them somewhere actually common
import { ElementSymbol, MaxAsa, DefaultMaxAsa, isPolymer, isNucleic, MoleculeType } from 'mol-model/structure/model/types';
import { VdwRadius } from 'mol-model/structure/model/properties/atomic/measures';
import { ParamDefinition as PD } from 'mol-util/param-definition'

// Chothia's amino acid atoms vdw radii
const trigonalCarbonVdw = 1.76;
const tetrahedralCarbonVdw = 1.87;
const trigonalNitrogenVdw = 1.65;
const tetrahedralNitrogenVdw = 1.50;
/** deviating radii from definition in types.ts */
const oxygenVdw = 1.40;
const sulfurVdw = 1.85;
// Chothia's nucleotide atom vdw radii
const nucCarbonVdw = 1.80;
const nucNitrogenVdw = 1.60;
const nucPhosphorusVdw = 1.40;
const missingAccessibleSurfaceAreaValue = -1.0;

interface AccessibleSurfaceAreaContext {
    unit: Unit.Atomic,
    params: PD.Values<AccessibleSurfaceAreaComputationParams>,
    spherePoints: Vec3[],
    cons: number,
    maxLookupRadius: number,
    atomRadius?: Float32Array,
    accessibleSurfaceArea?: Float32Array,
    relativeAccessibleSurfaceArea?: Float32Array,
    buried?: Uint8Array
}

/**
 * Adapts the BioJava implementation by Jose Duarte. That implementation is based on the publication by Shrake, A., and
 * J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin." JMB (1973).
 */
function computeAccessibleSurfaceArea(unit: Unit.Atomic, params?: PD.Values<AccessibleSurfaceAreaComputationParams>): AccessibleSurfaceArea {
    if (!params) params = PD.getDefaultValues(AccessibleSurfaceAreaComputationParams);

    // TODO non-polymer flag is currently useless as hetatms are located in different units - aim is not to color them, but to compute values correctly - relates to correct ASA computation for inter-chain contacts
    console.log(`computing accessible surface area for unit #${ unit.id + 1 } - ${ params.numberOfSpherePoints } points, ${ params.probeSize } probe size, ${ params.nonPolymer ? 'honoring' : 'ignoring'} non-polymer atoms`);

    const ctx = initialize(unit, params);
    assignRadiusForHeavyAtoms(ctx);
    computePerResidue(ctx);
    normalizeAccessibleSurfaceArea(ctx);

    return {
        atomRadius: ctx.atomRadius!,
        accessibleSurfaceArea: ctx.accessibleSurfaceArea!,
        relativeAccessibleSurfaceArea: ctx.relativeAccessibleSurfaceArea!,
        buried: ctx.buried!
    };
}

function normalizeAccessibleSurfaceArea(ctx: AccessibleSurfaceAreaContext) {
    const { residues, derived } = ctx.unit.model.atomicHierarchy;
    const { accessibleSurfaceArea, relativeAccessibleSurfaceArea } = ctx;

    for (let i = 0; i < residues.label_comp_id.rowCount; ++i) {
        // skip entities not part of a polymer chain
        if (!ctx.params.nonPolymer) {
            if (!isPolymer(derived.residue.moleculeType[i])) continue;
        }

        const maxAsa = (MaxAsa as any)[residues.label_comp_id.value(i)];
        const rasa = accessibleSurfaceArea![i] / (maxAsa === undefined ? DefaultMaxAsa : maxAsa);
        relativeAccessibleSurfaceArea![i] = rasa;
        ctx.buried![i] |= (rasa < ctx.params.buriedRasaThreshold ? SolventAccessibility.Flag.BURIED : SolventAccessibility.Flag.ACCESSIBLE)
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
        const radii1 = atomRadius![aI];
        if (radii1 === missingAccessibleSurfaceAreaValue) continue;

        // find suitable neighbor candidates by lookup
        const { indices, count } = lookup3d.find(x[aI], y[aI], z[aI], maxLookupRadius);
        position(aI, a1Pos);

        // refine list by actual criterion
        const cutoff = probeSize + probeSize + radii1;
        const filteredIndicies = []; // TODO might be better to use IntArray here and reuse it - how to find safe upper limit of possible neighborhood count - BioJava mentions 60 as relatively safe upper bound
        for (let ni = 0; ni < count; ni++) {
            const _bI = indices[ni];
            const bI = atoms[_bI];
            const radii2 = atomRadius![bI];
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
                const cutoff3 = (atomRadius![naI] + probeSize) * (atomRadius![naI] + probeSize);
                if (Vec3.squaredDistance(testPoint, a2Pos) < cutoff3) {
                    accessible = false;
                    break;
                }
            }

            if (accessible) ++accessiblePointCount;
        }

        const value = cons * accessiblePointCount * r * r;
        accessibleSurfaceArea![residueIndex[aI]] += value;
        // +30% computation by normalizing partial solutions
        // relativeAccessibleSurfaceArea[residueIndex[aI]] += value * (NormalizationFactors as any)[residueIndex[aI]];
    }
}

function assignRadiusForHeavyAtoms(ctx: AccessibleSurfaceAreaContext) {
    const atomCount = ctx.unit.elements.length;
    const { elements: atoms, residueIndex } = ctx.unit;
    const { residues } = ctx.unit.model.atomicHierarchy;
    const { moleculeType } = ctx.unit.model.atomicHierarchy.derived.residue;
    const { type_symbol, label_atom_id } = ctx.unit.model.atomicHierarchy.atoms;
    const { label_comp_id } = ctx.unit.model.atomicHierarchy.residues;

    const residueCount = residues.label_comp_id.rowCount;
    ctx.atomRadius = new Float32Array(atomCount - 1);
    ctx.accessibleSurfaceArea = new Float32Array(residueCount - 1);
    ctx.relativeAccessibleSurfaceArea = new Float32Array(residueCount - 1);
    ctx.buried = new Uint8Array(residueCount - 1);

    for (let _aI = 0; _aI < atomCount; ++_aI) {
        const aI =  atoms[_aI];
        const raI = residueIndex[aI];
        const aeI = getElementIdx(type_symbol.value(aI)!);

        // skip hydrogen atoms
        if (isHydrogen(aeI)) {
            ctx.atomRadius[aI] = missingAccessibleSurfaceAreaValue;
            continue;
        }

        // skip non-polymer groups
        if (!ctx.params.nonPolymer) {
            if (!isPolymer(moleculeType[raI])) {
                ctx.atomRadius[aI] = missingAccessibleSurfaceAreaValue;
                continue;
            }
        }

        const atomId = label_atom_id.value(aI);
        const element = type_symbol.value(aI);
        const resn = label_comp_id.value(raI)!;

        ctx.atomRadius[aI] = isNucleic(moleculeType[raI]) ? determineRadiusNucl(atomId, element, resn) : moleculeType[raI] === MoleculeType.protein ? determineRadiusAmino(atomId, element, resn) : VdwRadius(element);
    }
}

/**
 * Gets the van der Waals radius of the given atom following the values defined by Chothia (1976)
 * J.Mol.Biol.105,1-14. NOTE: the vdw values defined by the paper assume no Hydrogens and thus "inflates" slightly
 * the heavy atoms to account for Hydrogens.
 */
function determineRadiusAmino(atomId: string, element: ElementSymbol, compId: string): number {
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

function determineRadiusNucl(atomId: string, element: ElementSymbol, compId: string): number {
    switch (element) {
        case 'C':
        return nucCarbonVdw;
        case 'N':
        return nucNitrogenVdw;
        case 'P':
        return nucPhosphorusVdw;
        case 'O':
        return oxygenVdw;
    }
    return VdwRadius(element);
}

function initialize(unit: Unit.Atomic, params: PD.Values<AccessibleSurfaceAreaComputationParams>): AccessibleSurfaceAreaContext {
    return {
        unit: unit,
        params: params,
        spherePoints: generateSpherePoints(params.numberOfSpherePoints),
        cons: 4.0 * Math.PI / params.numberOfSpherePoints,
        maxLookupRadius: 1.4 + 1.4 + 1.87 + 1.87
    }
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