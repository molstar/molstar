import { Vec3 } from 'mol-math/linear-algebra';
import { AtomicHierarchy, AtomicConformation } from '../atomic';
import { ElementSymbol, VdwRadii, MaxAsa, DefaultMaxAsa } from '../../types';
import { max } from 'mol-data/int/impl/ordered-set';

const defaultNumberOfPoints = 960;
const defaultProbeSize = 1.4;
const trigonalCarbonVdw = 1.76;
const tetrahedralCarbonVdw = 1.87;
const trigonalNitrogenVdw = 1.65;
const tetrahedralNitrogenVdw = 1.50;
/** deviating radii from the definition in types.ts */
const oxygenVdw = 1.4;
const sulfurVdw = 1.85;

export function computeModelASA(hierarchy: AtomicHierarchy, conformation: AtomicConformation) {
    const numberOfSpherePoints = defaultNumberOfPoints;

    const ctx: ASAContext = {
        probeSize: defaultProbeSize,
        spherePoints: generateSpherePoints(numberOfSpherePoints),
        cons: 4.0 * Math.PI / numberOfSpherePoints,

        hierarchy,
        conformation
    }

    calculateASA(ctx);
}

const valueForIgnoredAtom = -1.0;

function calculateASA(ctx: ASAContext) {
    const { probeSize, spherePoints, cons } = ctx;
    const { atoms, residues, residueAtomSegments, derived } = ctx.hierarchy;
    const { type_symbol } = atoms;
    const { x, y, z } = ctx.conformation;

    const radii: number[] = [];
    // const atomAsa: number[] = [];
    const residueAsa: number[] = [];

    console.log(ctx.hierarchy);
    console.log(ctx.conformation);

    const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i]);
    const a1Pos = Vec3.zero();
    const a2Pos = Vec3.zero();

    // extract all heavy atoms
    // TODO can be more elegant by obtaining residue info once and then using offset to navigate to next residue
    for (let i = 0; i < type_symbol.rowCount; ++i) {
        // skip hydrogen atoms
        if (type_symbol.value(i) === 'H' || type_symbol.value(i) === 'D' || type_symbol.value(i) === 'T') {
            radii[radii.length] = valueForIgnoredAtom; // -1 is used to signal atoms not to be considered downstream
            continue;
        }

        // determine group this atom belongs to
        const groupIndex = residueAtomSegments.index[i];

        // skip entities not part of a peptide chain
        if (derived.residue.moleculeType[groupIndex] !== 4) {
            radii[radii.length] = valueForIgnoredAtom;
            continue;
        }

        const atomId = atoms.label_atom_id.value(i);
        const compId = residues.label_comp_id.value(groupIndex);
        const element = type_symbol.value(i);
        // assign radius to all heavy atoms - depends on element and bonding patterns
        radii[radii.length] = determineRadius(atomId, element, compId);
        // set ASA of corresponding group to 0
        residueAsa[groupIndex] = 0.0;
    }

    // calculate the individual ASA of each atom
    // TODO again might be more elegant to use offsets
    // TODO distance is symmetric, omit redudant calcuations
    for (let i = 0; i < radii.length; ++i) {
        const radius = radii[i];

        // skip invalid entries
        if (radius === valueForIgnoredAtom) {
            // atomAsa[atomAsa.length] = valueForIgnoredAtom;
            continue;
        }

        position(i, a1Pos);

        // collect all neighboring atoms
        const cutoff = probeSize + probeSize + radius;
        const neighborIndices: number[] = [];
        for (let k = 0; k < radii.length; ++k) {
            const radius2 = radii[k];
            if (i === k || radius2 === valueForIgnoredAtom)
                continue;

            position(k, a2Pos);

            if (Vec3.distance(a1Pos, a2Pos) < cutoff + radius2) {
                neighborIndices[neighborIndices.length] = k;
            }
        }

        const r = probeSize + radius;
        let accessiblePoints = 0;

        for (let k = 0; k < spherePoints.length; ++k) {
            const point = spherePoints[k];
            let accessible = true;
            const testPoint = [point[0] * r + a1Pos[0], point[1] * r + a1Pos[1], point[2] * r + a1Pos[2]] as Vec3;

            for (let j = 0; j < neighborIndices.length; ++j) {
                const naI = neighborIndices[j];
                // store position of neighboring atom in a2Pos
                position(naI, a2Pos);
                const neighboringAtomRadius = radii[naI] + probeSize;
                if (Vec3.squaredDistance(testPoint, a2Pos) < neighboringAtomRadius * neighboringAtomRadius) {
                    accessible = false;
                    break;
                }
            }

            if (accessible) {
                ++accessiblePoints;
            }
        }

        const value = cons * accessiblePoints * r * r;
        // atomAsa[atomAsa.length] = value;
        // sum up values for each residue
        residueAsa[residueAtomSegments.index[i]] += value;
    }

    console.log(residueAsa);

    // normalize by maximum value expected for a particular amino acid - needs lookup of max values
    for (let i = 0; i < residues.label_comp_id.rowCount; ++i) {
        // skip entities not part of a peptide chain
        if (derived.residue.moleculeType[i] !== 4)
            continue;

        const maxAsa = (MaxAsa as any)[residues.label_comp_id.value(i)];
        residueAsa[i] /= (maxAsa === undefined ? DefaultMaxAsa : maxAsa);
    }

    console.log(residueAsa);
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
    // TODO could potentially use logging or error thrown
    return (VdwRadii as any)[atomId];
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
        points[points.length] = [Math.cos(phi), y, Math.sin(phi) * r] as Vec3;
    }
    return points;
}

interface ASAContext {
    probeSize: number,
    spherePoints: Vec3[],
    cons: number,

    hierarchy: AtomicHierarchy,
    conformation: AtomicConformation
}