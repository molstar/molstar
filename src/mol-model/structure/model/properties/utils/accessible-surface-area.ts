import { Vec3 } from 'mol-math/linear-algebra';
import { AtomicHierarchy, AtomicConformation } from '../atomic';
import { ElementSymbol, VdwRadii } from '../../types';

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

function calculateASA(ctx: ASAContext) {
    const { atoms, residues, residueAtomSegments, derived } = ctx.hierarchy;
    const { type_symbol } = atoms;
    const radii: number[] = [];
    const asa: number[] = [];

    // 1. extract all heavy atoms
    // TODO can be more elegant by obtaining residue info once and then using offset to navigate to next residue
    for (let i = 0; i < type_symbol.rowCount; ++i) {
        // skip hydrogen atoms
        if (type_symbol.value(i) === 'H' || type_symbol.value(i) === 'D' || type_symbol.value(i) === 'T')
            continue;

        // determine group this atom belongs to
        const groupIndex = residueAtomSegments.index[i];

        // skip entities not part of a peptide chain
        if (derived.residue.moleculeType[groupIndex] !== 4)
            continue;

        const atomId = atoms.label_atom_id.value(i);
        const compId = residues.label_comp_id.value(groupIndex);
        const element = type_symbol.value(i);
        // 2. assign radius to all heavy atoms - depends on element and bonding patterns
        radii[radii.length] = determineRadius(atomId, element, compId);
    }

    // 3. for each residue
    // 3a. calculate the individual ASA of each atom
    // 3b. sum up
    // 3c. (optionally) normalize by maximum value expected for a particular amino acid - needs lookup of max values

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