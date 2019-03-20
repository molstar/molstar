import { Vec3 } from 'mol-math/linear-algebra';
import { AtomicHierarchy, AtomicConformation } from '../atomic';
import { ElementSymbol, VdwRadii } from '../../types';
import { skipUntil } from 'rxjs/operators';
import { atomicHet } from 'mol-model/structure/query/queries/internal';
import { compile } from 'mol-script/runtime/query/compiler';

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
    console.log(ctx.hierarchy)
    // 1. extract all heavy atoms
    // 2. assign radius to all heavy atoms - depends on element
    // 3. for each residue
    // 3a. calculate the individual ASA of each atom
    // 3b. sum up
    // 3c. (optionally) normalize by maximum value expected for a particular amino acid - needs lookup of max values

    // const { type_symbol } = ctx.hierarchy.atoms;
    // console.log(type_symbol.value(0));

    const { atoms, residues, residueAtomSegments, derived } = ctx.hierarchy;
    const { type_symbol } = atoms;
    const radii: number[] = [];

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
        radii[radii.length] = determineRadius(atomId, element, compId);
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
        if (atomId === 'C' || atomId === 'CE1' || atomId === 'CE2' || atomId === 'CE3' ||
            atomId === 'CH2' || atomId === 'CZ' || atomId === 'CZ2' || atomId === 'CZ3') {
            return trigonalCarbonVdw;
        } else if (atomId === 'CA' || atomId === 'CB' || atomId === 'CE' || atomId === 'CG1' ||
            atomId === 'CG2') {
            return tetrahedralCarbonVdw;
        } else if (compId === 'PHE' || compId === 'TRP' || compId === 'TYR' || compId === 'HIS' ||
            compId === 'ASP' || compId === 'ASN') {
            return trigonalCarbonVdw;
        } else if (compId === 'PRO' || compId === 'LYS' || compId === 'ARG' || compId === 'MET' ||
            compId === 'ILE' || compId === 'LEU') {
            return tetrahedralCarbonVdw;
        } else if (compId === 'GLU' || compId === 'GLN') {
            return atomId === 'CD' ? trigonalCarbonVdw : tetrahedralCarbonVdw;
        }
        default:
        return (VdwRadii as any)[atomId];
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