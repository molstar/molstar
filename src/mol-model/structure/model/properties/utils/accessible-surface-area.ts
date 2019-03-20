import { ResidueIndex } from 'mol-model/structure';
import { MoleculeType } from 'mol-model/structure/model/types';
import { GridLookup3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import { SortedArray } from 'mol-data/int';
import { ElementIndex } from 'mol-model/structure/model/indexing';
import { AtomicHierarchy, AtomicConformation } from '../atomic';

export function computeModelASA(hierarchy: AtomicHierarchy, conformation: AtomicConformation) {
    const { lookup3d, proteinResidues } = calcAtomicTraceLookup3D(hierarchy, conformation)
    const backboneIndices = calcBackboneAtomIndices(hierarchy, proteinResidues)

    const residueCount = proteinResidues.length
    const numberOfSpherePoints = 960

    const ctx: ASAContext = {
        probeSize: 1.4,
        spherePoints: generateSpherePoints(numberOfSpherePoints),
        cons: 4.0 * Math.PI / numberOfSpherePoints,

        hierarchy,
        proteinResidues
   }

    calculateASA(ctx)
}

function calculateASA(ctx: ASAContext) {

}

function generateSpherePoints(numberOfSpherePoints: number): Vec3[] {
    const points: Vec3[] = []
    const inc = Math.PI

    return points
}

interface ASAContext {
    probeSize: number,
    spherePoints: Vec3[],
    cons: number,

    hierarchy: AtomicHierarchy,
    proteinResidues: SortedArray<ResidueIndex>
}

interface BackboneAtomIndices {
    cIndices: ArrayLike<ElementIndex | -1>
    hIndices: ArrayLike<ElementIndex | -1>
    oIndices: ArrayLike<ElementIndex | -1>
    nIndices: ArrayLike<ElementIndex | -1>
}

function calcAtomicTraceLookup3D(hierarchy: AtomicHierarchy, conformation: AtomicConformation) {
    const { x, y, z } = conformation;
    const { moleculeType, traceElementIndex } = hierarchy.derived.residue
    const indices: number[] = []
    const _proteinResidues: number[] = []
    for (let i = 0, il = moleculeType.length; i < il; ++i) {
        if (moleculeType[i] === MoleculeType.protein) {
            indices[indices.length] = traceElementIndex[i]
            _proteinResidues[_proteinResidues.length] = i
        }
    }
    const lookup3d = GridLookup3D({ x, y, z, indices: SortedArray.ofSortedArray(indices) }, 4);
    const proteinResidues = SortedArray.ofSortedArray<ResidueIndex>(_proteinResidues)
    return { lookup3d, proteinResidues }
}


function calcBackboneAtomIndices(hierarchy: AtomicHierarchy, proteinResidues: SortedArray<ResidueIndex>): BackboneAtomIndices {
    const residueCount = proteinResidues.length
    const { index } = hierarchy

    const c = new Int32Array(residueCount)
    const h = new Int32Array(residueCount)
    const o = new Int32Array(residueCount)
    const n = new Int32Array(residueCount)

    for (let i = 0, il = residueCount; i < il; ++i) {
        const rI = proteinResidues[i]
        c[i] = index.findAtomOnResidue(rI, 'C')
        h[i] = index.findAtomOnResidue(rI, 'H')
        o[i] = index.findAtomOnResidue(rI, 'O')
        n[i] = index.findAtomOnResidue(rI, 'N')
    }

    return {
        cIndices: c as unknown as ArrayLike<ElementIndex | -1>,
        hIndices: h as unknown as ArrayLike<ElementIndex | -1>,
        oIndices: o as unknown as ArrayLike<ElementIndex | -1>,
        nIndices: n as unknown as ArrayLike<ElementIndex | -1>,
    }
}