/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SecondaryStructure } from 'mol-model/structure/model/properties/seconday-structure';
import { ResidueIndex } from 'mol-model/structure';
import { MoleculeType, SecondaryStructureType } from 'mol-model/structure/model/types';
import { Vec3 } from 'mol-math/linear-algebra';
import { GridLookup3D } from 'mol-math/geometry';
import { SortedArray } from 'mol-data/int';
import { IntAdjacencyGraph } from 'mol-math/graph';
import { BitFlags } from 'mol-util';
import { ElementIndex } from 'mol-model/structure/model/indexing';
import { AtomicHierarchy, AtomicConformation } from '../atomic';

export function computeSecondaryStructure(hierarchy: AtomicHierarchy, conformation: AtomicConformation): SecondaryStructure {
    // TODO use Zhang-Skolnik for CA alpha only parts or for coarse parts with per-residue elements
    return computeModelDSSP(hierarchy, conformation)
}

export function computeModelDSSP(hierarchy: AtomicHierarchy, conformation: AtomicConformation) {
    const { lookup3d, proteinResidues } = calcAtomicTraceLookup3D(hierarchy, conformation)
    const backboneIndices = calcBackboneAtomIndices(hierarchy, proteinResidues)
    const hbonds = calcBackboneHbonds(hierarchy, conformation, proteinResidues, backboneIndices, lookup3d)

    const residueCount = proteinResidues.length
    const flags = new Uint32Array(residueCount)

    const ctx: DSSPContext = {
        hierarchy,
        proteinResidues,
        flags,
        hbonds
    }

    assignBends(ctx)
    assignTurns(ctx)
    assignHelices(ctx)
    assignBridges(ctx)
    assignLadders(ctx)
    assignSheets(ctx)

    const assignment = getDSSPAssignment(flags)

    const type = new Uint32Array(hierarchy.residues._rowCount) as unknown as SecondaryStructureType[]
    for (let i = 0, il = proteinResidues.length; i < il; ++i) {
        type[proteinResidues[i]] = assignment[i]
    }
    
    const secondaryStructure: SecondaryStructure = {
        type,
        key: [], // TODO
        elements: [] // TODO
    }
    return secondaryStructure
}

interface DSSPContext {
    hierarchy: AtomicHierarchy
    proteinResidues: SortedArray<ResidueIndex>
    /** flags for each residue */
    flags: Uint32Array

    hbonds: DsspHbonds
}

interface DSSPType extends BitFlags<DSSPType.Flag> { }
namespace DSSPType {
    export const is: (t: DSSPType, f: Flag) => boolean = BitFlags.has
    export const create: (f: Flag) => DSSPType = BitFlags.create
    export const enum Flag {
        _ = 0x0,
        H = 0x1,
        B = 0x2,
        E = 0x4,
        G = 0x8,
        I = 0x10,
        S = 0x20,
        T = 0x40,
        T3 = 0x80,
        T4 = 0x100,
        T5 = 0x200,
    }
}

/** max distance between two C-alpha atoms to check for hbond */
const caMaxDist = 7.0;

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
    const lookup3d = GridLookup3D({ x, y, z, indices: SortedArray.ofSortedArray(indices) });
    const proteinResidues = SortedArray.ofSortedArray<ResidueIndex>(_proteinResidues)
    return { lookup3d, proteinResidues }
}

interface BackboneAtomIndices {
    cIndices: ArrayLike<ElementIndex | -1>
    hIndices: ArrayLike<ElementIndex | -1>
    oIndices: ArrayLike<ElementIndex | -1>
    nIndices: ArrayLike<ElementIndex | -1>
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

type DsspHbonds = IntAdjacencyGraph<{ readonly energies: ArrayLike<number> }>

function calcBackboneHbonds(hierarchy: AtomicHierarchy, conformation: AtomicConformation, proteinResidues: SortedArray<ResidueIndex>, backboneIndices: BackboneAtomIndices, lookup3d: GridLookup3D): DsspHbonds {
    const { cIndices, hIndices, nIndices, oIndices } = backboneIndices
    const { index } = hierarchy
    const { x, y, z } = conformation
    const { traceElementIndex } = hierarchy.derived.residue

    const residueCount = proteinResidues.length
    const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i])

    const oAtomResidues: number[] = [];
    const nAtomResidues: number[] = [];
    const energies: number[] = [];

    const oPos = Vec3.zero()
    const cPos = Vec3.zero()
    const caPos = Vec3.zero()
    const nPos = Vec3.zero()
    const hPos = Vec3.zero()

    const cPosPrev = Vec3.zero()
    const oPosPrev = Vec3.zero()

    for (let i = 0, il = proteinResidues.length; i < il; ++i) {
        const oPI = i
        const oRI = proteinResidues[i]

        const oAtom = oIndices[oPI]
        const cAtom = cIndices[oPI]
        const caAtom = traceElementIndex[oRI]

        // continue if residue is missing O or C atom
        if (oAtom === -1 || cAtom === -1) continue

        // ignore C-terminal residue as acceptor
        if (index.findAtomOnResidue(oRI, 'OXT') !== -1) continue

        position(oAtom, oPos)
        position(cAtom, cPos)
        position(caAtom, caPos)

        const { indices, count } = lookup3d.find(caPos[0], caPos[1], caPos[2], caMaxDist)

        for (let j = 0; j < count; ++j) {
            const nPI = indices[j]

            // ignore bonds within a residue or to prev or next residue, TODO take chain border into account
            if (nPI === oPI || nPI - 1 === oPI || nPI + 1 === oPI) continue

            const nAtom = nIndices[nPI]
            if (nAtom === -1) continue

            position(nAtom, nPos)

            const hAtom = hIndices[nPI]
            if (hAtom === -1) {
                // approximate calculation of H position, TODO factor out
                if (nPI === 0) continue
                const nPIprev = nPI - 1

                const oAtomPrev = oIndices[nPIprev]
                const cAtomPrev = cIndices[nPIprev]
                if (oAtomPrev === -1 || cAtomPrev === -1) continue

                position(oAtomPrev, oPosPrev)
                position(cAtomPrev, cPosPrev)

                Vec3.sub(hPos, cPosPrev, oPosPrev)
                const dist = Vec3.distance(oPosPrev, cPosPrev)
                Vec3.scaleAndAdd(hPos, nPos, hPos, 1 / dist)
            } else {
                position(hAtom, hPos)
            }

            const e = calcHbondEnergy(oPos, cPos, nPos, hPos)
            if (e > hbondEnergyCutoff) continue

            oAtomResidues[oAtomResidues.length] = oPI
            nAtomResidues[nAtomResidues.length] = nPI
            energies[energies.length] = e
        }
    }

    return buildHbondGraph(residueCount, oAtomResidues, nAtomResidues, energies);
}

function buildHbondGraph(residueCount: number, oAtomResidues: number[], nAtomResidues: number[], energies: number[]) {
    const builder = new IntAdjacencyGraph.DirectedEdgeBuilder(residueCount, oAtomResidues, nAtomResidues);
    const _energies = new Float32Array(builder.slotCount);

    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge();
        builder.assignProperty(_energies, energies[i]);
    }

    return builder.createGraph({ energies });
}

/** Original priority: H,B,E,G,I,T,S */
function getOriginalResidueFlag(f: DSSPType) {
    if (DSSPType.is(f, DSSPType.Flag.H)) return SecondaryStructureType.SecondaryStructureDssp.H
    if (DSSPType.is(f, DSSPType.Flag.B)) return SecondaryStructureType.SecondaryStructureDssp.B
    if (DSSPType.is(f, DSSPType.Flag.E)) return SecondaryStructureType.SecondaryStructureDssp.E
    if (DSSPType.is(f, DSSPType.Flag.G)) return SecondaryStructureType.SecondaryStructureDssp.G
    if (DSSPType.is(f, DSSPType.Flag.I)) return SecondaryStructureType.SecondaryStructureDssp.I
    if (DSSPType.is(f, DSSPType.Flag.T)) return SecondaryStructureType.SecondaryStructureDssp.T
    if (DSSPType.is(f, DSSPType.Flag.S)) return SecondaryStructureType.SecondaryStructureDssp.S
    return SecondaryStructureType.Flag.None
}

/** Version 2.1.0 priority: I,H,B,E,G,T,S */
function getUpdatedResidueFlag(f: DSSPType) {
    if (DSSPType.is(f, DSSPType.Flag.I)) return SecondaryStructureType.SecondaryStructureDssp.I
    if (DSSPType.is(f, DSSPType.Flag.H)) return SecondaryStructureType.SecondaryStructureDssp.H
    if (DSSPType.is(f, DSSPType.Flag.B)) return SecondaryStructureType.SecondaryStructureDssp.B
    if (DSSPType.is(f, DSSPType.Flag.E)) return SecondaryStructureType.SecondaryStructureDssp.E
    if (DSSPType.is(f, DSSPType.Flag.G)) return SecondaryStructureType.SecondaryStructureDssp.G
    if (DSSPType.is(f, DSSPType.Flag.T)) return SecondaryStructureType.SecondaryStructureDssp.T
    if (DSSPType.is(f, DSSPType.Flag.S)) return SecondaryStructureType.SecondaryStructureDssp.S
    return SecondaryStructureType.Flag.None
}

// function geFlagName(f: DSSPType) {
//     if (DSSPType.is(f, DSSPType.Flag.I)) return 'I'
//     if (DSSPType.is(f, DSSPType.Flag.H)) return 'H'
//     if (DSSPType.is(f, DSSPType.Flag.B)) return 'B'
//     if (DSSPType.is(f, DSSPType.Flag.E)) return 'E'
//     if (DSSPType.is(f, DSSPType.Flag.G)) return 'G'
//     if (DSSPType.is(f, DSSPType.Flag.T)) return 'T'
//     if (DSSPType.is(f, DSSPType.Flag.S)) return 'S'
//     return '-'
// }

function getDSSPAssignment(flags: Uint32Array, useOriginal = false) {
    const getResidueFlag = useOriginal ? getOriginalResidueFlag : getUpdatedResidueFlag
    const type = new Uint32Array(flags.length)
    for (let i = 0, il = flags.length; i < il; ++i) {
        const f = DSSPType.create(flags[i])
        // console.log(i, geFlagName(f))
        type[i] = getResidueFlag(f)
    }
    return type as unknown as ArrayLike<SecondaryStructureType>
}

/**
 * Constant for electrostatic energy in kcal/mol
 *      f  *  q1 *   q2
 * Q = -332 * 0.42 * 0.20
 *
 * f is the dimensional factor
 * 
 * q1 and q2 are partial charges which are placed on the C,O
 * (+q1,-q1) and N,H (-q2,+q2)
 */
const Q = -27.888

/** cutoff for hbonds in kcal/mol, must be lower to be consider as an hbond */
const hbondEnergyCutoff = -0.5

/**
 * E = Q * (1/r(ON) + l/r(CH) - l/r(OH) - l/r(CN))
 */
function calcHbondEnergy(oPos: Vec3, cPos: Vec3, nPos: Vec3, hPos: Vec3) {
    const distOH = Vec3.distance(oPos, hPos)
    const distCH = Vec3.distance(cPos, hPos)
    const distCN = Vec3.distance(cPos, nPos)
    const distON = Vec3.distance(oPos, nPos)

    const e1 = Q / distOH - Q / distCH
	const e2 = Q / distCN - Q / distON
    return e1 + e2
}

/**
 * The basic turn pattern is a single H bond of type (i, i + n).
 * We assign an n-turn at residue i if there is an H bond from CO(i) to NH(i + n),
 * i.e., “n-turn(i)=: Hbond(i, i + n), n = 3, 4, 5.”
 * 
 * Type: T
 */
function assignTurns(ctx: DSSPContext) {
    const { proteinResidues, hbonds, flags, hierarchy } = ctx
    const { chains, residueAtomSegments, chainAtomSegments } = hierarchy
    const { label_asym_id } = chains

    const turnFlag = [ 0, 0, 0, DSSPType.Flag.T3, DSSPType.Flag.T4, DSSPType.Flag.T5 ]

    for (let i = 0, il = proteinResidues.length; i < il; ++i) {
        const rI = proteinResidues[i]
        const cI = chainAtomSegments.index[residueAtomSegments.offsets[rI]]

        // TODO should take sequence gaps into account
        for (let k = 3; k <= 5; ++k) {
            if (i + k >= proteinResidues.length) continue

            const rN = proteinResidues[i + k]
            const cN = chainAtomSegments.index[residueAtomSegments.offsets[rN]]
            // check if on same chain
            if (!label_asym_id.areValuesEqual(cI, cN)) continue

            // check if hbond exists
            if (hbonds.getDirectedEdgeIndex(i, i + k) !== -1) {
                flags[i] |= turnFlag[k] | DSSPType.Flag.T
            }
        }
    }
}

/**
 * Two nonoverlapping stretches of three residues each, i - 1, i, i + 1 and j - 1, j, j + 1,
 * form either a parallel or antiparallel bridge, depending on which of 
 * two basic patterns is matched. We assign a bridge between residues i and j
 * if there are two H bonds characteristic of P-structure; in particular,
 * 
 * Parallel Bridge(i, j) =:
 *      [Hbond(i - 1, j) and Hbond(j, i + 1)] or
 *      [Hbond(j - 1, i) and Hbond(i, j + 1)]
 * 
 * Antiparallel Bridge(i, j) =:
 *      [Hbond(i, j) and Hbond(j, i)] or
 *      [Hbond(i - 1, j + 1) and Hbond(j - 1, i + l)]
 * 
 * Type: B
 */
function assignBridges(ctx: DSSPContext) {
    const { proteinResidues, hbonds, flags } = ctx

    const { offset, b } = hbonds
    let i: number, j: number

    for (let k = 0, kl = proteinResidues.length; k < kl; ++k) {
        for (let t = offset[k], _t = offset[k + 1]; t < _t; t++) {
            const l = b[t]
            if (k > l) continue
            
            // Parallel Bridge(i, j) =: [Hbond(i - 1, j) and Hbond(j, i + 1)]
            i = k + 1 // k is i - 1
            j = l
            if (i !== j && hbonds.getDirectedEdgeIndex(j, i + 1) !== -1) {
                flags[i] |= DSSPType.Flag.B
                flags[j] |= DSSPType.Flag.B
            }

            // Parallel Bridge(i, j) =: [Hbond(j - 1, i) and Hbond(i, j + 1)]
            i = k
            j = l - 1 // l is j + 1
            if (i !== j && hbonds.getDirectedEdgeIndex(j - 1, i) !== -1) {
                flags[i] |= DSSPType.Flag.B
                flags[j] |= DSSPType.Flag.B
            }

            // Antiparallel Bridge(i, j) =: [Hbond(i, j) and Hbond(j, i)]
            i = k
            j = l
            if (i !== j && hbonds.getDirectedEdgeIndex(j, i) !== -1) {
                flags[i] |= DSSPType.Flag.B
                flags[j] |= DSSPType.Flag.B
            }

            // Antiparallel Bridge(i, j) =: [Hbond(i - 1, j + 1) and Hbond(j - 1, i + l)]
            i = k + 1
            j = l - 1
            if (i !== j && hbonds.getDirectedEdgeIndex(j - 1, i + 1) !== -1) {
                flags[i] |= DSSPType.Flag.B
                flags[j] |= DSSPType.Flag.B
            }
        }
    }
}

/**
 * A minimal helix is defined by two consecutive n-turns.
 * For example, a 4-helix, of minimal length 4 from residues i to i + 3,
 * requires 4-turns at residues i - 1 and i,
 * 
 *      3-helix(i,i + 2)=: [3-turn(i - 1) and 3-turn(i)]
 *      4-helix(i,i + 3)=: [4-turn(i - 1) and 4-turn(i)]
 *      5-helix(i,i + 4)=: [5-turn(i - 1) and 5-turn(i)]
 * 
 * Type: G (n=3), H (n=4), I (n=5)
 */
function assignHelices(ctx: DSSPContext) {
    const { proteinResidues, flags } = ctx
    
    const turnFlag = [ 0, 0, 0, DSSPType.Flag.T3, DSSPType.Flag.T4, DSSPType.Flag.T5 ]
    const helixFlag = [ 0, 0, 0, DSSPType.Flag.G, DSSPType.Flag.H, DSSPType.Flag.I ]

    for (let i = 1, il = proteinResidues.length; i < il; ++i) {
        const fI = DSSPType.create(flags[i])
        const fI1 = DSSPType.create(flags[i - 1])

        for (let k = 3; k <= 5; ++k) {
            if (DSSPType.is(fI, turnFlag[k]) && DSSPType.is(fI1, turnFlag[k])) {
                for (let l = 0; l < k; ++l) {
                    flags[i + l] |= helixFlag[k]
                }
            }
        }
    }
}

/**
 * ladder=: set of one or more consecutive bridges of identical type
 * 
 * Type: E
 */
function assignLadders(ctx: DSSPContext) {
    // TODO
}

/**
 * sheet=: set of one or more ladders connected by shared residues
 * 
 * Type: E
 */
function assignSheets(ctx: DSSPContext) {
    // TODO
}

/**
 * Bend(i) =: [angle ((CW - Ca(i - 2)),(C"(i + 2) - C"(i))) > 70"]
 * 
 * Type: S
 */
function assignBends(ctx: DSSPContext) {
    // TODO
}