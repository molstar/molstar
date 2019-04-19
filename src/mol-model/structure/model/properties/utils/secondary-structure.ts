/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
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
import { ParamDefinition as PD } from 'mol-util/param-definition'

/**
 * TODO bugs to fix:
 * - some turns are not detected correctly: see e.g. pdb:1acj - maybe more than 2 hbonds require some residue to donate electrons
 * - some sheets are not extended correctly: see e.g. pdb:1acj
 * - validate new helix definition
 * - validate new ordering of secondary structure elements
 */

 /** max distance between two C-alpha atoms to check for hbond */
const caMaxDist = 9.0;

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
/** prevent extremely low hbond energies */
const hbondEnergyMinimal = -9.9

interface DSSPContext {
    params: Partial<PD.Values<SecondaryStructureComputationParams>>,
    getResidueFlag: (f: DSSPType) => SecondaryStructureType,
    getFlagName: (f: DSSPType) => String,

    hierarchy: AtomicHierarchy
    proteinResidues: SortedArray<ResidueIndex>
    /** flags for each residue */
    flags: Uint32Array
    hbonds: DsspHbonds,

    torsionAngles: { phi: Float32Array, psi: Float32Array },
    backboneIndices: BackboneAtomIndices,
    conformation: AtomicConformation,
    ladders: Ladder[],
    bridges: Bridge[]
}

interface Ladder {
    previousLadder: number,
    nextLadder: number,
    firstStart: number,
    secondStart: number,
    secondEnd: number,
    firstEnd: number,
    type: BridgeType
}

const enum BridgeType {
    PARALLEL = 0x0,
    ANTI_PARALLEL = 0x1
}

class Bridge {
    partner1: number;
    partner2: number;
    type: BridgeType;

    constructor(p1: number, p2: number, type: BridgeType) {
        this.partner1 = Math.min(p1, p2)
        this.partner2 = Math.max(p1, p2)
        this.type = type
    }
}

type DSSPType = BitFlags<DSSPType.Flag>
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
        T3S = 0x400, // marks 3-turn start
        T4S = 0x800,
        T5S = 0x1000
    }
}

export const SecondaryStructureComputationParams = {
    oldDefinition: PD.Boolean(true, { description: 'Whether to use the old DSSP convention for the annotation of turns and helices, causes them to be two residues shorter' }),
    oldOrdering: PD.Boolean(true, { description: 'Alpha-helices are preferred over 3-10 helices' })
}
export type SecondaryStructureComputationParams = typeof SecondaryStructureComputationParams

export function computeSecondaryStructure(hierarchy: AtomicHierarchy,
    conformation: AtomicConformation) {
    // TODO use Zhang-Skolnik for CA alpha only parts or for coarse parts with per-residue elements
    return computeModelDSSP(hierarchy, conformation)
}

export function computeModelDSSP(hierarchy: AtomicHierarchy,
    conformation: AtomicConformation,
    params: Partial<PD.Values<SecondaryStructureComputationParams>> = {}): SecondaryStructure {
    params = { ...PD.getDefaultValues(SecondaryStructureComputationParams), ...params };

    const { lookup3d, proteinResidues } = calcAtomicTraceLookup3D(hierarchy, conformation)
    const backboneIndices = calcBackboneAtomIndices(hierarchy, proteinResidues)
    const hbonds = calcBackboneHbonds(hierarchy, conformation, proteinResidues, backboneIndices, lookup3d)

    const residueCount = proteinResidues.length
    const flags = new Uint32Array(residueCount)

    // console.log(`calculating secondary structure elements using ${ params.oldDefinition ? 'old' : 'revised'} definition and ${ params.oldOrdering ? 'old' : 'revised'} ordering of secondary structure elements`)

    const torsionAngles = calculateDihedralAngles(hierarchy, conformation, proteinResidues, backboneIndices)

    const ladders: Ladder[] = []
    const bridges: Bridge[] = []

    const getResidueFlag = params.oldDefinition ? getOriginalResidueFlag : getUpdatedResidueFlag
    const getFlagName = params.oldOrdering ? getOriginalFlagName : getUpdatedFlagName

    const ctx: DSSPContext = {
        params,
        getResidueFlag,
        getFlagName,

        hierarchy,
        proteinResidues,
        flags,
        hbonds,

        torsionAngles,
        backboneIndices,
        conformation,
        ladders,
        bridges
    }

    assignTurns(ctx)
    assignHelices(ctx)
    assignBends(ctx)
    assignBridges(ctx)
    assignLadders(ctx)
    assignSheets(ctx)

    const assignment = getDSSPAssignment(flags, getResidueFlag)
    const type = new Uint32Array(hierarchy.residues._rowCount) as unknown as SecondaryStructureType[]
    const keys: number[] = []
    const elements: SecondaryStructure.Element[] = []

    for (let i = 0, il = proteinResidues.length; i < il; ++i) {
        const assign = assignment[i]
        type[proteinResidues[i]] = assign
        const flag = getResidueFlag(flags[i])
        // TODO is this expected behavior? elements will be strictly split depending on 'winning' flag
        if (elements.length === 0 /* would fail at very start */ || flag !== (elements[elements.length - 1] as SecondaryStructure.Helix | SecondaryStructure.Sheet | SecondaryStructure.Turn).flags /* flag changed */) {
            elements[elements.length] = createElement(mapToKind(assign), flags[i], getResidueFlag)
        }
        keys[i] = elements.length - 1
    }

    const secondaryStructure: SecondaryStructure = {
        type,
        key: keys,
        elements: elements
    }

    return secondaryStructure
}

function createElement(kind: string, flag: DSSPType.Flag, getResidueFlag: (f: DSSPType) => SecondaryStructureType): SecondaryStructure.Element {
    // TODO would be nice to add more detailed information
    if (kind === 'helix') {
        return {
            kind: 'helix',
            flags: getResidueFlag(flag)
        } as SecondaryStructure.Helix
    } else if (kind === 'sheet') {
        return {
            kind: 'sheet',
            flags: getResidueFlag(flag)
        } as SecondaryStructure.Sheet
    } else if (kind === 'turn' || kind === 'bend') {
        return {
            kind: 'turn',
            flags: getResidueFlag(flag)
        }
    } else {
        return {
            kind: 'none'
        }
    }
}

function mapToKind(assignment: SecondaryStructureType.Flag) {
    if (assignment === SecondaryStructureType.SecondaryStructureDssp.H || assignment === SecondaryStructureType.SecondaryStructureDssp.G || assignment === SecondaryStructureType.SecondaryStructureDssp.I) {
        return 'helix'
    } else if (assignment === SecondaryStructureType.SecondaryStructureDssp.B || assignment === SecondaryStructureType.SecondaryStructureDssp.E) {
        return 'sheet'
    } else if (assignment === SecondaryStructureType.SecondaryStructureDssp.T) {
        return 'turn'
    } else if (assignment === SecondaryStructureType.SecondaryStructureDssp.S) {
        return 'bend'
    } else {
        return 'none'
    }
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

/**
 * Bend(i) =: [angle ((CW - Ca(i - 2)),(C"(i + 2) - C"(i))) > 70"]
 *
 * Type: S
 */
function assignBends(ctx: DSSPContext) {
    const flags = ctx.flags
    const { x, y, z } = ctx.conformation
    const { traceElementIndex } = ctx.hierarchy.derived.residue

    const proteinResidues = ctx.proteinResidues
    const residueCount = proteinResidues.length

    const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i])

    const caPosPrev2 = Vec3()
    const caPos = Vec3()
    const caPosNext2 = Vec3()

    const nIndices = ctx.backboneIndices.nIndices
    const cPos = Vec3()
    const nPosNext = Vec3()

    const caMinus2 = Vec3()
    const caPlus2 = Vec3()

    f1: for (let i = 2; i < residueCount - 2; i++) {
        // check for peptide bond
        for (let k = 0; k < 4; k++) {
            let index = i + k - 2
            position(traceElementIndex[index], cPos)
            position(nIndices[index + 1], nPosNext)
            if (Vec3.squaredDistance(cPos, nPosNext) > 6.25 /* max squared peptide bond distance allowed */) {
                continue f1
            }
        }

        const oRIprev2 = proteinResidues[i - 2]
        const oRI = proteinResidues[i]
        const oRInext2 = proteinResidues[i + 2]

        const caAtomPrev2 = traceElementIndex[oRIprev2]
        const caAtom = traceElementIndex[oRI]
        const caAtomNext2 = traceElementIndex[oRInext2]

        position(caAtomPrev2, caPosPrev2)
        position(caAtom, caPos)
        position(caAtomNext2, caPosNext2)

        Vec3.sub(caMinus2, caPosPrev2, caPos)
        Vec3.sub(caPlus2, caPos, caPosNext2)

        const angle = Vec3.angle(caMinus2, caPlus2) * 360 / (2 * Math.PI)
        if (angle && angle > 70.00) {
            flags[i] |= DSSPType.Flag.S
        }
    }
}

function calculateDihedralAngles(hierarchy: AtomicHierarchy, conformation: AtomicConformation, proteinResidues: SortedArray<ResidueIndex>, backboneIndices: BackboneAtomIndices): { phi: Float32Array, psi: Float32Array } {
    const { cIndices, nIndices } = backboneIndices
    const { index } = hierarchy
    const { x, y, z } = conformation
    const { traceElementIndex } = hierarchy.derived.residue

    const residueCount = proteinResidues.length
    const position = (i: number, v: Vec3) => i === -1 ? Vec3.setNaN(v) : Vec3.set(v, x[i], y[i], z[i])

    let cPosPrev = Vec3(), caPosPrev = Vec3(), nPosPrev = Vec3()
    let cPos = Vec3(), caPos = Vec3(), nPos = Vec3()
    let cPosNext = Vec3(), caPosNext = Vec3(), nPosNext = Vec3()

    if (residueCount === 0) return { phi: new Float32Array(0), psi: new Float32Array(0) }

    const phi: Float32Array = new Float32Array(residueCount - 1)
    const psi: Float32Array = new Float32Array(residueCount - 1)

    position(-1, cPosPrev)
    position(-1, caPosPrev)
    position(-1, nPosPrev)

    position(cIndices[0], cPos)
    position(traceElementIndex[proteinResidues[0]], caPos)
    position(nIndices[0], nPos)

    position(cIndices[1], cPosNext)
    position(traceElementIndex[proteinResidues[1]], caPosNext)
    position(nIndices[1], nPosNext)

    for (let i = 0; i < residueCount - 1; ++i) {
        // ignore C-terminal residue as acceptor
        if (index.findAtomOnResidue(proteinResidues[i], 'OXT') !== -1) continue

        // returns NaN for missing atoms
        phi[i] = Vec3.dihedralAngle(cPosPrev, nPos, caPos, cPos)
        psi[i] = Vec3.dihedralAngle(nPos, caPos, cPos, nPosNext)

        cPosPrev = cPos, caPosPrev = caPos, nPosPrev = nPos
        cPos = cPosNext, caPos = caPosNext, nPos = nPosNext

        position(cIndices[i + 1], cPosNext)
        position(traceElementIndex[proteinResidues[i + 1]], caPosNext)
        position(nIndices[i + 1], nPosNext)
    }

    return { phi, psi };
}

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

    const oPos = Vec3()
    const cPos = Vec3()
    const caPos = Vec3()
    const nPos = Vec3()
    const hPos = Vec3()

    const cPosPrev = Vec3()
    const oPosPrev = Vec3()

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
    if (DSSPType.is(f, DSSPType.Flag.E)) return SecondaryStructureType.SecondaryStructureDssp.E
    if (DSSPType.is(f, DSSPType.Flag.B)) return SecondaryStructureType.SecondaryStructureDssp.B
    if (DSSPType.is(f, DSSPType.Flag.G)) return SecondaryStructureType.SecondaryStructureDssp.G
    if (DSSPType.is(f, DSSPType.Flag.I)) return SecondaryStructureType.SecondaryStructureDssp.I
    if (DSSPType.is(f, DSSPType.Flag.T)) return SecondaryStructureType.SecondaryStructureDssp.T
    if (DSSPType.is(f, DSSPType.Flag.S)) return SecondaryStructureType.SecondaryStructureDssp.S
    return SecondaryStructureType.Flag.None
}

function getOriginalFlagName(f: DSSPType) {
    if (DSSPType.is(f, DSSPType.Flag.H)) return 'H'
    if (DSSPType.is(f, DSSPType.Flag.E)) return 'E'
    if (DSSPType.is(f, DSSPType.Flag.B)) return 'B'
    if (DSSPType.is(f, DSSPType.Flag.G)) return 'G'
    if (DSSPType.is(f, DSSPType.Flag.I)) return 'I'
    if (DSSPType.is(f, DSSPType.Flag.T)) return 'T'
    if (DSSPType.is(f, DSSPType.Flag.S)) return 'S'
    return '-'
}

/** Version 2.1.0 priority: I,H,B,E,G,T,S */
function getUpdatedResidueFlag(f: DSSPType) {
    if (DSSPType.is(f, DSSPType.Flag.I)) return SecondaryStructureType.SecondaryStructureDssp.I
    if (DSSPType.is(f, DSSPType.Flag.H)) return SecondaryStructureType.SecondaryStructureDssp.H
    if (DSSPType.is(f, DSSPType.Flag.E)) return SecondaryStructureType.SecondaryStructureDssp.E
    if (DSSPType.is(f, DSSPType.Flag.B)) return SecondaryStructureType.SecondaryStructureDssp.B
    if (DSSPType.is(f, DSSPType.Flag.G)) return SecondaryStructureType.SecondaryStructureDssp.G
    if (DSSPType.is(f, DSSPType.Flag.T)) return SecondaryStructureType.SecondaryStructureDssp.T
    if (DSSPType.is(f, DSSPType.Flag.S)) return SecondaryStructureType.SecondaryStructureDssp.S
    return SecondaryStructureType.Flag.None
}

function getUpdatedFlagName(f: DSSPType) {
    if (DSSPType.is(f, DSSPType.Flag.I)) return 'I'
    if (DSSPType.is(f, DSSPType.Flag.H)) return 'H'
    if (DSSPType.is(f, DSSPType.Flag.E)) return 'E'
    if (DSSPType.is(f, DSSPType.Flag.B)) return 'B'
    if (DSSPType.is(f, DSSPType.Flag.G)) return 'G'
    if (DSSPType.is(f, DSSPType.Flag.T)) return 'T'
    if (DSSPType.is(f, DSSPType.Flag.S)) return 'S'
    return '-'
}

function getDSSPAssignment(flags: Uint32Array, getResidueFlag: (f: DSSPType) => SecondaryStructureType) {
    const type = new Uint32Array(flags.length)
    for (let i = 0, il = flags.length; i < il; ++i) {
        const f = DSSPType.create(flags[i])
        type[i] = getResidueFlag(f)
    }

    return type as unknown as ArrayLike<SecondaryStructureType>
}

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
    const e = e1 + e2

    // cap lowest possible energy
    if (e < hbondEnergyMinimal)
        return hbondEnergyMinimal

    return e
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

    const turnFlag = [DSSPType.Flag.T3S, DSSPType.Flag.T4S, DSSPType.Flag.T5S, DSSPType.Flag.T3, DSSPType.Flag.T4, DSSPType.Flag.T5]

    for (let idx = 0; idx < 3; idx++) {
        for (let i = 0, il = proteinResidues.length - 1; i < il; ++i) {
            const rI = proteinResidues[i]
            const cI = chainAtomSegments.index[residueAtomSegments.offsets[rI]]

            // TODO should take sequence gaps into account
            const rN = proteinResidues[i + idx + 3]
            const cN = chainAtomSegments.index[residueAtomSegments.offsets[rN]]
            // check if on same chain
            if (!label_asym_id.areValuesEqual(cI, cN)) continue

            // check if hbond exists
            if (hbonds.getDirectedEdgeIndex(i, i + idx + 3) !== -1) {
                flags[i] |= turnFlag[idx + 3] | turnFlag[idx]
                if (ctx.params.oldDefinition) {
                    for (let k = 1; k < idx + 3; ++k) {
                        flags[i + k] |= turnFlag[idx + 3] | DSSPType.Flag.T
                    }
                } else {
                    for (let k = 0; k <= idx + 3; ++k) {
                        flags[i + k] |= turnFlag[idx + 3] | DSSPType.Flag.T
                    }
                }
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
    const { proteinResidues, hbonds, flags, bridges } = ctx

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
                // TODO move to constructor, actually omit object all together
                bridges[bridges.length] = new Bridge(i, j, BridgeType.PARALLEL)
            }

            // Parallel Bridge(i, j) =: [Hbond(j - 1, i) and Hbond(i, j + 1)]
            i = k
            j = l - 1 // l is j + 1
            if (i !== j && hbonds.getDirectedEdgeIndex(j - 1, i) !== -1) {
                flags[i] |= DSSPType.Flag.B
                flags[j] |= DSSPType.Flag.B
                bridges[bridges.length] = new Bridge(j, i, BridgeType.PARALLEL)
            }

            // Antiparallel Bridge(i, j) =: [Hbond(i, j) and Hbond(j, i)]
            i = k
            j = l
            if (i !== j && hbonds.getDirectedEdgeIndex(j, i) !== -1) {
                flags[i] |= DSSPType.Flag.B
                flags[j] |= DSSPType.Flag.B
                bridges[bridges.length] = new Bridge(j, i, BridgeType.ANTI_PARALLEL)
            }

            // Antiparallel Bridge(i, j) =: [Hbond(i - 1, j + 1) and Hbond(j - 1, i + l)]
            i = k + 1
            j = l - 1
            if (i !== j && hbonds.getDirectedEdgeIndex(j - 1, i + 1) !== -1) {
                flags[i] |= DSSPType.Flag.B
                flags[j] |= DSSPType.Flag.B
                bridges[bridges.length] = new Bridge(j, i, BridgeType.ANTI_PARALLEL)
            }
        }
    }

    bridges.sort((a, b) => a.partner1 > b.partner1 ? 1 : a.partner1 < b.partner1 ? -1 : 0)
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

    const turnFlag = [DSSPType.Flag.T3S, DSSPType.Flag.T4S, DSSPType.Flag.T5S, DSSPType.Flag.T3, DSSPType.Flag.T4, DSSPType.Flag.T5]
    const helixFlag = [0, 0, 0, DSSPType.Flag.G, DSSPType.Flag.H, DSSPType.Flag.I]

    const helixCheckOrder = ctx.params.oldOrdering ? [4, 3, 5] : [3, 4, 5]
    for (let ni = 0; ni < helixCheckOrder.length; ni++) {
        const n = helixCheckOrder[ni]

        for (let i = 1, il = proteinResidues.length - n; i < il; i++) {
            const fI = DSSPType.create(flags[i])
            const fI1 = DSSPType.create(flags[i - 1])
            const fI2 = DSSPType.create(flags[i + 1])

            // TODO rework to elegant solution which will not break instantly
            if (ctx.params.oldOrdering) {
                if ((n === 3 && (DSSPType.is(fI, DSSPType.Flag.H) || DSSPType.is(fI2, DSSPType.Flag.H)) || // for 3-10 yield to alpha helix
                    (n === 5 && ((DSSPType.is(fI, DSSPType.Flag.H) || DSSPType.is(fI, DSSPType.Flag.G)) || (DSSPType.is(fI2, DSSPType.Flag.H) || DSSPType.is(fI2, DSSPType.Flag.G)))))) { // for pi yield to all other helices
                    continue
                }
            } else {
                if ((n === 4 && (DSSPType.is(fI, DSSPType.Flag.G) || DSSPType.is(fI2, DSSPType.Flag.G)) || // for alpha helix yield to 3-10
                    (n === 5 && ((DSSPType.is(fI, DSSPType.Flag.H) || DSSPType.is(fI, DSSPType.Flag.G)) || (DSSPType.is(fI2, DSSPType.Flag.H) || DSSPType.is(fI2, DSSPType.Flag.G)))))) { // for pi yield to all other helices
                    continue
                }
            }

            if (DSSPType.is(fI, turnFlag[n]) && DSSPType.is(fI, turnFlag[n - 3]) && // check fI for turn start of proper type
                DSSPType.is(fI1, turnFlag[n]) && DSSPType.is(fI1, turnFlag[n - 3])) { // check fI1 accordingly
                if (ctx.params.oldDefinition) {
                    for (let k = 0; k < n; k++) {
                        flags[i + k] |= helixFlag[n]
                    }
                } else {
                    for (let k = -1; k <= n; k++) {
                        flags[i + k] |= helixFlag[n]
                    }
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
    const { bridges, ladders } = ctx

    // create ladders
    for (let bridgeIndex = 0; bridgeIndex < bridges.length; bridgeIndex++) {
        const bridge = bridges[bridgeIndex]
        let found = false
        for (let ladderIndex = 0; ladderIndex < ladders.length; ladderIndex++) {
            const ladder = ladders[ladderIndex]
            if (shouldExtendLadder(ladder, bridge)) {
                found = true
                ladder.firstEnd++
                if (bridge.type === BridgeType.PARALLEL) {
                    ladder.secondEnd++
                } else {
                    ladder.secondStart--
                }
            }
        }

        // no suitable assignment: create new ladder with single bridge as content
        if (!found) {
            ladders[ladders.length] = {
                previousLadder: 0,
                nextLadder: 0,
                firstStart: bridge.partner1,
                firstEnd: bridge.partner1,
                secondStart: bridge.partner2,
                secondEnd: bridge.partner2,
                type: bridge.type
            }
        }
    }

    // connect ladders
    for (let ladderIndex1 = 0; ladderIndex1 < ladders.length; ladderIndex1++) {
        const ladder1 = ladders[ladderIndex1]
        for (let ladderIndex2 = ladderIndex1; ladderIndex2 < ladders.length; ladderIndex2++) {
            const ladder2 = ladders[ladderIndex2]
            if (resemblesBulge(ladder1, ladder2)) {
                ladder1.nextLadder = ladderIndex2
                ladder2.previousLadder = ladderIndex1
            }
        }
    }
}

/**
 * For beta structures, we define: a bulge-linked ladder consists of two ladders or bridges of the same type
 * connected by at most one extra residue of one strand and at most four extra residues  on the other strand,
 * all residues in bulge-linked ladders are marked E, including any extra residues.
 */
function resemblesBulge(ladder1: Ladder, ladder2: Ladder) {
    if (!(ladder1.type === ladder2.type && ladder2.firstStart - ladder1.firstEnd < 6 &&
        ladder1.firstStart < ladder2.firstStart && ladder2.nextLadder === 0)) return false

    if (ladder1.type === BridgeType.PARALLEL) {
        return bulgeCriterion2(ladder1, ladder2)
    } else {
        return bulgeCriterion2(ladder2, ladder1)
    }
}

function bulgeCriterion2(ladder1: Ladder, ladder2: Ladder) {
    return ladder2.secondStart - ladder1.secondEnd > 0 && ((ladder2.secondStart - ladder1.secondEnd < 6 &&
        ladder2.firstStart - ladder1.firstEnd < 3) || ladder2.secondStart - ladder1.secondEnd < 3)
}

function shouldExtendLadder(ladder: Ladder, bridge: Bridge): boolean {
    // in order to extend ladders, same type must be present
    if (bridge.type !== ladder.type) return false

    // only extend if residue 1 is sequence neighbor with regard to ladder
    if (bridge.partner1 !== ladder.firstEnd + 1) return false

    if (bridge.type === BridgeType.PARALLEL) {
        if (bridge.partner2 === ladder.secondEnd + 1) {
            return true
        }
    } else {
        if (bridge.partner2 === ladder.secondStart - 1) {
            return true
        }
    }

    return false
}

function isHelixType(f: DSSPType) {
    return DSSPType.is(f, DSSPType.Flag.G) || DSSPType.is(f, DSSPType.Flag.H) || DSSPType.is(f, DSSPType.Flag.I)
}

/**
 * sheet=: set of one or more ladders connected by shared residues
 *
 * Type: E
 */
function assignSheets(ctx: DSSPContext) {
    const { ladders, flags } = ctx
    for (let ladderIndex = 0; ladderIndex < ladders.length; ladderIndex++) {
        const ladder = ladders[ladderIndex]
        for (let lcount = ladder.firstStart; lcount <= ladder.firstEnd; lcount++) {
            const diff = ladder.firstStart - lcount
            const l2count = ladder.secondStart - diff

            if (ladder.firstStart !== ladder.firstEnd) {
                flags[lcount] |= DSSPType.Flag.E
                flags[l2count] |= DSSPType.Flag.E
            } else {
                if (!isHelixType(flags[lcount]) && DSSPType.is(flags[lcount], DSSPType.Flag.E)) {
                    flags[lcount] |= DSSPType.Flag.B
                }
                if (!isHelixType(flags[l2count]) && DSSPType.is(flags[l2count], DSSPType.Flag.E)) {
                    flags[l2count] |= DSSPType.Flag.B
                }
            }
        }

        if (ladder.nextLadder === 0) continue

        const conladder = ladders[ladder.nextLadder]
        for (let lcount = ladder.firstStart; lcount <= conladder.firstEnd; lcount++) {
            flags[lcount] |= DSSPType.Flag.E
        }
        if (ladder.type === BridgeType.PARALLEL) {
            for (let lcount = ladder.secondStart; lcount <= conladder.secondEnd; lcount++) {
                flags[lcount] |= DSSPType.Flag.E
            }
        } else {
            for (let lcount = conladder.secondEnd; lcount <= ladder.secondStart; lcount++) {
                flags[lcount] |= DSSPType.Flag.E
            }
        }
    }
}