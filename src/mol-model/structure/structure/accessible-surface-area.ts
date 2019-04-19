/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import Structure from './structure';
import { Task, RuntimeContext } from 'mol-task';
import { BitFlags } from 'mol-util';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { Vec3 } from 'mol-math/linear-algebra';
import { isPolymer, ElementSymbol, isNucleic, MoleculeType } from '../model/types';
import { VdwRadius } from '../model/properties/atomic';
import { isHydrogen, getElementIdx } from './unit/links/common';

namespace AccessibleSurfaceArea {
    // Chothia's amino acid and nucleotide atom vdw radii
    export const VdWLookup = [
        -1.0, // 0: missing
        1.76, // 1: trigonal C
        1.87, // 2: tetrahedral C
        1.65, // 3: trigonal N
        1.50, // 4: tetrahedral N
        1.40, // 5: O
        1.85, // 6: S
        1.80, // 7: C (nucleic)
        1.60, // 8: N (nucleic)
        1.40  // 9: P (nucleic)
    ] // can still be appended on-the-fly for rare elements like selenium

    /**
     * Adapts the BioJava implementation by Jose Duarte. That implementation is based on the publication by Shrake, A., and
     * J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin." JMB (1973).
     */
    export function compute(structure: Structure,
        params: Partial<PD.Values<AccessibleSurfaceAreaComputationParams>> = {}) {
        params = { ...PD.getDefaultValues(AccessibleSurfaceAreaComputationParams), ...params };
        return Task.create('Compute Accessible Surface Area', async rtctx => {
            return await _compute(rtctx, structure, params);
        }).run();
    }

    async function _compute(rtctx: RuntimeContext, structure: Structure, params: Partial<PD.Values<AccessibleSurfaceAreaComputationParams>> = {}): Promise<AccessibleSurfaceArea> {
        const ctx = initialize(rtctx, structure, params);

        assignRadiusForHeavyAtoms(ctx);
        computePerResidue(ctx);
        normalizeAccessibleSurfaceArea(ctx);

        return {
            accessibleSurfaceArea: ctx.accessibleSurfaceArea!,
            relativeAccessibleSurfaceArea: ctx.relativeAccessibleSurfaceArea!,
            buried: (index: number) => ctx.relativeAccessibleSurfaceArea![index] < 0.16 // TODO this doesnt seem super elegant
        };
    }

    interface AccessibleSurfaceAreaContext {
        rtctx: RuntimeContext,
        structure: Structure,
        params: Partial<PD.Values<AccessibleSurfaceAreaComputationParams>>,
        spherePoints: Vec3[],
        cons: number,
        maxLookupRadius: number,
        atomRadius?: Int8Array,
        accessibleSurfaceArea?: Float32Array,
        relativeAccessibleSurfaceArea?: Float32Array
    }

    function normalizeAccessibleSurfaceArea(ctx: AccessibleSurfaceAreaContext) {
        const { accessibleSurfaceArea, relativeAccessibleSurfaceArea, structure } = ctx;
        const { residues, derived } = structure.model.atomicHierarchy;

        for (let i = 0; i < residues.label_comp_id.rowCount; ++i) {
            // skip entities not part of a polymer chain
            if (!ctx.params.nonPolymer) {
                if (!isPolymer(derived.residue.moleculeType[i])) continue;
            }

            const maxAsa = (MaxAsa as any)[residues.label_comp_id.value(i)];
            const rasa = accessibleSurfaceArea![i] / (maxAsa === undefined ? DefaultMaxAsa : maxAsa);
            relativeAccessibleSurfaceArea![i] = rasa;
        }
    }

    async function computePerResidue(ctx: AccessibleSurfaceAreaContext) {
        const { structure, atomRadius, accessibleSurfaceArea, spherePoints, cons, params, maxLookupRadius } = ctx;
        const { probeSize } = params;
        const { model, elementCount: atomCount } = structure;
        const { x, y, z } = model.atomicConformation;
        const { residueAtomSegments } = model.atomicHierarchy;
        const { lookup3d } = structure;

        const position = (i: number, v: Vec3) => Vec3.set(v, x[i], y[i], z[i]);
        const aPos = Vec3.zero();
        const bPos = Vec3.zero();
        let testPoint = Vec3.zero();

        for (let aI = 0; aI < atomCount; ++aI) {
            if (aI % 10000 === 0) {
                if (ctx.rtctx.shouldUpdate) {
                    await ctx.rtctx.update({ message: 'calculating accessible surface area', current: aI, max: atomCount });
                }
            }

            const radius1 = VdWLookup[atomRadius![aI]];
            if (radius1 === VdWLookup[0]) continue;

            // pre-filter by lookup3d
            // 36275 ms - lookup ~3000 ms
            const { count, units, indices, squaredDistances } = lookup3d.find(x[aI], y[aI], z[aI], maxLookupRadius);

            // collect neighbors for each atom
            const cutoff1 = probeSize! + probeSize! + radius1;
            const neighbors = [];
            for (let iI = 0; iI < count; ++iI) {
                const bI = units[iI].elements[indices[iI]];
                const radius2 = VdWLookup[atomRadius![bI]];
                if (aI === bI || radius2 === VdWLookup[0]) continue;

                const cutoff2 = (cutoff1 + radius2) * (cutoff1 + radius2);
                if (squaredDistances[iI] < cutoff2) {
                    neighbors[neighbors.length] = bI;
                }
            }

            // for all neighbors: test all sphere points
            position(aI, aPos);
            const scalar = probeSize! + radius1;
            let accessiblePointCount = 0;
            for (let sI = 0; sI < spherePoints.length; ++sI) {
                const spherePoint = spherePoints[sI];
                testPoint = [spherePoint[0] * scalar + aPos[0], spherePoint[1] * scalar + aPos[1], spherePoint[2] * scalar + aPos[2]] as Vec3;
                let accessible = true;

                for (let _nI = 0; _nI < neighbors.length; ++_nI) {
                    const nI = neighbors[_nI];
                    position(nI, bPos);
                    const radius3 = VdWLookup[atomRadius![nI]];
                    const cutoff3 = (radius3 + probeSize!) * (radius3 + probeSize!);
                    if (Vec3.squaredDistance(testPoint, bPos) < cutoff3) {
                        accessible = false;
                        break;
                    }
                }

                if (accessible) ++accessiblePointCount;
            }

            accessibleSurfaceArea![residueAtomSegments.index[aI]] += cons * accessiblePointCount * scalar * scalar;
        }
    }

    function assignRadiusForHeavyAtoms(ctx: AccessibleSurfaceAreaContext) {
        const { structure } = ctx;
        const { model, elementCount: atomCount } = structure;
        const { atoms: atomInfo, derived, residues, residueAtomSegments } = model.atomicHierarchy;
        const { label_comp_id } = residues;
        const { moleculeType } = derived.residue;
        const { type_symbol, label_atom_id } = atomInfo;
        const residueCount = moleculeType.length;

        // with atom and residue count at hand: initialize arrays
        ctx.atomRadius = new Int8Array(atomCount - 1);
        ctx.accessibleSurfaceArea = new Float32Array(residueCount - 1);
        ctx.relativeAccessibleSurfaceArea = new Float32Array(residueCount - 1);

        for (let aI = 0; aI < atomCount; ++aI) {
            const rI = residueAtomSegments.index[aI];
            const element = type_symbol.value(aI);
            const elementIdx = getElementIdx(element);
            // skip hydrogen atoms
            if (isHydrogen(elementIdx)) {
                ctx.atomRadius[aI] = VdWLookup[0];
                continue;
            }

            const residueType = moleculeType[rI];
            // skip non-polymer groups
            if (!ctx.params.nonPolymer) {
                if (!isPolymer(residueType)) {
                    ctx.atomRadius[aI] = VdWLookup[0];
                    continue;
                }
            }

            const atomId = label_atom_id.value(aI);
            let compId = label_comp_id.value(rI);

            // handle modified residues
            const parentId = model.properties.modifiedResidues.parentId.get(compId);
            if (parentId !== void 0) compId = parentId;

            if (isNucleic(residueType)) {
                 ctx.atomRadius[aI] = determineRadiusNucl(atomId, element, compId);
            } else if (residueType === MoleculeType.protein) {
                ctx.atomRadius[aI] = determineRadiusAmino(atomId, element, compId);
            } else {
                ctx.atomRadius[aI] = handleNonStandardCase(element);
            }
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
            return 5;
            case 'S':
            return 6;
            case 'N':
            return atomId === 'NZ' ? 4 : 3;
            case 'C':
            switch (atomId) {
                case 'C': case 'CE1': case'CE2': case 'CE3': case 'CH2': case 'CZ': case 'CZ2': case 'CZ3':
                return 1;
                case 'CA': case 'CB': case 'CE': case 'CG1': case 'CG2':
                return 2;
                default:
                switch (compId) {
                    case 'PHE': case 'TRP': case 'TYR': case 'HIS': case 'ASP': case 'ASN':
                    return 1;
                    case 'PRO': case 'LYS': case 'ARG': case 'MET': case 'ILE': case 'LEU':
                    return 2;
                    case 'GLU': case 'GLN':
                    return atomId === 'CD' ? 1 : 2;
                }
            }
        }
        return handleNonStandardCase(element);
    }

    function determineRadiusNucl(atomId: string, element: ElementSymbol, compId: string): number {
        switch (element) {
            case 'C':
            return 7;
            case 'N':
            return 8;
            case 'P':
            return 9;
            case 'O':
            return 5;
        }
        return handleNonStandardCase(element);
    }

    function handleNonStandardCase(element: ElementSymbol): number {
        const radius = VdwRadius(element);
        let index = VdWLookup.indexOf(radius);
        if (index === -1) {
            // add novel value to lookup array
            index = VdWLookup.length;
            VdWLookup[index] = radius;
        }
        return index;
    }

    function initialize(rtctx: RuntimeContext, structure: Structure, params: Partial<PD.Values<AccessibleSurfaceAreaComputationParams>>): AccessibleSurfaceAreaContext {
        return {
            rtctx: rtctx,
            structure: structure,
            params: params,
            spherePoints: generateSpherePoints(params.numberOfSpherePoints!),
            cons: 4.0 * Math.PI / params.numberOfSpherePoints!,
            maxLookupRadius: 2 * params.probeSize! + 2 * VdWLookup[2] // 2x probe size + 2x largest VdW
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

    /** Maximum accessible surface area observed for amino acids. Taken from: http://dx.doi.org/10.1371/journal.pone.0080635 */
    export const MaxAsa = {
        'ALA': 121.0,
        'ARG': 265.0,
        'ASN': 187.0,
        'ASP': 187.0,
        'CYS': 148.0,
        'GLU': 214.0,
        'GLN': 214.0,
        'GLY': 97.0,
        'HIS': 216.0,
        'ILE': 195.0,
        'LEU': 191.0,
        'LYS': 230.0,
        'MET': 203.0,
        'PHE': 228.0,
        'PRO': 154.0,
        'SER': 143.0,
        'THR': 163.0,
        'TRP': 264.0,
        'TYR': 255.0,
        'VAL': 165.0
    }
    export const DefaultMaxAsa = 121.0

    export const AccessibleSurfaceAreaComputationParams = {
        numberOfSpherePoints: PD.Numeric(92, {}, { description: 'number of sphere points to sample per atom: 92 (original paper), 960 (BioJava), 3000 (EPPIC) - see Shrake A, Rupley JA: Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 1973.' }),
        probeSize: PD.Numeric(1.4, {}, { description: 'corresponds to the size of a water molecule: 1.4 (original paper), 1.5 (occassionally used)' }),
        buriedRasaThreshold: PD.Numeric(0.16, { min: 0.0, max: 1.0 }, { description: 'below this cutoff of relative accessible surface area a residue will be considered buried - see: Rost B, Sander C: Conservation and prediction of solvent accessibility in protein families. Proteins 1994.' }),
        nonPolymer: PD.Boolean(false, { description: 'Include non-polymer atoms in computation.' })
    }
    export type AccessibleSurfaceAreaComputationParams = typeof AccessibleSurfaceAreaComputationParams

    export namespace SolventAccessibility {
        export const is: (t: number, f: Flag) => boolean = BitFlags.has
        export const create: (f: Flag) => number = BitFlags.create
        export const enum Flag {
            _ = 0x0,
            BURIED = 0x1,
            ACCESSIBLE = 0x2
        }
    }
}

interface AccessibleSurfaceArea {
    readonly accessibleSurfaceArea?: ArrayLike<number>
    readonly relativeAccessibleSurfaceArea?: ArrayLike<number>
    buried(index: number): boolean
}

export { AccessibleSurfaceArea }