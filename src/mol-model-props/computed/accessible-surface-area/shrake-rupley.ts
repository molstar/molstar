/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Task, RuntimeContext } from '../../../mol-task';
import { BitFlags } from '../../../mol-util';
import { ParamDefinition as PD } from '../../../mol-util/param-definition'
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Structure } from '../../../mol-model/structure';
import { assignRadiusForHeavyAtoms } from './shrake-rupley/radii';
import { ShrakeRupleyContext, VdWLookup } from './shrake-rupley/common';
import { computePerResidue } from './shrake-rupley/per-residue';
import { normalizeAccessibleSurfaceArea } from './shrake-rupley/normalize';

export const ShrakeRupleyComputationParams = {
    numberOfSpherePoints: PD.Numeric(92, {}, { description: 'number of sphere points to sample per atom: 92 (original paper), 960 (BioJava), 3000 (EPPIC) - see Shrake A, Rupley JA: Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 1973.' }),
    probeSize: PD.Numeric(1.4, {}, { description: 'corresponds to the size of a water molecule: 1.4 (original paper), 1.5 (occassionally used)' }),
    buriedRasaThreshold: PD.Numeric(0.16, { min: 0.0, max: 1.0 }, { description: 'below this cutoff of relative accessible surface area a residue will be considered buried - see: Rost B, Sander C: Conservation and prediction of solvent accessibility in protein families. Proteins 1994.' }),
    nonPolymer: PD.Boolean(false, { description: 'Include non-polymer atoms in computation.' })
}
export type ShrakeRupleyComputationParams = typeof ShrakeRupleyComputationParams

namespace AccessibleSurfaceArea {
    /**
     * Adapts the BioJava implementation by Jose Duarte. That implementation is based on the publication by Shrake, A., and
     * J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin." JMB (1973).
     */
    export function compute(structure: Structure,
        params: Partial<PD.Values<ShrakeRupleyComputationParams>> = {}) {
        params = { ...PD.getDefaultValues(ShrakeRupleyComputationParams), ...params };
        return Task.create('Compute Accessible Surface Area', async rtctx => {
            return await _compute(rtctx, structure, params);
        }).run();
    }

    async function _compute(rtctx: RuntimeContext, structure: Structure, params: Partial<PD.Values<ShrakeRupleyComputationParams>> = {}): Promise<AccessibleSurfaceArea> {
        const ctx = initialize(rtctx, structure, params);

        assignRadiusForHeavyAtoms(ctx);
        computePerResidue(ctx);
        normalizeAccessibleSurfaceArea(ctx);

        return {
            accessibleSurfaceArea: ctx.accessibleSurfaceArea,
            relativeAccessibleSurfaceArea: ctx.relativeAccessibleSurfaceArea,
            buried: (index: number) => ctx.relativeAccessibleSurfaceArea[index] < 0.16
        };
    }

    function initialize(rtctx: RuntimeContext, structure: Structure, params: Partial<PD.Values<ShrakeRupleyComputationParams>>): ShrakeRupleyContext {
        const { elementCount, polymerResidueCount } = structure;

        return {
            rtctx: rtctx,
            structure: structure,
            probeSize: params.probeSize!,
            nonPolymer: params.nonPolymer!,
            spherePoints: generateSpherePoints(params.numberOfSpherePoints!),
            cons: 4.0 * Math.PI / params.numberOfSpherePoints!,
            maxLookupRadius: 2 * params.probeSize! + 2 * VdWLookup[2], // 2x probe size + 2x largest VdW
            atomRadius: new Int8Array(elementCount),
            accessibleSurfaceArea: new Float32Array(polymerResidueCount),
            relativeAccessibleSurfaceArea: new Float32Array(polymerResidueCount),
            updateChunk: 25000
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
    readonly accessibleSurfaceArea: ArrayLike<number>
    readonly relativeAccessibleSurfaceArea: ArrayLike<number>
    buried(index: number): boolean
}

export { AccessibleSurfaceArea }