/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from '../../../mol-task';
// import { BitFlags } from '../../../mol-util';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { Structure, StructureElement, StructureProperties } from '../../../mol-model/structure';
import { assignRadiusForHeavyAtoms } from './shrake-rupley/radii';
import { ShrakeRupleyContext, VdWLookup, MaxAsa, DefaultMaxAsa } from './shrake-rupley/common';
import { computeArea } from './shrake-rupley/area';

export const ShrakeRupleyComputationParams = {
    numberOfSpherePoints: PD.Numeric(92, { min: 12, max: 360, step: 1 }, { description: 'Number of sphere points to sample per atom: 92 (original paper), 960 (BioJava), 3000 (EPPIC) - see Shrake A, Rupley JA: Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 1973.' }),
    probeSize: PD.Numeric(1.4, { min: 0.1, max: 4, step: 0.01 }, { description: 'Corresponds to the size of a water molecule: 1.4 (original paper), 1.5 (occassionally used)' }),
    // buriedRasaThreshold: PD.Numeric(0.16, { min: 0.0, max: 1.0 }, { description: 'below this cutoff of relative accessible surface area a residue will be considered buried - see: Rost B, Sander C: Conservation and prediction of solvent accessibility in protein families. Proteins 1994.' }),
    nonPolymer: PD.Boolean(false, { description: 'Include non-polymer atoms as occluders.' })
};
export type ShrakeRupleyComputationParams = typeof ShrakeRupleyComputationParams
export type ShrakeRupleyComputationProps = PD.Values<ShrakeRupleyComputationParams>

// TODO
// - add back buried and relative asa

export { AccessibleSurfaceArea };

interface AccessibleSurfaceArea {
    readonly serialResidueIndex: ArrayLike<number>
    readonly area: ArrayLike<number>
}

namespace AccessibleSurfaceArea {
    /**
     * Adapts the BioJava implementation by Jose Duarte. That implementation is based on the publication by Shrake, A., and
     * J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin." JMB (1973).
     */
    export function compute(structure: Structure, props: Partial<ShrakeRupleyComputationProps> = {}) {
        const p = { ...PD.getDefaultValues(ShrakeRupleyComputationParams), ...props };
        return Task.create('Compute Accessible Surface Area', async runtime => {
            return await calculate(runtime, structure, p);
        });
    }

    async function calculate(runtime: RuntimeContext, structure: Structure, props: ShrakeRupleyComputationProps): Promise<AccessibleSurfaceArea> {
        const ctx = initialize(structure, props);

        assignRadiusForHeavyAtoms(ctx);
        await computeArea(runtime, ctx);

        const { area, serialResidueIndex } = ctx;
        return { area, serialResidueIndex };
    }

    function initialize(structure: Structure, props: ShrakeRupleyComputationProps): ShrakeRupleyContext {
        const { elementCount, atomicResidueCount } = structure;
        const { probeSize, nonPolymer, numberOfSpherePoints } = props;

        return {
            structure,
            probeSize,
            nonPolymer,
            spherePoints: generateSpherePoints(numberOfSpherePoints),
            scalingConstant: 4.0 * Math.PI / numberOfSpherePoints,
            maxLookupRadius: 2 * props.probeSize + 2 * VdWLookup[2], // 2x probe size + 2x largest VdW
            atomRadiusType: new Int8Array(elementCount),
            serialResidueIndex: new Int32Array(elementCount),
            area: new Float32Array(atomicResidueCount)
        };
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
            points[points.length] = Vec3.create(Math.cos(phi) * r, y, Math.sin(phi) * r);
        }
        return points;
    }

    export const enum Flag {
        NA = 0x0,
        Buried = 0x1,
        Accessible = 0x2
    }

    /** Get relative area for a given component id */
    export function normalize(compId: string, asa: number) {
        const maxAsa = MaxAsa[compId] || DefaultMaxAsa;
        return asa / maxAsa;
    }

    export function getValue(location: StructureElement.Location, accessibleSurfaceArea: AccessibleSurfaceArea) {
        const { getSerialIndex } = location.structure.root.serialMapping;
        const { area, serialResidueIndex } = accessibleSurfaceArea;
        const rSI = serialResidueIndex[getSerialIndex(location.unit, location.element)];
        if (rSI === -1) return -1;
        return area[rSI];
    }

    export function getNormalizedValue(location: StructureElement.Location, accessibleSurfaceArea: AccessibleSurfaceArea) {
        const value = getValue(location, accessibleSurfaceArea);
        return value === -1 ? -1 : normalize(StructureProperties.atom.label_comp_id(location), value);
    }

    export function getFlag(location: StructureElement.Location, accessibleSurfaceArea: AccessibleSurfaceArea) {
        const value = getNormalizedValue(location, accessibleSurfaceArea);
        return value === -1 ? Flag.NA :
            value < 0.16 ? Flag.Buried :
                Flag.Accessible;
    }
}