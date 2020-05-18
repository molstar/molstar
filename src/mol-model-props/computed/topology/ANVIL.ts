/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from '../../../mol-task';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, StructureElement, StructureProperties } from '../../../mol-model/structure';
import { getElementMoleculeType } from '../../../mol-model/structure/util';
import { MoleculeType } from '../../../mol-model/structure/model/types';
import { CentroidHelper } from '../../../mol-math/geometry/centroid-helper';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { AccessibleSurfaceArea } from '../accessible-surface-area/shrake-rupley';
import { AccessibleSurfaceAreaProvider } from '../accessible-surface-area';

export const ANVILParams = {
    numberOfSpherePoints: PD.Numeric(350),
    stepSize: PD.Numeric(1),
    minThickness: PD.Numeric(20, { min: 10, max: 30, step: 1}, { description: 'Minimum membrane thickness used during refinement' }),
    maxThickness: PD.Numeric(40, { min: 30, max: 50, step: 1}, { description: 'Maximum membrane thickness used during refinement' }),
    afilter: PD.Numeric(40),
    membranePointDensity: PD.Numeric(2, { min: 0.1, max: 10, step: 0.1 }, { description: 'Distance betwween points representing membrane layer'})
};
export type ANVILParams = typeof ANVILParams
export type ANVILProps = PD.Values<ANVILParams>

export { Topology };

interface Topology {
    readonly serialResidueIndex: ArrayLike<number>
    readonly exposure: ArrayLike<Topology>
}

namespace Topology {
    /**
     * Implements:
     * Membrane positioning for high- and low-resolution protein structures through a binary classification approach
     * Guillaume Postic, Yassine Ghouzam, Vincent Guiraud, and Jean-Christophe Gelly
     * Protein Engineering, Design & Selection, 2015, 1–5
     * doi: 10.1093/protein/gzv063
     */
    export function compute(structure: Structure, props: Partial<ANVILProps> = {}) {
        const p = { ...PD.getDefaultValues(ANVILParams), ...props };
        return Task.create('Compute Membrane Topology', async runtime => {
            return await calculate(runtime, structure, p);
        });
    }

    const l = StructureElement.Location.create(void 0);
    const centroidHelper = new CentroidHelper();
    const vec = Vec3();
    export async function calculate(runtime: RuntimeContext, structure: Structure, params: ANVILProps): Promise<Topology> {
        const { label_atom_id, x, y, z } = StructureProperties.atom;
        const elementCount = structure.elementCount;
        centroidHelper.reset();
        l.structure = structure;

        const offsets = new Int32Array(elementCount);
        const exposed: Array<Boolean> = new Array(elementCount);

        // ensure ASA
        const accessibleSurfaceArea = structure && AccessibleSurfaceAreaProvider.get(structure);
        const asa = accessibleSurfaceArea.value!;
        
        let m = 0;
        for (let i = 0, il = structure.units.length; i < il; ++i) {
            const unit = structure.units[i];
            const { elements } = unit;
            l.unit = unit;

            for (let j = 0, jl = elements.length; j < jl; ++j) {
                const eI = elements[j];
                l.element = eI;

                // consider only amino acids
                if (getElementMoleculeType(unit, eI) !== MoleculeType.Protein) {
                    continue;
                }

                // only CA is considered for downstream operations
                if (label_atom_id(l) !== 'CA') {
                    continue;
                }

                // while iterating use first pass to compute centroid
                Vec3.set(vec, x(l), y(l), z(l));
                // console.log(vec);
                centroidHelper.includeStep(vec);

                // keep track of offsets and exposed state to reuse
                offsets[m] = l.element;
                exposed[m] = AccessibleSurfaceArea.getValue(l, asa) > params.afilter;

                m++;
            }
        }

        centroidHelper.finishedIncludeStep();
        console.log(centroidHelper.center);

        for (let k = 0; k < m; k++) {
            setLocation(l, structure, offsets[k]);
            Vec3.set(vec, x(l), y(l), z(l));
            // console.log(vec);
            centroidHelper.radiusStep(vec);
        }

        console.log(1.2 * Math.sqrt(centroidHelper.radiusSq));
        
        return {
            exposure: new Array(),
            serialResidueIndex: new Array()
        };
    }

    function setLocation(l: StructureElement.Location, structure: Structure, serialIndex: number) {
        l.structure = structure;
        l.unit = structure.units[structure.serialMapping.unitIndices[serialIndex]];
        l.element = structure.serialMapping.elementIndices[serialIndex];
        return l;
    }

    export const enum Flag {
        NA = 0x0,
        Membrane = 0x1,
        NotMembrane = 0x2
    }

    export function getValue(location: StructureElement.Location, topology: Topology) {
        const { getSerialIndex } = location.structure.root.serialMapping;
        const { exposure, serialResidueIndex } = topology;
        const rSI = serialResidueIndex[getSerialIndex(location.unit, location.element)];
        if (rSI === -1) return -1;
        return exposure[rSI];
    }

    export function getFlag(location: StructureElement.Location, topology: Topology) {
        const value = getValue(location, topology);
        return value === -1 ? Flag.NA :
            value ? Flag.NotMembrane :
                Flag.Membrane;
    }
}