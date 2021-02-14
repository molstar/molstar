/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureElement, StructureProperties } from '../../mol-model/structure';
import { Task, RuntimeContext } from '../../mol-task';
import { CentroidHelper } from '../../mol-math/geometry/centroid-helper';
import { AccessibleSurfaceAreaProvider } from '../../mol-model-props/computed/accessible-surface-area';
import { Vec3 } from '../../mol-math/linear-algebra';
import { getElementMoleculeType } from '../../mol-model/structure/util';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { AccessibleSurfaceArea } from '../../mol-model-props/computed/accessible-surface-area/shrake-rupley';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { MembraneOrientation } from './prop';

interface ANVILContext {
    structure: Structure,

    numberOfSpherePoints: number,
    stepSize: number,
    minThickness: number,
    maxThickness: number,
    asaCutoff: number,

    offsets: ArrayLike<number>,
    exposed: ArrayLike<boolean>,
    centroid: Vec3,
    extent: number
};

export const ANVILParams = {
    numberOfSpherePoints: PD.Numeric(120, { min: 35, max: 700, step: 1 }, { description: 'Number of spheres/directions to test for membrane placement. Original value is 350.' }),
    stepSize: PD.Numeric(1, { min: 0.25, max: 4, step: 0.25 }, { description: 'Thickness of membrane slices that will be tested' }),
    minThickness: PD.Numeric(20, { min: 10, max: 30, step: 1}, { description: 'Minimum membrane thickness used during refinement' }),
    maxThickness: PD.Numeric(40, { min: 30, max: 50, step: 1}, { description: 'Maximum membrane thickness used during refinement' }),
    asaCutoff: PD.Numeric(40, { min: 10, max: 100, step: 1 }, { description: 'Absolute ASA cutoff above which residues will be considered' })
};
export type ANVILParams = typeof ANVILParams
export type ANVILProps = PD.Values<ANVILParams>

/**
 * Implements:
 * Membrane positioning for high- and low-resolution protein structures through a binary classification approach
 * Guillaume Postic, Yassine Ghouzam, Vincent Guiraud, and Jean-Christophe Gelly
 * Protein Engineering, Design & Selection, 2015, 1–5
 * doi: 10.1093/protein/gzv063
 */
export function computeANVIL(structure: Structure, props: ANVILProps) {
    return Task.create('Compute Membrane Orientation', async runtime => {
        return await calculate(runtime, structure, props);
    });
}


const centroidHelper = new CentroidHelper();
function initialize(structure: Structure, props: ANVILProps): ANVILContext {
    const l = StructureElement.Location.create(structure);
    const { label_atom_id, x, y, z } = StructureProperties.atom;
    const elementCount = structure.polymerResidueCount;
    centroidHelper.reset();

    let offsets = new Int32Array(elementCount);
    let exposed = new Array<boolean>(elementCount);

    const accessibleSurfaceArea = structure && AccessibleSurfaceAreaProvider.get(structure);
    const asa = accessibleSurfaceArea.value!;

    const vec = Vec3();
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
            centroidHelper.includeStep(vec);

            // keep track of offsets and exposed state to reuse
            offsets[m] = structure.serialMapping.getSerialIndex(l.unit, l.element);
            exposed[m] = AccessibleSurfaceArea.getValue(l, asa) > props.asaCutoff;

            m++;
        }
    }

    // omit potentially empty tail1
    offsets = offsets.slice(0, m);
    exposed = exposed.slice(0, m);

    // calculate centroid and extent
    centroidHelper.finishedIncludeStep();
    const centroid = centroidHelper.center;
    for (let k = 0, kl = offsets.length; k < kl; k++) {
        setLocation(l, structure, offsets[k]);
        Vec3.set(vec, x(l), y(l), z(l));
        centroidHelper.radiusStep(vec);
    }
    const extent = 1.2 * Math.sqrt(centroidHelper.radiusSq);

    return {
        ...props,
        structure: structure,

        offsets: offsets,
        exposed: exposed,
        centroid: centroid,
        extent: extent
    };
}

export async function calculate(runtime: RuntimeContext, structure: Structure, params: ANVILProps): Promise<MembraneOrientation> {
    const { label_comp_id } = StructureProperties.atom;

    const ctx = initialize(structure, params);
    const initialHphobHphil = HphobHphil.filtered(ctx, label_comp_id);

    const initialMembrane = findMembrane(ctx, generateSpherePoints(ctx, ctx.numberOfSpherePoints), initialHphobHphil, label_comp_id);
    const alternativeMembrane = findMembrane(ctx, findProximateAxes(ctx, initialMembrane), initialHphobHphil, label_comp_id);

    const membrane = initialMembrane.qmax! > alternativeMembrane.qmax! ? initialMembrane : alternativeMembrane;

    return {
        planePoint1: membrane.planePoint1,
        planePoint2: membrane.planePoint2,
        normalVector: membrane.normalVector!,
        radius: ctx.extent,
        centroid: ctx.centroid
    };
}

interface MembraneCandidate {
    planePoint1: Vec3,
    planePoint2: Vec3,
    stats: HphobHphil,
    normalVector?: Vec3,
    spherePoint?: Vec3,
    qmax?: number
}

namespace MembraneCandidate {
    export function initial(c1: Vec3, c2: Vec3, stats: HphobHphil): MembraneCandidate {
        return {
            planePoint1: c1,
            planePoint2: c2,
            stats: stats
        };
    }

    export function scored(spherePoint: Vec3, c1: Vec3, c2: Vec3, stats: HphobHphil, qmax: number, centroid: Vec3): MembraneCandidate {
        const diam_vect = Vec3();
        Vec3.sub(diam_vect, centroid, spherePoint);
        return {
            planePoint1: c1,
            planePoint2: c2,
            stats: stats,
            normalVector: diam_vect,
            spherePoint: spherePoint,
            qmax: qmax
        };
    }
}

function findMembrane(ctx: ANVILContext, spherePoints: Vec3[], initialStats: HphobHphil, label_comp_id: StructureElement.Property<string>): MembraneCandidate {
    const { centroid, stepSize, minThickness, maxThickness } = ctx;
    // best performing membrane
    let membrane: MembraneCandidate;
    // score of the best performing membrane
    let qmax = 0;

    // construct slices of thickness 1.0 along the axis connecting the centroid and the spherePoint
    const diam = Vec3();
    for (let i = 0, il = spherePoints.length; i < il; i++) {
        const spherePoint = spherePoints[i];
        Vec3.sub(diam, centroid, spherePoint);
        Vec3.scale(diam, diam, 2);
        const diamNorm = Vec3.magnitude(diam);
        const qvartemp = [];

        for (let i = 0, il = diamNorm - stepSize; i < il; i += stepSize) {
            const c1 = Vec3();
            const c2 = Vec3();
            Vec3.scaleAndAdd(c1, spherePoint, diam, i / diamNorm);
            Vec3.scaleAndAdd(c2, spherePoint, diam, (i + stepSize) / diamNorm);

            // evaluate how well this membrane slice embeddeds the peculiar residues
            const stats = HphobHphil.filtered(ctx, label_comp_id, (testPoint: Vec3) => isInMembranePlane(testPoint, diam, c1, c2));
            qvartemp.push(MembraneCandidate.initial(c1, c2, stats));
        }

        let jmax = (minThickness / stepSize) - 1;

        for (let width = 0, widthl = maxThickness; width < widthl;) {
            const imax = qvartemp.length - 1 - jmax;

            for (let i = 0, il = imax; i < il; i++) {
                const c1 = qvartemp[i].planePoint1;
                const c2 = qvartemp[i + jmax].planePoint2;

                let hphob = 0;
                let hphil = 0;
                let total = 0;
                for (let j = 0; j < jmax; j++) {
                    const ij = qvartemp[i + j];
                    if (j === 0 || j === jmax - 1) {
                        hphob += 0.5 * ij.stats.hphob;
                        hphil += 0.5 * ij.stats.hphil;
                    } else {
                        hphob += ij.stats.hphob;
                        hphil += ij.stats.hphil;
                    }
                    total += ij.stats.total;
                }

                const stats = HphobHphil.of(hphob, hphil, total);

                if (hphob !== 0) {
                    const qvaltest = qValue(stats, initialStats);
                    if (qvaltest > qmax) {
                        qmax = qvaltest;
                        membrane = MembraneCandidate.scored(spherePoint, c1, c2, HphobHphil.of(hphob, hphil, total), qmax, centroid);
                    }
                }
            }
            jmax++;
            width = (jmax + 1) * stepSize;
        }
    }

    return membrane!;
}

function qValue(currentStats: HphobHphil, initialStats: HphobHphil): number {
    if(initialStats.hphob < 1) {
        initialStats.hphob = 0.1;
    }

    if(initialStats.hphil < 1) {
        initialStats.hphil += 1;
    }

    const part_tot = currentStats.hphob + currentStats.hphil;
    return (currentStats.hphob * (initialStats.hphil - currentStats.hphil) - currentStats.hphil * (initialStats.hphob - currentStats.hphob)) /
            Math.sqrt(part_tot * initialStats.hphob * initialStats.hphil * (initialStats.hphob + initialStats.hphil - part_tot));
}

export function isInMembranePlane(testPoint: Vec3, normalVector: Vec3, planePoint1: Vec3, planePoint2: Vec3): boolean {
    const d1 = -Vec3.dot(normalVector, planePoint1);
    const d2 = -Vec3.dot(normalVector, planePoint2);
    const d = -Vec3.dot(normalVector, testPoint);
    return d > Math.min(d1, d2) && d < Math.max(d1, d2);
}

// generates a defined number of points on a sphere with radius = extent around the specified centroid
function generateSpherePoints(ctx: ANVILContext, numberOfSpherePoints: number): Vec3[] {
    const { centroid, extent } = ctx;
    const points = [];
    let oldPhi = 0, h, theta, phi;
    for(let k = 1, kl = numberOfSpherePoints + 1; k < kl; k++) {
        h = -1 + 2 * (k - 1) / (numberOfSpherePoints - 1);
        theta = Math.acos(h);
        phi = (k === 1 || k === numberOfSpherePoints) ? 0 : (oldPhi + 3.6 / Math.sqrt(numberOfSpherePoints * (1 - h * h))) % (2 * Math.PI);

        const point = Vec3.create(
            extent * Math.sin(phi) * Math.sin(theta) + centroid[0],
            extent * Math.cos(theta) + centroid[1],
            extent * Math.cos(phi) * Math.sin(theta) + centroid[2]
        );
        points[k - 1] = point;
        oldPhi = phi;
    }

    return points;
}

// generates sphere points close to that of the initial membrane
function findProximateAxes(ctx: ANVILContext, membrane: MembraneCandidate): Vec3[] {
    const { numberOfSpherePoints, extent } = ctx;
    const points = generateSpherePoints(ctx, 30000);
    let j = 4;
    let sphere_pts2: Vec3[] = [];
    while (sphere_pts2.length < numberOfSpherePoints) {
        const d = 2 * extent / numberOfSpherePoints + j;
        const dsq = d * d;
        sphere_pts2 = [];
        for (let i = 0, il = points.length; i < il; i++) {
            if (Vec3.squaredDistance(points[i], membrane.spherePoint!) < dsq) {
                sphere_pts2.push(points[i]);
            }
        }
        j += 0.2;
    }
    return sphere_pts2;
}

interface HphobHphil {
    hphob: number,
    hphil: number,
    total: number
}

namespace HphobHphil {
    export function of(hphob: number, hphil: number, total?: number) {
        return {
            hphob: hphob,
            hphil: hphil,
            total: !!total ? total : hphob + hphil
        };
    }

    const testPoint = Vec3();
    export function filtered(ctx: ANVILContext, label_comp_id: StructureElement.Property<string>, filter?: (test: Vec3) => boolean): HphobHphil {
        const { offsets, exposed, structure } = ctx;
        const l = StructureElement.Location.create(structure);
        const { x, y, z } = StructureProperties.atom;
        let hphob = 0;
        let hphil = 0;
        for (let k = 0, kl = offsets.length; k < kl; k++) {
            // ignore buried residues
            if (!exposed[k]) {
                continue;
            }

            setLocation(l, structure, offsets[k]);
            Vec3.set(testPoint, x(l), y(l), z(l));

            // testPoints have to be in putative membrane layer
            if (filter && !filter(testPoint)) {
                continue;
            }

            if (isHydrophobic(label_comp_id(l))) {
                hphob++;
            } else {
                hphil++;
            }
        }
        return of(hphob, hphil);
    }
}

// ANVIL-specific (not general) definition of membrane-favoring amino acids
const HYDROPHOBIC_AMINO_ACIDS = new Set(['ALA', 'CYS', 'GLY', 'HIS', 'ILE', 'LEU', 'MET', 'PHE', 'SER', 'THR', 'VAL']);
export function isHydrophobic(label_comp_id: string): boolean {
    return HYDROPHOBIC_AMINO_ACIDS.has(label_comp_id);
}

function setLocation(l: StructureElement.Location, structure: Structure, serialIndex: number) {
    l.structure = structure;
    l.unit = structure.units[structure.serialMapping.unitIndices[serialIndex]];
    l.element = structure.serialMapping.elementIndices[serialIndex];
    return l;
}