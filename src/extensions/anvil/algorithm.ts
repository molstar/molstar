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

const LARGE_CA_THRESHOLD = 5000;
const UPDATE_INTERVAL = 10;

interface ANVILContext {
    structure: Structure,

    numberOfSpherePoints: number,
    stepSize: number,
    minThickness: number,
    maxThickness: number,
    asaCutoff: number,
    adjust: number,

    offsets: ArrayLike<number>,
    exposed: ArrayLike<number>,
    centroid: Vec3,
    extent: number
};

export const ANVILParams = {
    numberOfSpherePoints: PD.Numeric(140, { min: 35, max: 700, step: 1 }, { description: 'Number of spheres/directions to test for membrane placement. Original value is 350.' }),
    stepSize: PD.Numeric(1, { min: 0.25, max: 4, step: 0.25 }, { description: 'Thickness of membrane slices that will be tested' }),
    minThickness: PD.Numeric(20, { min: 10, max: 30, step: 1}, { description: 'Minimum membrane thickness used during refinement' }),
    maxThickness: PD.Numeric(40, { min: 30, max: 50, step: 1}, { description: 'Maximum membrane thickness used during refinement' }),
    asaCutoff: PD.Numeric(40, { min: 10, max: 100, step: 1 }, { description: 'Relative ASA cutoff above which residues will be considered' }),
    adjust: PD.Numeric(14, { min: 0, max: 30, step: 1 }, { description: 'Minimum length of membrane-spanning regions (original values: 14 for alpha-helices and 5 for beta sheets). Set to 0 to not optimize membrane thickness.' })
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
    const { label_atom_id, label_comp_id, x, y, z } = StructureProperties.atom;
    const elementCount = structure.polymerResidueCount;
    centroidHelper.reset();

    let offsets = new Array<number>(elementCount);
    let exposed = new Array<number>(elementCount);

    const accessibleSurfaceArea = structure && AccessibleSurfaceAreaProvider.get(structure);
    const asa = accessibleSurfaceArea.value!;
    const asaCutoff = props.asaCutoff / 100;

    const vec = Vec3();
    let m = 0;
    let n = 0;
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
            if (label_atom_id(l) !== 'CA' && label_atom_id(l) !== 'BB') {
                continue;
            }

            // original ANVIL only considers canonical amino acids
            if (!MaxAsa[label_comp_id(l)]) {
                continue;
            }

            // while iterating use first pass to compute centroid
            Vec3.set(vec, x(l), y(l), z(l));
            centroidHelper.includeStep(vec);

            // keep track of offsets and exposed state to reuse
            offsets[m++] = structure.serialMapping.getSerialIndex(l.unit, l.element);
            if (AccessibleSurfaceArea.getValue(l, asa) / MaxAsa[label_comp_id(l)] > asaCutoff) {
                exposed[n++] = structure.serialMapping.getSerialIndex(l.unit, l.element);
            }
        }
    }

    // omit potentially empty tail
    offsets = offsets.slice(0, m);
    exposed = exposed.slice(0, n);

    // calculate centroid and extent
    centroidHelper.finishedIncludeStep();
    const centroid = Vec3.clone(centroidHelper.center);
    for (let k = 0, kl = offsets.length; k < kl; k++) {
        setLocation(l, structure, offsets[k]);
        Vec3.set(vec, x(l), y(l), z(l));
        centroidHelper.radiusStep(vec);
    }
    const extent = 1.2 * Math.sqrt(centroidHelper.radiusSq);

    return {
        ...props,
        structure,

        offsets,
        exposed,
        centroid,
        extent
    };
}

export async function calculate(runtime: RuntimeContext, structure: Structure, params: ANVILProps): Promise<MembraneOrientation> {
    const ctx = initialize(structure, params);
    const initialHphobHphil = HphobHphil.filtered(ctx);

    const initialMembrane = (await findMembrane(runtime, 'Placing initial membrane...', ctx, generateSpherePoints(ctx, ctx.numberOfSpherePoints), initialHphobHphil))!;
    const refinedMembrane = (await findMembrane(runtime, 'Refining membrane placement...', ctx, findProximateAxes(ctx, initialMembrane), initialHphobHphil))!;
    let membrane = initialMembrane.qmax! > refinedMembrane.qmax! ? initialMembrane : refinedMembrane;

    if (ctx.adjust && ctx.offsets.length < LARGE_CA_THRESHOLD) {
        membrane = await adjustThickness(runtime, 'Adjusting membrane thickness...', ctx, membrane, initialHphobHphil);
    }

    const normalVector = Vec3.zero();
    const center =  Vec3.zero();
    const jitter = [0.001, 0.001, 0.001] as Vec3;
    Vec3.sub(normalVector, membrane.planePoint1, membrane.planePoint2);
    Vec3.normalize(normalVector, normalVector);

    // prevent degenerate matrices (e.g., 5cbg - assembly 1 which is oriented along the y-axis)
    if (Math.abs(normalVector[0]) === 1 || Math.abs(normalVector[1]) === 1 || Math.abs(normalVector[2]) === 1) {
        Vec3.add(normalVector, normalVector, jitter);
    }

    Vec3.add(center, membrane.planePoint1, membrane.planePoint2);
    Vec3.scale(center, center, 0.5);
    Vec3.add(center, center, jitter);
    const extent = adjustExtent(ctx, membrane, center);

    return {
        planePoint1: membrane.planePoint1,
        planePoint2: membrane.planePoint2,
        normalVector,
        centroid: center,
        radius: extent
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
            stats
        };
    }

    export function scored(spherePoint: Vec3, planePoint1: Vec3, planePoint2: Vec3, stats: HphobHphil, qmax: number, centroid: Vec3): MembraneCandidate {
        const normalVector = Vec3();
        Vec3.sub(normalVector, centroid, spherePoint);
        return {
            planePoint1,
            planePoint2,
            stats,
            normalVector,
            spherePoint,
            qmax
        };
    }
}

async function findMembrane(runtime: RuntimeContext, message: string | undefined, ctx: ANVILContext, spherePoints: Vec3[], initialStats: HphobHphil): Promise<MembraneCandidate | undefined> {
    const { centroid, stepSize, minThickness, maxThickness } = ctx;
    // best performing membrane
    let membrane: MembraneCandidate | undefined;
    // score of the best performing membrane
    let qmax = 0;

    // construct slices of thickness 1.0 along the axis connecting the centroid and the spherePoint
    const diam = Vec3();
    for (let n = 0, nl = spherePoints.length; n < nl; n++) {
        if (runtime?.shouldUpdate && message && (n + 1) % UPDATE_INTERVAL === 0) {
            await runtime.update({ message, current: (n + 1), max: nl });
        }

        const spherePoint = spherePoints[n];
        Vec3.sub(diam, centroid, spherePoint);
        Vec3.scale(diam, diam, 2);
        const diamNorm = Vec3.magnitude(diam);
        const qvartemp = [];

        for (let i = 0, il = diamNorm - stepSize; i < il; i += stepSize) {
            const c1 = Vec3();
            const c2 = Vec3();
            Vec3.scaleAndAdd(c1, spherePoint, diam, i / diamNorm);
            Vec3.scaleAndAdd(c2, spherePoint, diam, (i + stepSize) / diamNorm);
            const d1 = -Vec3.dot(diam, c1);
            const d2 = -Vec3.dot(diam, c2);
            const dMin = Math.min(d1, d2);
            const dMax = Math.max(d1, d2);

            // evaluate how well this membrane slice embeddeds the peculiar residues
            const stats = HphobHphil.filtered(ctx, { normalVector: diam, dMin, dMax });
            qvartemp.push(MembraneCandidate.initial(c1, c2, stats));
        }

        let jmax = Math.floor((minThickness / stepSize) - 1);

        for (let width = 0, widthl = maxThickness; width <= widthl;) {
            for (let i = 0, il = qvartemp.length - 1 - jmax; i < il; i++) {
                let hphob = 0;
                let hphil = 0;
                for (let j = 0; j < jmax; j++) {
                    const ij = qvartemp[i + j];
                    if (j === 0 || j === jmax - 1) {
                        hphob += Math.floor(0.5 * ij.stats.hphob);
                        hphil += 0.5 * ij.stats.hphil;
                    } else {
                        hphob += ij.stats.hphob;
                        hphil += ij.stats.hphil;
                    }
                }

                if (hphob !== 0) {
                    const stats = HphobHphil.of(hphob, hphil);
                    const qvaltest = qValue(stats, initialStats);
                    if (qvaltest >= qmax) {
                        qmax = qvaltest;
                        membrane = MembraneCandidate.scored(spherePoint, qvartemp[i].planePoint1, qvartemp[i + jmax].planePoint2, stats, qmax, centroid);
                    }
                }
            }
            jmax++;
            width = (jmax + 1) * stepSize;
        }
    }

    return membrane;
}

/** Adjust membrane thickness by maximizing the number of membrane segments. */
async function adjustThickness(runtime: RuntimeContext, message: string | undefined, ctx: ANVILContext, membrane: MembraneCandidate, initialHphobHphil: HphobHphil): Promise<MembraneCandidate> {
    const { minThickness } = ctx;
    const step = 0.3;
    let maxThickness = Vec3.distance(membrane.planePoint1, membrane.planePoint2);

    let maxNos = membraneSegments(ctx, membrane).length;
    let optimalThickness = membrane;

    let n = 0;
    const nl = Math.ceil((maxThickness - minThickness) / step);
    while (maxThickness > minThickness) {
        n++;
        if (runtime?.shouldUpdate && message && n % UPDATE_INTERVAL === 0) {
            await runtime.update({ message, current: n, max: nl });
        }

        const p = {
            ...ctx,
            maxThickness,
            stepSize: step
        };
        const temp = await findMembrane(runtime, void 0, p, [membrane.spherePoint!], initialHphobHphil);
        if (temp) {
            const nos = membraneSegments(ctx, temp).length;
            if (nos > maxNos) {
                maxNos = nos;
                optimalThickness = temp;
            }
        }
        maxThickness -= step;
    }

    return optimalThickness;
}

/** Report auth_seq_ids for all transmembrane segments. Will reject segments that are shorter than the adjust parameter specifies. Missing residues are considered in-membrane. */
function membraneSegments(ctx: ANVILContext, membrane: MembraneCandidate): ArrayLike<{ start: number, end: number }> {
    const { offsets, structure, adjust } = ctx;
    const { normalVector, planePoint1, planePoint2 } = membrane;
    const l = StructureElement.Location.create(structure);
    const testPoint = Vec3();
    const { auth_asym_id } = StructureProperties.chain;
    const { auth_seq_id } = StructureProperties.residue;

    const d1 = -Vec3.dot(normalVector!, planePoint1);
    const d2 = -Vec3.dot(normalVector!, planePoint2);
    const dMin = Math.min(d1, d2);
    const dMax = Math.max(d1, d2);

    const inMembrane: { [k: string]: Set<number> } = Object.create(null);
    const outMembrane: { [k: string]: Set<number> } = Object.create(null);
    const segments: Array<{ start: number, end: number }> = [];
    setLocation(l, structure, offsets[0]);
    let authAsymId;
    let lastAuthAsymId = null;
    let authSeqId;
    let lastAuthSeqId = auth_seq_id(l) - 1;
    let startOffset = 0;
    let endOffset = 0;

    // collect all residues in membrane layer
    for (let k = 0, kl = offsets.length; k < kl; k++) {
        setLocation(l, structure, offsets[k]);
        authAsymId = auth_asym_id(l);
        if (authAsymId !== lastAuthAsymId) {
            if (!inMembrane[authAsymId]) inMembrane[authAsymId] = new Set<number>();
            if (!outMembrane[authAsymId]) outMembrane[authAsymId] = new Set<number>();
            lastAuthAsymId = authAsymId;
        }

        authSeqId = auth_seq_id(l);
        Vec3.set(testPoint, l.unit.conformation.x(l.element), l.unit.conformation.y(l.element), l.unit.conformation.z(l.element));
        if (_isInMembranePlane(testPoint, normalVector!, dMin, dMax)) {
            inMembrane[authAsymId].add(authSeqId);
        } else {
            outMembrane[authAsymId].add(authSeqId);
        }
    }

    for (let k = 0, kl = offsets.length; k < kl; k++) {
        setLocation(l, structure, offsets[k]);
        authAsymId = auth_asym_id(l);
        authSeqId = auth_seq_id(l);
        if (inMembrane[authAsymId].has(authSeqId)) {
            // chain change
            if (authAsymId !== lastAuthAsymId) {
                segments.push({ start: startOffset, end: endOffset });
                lastAuthAsymId = authAsymId;
                startOffset = k;
                endOffset = k;
            }

            // sequence gaps
            if (authSeqId !== lastAuthSeqId + 1) {
                if (outMembrane[authAsymId].has(lastAuthSeqId + 1)) {
                    segments.push({ start: startOffset, end: endOffset });
                    startOffset = k;
                }
                lastAuthSeqId = authSeqId;
                endOffset = k;
            } else  {
                lastAuthSeqId++;
                endOffset++;
            }
        }
    }
    segments.push({ start: startOffset, end: endOffset });

    let startAuth;
    let endAuth;
    const refinedSegments: Array<{ start: number, end: number }> = [];
    for (let k = 0, kl = segments.length; k < kl; k++) {
        const { start, end } = segments[k];
        if (start === 0 || end === offsets.length - 1) continue;

        // evaluate residues 1 pos outside of membrane
        setLocation(l, structure, offsets[start - 1]);
        Vec3.set(testPoint, l.unit.conformation.x(l.element), l.unit.conformation.y(l.element), l.unit.conformation.z(l.element));
        const d3 = -Vec3.dot(normalVector!, testPoint);

        setLocation(l, structure, offsets[end + 1]);
        Vec3.set(testPoint, l.unit.conformation.x(l.element), l.unit.conformation.y(l.element), l.unit.conformation.z(l.element));
        const d4 = -Vec3.dot(normalVector!, testPoint);

        if (Math.min(d3, d4) < dMin && Math.max(d3, d4) > dMax) {
            // reject this refinement
            setLocation(l, structure, offsets[start]);
            startAuth = auth_seq_id(l);
            setLocation(l, structure, offsets[end]);
            endAuth = auth_seq_id(l);
            if (Math.abs(startAuth - endAuth) + 1 < adjust) {
                return [];
            }
            refinedSegments.push(segments[k]);
        }
    }

    return refinedSegments;
}

/** Filter for membrane residues and calculate the final extent of the membrane layer */
function adjustExtent(ctx: ANVILContext, membrane: MembraneCandidate, centroid: Vec3): number {
    const { offsets, structure } = ctx;
    const { normalVector, planePoint1, planePoint2 } = membrane;
    const l = StructureElement.Location.create(structure);
    const testPoint = Vec3();
    const { x, y, z } = StructureProperties.atom;

    const d1 = -Vec3.dot(normalVector!, planePoint1);
    const d2 = -Vec3.dot(normalVector!, planePoint2);
    const dMin = Math.min(d1, d2);
    const dMax = Math.max(d1, d2);
    let extent = 0;

    for (let k = 0, kl = offsets.length; k < kl; k++) {
        setLocation(l, structure, offsets[k]);
        Vec3.set(testPoint, x(l), y(l), z(l));
        if (_isInMembranePlane(testPoint, normalVector!, dMin, dMax)) {
            const dsq = Vec3.squaredDistance(testPoint, centroid);
            if (dsq > extent) extent = dsq;
        }
    }

    return Math.sqrt(extent);
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
    return _isInMembranePlane(testPoint, normalVector, Math.min(d1, d2), Math.max(d1, d2));
}

function _isInMembranePlane(testPoint: Vec3, normalVector: Vec3, min: number, max: number): boolean {
    const d = -Vec3.dot(normalVector, testPoint);
    return d > min && d < max;
}

/** Generates a defined number of points on a sphere with radius = extent around the specified centroid */
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

/** Generates sphere points close to that of the initial membrane */
function findProximateAxes(ctx: ANVILContext, membrane: MembraneCandidate): Vec3[] {
    const { numberOfSpherePoints, extent } = ctx;
    const points = generateSpherePoints(ctx, 30000);
    let j = 4;
    let sphere_pts2: Vec3[] = [];
    const s = 2 * extent / numberOfSpherePoints;
    while (sphere_pts2.length < numberOfSpherePoints) {
        const dsq = (s + j) * (s + j);
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
    hphil: number
}

namespace HphobHphil {
    export function of(hphob: number, hphil: number) {
        return {
            hphob: hphob,
            hphil: hphil
        };
    }

    const testPoint = Vec3();
    export function filtered(ctx: ANVILContext, filter?: { normalVector: Vec3, dMin: number, dMax: number }): HphobHphil {
        const { exposed, structure } = ctx;
        const { label_comp_id } = StructureProperties.atom;
        const l = StructureElement.Location.create(structure);
        let hphob = 0;
        let hphil = 0;
        for (let k = 0, kl = exposed.length; k < kl; k++) {
            setLocation(l, structure, exposed[k]);

            // testPoints have to be in putative membrane layer
            if (filter) {
                Vec3.set(testPoint, l.unit.conformation.x(l.element), l.unit.conformation.y(l.element), l.unit.conformation.z(l.element));
                if (!_isInMembranePlane(testPoint, filter.normalVector, filter.dMin, filter.dMax)) {
                    continue;
                }
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

/** ANVIL-specific (not general) definition of membrane-favoring amino acids */
const HYDROPHOBIC_AMINO_ACIDS = new Set(['ALA', 'CYS', 'GLY', 'HIS', 'ILE', 'LEU', 'MET', 'PHE', 'SER', 'TRP', 'VAL']);
/** Returns true if ANVIL considers this as amino acid that favors being embedded in a membrane */
export function isHydrophobic(label_comp_id: string): boolean {
    return HYDROPHOBIC_AMINO_ACIDS.has(label_comp_id);
}

/** Accessible surface area used for normalization. ANVIL uses 'Total-Side REL' values from NACCESS, from: Hubbard, S. J., & Thornton, J. M. (1993). naccess. Computer Program, Department of Biochemistry and Molecular Biology, University College London, 2(1). */
export const MaxAsa: { [k: string]: number } = {
    'ALA': 69.41,
    'ARG': 201.25,
    'ASN': 106.24,
    'ASP': 102.69,
    'CYS': 96.75,
    'GLU': 134.74,
    'GLN': 140.99,
    'GLY': 32.33,
    'HIS': 147.08,
    'ILE': 137.96,
    'LEU': 141.12,
    'LYS': 163.30,
    'MET': 156.64,
    'PHE': 164.11,
    'PRO': 119.90,
    'SER': 78.11,
    'THR': 101.70,
    'TRP': 211.26,
    'TYR': 177.38,
    'VAL': 114.28
};

function setLocation(l: StructureElement.Location, structure: Structure, serialIndex: number) {
    l.structure = structure;
    l.unit = structure.units[structure.serialMapping.unitIndices[serialIndex]];
    l.element = structure.serialMapping.elementIndices[serialIndex];
    return l;
}