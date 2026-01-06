import { Vec3 } from '../../mol-math/linear-algebra';

/**
 * Type representing alignment data for molecular structures.
 * It can either be explicit atom coordinates and atom types,
 * or a translation vector and rotation matrix.
 */
export type AlignmentData =
    | { x: number[]; y: number[]; z: number[]; atomType: string[] }
    | { centroidReference: Vec3; centroid: Vec3; rotMat: number[] };