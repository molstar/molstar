/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Diego del Alamo <diego.delalamo@gmail.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * TM-align: Structure-based protein alignment algorithm
 *
 * References:
 * Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)
 * "TM-align: a protein structure alignment algorithm based on the TM-score"
 */

import { Mat4 } from './mat4';
import { Vec3 } from './vec3';
import { MinimizeRmsd } from './minimize-rmsd';
import type { UintArray } from '../../../mol-util/type-helpers';

export { TMAlign };

namespace TMAlign {
    export interface Result {
        /** Transformation matrix for structure B to superpose onto A */
        bTransform: Mat4;
        /** TM-score normalized by length of structure A */
        tmScoreA: number;
        /** TM-score normalized by length of structure B */
        tmScoreB: number;
        /** RMSD of aligned residues */
        rmsd: number;
        /** Number of aligned residue pairs */
        alignedLength: number;
        /** Sequence identity of aligned residues (if sequences provided) */
        sequenceIdentity: number;
        /** Alignment mapping: alignmentA[i] aligns with alignmentB[i] */
        alignmentA: number[];
        alignmentB: number[];
    }

    /** Reuse MinimizeRmsd.Positions type for consistency */
    export type Positions = MinimizeRmsd.Positions;
    export const Positions = MinimizeRmsd.Positions;

    export interface Input {
        /** Coordinates of structure A (reference) */
        a: Positions;
        /** Coordinates of structure B (mobile) */
        b: Positions;
        /** Optional: sequence of structure A for identity calculation */
        seqA?: string;
        /** Optional: sequence of structure B for identity calculation */
        seqB?: string;
    }

    /**
     * Compute TM-align between two structures
     */
    export function compute(input: Input): Result {
        const { a, b } = input;
        const lenA = a.x.length;
        const lenB = b.x.length;

        if (lenA === 0 || lenB === 0) {
            return {
                bTransform: Mat4.identity(),
                tmScoreA: 0,
                tmScoreB: 0,
                rmsd: 0,
                alignedLength: 0,
                sequenceIdentity: 0,
                alignmentA: [],
                alignmentB: []
            };
        }

        // Convert to 2D arrays for internal computation
        const xa = positionsToArray(a);
        const ya = positionsToArray(b);

        // Calculate d0 parameters for both normalizations
        const d0A = calculateD0(lenA);
        const d0B = calculateD0(lenB);
        const minLen = Math.min(lenA, lenB);
        const d0Search = Math.min(Math.max(Math.min(d0A, d0B), 4.5), 8);
        const scoreD8 = 1.5 * Math.pow(minLen, 0.3) + 3.5;

        // Run the TM-align algorithm
        const state = new TMAlignState(xa, ya, lenA, lenB, d0A, d0Search, scoreD8);
        state.align();

        // Get the final alignment and transformation
        const { alignmentA, alignmentB, transform } = state.getBestAlignment();

        // Calculate final scores
        const tmScoreA = calculateTMScore(xa, ya, alignmentA, alignmentB, transform, d0A, lenA);
        const tmScoreB = calculateTMScore(xa, ya, alignmentA, alignmentB, transform, d0B, lenB);
        const rmsd = calculateRMSD(xa, ya, alignmentA, alignmentB, transform);

        // Calculate sequence identity if sequences provided
        let sequenceIdentity = 0;
        if (input.seqA && input.seqB) {
            let identical = 0;
            for (let i = 0; i < alignmentA.length; i++) {
                if (input.seqA[alignmentA[i]] === input.seqB[alignmentB[i]]) {
                    identical++;
                }
            }
            sequenceIdentity = alignmentA.length > 0 ? identical / alignmentA.length : 0;
        }

        return {
            bTransform: transform,
            tmScoreA,
            tmScoreB,
            rmsd,
            alignedLength: alignmentA.length,
            sequenceIdentity,
            alignmentA,
            alignmentB
        };
    }

    /**
     * Calculate the d0 normalization parameter
     * d0 = 1.24 * (L - 15)^(1/3) - 1.8 for L > 21
     * d0 = 0.5 for L <= 21
     */
    export function calculateD0(length: number): number {
        if (length <= 21) return 0.5;
        const d0 = 1.24 * Math.pow(length - 15, 1 / 3) - 1.8;
        return Math.max(d0, 0.5);
    }

    /**
     * Calculate TM-score for a given alignment and transformation
     */
    function calculateTMScore(
        xa: Vec3[],
        ya: Vec3[],
        alignA: number[],
        alignB: number[],
        transform: Mat4,
        d0: number,
        normLength: number
    ): number {
        if (alignA.length === 0) return 0;

        const d02 = d0 * d0;
        let score = 0;
        const v = Vec3();

        for (let i = 0; i < alignA.length; i++) {
            const ai = alignA[i];
            const bi = alignB[i];

            Vec3.transformMat4(v, ya[bi], transform);
            const distSq = Vec3.squaredDistance(xa[ai], v);
            score += 1.0 / (1.0 + distSq / d02);
        }

        return score / normLength;
    }

    /**
     * Calculate RMSD for aligned pairs after transformation
     */
    function calculateRMSD(
        xa: Vec3[],
        ya: Vec3[],
        alignA: number[],
        alignB: number[],
        transform: Mat4
    ): number {
        if (alignA.length === 0) return 0;

        let sumSq = 0;
        const v = Vec3();

        for (let i = 0; i < alignA.length; i++) {
            Vec3.transformMat4(v, ya[alignB[i]], transform);
            sumSq += Vec3.squaredDistance(xa[alignA[i]], v);
        }

        return Math.sqrt(sumSq / alignA.length);
    }

    /** Convert Positions to 2D array for internal computation */
    function positionsToArray(pos: Positions): Vec3[] {
        const n = pos.x.length;
        const result: Vec3[] = new Array(n);
        for (let i = 0; i < n; i++) {
            result[i] = <any>[pos.x[i], pos.y[i], pos.z[i]];
        }
        return result;
    }
}

/**
 * Internal state class for TM-align computation
 */
class TMAlignState {
    private xa: Vec3[];
    private ya: Vec3[];
    private lenA: number;
    private lenB: number;
    private d0A: number;
    private d0Search: number;
    private scoreD8: number;

    private bestScore: number = -1;
    private bestAlignmentA: number[] = [];
    private bestAlignmentB: number[] = [];
    private bestTransform: Mat4 = Mat4.identity();

    // Working arrays for DP
    private dpPath: boolean[][];
    private dpVal: number[][];
    private j2i: number[];

    // Working pos vectors
    private posA: MinimizeRmsd.Positions;
    private posB: MinimizeRmsd.Positions;

    // Transformed coordinates of structure B
    private yt: Vec3[];

    // Working arrays for alignment indices
    private tmpAlignA: Uint16Array;
    private tmpAlignB: Uint16Array;

    constructor(
        xa: number[][],
        ya: number[][],
        lenA: number,
        lenB: number,
        d0A: number,
        d0Search: number,
        scoreD8: number
    ) {
        this.xa = <any> xa;
        this.ya = <any> ya;
        this.lenA = lenA;
        this.lenB = lenB;
        this.d0A = d0A;
        this.d0Search = d0Search;
        this.scoreD8 = scoreD8;

        // Initialize DP arrays
        this.dpPath = new Array(lenA + 1);
        this.dpVal = new Array(lenA + 1);
        for (let i = 0; i <= lenA; i++) {
            this.dpPath[i] = new Array(lenB + 1).fill(false);
            this.dpVal[i] = new Array(lenB + 1).fill(0.0);
        }
        this.j2i = new Array(lenB).fill(-1);

        this.posA = MinimizeRmsd.Positions.empty(lenA);
        this.posB = MinimizeRmsd.Positions.empty(lenB);

        this.yt = Array.from({ length: lenB }, () => [0.1, 0.0, 0.0] as unknown as Vec3);

        this.tmpAlignA = new Uint16Array(lenA);
        this.tmpAlignB = new Uint16Array(lenB);
    }

    /**
     * Main alignment procedure
     */
    align(): void {
        const { lenA, lenB } = this;
        const minLen = Math.min(lenA, lenB);

        // Strategy 1: Initial global DP alignment
        this.getInitialAlignment();
        this.refineAlignment();

        // Strategy 2: Gapless threading with various offsets
        this.tryGaplessThreading();

        // Strategy 3: Fragment-based initialization
        // Try small fragments (5-12) with medium search (thorough was too slow)
        for (let fragLen = 5; fragLen <= Math.min(12, minLen); fragLen++) {
            this.tryFragmentInitializationMedium(fragLen);
        }

        // Try medium fragments (16-30) with medium search
        for (let fragLen = 16; fragLen <= Math.min(30, minLen); fragLen += 2) {
            this.tryFragmentInitializationMedium(fragLen);
        }

        // Try larger fragments with strategic sizes
        const fragLengths = [
            Math.floor(minLen * 0.15),
            Math.floor(minLen * 0.25),
            Math.floor(minLen * 0.4),
            Math.floor(minLen * 0.6)
        ].filter(f => f > 30);

        for (const fragLen of fragLengths) {
            this.tryFragmentInitializationFast(fragLen);
        }

        // Final refinement passes with different cutoffs
        this.refineAlignment();

        // TM-align uses multiple refinement passes with different d0 values
        for (let d = Math.max(1, this.d0A - 3); d <= this.d0A + 3; d++) {
            this.refineWithCutoff(d);
        }

        // Final pass with tighter cutoff
        this.refineWithCutoff(this.d0A * 0.8);

        // TM-align performs a final optimization to maximize TM-score
        // by iteratively removing poorly aligned pairs
        this.optimizeFinalAlignment();
    }

    /**
     * Optimize final alignment by removing poorly aligned pairs
     * to maximize TM-score
     */
    private optimizeFinalAlignment(): void {
        const { xa, ya, d0A, lenA, tmpAlignA: filteredA, tmpAlignB: filteredB } = this;

        // Try different distance cutoffs to find optimal alignment
        const cutoffs = [8.0, 7.0, 6.0, 5.5, 5.0, 4.5, 4.0];
        const v = Vec3();
        let n = 0;

        for (const cutoff of cutoffs) {
            const cutoffSq = cutoff * cutoff;

            // Filter pairs by distance
            n = 0;

            for (let i = 0; i < this.bestAlignmentA.length; i++) {
                const ai = this.bestAlignmentA[i];
                const bi = this.bestAlignmentB[i];

                Vec3.transformMat4(v, ya[bi], this.bestTransform);
                const distSq = Vec3.squaredDistance(v, xa[ai]);

                if (distSq <= cutoffSq) {
                    filteredA[n] = ai;
                    filteredB[n] = bi;
                    n++;
                }
            }

            if (n < 3) break;

            // Compute Kabsch on filtered pairs
            const transform = this.kabsch(filteredA.subarray(0, n), filteredB.subarray(0, n));

            // Score using FULL alignment (not filtered) with new transform
            const score = this.scoreTM(
                this.bestAlignmentA,
                this.bestAlignmentB,
                transform,
                d0A,
                lenA
            );

            if (score > this.bestScore) {
                this.bestScore = score;
                Mat4.copy(this.bestTransform, transform);
            }
        }

        // Recompute final score
        this.bestScore = this.scoreTM(
            this.bestAlignmentA,
            this.bestAlignmentB,
            this.bestTransform,
            d0A,
            lenA
        );
    }

    /**
     * Calculate TM-score without distance cutoff (for final result)
     */
    private scoreTM(
        alignA: number[],
        alignB: number[],
        transform: Mat4,
        d0: number,
        normLen: number
    ): number {
        const { xa, ya } = this;
        const d02 = d0 * d0;
        let score = 0;
        const v = Vec3();

        for (let i = 0; i < alignA.length; i++) {
            const ai = alignA[i];
            const bi = alignB[i];

            Vec3.transformMat4(v, ya[bi], transform);

            const distSq = Vec3.squaredDistance(v, xa[ai]);

            score += 1.0 / (1.0 + distSq / d02);
        }

        return score / normLen;
    }

    /**
     * Refine with a specific distance cutoff
     */
    private refineWithCutoff(cutoff: number): void {
        if (cutoff < 1) return;
        const { xa, ya, yt, lenA, d0A, d0Search, scoreD8, tmpAlignA, tmpAlignB } = this;
        const cutoffSq = cutoff * cutoff;
        const v = Vec3();
        let n = 0;

        for (let iter = 0; iter < 10; iter++) {
            if (this.bestAlignmentA.length < 3) break;

            // Use trimmed Kabsch with this cutoff
            n = 0;

            for (let i = 0; i < this.bestAlignmentA.length; i++) {
                const ai = this.bestAlignmentA[i];
                const bi = this.bestAlignmentB[i];

                Vec3.transformMat4(v, ya[bi], this.bestTransform);

                const distSq = Vec3.squaredDistance(v, xa[ai]);
                if (distSq <= cutoffSq) {
                    tmpAlignA[n] = ai;
                    tmpAlignB[n] = bi;
                    n++;
                }
            }

            if (n < 3) break;

            const trimmedA = tmpAlignA.subarray(0, n);
            const trimmedB = tmpAlignB.subarray(0, n);
            const transform = this.kabsch(trimmedA, trimmedB);

            // Transform ya
            for (let i = 0; i < this.lenB; i++) {
                Vec3.transformMat4(yt[i], ya[i], transform);
            }

            // Run DP with transformed coordinates
            const d = d0Search + 1;
            this.nwdpStructure(xa, yt, d * d, -0.6);

            // Extract new alignment
            n = 0;
            for (let j = 0; j < this.lenB; j++) {
                if (this.j2i[j] >= 0) {
                    tmpAlignA[n] = this.j2i[j];
                    tmpAlignB[n] = j;
                    n++;
                }
            }

            if (n < 3) break;

            const newAlignA = tmpAlignA.subarray(0, n);
            const newAlignB = tmpAlignB.subarray(0, n);

            // Compute final Kabsch and score
            const newTransform = this.trimmedKabsch(newAlignA, newAlignB, transform, scoreD8);
            const newScore = this.scoreTMWithCutoff(newAlignA, newAlignB, newTransform, d0A, lenA, scoreD8);

            if (newScore > this.bestScore) {
                this.bestScore = newScore;
                this.bestAlignmentA = Array.from(newAlignA);
                this.bestAlignmentB = Array.from(newAlignB);
                this.bestTransform = newTransform;
            } else {
                break;
            }
        }
    }

    /**
     * Try gapless threading - align residue i of A with residue i+offset of B
     * This is O(n) per offset and provides good initial seeds
     */
    private tryGaplessThreading(): void {
        const { lenA, lenB, d0A, scoreD8, tmpAlignA, tmpAlignB } = this;
        let n = 0;

        // Try various offsets
        for (let offset = -lenB + 4; offset <= lenA - 4; offset++) {
            n = 0;
            // TODO: use an increment array and compute the boundaries depending on lenA, lenB and offset

            // Build gapless alignment with this offset
            for (let i = 0; i < lenA; i++) {
                const j = i - offset;
                if (j >= 0 && j < lenB) {
                    tmpAlignA[n] = i;
                    tmpAlignB[n] = j;
                    n++;
                }
            }

            if (n < 4) continue;

            const alignA = tmpAlignA.subarray(0, n);
            const alignB = tmpAlignB.subarray(0, n);

            // Compute Kabsch
            const transform = this.kabsch(alignA, alignB);

            // Score
            const score = this.scoreTMWithCutoff(alignA, alignB, transform, d0A, lenA, scoreD8);

            if (score > this.bestScore * 0.8) {
                // Promising seed - refine it
                this.extendFromSeed(transform);
            }
        }
    }

    /**
     * Medium-thoroughness fragment-based initialization
     */
    private tryFragmentInitializationMedium(fragLen: number): void {
        const { lenA, lenB, d0A, scoreD8, tmpAlignA, tmpAlignB } = this;
        const maxStartA = lenA - fragLen;
        const maxStartB = lenB - fragLen;
        // Use fragLen/2 as step for more thorough search
        const step = Math.max(2, Math.floor(fragLen / 2));

        const incrSequence = lenA < lenB ? tmpAlignB : tmpAlignA;
        for (let i = 0; i < incrSequence.length; i++) {
            incrSequence[i] = i;
        }

        for (let startA = 0; startA <= maxStartA; startA += step) {
            const fragA = incrSequence.subarray(startA, startA + fragLen);

            for (let startB = 0; startB <= maxStartB; startB += step) {
                const fragB = incrSequence.subarray(startB, startB + fragLen);

                // Superpose fragments
                const transform = this.kabsch(fragA, fragB);

                // Quick score check
                const quickScore = this.scoreTMWithCutoff(fragA, fragB, transform, d0A, lenA, scoreD8);

                if (quickScore > this.bestScore * 0.3) {
                    this.extendFromSeed(transform);
                }
            }
        }
    }

    /**
     * Fast fragment-based initialization with large steps
     */
    private tryFragmentInitializationFast(fragLen: number): void {
        const { lenA, lenB, d0A, scoreD8, tmpAlignA, tmpAlignB } = this;
        const maxStartA = lenA - fragLen;
        const maxStartB = lenB - fragLen;
        // Use fragment length as step - only try diagonal and near-diagonal positions
        const step = Math.max(fragLen, 10);

        const incrSequence = lenA < lenB ? tmpAlignB : tmpAlignA;
        for (let i = 0; i < incrSequence.length; i++) {
            incrSequence[i] = i;
        }

        for (let startA = 0; startA <= maxStartA; startA += step) {
            const fragA = incrSequence.subarray(startA, startA + fragLen);

            for (let startB = 0; startB <= maxStartB; startB += step) {
                const fragB = incrSequence.subarray(startB, startB + fragLen);

                // Superpose fragments
                const transform = this.kabsch(fragA, fragB);

                // Quick score check before full refinement
                const quickScore = this.scoreTMWithCutoff(fragA, fragB, transform, d0A, lenA, scoreD8);

                if (quickScore > this.bestScore * 0.5) {
                    // Promising seed - refine it
                    this.extendFromSeed(transform);
                }
            }
        }
    }

    /**
     * Get initial alignment using length-independent approach
     */
    private getInitialAlignment(): void {
        const { xa, ya, lenB, d0Search, j2i } = this;
        const d02 = d0Search * d0Search;

        this.nwdpStructure(xa, ya, d02, -0.6);

        // Extract alignment from j2i
        this.bestAlignmentA = [];
        this.bestAlignmentB = [];
        for (let jj = 0; jj < lenB; jj++) {
            if (j2i[jj] >= 0) {
                this.bestAlignmentA.push(j2i[jj]);
                this.bestAlignmentB.push(jj);
            }
        }

        if (this.bestAlignmentA.length >= 3) {
            this.bestTransform = this.kabsch(this.bestAlignmentA, this.bestAlignmentB);
            this.bestScore = this.scoreTMWithCutoff(
                this.bestAlignmentA,
                this.bestAlignmentB,
                this.bestTransform,
                this.d0A,
                this.lenA,
                this.scoreD8
            );
        }
    }

    /**
     * Extend alignment from a seed transformation with iterative refinement
     */
    private extendFromSeed(initialTransform: Mat4): void {
        const { xa, ya, yt, lenA, lenB, d0A, d0Search, scoreD8, tmpAlignA, tmpAlignB } = this;

        let currentTransform = initialTransform;
        let currentAlignA: number[] = [];
        let currentAlignB: number[] = [];
        let prevScore = -1;
        let n = 0;

        // Iteratively refine this seed
        for (let iter = 0; iter < 20; iter++) {
            // Transform ya by current transform
            for (let i = 0; i < lenB; i++) {
                Vec3.transformMat4(yt[i], ya[i], currentTransform);
            }

            // Run structure-based DP
            const d = d0Search + 1;
            this.nwdpStructure(xa, yt, d * d, -0.6);

            // Extract alignment from j2i
            n = 0;
            for (let j = 0; j < lenB; j++) {
                if (this.j2i[j] >= 0) {
                    tmpAlignA[n] = this.j2i[j];
                    tmpAlignB[n] = j;
                    n++;
                }
            }

            if (n < 3) break;
            const alignA = tmpAlignA.subarray(0, n);
            const alignB = tmpAlignB.subarray(0, n);

            // Check convergence
            if (this.alignmentsEqual(alignA, alignB, currentAlignA, currentAlignB)) {
                break;
            }

            // Use trimmed Kabsch - only use well-aligned pairs for superposition
            const transform = this.trimmedKabsch(alignA, alignB, currentTransform, scoreD8);

            // Score this alignment
            const score = this.scoreTMWithCutoff(alignA, alignB, transform, d0A, lenA, scoreD8);

            // Check if score improved
            if (score <= prevScore) break;
            prevScore = score;

            // Update current state
            currentAlignA = Array.from(alignA);
            currentAlignB = Array.from(alignB);
            currentTransform = transform;

            // Update global best if this is better
            if (score > this.bestScore) {
                this.bestScore = score;
                this.bestAlignmentA = currentAlignA;
                this.bestAlignmentB = currentAlignB;
                Mat4.copy(this.bestTransform, transform);
            }
        }
    }

    /**
     * Kabsch using only well-aligned pairs (distance < cutoff)
     */
    private trimmedKabsch(alignA: UintArray, alignB: UintArray, currentTransform: Mat4, cutoff: number): Mat4 {
        const { xa, ya } = this;
        const cutoffSq = cutoff * cutoff;

        // Find well-aligned pairs
        const trimmedAlignA: number[] = [];
        const trimmedAlignB: number[] = [];
        const v = Vec3();

        for (let i = 0; i < alignA.length; i++) {
            const ai = alignA[i];
            const bi = alignB[i];

            Vec3.transformMat4(v, ya[bi], currentTransform);
            const distSq = Vec3.squaredDistance(v, xa[ai]);

            if (distSq <= cutoffSq) {
                trimmedAlignA.push(ai);
                trimmedAlignB.push(bi);
            }
        }

        // Need at least 3 pairs for Kabsch
        if (trimmedAlignA.length < 3) {
            // Fall back to using all pairs
            return this.kabsch(alignA, alignB);
        }

        return this.kabsch(trimmedAlignA, trimmedAlignB);
    }

    /**
     * Refine the current best alignment iteratively
     */
    private refineAlignment(): void {
        const { xa, ya, yt, lenA, lenB, d0A, d0Search, scoreD8 } = this;
        const maxIterations = 20;

        for (let iter = 0; iter < maxIterations; iter++) {
            if (this.bestAlignmentA.length < 3) break;

            // Use trimmed Kabsch for refinement
            const transform = this.trimmedKabsch(
                this.bestAlignmentA,
                this.bestAlignmentB,
                this.bestTransform,
                scoreD8
            );

            // Transform ya
            for (let i = 0; i < lenB; i++) {
                Vec3.transformMat4(yt[i], ya[i], transform);
            }

            // Run DP with transformed coordinates
            const d = d0Search + 1;
            this.nwdpStructure(xa, yt, d * d, -0.6);

            // Extract new alignment
            const newAlignA: number[] = [];
            const newAlignB: number[] = [];
            for (let j = 0; j < this.lenB; j++) {
                if (this.j2i[j] >= 0) {
                    newAlignA.push(this.j2i[j]);
                    newAlignB.push(j);
                }
            }

            if (newAlignA.length < 3) break;

            // Check convergence
            if (this.alignmentsEqual(newAlignA, newAlignB, this.bestAlignmentA, this.bestAlignmentB)) {
                break;
            }

            // Compute new Kabsch and score
            const newTransform = this.kabsch(newAlignA, newAlignB);
            const newScore = this.scoreTMWithCutoff(newAlignA, newAlignB, newTransform, d0A, lenA, scoreD8);

            if (newScore > this.bestScore) {
                this.bestScore = newScore;
                this.bestAlignmentA = newAlignA;
                this.bestAlignmentB = newAlignB;
                Mat4.copy(this.bestTransform, newTransform);
            } else {
                break;
            }
        }
    }

    /**
     * Check if two alignments are equal
     */
    private alignmentsEqual(a1: UintArray, b1: UintArray, a2: UintArray, b2: UintArray): boolean {
        if (a1.length !== a2.length) return false;
        for (let i = 0; i < a1.length; i++) {
            if (a1[i] !== a2[i] || b1[i] !== b2[i]) return false;
        }
        return true;
    }

    /**
     * Calculate TM-score with distance cutoff (score_d8)
     */
    private scoreTMWithCutoff(
        alignA: UintArray,
        alignB: UintArray,
        transform: Mat4,
        d0: number,
        normLen: number,
        scoreD8: number
    ): number {
        const { xa, ya } = this;
        const d02 = d0 * d0;
        const scoreD8Sq = scoreD8 * scoreD8;
        let score = 0;
        const v = Vec3();

        for (let i = 0; i < alignA.length; i++) {
            const ai = alignA[i];
            const bi = alignB[i];

            Vec3.transformMat4(v, ya[bi], transform);
            const distSq = Vec3.squaredDistance(v, xa[ai]);

            if (distSq <= scoreD8Sq) {
                score += 1.0 / (1.0 + distSq / d02);
            }
        }

        return score / normLen;
    }

    /**
     * Structure-based DP (coordinates already transformed)
     */
    private nwdpStructure(xa: Vec3[], yt: Vec3[], d02: number, gapOpen: number): void {
        const { lenA, lenB, dpPath, dpVal, j2i } = this;

        // Initialize
        for (let i = 0; i <= lenA; i++) {
            dpVal[i][0] = 0.0;
            dpPath[i][0] = false;
        }
        for (let j = 0; j <= lenB; j++) {
            dpVal[0][j] = 0.0;
            dpPath[0][j] = false;
            j2i[j] = -1;
        }

        // Fill DP matrix
        for (let i = 1; i <= lenA; i++) {
            for (let j = 1; j <= lenB; j++) {
                const distSq = Vec3.squaredDistance(xa[i - 1], yt[j - 1]);

                const d = dpVal[i - 1][j - 1] + 1.0 / (1.0 + distSq / d02);
                let h = dpVal[i - 1][j];
                if (dpPath[i - 1][j]) h += gapOpen;
                let v = dpVal[i][j - 1];
                if (dpPath[i][j - 1]) v += gapOpen;

                if (d >= h && d >= v) {
                    dpPath[i][j] = true;
                    dpVal[i][j] = d;
                } else {
                    dpPath[i][j] = false;
                    dpVal[i][j] = v >= h ? v : h;
                }
            }
        }

        // Traceback
        let i = lenA;
        let j = lenB;
        while (i > 0 && j > 0) {
            if (dpPath[i][j]) {
                j2i[j - 1] = i - 1;
                i--;
                j--;
            } else {
                let h = dpVal[i - 1][j];
                if (dpPath[i - 1][j]) h += gapOpen;
                let v = dpVal[i][j - 1];
                if (dpPath[i][j - 1]) v += gapOpen;
                if (v >= h) j--;
                else i--;
            }
        }
    }

    /**
     * Kabsch superposition using MinimizeRmsd
     */
    private kabsch(alignA: UintArray, alignB: UintArray): Mat4 {
        const n = alignA.length;
        if (n < 3) return Mat4.identity();

        const { xa, ya, posA, posB } = this;

        let xai = xa[0];
        let ybi = ya[0];

        for (let i = 0; i < n; i++) {
            xai = xa[alignA[i]];
            ybi = ya[alignB[i]];
            (<any>posA.x)[i] = xai[0];
            (<any>posA.y)[i] = xai[1];
            (<any>posA.z)[i] = xai[2];
            (<any>posB.x)[i] = ybi[0];
            (<any>posB.y)[i] = ybi[1];
            (<any>posB.z)[i] = ybi[2];
        }

        const result = MinimizeRmsd.compute({ a: posA, b: posB, length: n });
        return result.bTransform;
    }

    /**
     * Get the best alignment result
     */
    getBestAlignment(): { alignmentA: number[]; alignmentB: number[]; transform: Mat4 } {
        return {
            alignmentA: this.bestAlignmentA,
            alignmentB: this.bestAlignmentB,
            transform: this.bestTransform
        };
    }
}
