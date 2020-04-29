/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SubstitutionMatrix, SubstitutionMatrices, SubstitutionMatrixData } from './substitution-matrix';

const DefaultAlignmentOptions = {
    gapPenalty: -11,
    gapExtensionPenalty: -1,
    substMatrix: 'blosum62' as SubstitutionMatrix
};
export type AlignmentOptions = typeof DefaultAlignmentOptions;

export function align(seq1: string, seq2: string, options: Partial<AlignmentOptions> = {}) {
    const o = { ...DefaultAlignmentOptions, ...options };
    const alignment = new Alignment(seq1, seq2, o);
    alignment.calculate();
    alignment.trace();
    return {
        ali1: alignment.ali1,
        ali2: alignment.ali2,
        score: alignment.score,
    };
}

class Alignment {
    gapPenalty: number; gapExtensionPenalty: number
    substMatrix: SubstitutionMatrixData

    n: number; m: number
    S: number[][]; V: number[][]; H: number[][]
    ali1: string; ali2: string;
    score: number

    constructor (readonly seq1: string, readonly seq2: string, options: AlignmentOptions) {
        this.gapPenalty = options.gapPenalty;
        this.gapExtensionPenalty = options.gapExtensionPenalty;
        this.substMatrix = SubstitutionMatrices[options.substMatrix];

        this.n = this.seq1.length;
        this.m = this.seq2.length;
    }

    private initMatrices () {
        this.S = [];
        this.V = [];
        this.H = [];

        for (let i = 0; i <= this.n; ++i) {
            this.S[i] = [];
            this.V[i] = [];
            this.H[i] = [];

            for (let j = 0; j <= this.m; ++j) {
                this.S[i][j] = 0;
                this.V[i][j] = 0;
                this.H[i][j] = 0;
            }
        }

        for (let i = 0; i <= this.n; ++i) {
            this.S[i][0] = this.gap(0);
            this.H[i][0] = -Infinity;
        }

        for (let j = 0; j <= this.m; ++j) {
            this.S[0][j] = this.gap(0);
            this.V[0][j] = -Infinity;
        }

        this.S[0][0] = 0;
    }

    private gap (len: number) {
        return this.gapPenalty + len * this.gapExtensionPenalty;
    }

    private makeScoreFn () {
        const seq1 = this.seq1;
        const seq2 = this.seq2;

        const substMatrix = this.substMatrix;

        if (substMatrix) {
            return function score (i: number, j: number) {
                const c1 = seq1[i];
                const c2 = seq2[j];

                try {
                    return substMatrix[c1][c2];
                } catch (e) {
                    return -4;
                }
            };
        } else {
            return function scoreNoSubstMat (i: number, j: number) {
                const c1 = seq1[i];
                const c2 = seq2[j];

                return c1 === c2 ? 5 : -3;
            };
        }
    }

    calculate () {
        this.initMatrices();

        const gap0 = this.gap(0);
        const scoreFn = this.makeScoreFn();
        const { V, H, S, n, m, gapExtensionPenalty } = this;

        let Vi1, Si1, Vi, Hi, Si;
        for (let i = 1; i <= n; ++i) {
            Si1 = S[i - 1], Vi1 = V[i - 1];
            Vi = V[i], Hi = H[i], Si = S[i];

            for (let j = 1; j <= m; ++j) {
                Vi[j] = Math.max(
                    Si1[j] + gap0,
                    Vi1[j] + gapExtensionPenalty
                );

                Hi[j] = Math.max(
                    Si[j - 1] + gap0,
                    Hi[j - 1] + gapExtensionPenalty
                );

                Si[j] = Math.max(
                    Si1[j - 1] + scoreFn(i - 1, j - 1), // match
                    Vi[j], // del
                    Hi[j]  // ins
                );
            }
        }
    }

    trace () {
        this.ali1 = '';
        this.ali2 = '';

        const scoreFn = this.makeScoreFn();

        let i = this.n;
        let j = this.m;
        let mat: 'S' | 'V' | 'H';

        if (this.S[i][j] >= this.V[i][j]) {
            mat = 'S';
            this.score = this.S[i][j];
        } else if (this.V[i][j] >= this.H[i][j]) {
            mat = 'V';
            this.score = this.V[i][j];
        } else {
            mat = 'H';
            this.score = this.H[i][j];
        }

        while (i > 0 && j > 0) {
            if (mat === 'S') {
                if (this.S[i][j] === this.S[i - 1][j - 1] + scoreFn(i - 1, j - 1)) {
                    this.ali1 = this.seq1[i - 1] + this.ali1;
                    this.ali2 = this.seq2[j - 1] + this.ali2;
                    --i;
                    --j;
                    mat = 'S';
                } else if (this.S[i][j] === this.V[i][j]) {
                    mat = 'V';
                } else if (this.S[i][j] === this.H[i][j]) {
                    mat = 'H';
                } else {
                    --i;
                    --j;
                }
            } else if (mat === 'V') {
                if (this.V[i][j] === this.V[i - 1][j] + this.gapExtensionPenalty) {
                    this.ali1 = this.seq1[i - 1] + this.ali1;
                    this.ali2 = '-' + this.ali2;
                    --i;
                    mat = 'V';
                } else if (this.V[i][j] === this.S[i - 1][j] + this.gap(0)) {
                    this.ali1 = this.seq1[i - 1] + this.ali1;
                    this.ali2 = '-' + this.ali2;
                    --i;
                    mat = 'S';
                } else {
                    --i;
                }
            } else if (mat === 'H') {
                if (this.H[i][j] === this.H[i][j - 1] + this.gapExtensionPenalty) {
                    this.ali1 = '-' + this.ali1;
                    this.ali2 = this.seq2[j - 1] + this.ali2;
                    --j;
                    mat = 'H';
                } else if (this.H[i][j] === this.S[i][j - 1] + this.gap(0)) {
                    this.ali1 = '-' + this.ali1;
                    this.ali2 = this.seq2[j - 1] + this.ali2;
                    --j;
                    mat = 'S';
                } else {
                    --j;
                }
            }
        }

        while (i > 0) {
            this.ali1 = this.seq1[i - 1] + this.ali1;
            this.ali2 = '-' + this.ali2;
            --i;
        }

        while (j > 0) {
            this.ali1 = '-' + this.ali1;
            this.ali2 = this.seq2[j - 1] + this.ali2;
            --j;
        }
    }
}