/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SubstitutionMatrix, SubstitutionMatrices, SubstitutionMatrixData } from './substitution-matrix';

const DefaultAlignmentOptions = {
    gapPenalty: -11,
    gapExtensionPenalty: -1,
    substMatrix: 'default' as SubstitutionMatrix | 'default'
};
export type AlignmentOptions = typeof DefaultAlignmentOptions;

export function align(seqA: ArrayLike<string>, seqB: ArrayLike<string>, options: Partial<AlignmentOptions> = {}) {
    const o = { ...DefaultAlignmentOptions, ...options };
    const alignment = new Alignment(seqA, seqB, o);
    alignment.calculate();
    return alignment.trace();
}

class Alignment {
    readonly gapPenalty: number; readonly gapExtensionPenalty: number
    readonly substMatrix: SubstitutionMatrixData | undefined

    readonly n: number; readonly m: number
    readonly S: number[][] = []; readonly V: number[][] = []; readonly H: number[][] = []

    constructor (readonly seqA: ArrayLike<string>, readonly seqB: ArrayLike<string>, options: AlignmentOptions) {
        this.gapPenalty = options.gapPenalty;
        this.gapExtensionPenalty = options.gapExtensionPenalty;
        this.substMatrix = options.substMatrix === 'default' ? undefined : SubstitutionMatrices[options.substMatrix];

        this.n = this.seqA.length;
        this.m = this.seqB.length;
    }

    private initMatrices () {
        const { n, m, gapPenalty, S, V, H } = this;

        for (let i = 0; i <= n; ++i) {
            S[i] = [], V[i] = [], H[i] = [];

            for (let j = 0; j <= m; ++j) {
                S[i][j] = 0, V[i][j] = 0, H[i][j] = 0;
            }
        }

        for (let i = 0; i <= n; ++i) {
            S[i][0] = gapPenalty;
            H[i][0] = -Infinity;
        }

        for (let j = 0; j <= m; ++j) {
            S[0][j] = gapPenalty;
            V[0][j] = -Infinity;
        }

        S[0][0] = 0;
    }

    private makeScoreFn () {
        const { seqA, seqB, substMatrix } = this;

        if (substMatrix) {
            return function score (i: number, j: number) {
                const cA = seqA[i];
                const cB = seqB[j];
                return substMatrix[cA]?.[cB] ?? -4;
            };
        } else {
            return function scoreNoSubstMat (i: number, j: number) {
                const cA = seqA[i];
                const cB = seqB[j];
                return cA === cB ? 5 : -3;
            };
        }
    }

    calculate () {
        this.initMatrices();

        const scoreFn = this.makeScoreFn();
        const { V, H, S, n, m, gapExtensionPenalty, gapPenalty } = this;

        let Vi1, Si1, Vi, Hi, Si;
        for (let i = 1; i <= n; ++i) {
            Si1 = S[i - 1], Vi1 = V[i - 1];
            Vi = V[i], Hi = H[i], Si = S[i];

            for (let j = 1; j <= m; ++j) {
                Vi[j] = Math.max(
                    Si1[j] + gapPenalty,
                    Vi1[j] + gapExtensionPenalty
                );

                Hi[j] = Math.max(
                    Si[j - 1] + gapPenalty,
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

    trace (): { aliA: ArrayLike<string>, aliB: ArrayLike<string>, score: number } {
        const scoreFn = this.makeScoreFn();
        const { V, H, S, seqA, seqB, gapExtensionPenalty, gapPenalty } = this;

        let i = this.n;
        let j = this.m;
        let mat: 'S' | 'V' | 'H';
        let score: number;

        let aliA = '';
        let aliB = '';

        if (S[i][j] >= V[i][j]) {
            mat = 'S';
            score = S[i][j];
        } else if (V[i][j] >= H[i][j]) {
            mat = 'V';
            score = V[i][j];
        } else {
            mat = 'H';
            score = H[i][j];
        }

        while (i > 0 && j > 0) {
            if (mat === 'S') {
                if (S[i][j] === S[i - 1][j - 1] + scoreFn(i - 1, j - 1)) {
                    aliA = seqA[i - 1] + aliA;
                    aliB = seqB[j - 1] + aliB;
                    --i;
                    --j;
                    mat = 'S';
                } else if (S[i][j] === V[i][j]) {
                    mat = 'V';
                } else if (S[i][j] === H[i][j]) {
                    mat = 'H';
                } else {
                    --i;
                    --j;
                }
            } else if (mat === 'V') {
                if (V[i][j] === V[i - 1][j] + gapExtensionPenalty) {
                    aliA = seqA[i - 1] + aliA;
                    aliB = '-' + aliB;
                    --i;
                    mat = 'V';
                } else if (V[i][j] === S[i - 1][j] + gapPenalty) {
                    aliA = seqA[i - 1] + aliA;
                    aliB = '-' + aliB;
                    --i;
                    mat = 'S';
                } else {
                    --i;
                }
            } else if (mat === 'H') {
                if (H[i][j] === H[i][j - 1] + gapExtensionPenalty) {
                    aliA = '-' + aliA;
                    aliB = seqB[j - 1] + aliB;
                    --j;
                    mat = 'H';
                } else if (H[i][j] === S[i][j - 1] + gapPenalty) {
                    aliA = '-' + aliA;
                    aliB = seqB[j - 1] + aliB;
                    --j;
                    mat = 'S';
                } else {
                    --j;
                }
            }
        }

        while (i > 0) {
            aliA = seqA[i - 1] + aliA;
            aliB = '-' + aliB;
            --i;
        }

        while (j > 0) {
            aliA = '-' + aliA;
            aliB = seqB[j - 1] + aliB;
            --j;
        }

        return { aliA, aliB, score };
    }
}