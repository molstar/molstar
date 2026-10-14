/**
 * @file estring.ts
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description E-strings (Edit strings) are operators that act on a string to
 * add insertion characters.
 * The follofing example is taken from MUSCLE (a) presentation:
 * Sequence: MQTIF
 * Estring: <3, -1, 2>
 * Transformation: MQT-IF
 * Positive integers are the counts of letters conserved in a fragment, negative
 * numbers are the counts of insertions.
 * Note that the sum of the positive numbers is the size of the initial string,
 * and the sum of the absolute values is the final size of the string.
 * The alignment of two sequences can be described by the operators that
 * transform the initial strings in 2 strings of same size.
 *
 * (a) Edgar, R.C. MUSCLE: a multiple sequence alignment method with reduced
 * time and space complexity. BMC Bioinformatics 5, 113 (2004).
 * https://doi.org/10.1186/1471-2105-5-113
 */

import { DEFAULT_GAP_CHAR } from "../align/params.ts";
import type { TSequence } from "../sequence/sequence.ts";

type TAlignedSeqFormat = {
  gapchar: string;
};

export function applyEstring(pSeqArr: TSequence[], pEstring: number[][], opt?: Partial<TAlignedSeqFormat>) {
  let lGapchar = opt?.gapchar ?? DEFAULT_GAP_CHAR;
  let lResp: string[] = [];
  for (let i = 0; i < pSeqArr.length; i++) {
    lResp.push(estringTransform(pSeqArr[i].rawSeq, pEstring[i], { gapchar: lGapchar }));
  }
  return lResp;
}

/**
 * Edits a string such as a sequence, by inserting gap symbols based on the
 * content of and Edit string.
 *
 * @export
 * @param {string} text
 * @param {number[]} estring
 * @returns {string}    edited string
 */
export function estringTransform(text: string, estring: number[], opt: TAlignedSeqFormat) {
  const lPieces: string[] = [];
  const GAPCHAR = opt.gapchar;
  let lPos = 0; // cursor in text
  for (let i = 0; i < estring.length; i++) {
    let n = estring[i];
    if (n < 0) {
      lPieces.push(GAPCHAR.repeat(-n));
    } else {
      lPieces.push(text.substr(lPos, n));
      lPos += n;
    }
  }
  return lPieces.join("");
}

/**
 * returns an edit string equivalent to the composition of the two edit strings
 * passed as an argument
 *
 * @export
 * @param {number[]} estringB
 * @param {number[]} estringA
 * @returns {number[]}
 */
export function estringProduct(estringB: number[], estringA: number[]): number[] {
  const lSteps: number[] = [0];
  let lPosA = 0;
  const lEstringA = estringA.slice(); // make copy, will be modified

  // add new value in lSteps. Collapses tail value if they have same sign
  function pushOperator(val: number) {
    if ((val ^ lSteps[lSteps.length - 1]) >= 0) {
      // same sign
      lSteps[lSteps.length - 1] += val;
    } else {
      // opposite sign: new item in array
      lSteps.push(val);
    }
  }

  for (let i = 0; i < estringB.length; i++) {
    let n = estringB[i];
    if (n < 0) {
      pushOperator(n);
    } else {
      while (n > 0 && lPosA < lEstringA.length) {
        let lValA = lEstringA[lPosA];
        if (lValA < 0) {
          if (n < -lValA) {
            lEstringA[lPosA] += n;
            pushOperator(-n);
            n = 0;
          } else {
            pushOperator(lValA);
            n += lValA;
            lPosA++;
          }
        } else {
          if (n < lValA) {
            lEstringA[lPosA] -= n;
            pushOperator(n);
            n = 0;
          } else {
            pushOperator(lValA);
            n -= lValA;
            lPosA++;
          }
        }
      }
    }
  }

  return lSteps;
}

/**
 * estringA and estringB are estrings that are applied independently to the
 * same sequence. For example a sequence X is aligned to a sequence A using the
 * estringA operator, and is aligned to a sequenceB using the estringB operator.
 * estringMerge computes an estring that adds to sequence X the necessary gaps
 * to align it both with sequence A and with sequence B. The main goal of this
 * computation is to avoid introducing gaps twice when they interrupt the same
 * segments of X while aligning with A and aligning with B.
 * X: MQTIF
 * A: MQQTIIF       B: MQVTIFE      A: MQQTIIF-
 * X: MQ-TI-F       X: MT-TIF-      X: MG-TI-F- <2,-1,2,-1,1,-1>
 *  <2,-1,2,-1,1>   <2,-1,3,-1>     B: MQVTI-FE
 *
 * @export
 * @param {number[]} estringA
 * @param {number[]} estringB
 * @returns
 */
export function estringMerge(estringA: number[], estringB: number[]) {
  const lEstring: number[] = [];

  let i = 0;
  let j = 0;
  let lValA = estringA[0];
  let lValB = estringB[0];

  while (i < estringA.length || j < estringB.length) {
    if (lValA === lValB) {
      lEstring.push(lValB);
      lValA = estringA[++i];
      lValB = estringB[++j];
      continue;
    }

    if (Math.sign(lValA) === Math.sign(lValB)) {
      if (lValA > 0) {
        // positive signs: add length of smallest segment (next
        // iteration will be gaps insertions)

        if (lValA > lValB) {
          lEstring.push(lValB);
          lValA -= lValB;
          lValB = estringB[++j];
        } else {
          lEstring.push(lValA);
          lValB -= lValA;
          lValA = estringA[++i];
        }
      } else {
        // gap insertions at the same position. Insert the maximum
        // number of gaps

        lEstring.push(Math.min(lValA, lValB));
        lValA = estringA[++i];
        lValB = estringB[++j];
      }

      continue;
    }

    // opposite sign: insert gaps
    if (lValA < 0 || lValB === undefined) {
      lEstring.push(lValA);
      lValA = estringA[++i];
      continue;
    }

    lEstring.push(lValB);
    lValB = estringB[++j];
  }

  return lEstring;
}

/**
 * Computes estringC that can be combined with estringA in estringProduct() to
 * get estringB.
 * The use case is to apply to sequence A, a transformation to align it with
 * the profile sequence B/sequence C, knowing the estring that leads to the
 * final shape of sequence B and the estring that aligns sequence B to sequence
 * A.
 * Important: we assume that estringB contains estringA. This entails that there
 * is no continuous segment in estringB that can be longer than its homologous
 * continuous segment in estringA (there is no gap between residues from
 * applying estringA that is not found when applying estringB)
 * X: MQTIF
 * A: MQQTIIF       B: MQVTIFE      A: MQQTIIF-
 * X: MQ-TI-F       X: MT-TIF-      X: MQ-TI-F- <2,-1,2,-1,1,-1>
 *  <2,-1,2,-1,1>   <2,-1,3,-1>     B: MQVTI-FE
 *
 * difference <2,-1,2,-1,1,-1> - <2,-1,3,-1> = <5,-1,2>
 * difference <2,-1,2,-1,1,-1> - <2,-1,2,-1,1> = <7,-1>
 * @export
 * @param {number[]} estringB
 * @param {number[]} estringA
 * @returns
 */
export function estringDifference(estringB: number[], estringA: number[]) {
  const lEstring: number[] = [];
  let i = 0;
  let j = 0;
  let lValA = estringA[0];
  let lValB = estringB[0];
  let lCurVal = 0;

  while (i < estringA.length || j < estringB.length) {
    lCurVal = lEstring[lEstring.length - 1];
    if (lValA === lValB) {
      // Common segment. Collapse with current segment or create new one.
      if (lCurVal > 0) {
        lEstring[lEstring.length - 1] = lCurVal + Math.abs(lValA);
      } else {
        lEstring.push(Math.abs(lValA));
      }
      lValA = estringA[++i];
      lValB = estringB[++j];
      continue;
    } else if (Math.sign(lValA) === Math.sign(lValB)) {
      if (lValA > 0) {
        // positive sign. lValB can only be smaller than lValA
        // (otherwise estringA contains gaps not accounted for in
        // estringB)
        if (lValB > lValA) return undefined; // error

        if (lCurVal > 0) {
          lEstring[lEstring.length - 1] = lCurVal + lValB;
        } else {
          lEstring.push(lValB);
        }
        lValA -= lValB;
        lValB = estringB[++j];
        continue;
      }

      // negative sign. we keep as much as possible from string A

      if (lValA > lValB) {
        // take all lValA
        if (lCurVal > 0) {
          lEstring[lEstring.length - 1] = lCurVal - lValA;
        } else {
          lEstring.push(Math.abs(lValA));
        }
        lValB -= lValA;
        lValA = estringA[++i];
        continue;
      } else {
        // this shan't happen: estringB can't introduce less gaps than
        // estringA does.
        return undefined;
      }
    }

    // opposite signs: if valB is negative, it introduces a gap that must
    // be kept.
    // if valB is positive, then estringA introduces a gap that is not in
    // estringB, which is not possible.
    else if (lValB < 0) {
      if (lCurVal < 0) {
        lEstring[lEstring.length - 1] = lCurVal + lValB;
      } else {
        lEstring.push(lValB);
      }
      lValB = estringB[++j];
      continue;
    } else {
      return undefined;
    }
  }

  return lEstring;
}

/**
 * Edit strings catenation.
 *
 * @export
 * @param {number[]} estringA
 * @param {number[]} estringB
 * @returns {number[]}
 */
export function estringCat(estringA: number[], estringB: number[]) {
  if (Math.sign(estringA[estringA.length - 1]) === Math.sign(estringB[0])) {
    let lRet = [...estringA];
    lRet[lRet.length - 1] += estringB[0];
    lRet.push(...estringB.slice(1)); // avoid changing reference
    return lRet;
  }
  return [...estringA, ...estringB];
}

// STORIES-UNTOLD(B1): const enum as an `as const` object and a type (erasable TypeScript).
export const EPATH_2_STRING = {
  NO_REVERSE: 1,
} as const;
export type EPATH_2_STRING = (typeof EPATH_2_STRING)[keyof typeof EPATH_2_STRING];

export function epath2estring(epath: number[], opt = 0): number[] {
  const lEstring: number[] = [];

  let lAcc = 0;
  let lVal = 0;
  for (let i = 0; i < epath.length; i++) {
    lVal = epath[i];
    if ((lVal ^ lAcc) >= 0) lAcc += lVal;
    else {
      lEstring.push(lAcc);
      lAcc = lVal;
    }
  }
  lEstring.push(lAcc);
  if (lEstring[0] === 0) lEstring.shift();

  // Usually, the path comes from a traceback matrix and must be reversed

  if (!(opt & EPATH_2_STRING.NO_REVERSE)) lEstring.reverse();

  return lEstring;
}

export function estringLength(estring: number[]): number {
  let lSize = 0;
  let lVal = 0;
  for (let i = 0; i < estring.length; i++) {
    lVal = estring[i];
    lSize += lVal < 0 ? -lVal : lVal;
  }
  return lSize;
}

export function estringCountPositive(estring: number[]) {
  let lCount = 0;
  let lVal = 0;
  for (let i = 0; i < estring.length; i++) {
    lVal = estring[i];
    if (lVal < 0) continue;
    lCount += lVal;
  }
  return lCount;
}

/**
 * Convert an Edit string to a vector of indices. Negative estring values (in-
 * dels) are encoded as -1. Positive values are encoded as an increasing sequence
 * @example
 * const estring = [2, -3, 2];
 * estringToIdx(estring); // => [0, 1, -1, -1, -1, 2, 3]
 * @param {number[]} estring
 * @returns {number[]}
 */
export function estringToIdx(estring: number[]) {
  const lIdx: number[] = [];
  let idx = 0;
  let val = 0;
  for (let i = 0; i < estring.length; i++) {
    val = estring[i];
    if (val < 0) {
      for (let j = 0; j < -val; j++) {
        lIdx.push(-1);
      }
      continue;
    }

    for (let j = 0; j < val; j++) {
      lIdx.push(idx);
      idx++;
    }
  }
  return lIdx;
}
