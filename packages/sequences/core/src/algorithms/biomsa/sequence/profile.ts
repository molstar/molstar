/**
 * @file profile
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020
 */

import { ALIGNOPT } from "../align/align.ts";
import { SEQUENCE_TYPE, type TAlignmentParam } from "../align/params.ts";
import BitArray from "../utils/bitarray.ts";
import { estringLength, estringToIdx } from "../utils/estring.ts";
import { aaToNum, nucToNum, type TSequence } from "./sequence.ts";

/*
 * ProfPos definition ()
 */
export class ProfPos {
  m_bAllGaps: BitArray;
  /** Residues present at this position unordered */
  m_uSortOrder: Uint8Array;

  /** Residues counts */
  m_fcCounts: Uint8Array;

  /** Weighted contribution per residue */
  m_wCounts: Float32Array;
  // m_LL = 1;
  // m_LG = 0;
  // m_GL = 0;
  // m_GG = 0;

  /** Substitution score computed at this position */
  m_AAScores: Float32Array;

  /** Residue diversity (count per residue) */
  m_uResidueGroup: Uint8Array;

  /** Occupancy (non-indels) */
  m_fOcc: Float32Array;

  /** Frequency of opening gaps */
  m_fcStartOcc: Float32Array;

  /** Frequency of closing gaps */
  m_fcEndOcc: Float32Array;

  /** Gap Open penalty computed for this position */
  m_ScoreGapOpen: Float32Array;

  /** Gap Close penalty computed for this position */
  m_ScoreGapClose: Float32Array;

  /** Sequence count in this profile */
  nbSeq = 0;
  alphaSize: number;
  length = 0;
  weight = 0;

  constructor(pSize: number, pAlphaSize: number, pNbSeq: number, pWeight = 0) {
    this.length = pSize;
    this.m_bAllGaps = new BitArray(pSize);
    this.m_uSortOrder = new Uint8Array(pAlphaSize * pSize);
    this.m_fcCounts = new Uint8Array(pAlphaSize * pSize);
    this.m_wCounts = new Float32Array(pAlphaSize * pSize);
    this.m_AAScores = new Float32Array(pAlphaSize * pSize);
    this.m_uResidueGroup = new Uint8Array(pSize);
    this.m_fOcc = new Float32Array(pSize);
    this.m_fcStartOcc = new Float32Array(pSize);
    this.m_fcEndOcc = new Float32Array(pSize);
    this.m_ScoreGapOpen = new Float32Array(pSize);
    this.m_ScoreGapClose = new Float32Array(pSize);
    this.nbSeq = pNbSeq;
    this.alphaSize = pAlphaSize;
    this.weight = pWeight;
  }

  getProxy(id = 0) {
    return new ProfProxy(this, id);
  }
}

class ProfProxy {
  prof: ProfPos;
  idx: number = 0;
  alphaSize = 4;
  offset = 0;

  constructor(pProf: ProfPos, pIdx = 0) {
    this.prof = pProf;
    this.idx = pIdx;
    this.alphaSize = pProf.alphaSize;
    this.offset = this.idx * this.alphaSize;
  }

  setProxy(pIdx: number) {
    this.idx = pIdx;
    this.offset = pIdx * this.alphaSize;
  }

  get m_bAllGaps() {
    return this.prof.m_bAllGaps.get(this.idx);
  }
  set m_bAllGaps(val: boolean) {
    if (val) this.prof.m_bAllGaps.set(this.idx);
    else this.prof.m_bAllGaps.clear(this.idx);
  }

  get m_uSortOrder() {
    return this.prof.m_uSortOrder.subarray(this.offset, this.offset + this.alphaSize);
  }
  set m_uSortOrder(val: Uint8Array) {
    this.prof.m_uSortOrder.set(val, this.offset);
  }

  get m_uResidueGroup() {
    return this.prof.m_uResidueGroup[this.idx];
  }
  set m_uResidueGroup(val: number) {
    this.prof.m_uResidueGroup[this.idx] = val;
  }

  get m_fcCounts() {
    return this.prof.m_fcCounts.subarray(this.offset, this.offset + this.alphaSize);
  }
  set m_fcCounts(val: Uint8Array) {
    this.prof.m_fcCounts.set(val, this.offset);
  }

  get m_wCounts() {
    return this.prof.m_wCounts.subarray(this.offset, this.offset + this.alphaSize);
  }
  set m_wCounts(val: Float32Array) {
    this.prof.m_wCounts.set(val, this.offset);
  }

  get m_AAScores() {
    return this.prof.m_AAScores.subarray(this.offset, this.offset + this.alphaSize);
  }
  set m_AAScores(val: Float32Array) {
    this.prof.m_AAScores.set(val, this.offset);
  }

  get m_fOcc() {
    return this.prof.m_fOcc[this.idx];
  }
  set m_fOcc(val: number) {
    this.prof.m_fOcc[this.idx] = val;
  }

  get m_fcStartOcc() {
    return this.prof.m_fcStartOcc[this.idx];
  }
  set m_fcStartOcc(val: number) {
    this.prof.m_fcStartOcc[this.idx] = val;
  }

  get m_fcEndOcc() {
    return this.prof.m_fcEndOcc[this.idx];
  }
  set m_fcEndOcc(val: number) {
    this.prof.m_fcEndOcc[this.idx] = val;
  }

  get m_ScoreGapOpen() {
    return this.prof.m_ScoreGapOpen[this.idx];
  }
  set m_ScoreGapOpen(val: number) {
    this.prof.m_ScoreGapOpen[this.idx] = val;
  }

  get m_ScoreGapClose() {
    return this.prof.m_ScoreGapClose[this.idx];
  }
  set m_ScoreGapClose(val: number) {
    this.prof.m_ScoreGapClose[this.idx] = val;
  }
}

/**
 * Compute Profile from MSA
 * A profile stores at each column of the alignment a statistical representation
 * of the alignment properties suitable for computing alignment scores against
 * a sequence or another profile.
 * The data structure modeling the profile is a store containing numerical
 * values for the different properties (residue counts, etc...)
 * @param {string[]} pMSA: aligned sequences as an array of strings with `-` as
 *                          indels.
 * @param {number}  pGapO: gap opening penalty
 * @param {number[]} pWeights : weight of each sequence
 */
export function profileFromMSA(pMSA: string[], pWeights: number[], params: TAlignmentParam, opt: number = 0) {
  /** alignment length (columns count) */
  const l = pMSA[0].length;

  let aa = "", // amino acid name (1 letter code)
    na = 0, // amino acid # (in lexical order)
    tabCodeResNum: Uint8Array;
  if (params.type === SEQUENCE_TYPE.NUCLEIC) {
    tabCodeResNum = nucToNum;
  } else {
    tabCodeResNum = aaToNum;
  }

  //uHydrophobicRunLength = 0, // hydrophobic window (to avoid introducing gaps within)
  let w = 0, //sequence weight
    totalW = pWeights.reduce((a, b) => a + b, 0);

  const lProf = new ProfPos(l, params.abSize, pMSA.length, totalW);

  const lGapO = params.gapOP;
  let lProx = lProf.getProxy();

  for (let col = 0; col < l; col++) {
    lProx.setProxy(col);
    let fGap = 0;

    for (let numSeq = 0, max = pMSA.length; numSeq < max; numSeq++) {
      w = pWeights[numSeq];
      aa = pMSA[numSeq][col];
      if (aa === "-") {
        fGap += w;

        if (col === 0 || pMSA[numSeq][col - 1] !== "-") {
          lProx.m_fcStartOcc += w;
        }

        if (col === l - 1 || pMSA[numSeq][col + 1] !== "-") {
          lProx.m_fcEndOcc += w;
        }
      } else {
        na = tabCodeResNum[aa.charCodeAt(0)];
        lProx.m_fcCounts[na]++;
        lProx.m_wCounts[na] += w;
        if (lProx.m_fcCounts[na] === 1) {
          lProx.m_uSortOrder[lProx.m_uResidueGroup] = na;
          lProx.m_uResidueGroup++;
        }
      }
    }
    lProx.m_fOcc = totalW - fGap;
    lProx.m_ScoreGapOpen = (lGapO / 2) * (1 - lProx.m_fcStartOcc);
    lProx.m_ScoreGapClose = (lGapO / 2) * (1 - lProx.m_fcEndOcc);

    // Cache substitution scores at this position

    let imax = lProx.m_uResidueGroup;
    let lSortOrder = lProx.m_uSortOrder;
    let lScore = 0,
      lResProfNb = 0;
    for (let aaNum = 0; aaNum < params.abSize; aaNum++) {
      lScore = 0;
      for (let i = 0; i < imax; i++) {
        lResProfNb = lSortOrder[i];
        lScore += lProx.m_wCounts[lResProfNb] * params.scoringMatrix[aaNum][lResProfNb];
      }
      lProx.m_AAScores[aaNum] = lScore;
    }
  }

  //  Corrections for end/begin gap penalties
  if (opt & ALIGNOPT.DISABLE_FAVOR_START_GAP) {
    lProf.m_ScoreGapOpen[0] /= 2;
    lProf.m_ScoreGapOpen[l - 1] /= 2;
  }
  if (opt & ALIGNOPT.DISABLE_FAVOR_END_GAP) {
    lProf.m_ScoreGapClose[0] /= 2;
    lProf.m_ScoreGapClose[l - 1] /= 2;
  }

  return lProf;
}

/**
 * Utility function to make a profile object from a single sequence
 * @param {TSequence} pSeq: Sequence object
 * @param {number} pWeight: Sequence weight
 * @param {TAlignmentParam} pParams: Alignment parameters. Defines penalties.
 * @param {number} opt: option as a bitmask
 */
export function seqToProf(pSeq: TSequence, pWeight: number, params: TAlignmentParam, opt = 0) {
  const lEncSeq = pSeq.encodedSeq;
  /** alignment length (columns count) */
  const l = lEncSeq.length;

  const lProf = new ProfPos(l, params.abSize, 1, pWeight);

  let w = pWeight; //sequence weight

  const lGapO = params.gapOP;
  let lProx = lProf.getProxy();
  let lResNb = 0;

  for (let col = 0; col < l; col++) {
    lProx.setProxy(col);

    lResNb = lEncSeq[col];
    lProx.m_fcCounts[lResNb] = 1;
    lProx.m_wCounts[lResNb] = w;
    lProx.m_uSortOrder[lProx.m_uResidueGroup] = lResNb;
    lProx.m_uResidueGroup = 1;
    lProx.m_fOcc = w;
    lProx.m_ScoreGapOpen = (lGapO / 2) * w;
    lProx.m_ScoreGapClose = (lGapO / 2) * w;

    // Cache substitution scores at this position

    for (let aaNum = 0; aaNum < params.abSize; aaNum++) {
      lProx.m_AAScores[aaNum] = w * params.scoringMatrix[aaNum][lResNb];
    }
  }

  //  Corrections for end/begin gap penalties
  if (!(opt & ALIGNOPT.DISABLE_FAVOR_START_GAP)) {
    lProf.m_ScoreGapOpen[0] /= 2;
  }
  if (!(opt & ALIGNOPT.DISABLE_FAVOR_END_GAP)) {
    lProf.m_ScoreGapOpen[l - 1] /= 2;
  }

  return lProf;
}

/**
 * Compute a profile from two profiles and two edit strings
 *
 * @export
 * @param {ProfPos} pProfA
 * @param {ProfPos} pProfB
 * @param {number[]} pESA
 * @param {number[]} pESB
 * @param {TAlignmentParam} params
 * @param {number} [opt=0]
 * @returns
 */
export function mergeProfiles(
  pProfA: ProfPos,
  pProfB: ProfPos,
  pESA: number[],
  pESB: number[],
  params: TAlignmentParam,
  opt: number = 0,
) {
  /** profile length (alignment length, columns count) */
  const l = estringLength(pESA);
  const n = pProfA.nbSeq + pProfB.nbSeq;
  const w = pProfA.weight + pProfB.weight;

  const lProf = new ProfPos(l, params.abSize, n, w);

  const lGapO = params.gapOP;
  let lProx = lProf.getProxy();
  let lProxA = pProfA.getProxy();
  let lProxB = pProfB.getProxy();

  // indices to map columns in the new profile to the columns in the original
  // profiles.
  let lAidx = estringToIdx(pESA);
  let lBidx = estringToIdx(pESB);
  let lAi = 0;
  let lBi = 0;

  let lAiPrev = 1; // used to detect gap opening/closing
  let lBiPrev = 1;
  let lAisGapOpening = false;
  let lAisGapClosing = false;
  let lBisGapOpening = false;
  let lBisGapClosing = false;

  let lNbRes = 0; // nb of distinct residues at this column

  for (let col = 0; col < l; col++) {
    lProx.setProxy(col);

    lAi = lAidx[col];
    lProxA.setProxy(lAi);
    if (lAi === -1) {
      lAisGapOpening = lAiPrev >= 0;
      lAisGapClosing =
        col === l - 1 || // end of profile
        lAidx[col + 1] > 0; // next is not part of this indel
    } else {
      lAisGapOpening = false;
      lAisGapClosing = false;
    }
    lAiPrev = lAi;

    lBi = lBidx[col];
    lProxB.setProxy(lBi);
    if (lBi === -1) {
      lBisGapOpening = lBiPrev >= 0;
      lBisGapClosing = col == l - 1 || lBidx[col + 1] > 0;
    } else {
      lBisGapOpening = false;
      lBisGapClosing = false;
    }
    lBiPrev = lBi;

    lNbRes = 0;

    // First case: no gap introduced at this position

    if (lAi !== -1 && lBi !== -1) {
      lProx.m_fcStartOcc = lProxA.m_fcStartOcc + lProxB.m_fcStartOcc;
      lProx.m_fcEndOcc = lProxA.m_fcEndOcc + lProxB.m_fcEndOcc;
      lProx.m_fOcc = lProxA.m_fOcc + lProxB.m_fOcc;
      lProx.m_fcCounts.set(lProxA.m_fcCounts);
      lProx.m_wCounts.set(lProxA.m_wCounts);
      for (let i = 0; i < pProfB.alphaSize; i++) {
        let v = lProxB.m_fcCounts[i];
        if (v) {
          lProx.m_fcCounts[i] += v;
          lProx.m_wCounts[i] += lProxB.m_wCounts[i];
        }
        if (lProx.m_fcCounts[i]) {
          lProx.m_uSortOrder[lNbRes] = i;
          lNbRes++;
        }
        lProx.m_AAScores[i] = lProxA.m_AAScores[i] + lProxB.m_AAScores[i];
      }
      lProx.m_uResidueGroup = lNbRes;
    }

    // Gap introduced in profile A only
    else if (lAi == -1) {
      lProx.m_fcStartOcc = lProxB.m_fcStartOcc + (lAisGapOpening ? pProfA.weight : 0);
      lProx.m_fcEndOcc = lProxB.m_fcEndOcc + (lAisGapClosing ? pProfA.weight : 0);
      lProx.m_fOcc = lProxB.m_fOcc;
      lProx.m_fcCounts.set(lProxB.m_fcCounts);
      lProx.m_wCounts.set(lProxB.m_wCounts);
      lProx.m_uSortOrder.set(lProxB.m_uSortOrder);
      lProx.m_AAScores.set(lProxB.m_AAScores);
      lProx.m_uResidueGroup = lProxB.m_uResidueGroup;
    }

    // Gap introduced in profile B only
    else if (lBi == -1) {
      lProx.m_fcStartOcc = lProxA.m_fcStartOcc + (lBisGapOpening ? pProfB.weight : 0);
      lProx.m_fcEndOcc = lProxA.m_fcEndOcc + (lBisGapClosing ? pProfB.weight : 0);
      lProx.m_fOcc = lProxA.m_fOcc;
      lProx.m_fcCounts.set(lProxA.m_fcCounts);
      lProx.m_wCounts.set(lProxA.m_wCounts);
      lProx.m_uSortOrder.set(lProxA.m_uSortOrder);
      lProx.m_AAScores.set(lProxA.m_AAScores);
      lProx.m_uResidueGroup = lProxA.m_uResidueGroup;
    }

    // Open gap penalty is half the global gap open penalty,
    // weighted, reduced by the number of gaps already there
    // and opened at the same position.

    lProx.m_ScoreGapOpen = (lGapO / 2) * (1 - lProx.m_fcStartOcc) * lProx.m_fOcc;
    lProx.m_ScoreGapClose = (lGapO / 2) * (1 - lProx.m_fcEndOcc) * lProx.m_fOcc;
  }

  //  Corrections to end/begin gap penalties
  if (!(opt & ALIGNOPT.DISABLE_FAVOR_START_GAP)) {
    lProf.m_ScoreGapOpen[0] /= 2;
  }
  if (!(opt & ALIGNOPT.DISABLE_FAVOR_END_GAP)) {
    lProf.m_ScoreGapOpen[l - 1] /= 2;
  }

  return lProf;
}
