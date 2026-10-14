/**
 * @file noalign
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description This module contains code to make a fast alignment computation
 * without doing the full Dynamic Programming routine as it is overexpensive
 * for long sequences, running in O(n2) time and space. For example aligning
 * a pair of mitochondrial DNA sequences (around 16kpb) requires 128MB for the
 * storage of the trace-back (TB) matrix.
 * It is noticeable that most of the TB is unused, and most of the DP
 * computations are useless when done in regions of low homology between
 * sequences.
 * To speed-up the process for long sequences, we'll try to restrict the DP
 * to the smallest regions possible. This entails to detect the regions where
 * the sequences are identical between sequences to anchor alignments from there
 * without doing any DP recursion.
 * Here is the proposed algorithm for aligning a pair of long genomic sequences
 * SeqA and SeqB:
 * -1. Extract minimizers from SeqA and SeqB and store them in hashes with their
 *     position and the actual sequence they minimize.
 *     Minimizers (Roberts 2004) are efficient constructs based on K-mers that
 *     store in a compressed way a longer sequence for allowing fast comparisons
 *     A K-mer is a substring of fixed length (for example 8 nucleotides).
 *     Storing and comparing all K-mers in a sequence is doable but takes time
 *     and space. A Minimizer is the K-mer with the minimal value in a larger
 *     window (for example a window of 30nt). The minimal value can be the
 *     lowest lexicographic order for example.
 *     If two sequences share an identical window of 30nt, they also share
 *     the same minimizer. Minimizers can be an order of magnitude less
 *     numerous than the exhaustive K-mers enumeration, and they allow to
 *     compare larger segments in one iteration.
 * -2. Find common minimizers and collect the ones minimizing common substrings
 *     At the end of this process a collection of ranges each associated with a
 *     diagonal has been established.
 * -3. Merge/filter ranges. Contiguous ranges on the same diagonal can be merged
 *     A range that is overlapped on its both sides by other ranges can be
 *     removed from the collection. REF to add here
 * -4. Find an optimal order of the remaining diagonals (similarly to a
 *     Longest Increasing Sequence problem).
 * -5. Extend diagonals with exact and approximate matches
 * -6. Solve alignments between diagonals using the DP method (should be very
 *     fast then as the remaining DP matrices should be a fraction of the
 *     initial one)
 *
 *
 * (a) Roberts M. et al. BIOINFORMATICS Vol. 20 no. 18 2004, pages 3363–3369
 * doi:10.1093/bioinformatics/bth408 - Reducing storage requirements for
 * biological sequence comparison
 *
 */

import type { TSequence } from "../sequence/sequence.js";
import { DEQueue } from "../utils/queue.js";
import { DEBUG, type TAlignmentParam } from "./params.js";
import Log from "../utils/logger.js";
import { hammingWeight } from "../utils/bitarray.js";
import { ALIGNOPT, pairwiseAlignment } from "./align.js";
import {
  epath2estring,
  EPATH_2_STRING,
  estringCat,
  estringCountPositive,
  estringDifference,
  estringLength,
  estringMerge,
  estringProduct,
} from "../utils/estring.js";

type TRange = {
  /** diagonal number. 0 is the main diagonal. negative numbers are diagonals
   * shifted towards sequence B (gaps inserted in B). Positive numbers are
   * diagonals shifted towards sequence A (gaps inserted in A)
   */
  diagId: number;
  /** X axis (sequence A is the reference for this axis) coordinate for the
   * start (top tip) of this diagonal
   */
  begin: number;
  /** X axis coordinate for the end (bottom tip) of this diagonal */
  end: number;
};

/** A missing segment is the zone defined between the end points of two
 * successive diagonals */
type TMissingSegment = {
  beginDiagId: number;
  endDiagId: number;
  begin: number;
  end: number;
};

type TMinzStore = {
  /** Vector of minimizing kmers found in sequence */
  kmer: Uint16Array;
  // STORIES-UNTOLD(B8): the three position vectors are Int32Array; Uint16Array wrapped past 65,535 letters.
  /** Start position in sequence of kmer at same index in vector `kmer` */
  kmerPos: Int32Array;
  /** Start position of the window minimized by kmer at same index in vector `kmer` */
  winPos: Int32Array;
  /** End position of the window minimized by kmer at same index in vector `kmer` */
  winPosEnd: Int32Array;
  /** Number of minimizers found in store */
  count: number;
};

const KSIZE = 8; // fits in 16bits
const WSIZE = 24; // empirical
const EXTENSION_THRESHOLD = KSIZE / 4; // don't extend when more than 25% difference

/**
 * Given a sequence, computes its kmers and extracts the subset of
 * kmers that are minimizers (kmer with minimal value) in a sliding window.
 *
 * @export
 * @param {TSequence} seq Sequence object
 *
 * @see https://doi.org/10.1093/bioinformatics/btaa435 Chirag Jain, Arang Rhie,
 * Haowen Zhang, Claudia Chu, Brian P Walenz, Sergey Koren, Adam M Phillippy,
 * Weighted minimizer sampling improves long read mapping, Bioinformatics,
 * Volume 36, Issue Supplement_1, July 2020, Pages i111–i118
 */
export function extractMinimizers(seq: TSequence, ksize: number, wsize: number): TMinzComp {
  const lMinzMap: Map<number, number[]> = new Map();
  const lQueue = new DEQueue<number>(wsize);
  const lKSIZE = ksize | 0;

  // STORIES-UNTOLD(B9): a sequence shorter than k has no k-mers; an empty store sends the pair to the NW fallback
  // (`new Uint16Array(length - k + 1)` below had a negative length).
  if (seq.encodedSeq.length < lKSIZE) {
    const empty = new Int32Array(0);
    return [
      lMinzMap,
      { kmer: new Uint16Array(0), kmerPos: empty, winPos: empty, winPosEnd: empty, count: 0 },
      new Uint16Array(0),
    ];
  }

  // Create an array of 8-kmers from sequence
  // 8-kmers, with 2 bits per letter in a 4 letter alphabet require 16
  // bits. Let's compute it using a sliding window (in O(n) )

  const lKarr = new Uint16Array(seq.encodedSeq.length - lKSIZE + 1);
  let kval = 0 | 0;
  let lIdx = 0;
  for (let i = 0; i < lKSIZE; i++) {
    kval |= seq.encodedSeq[lKSIZE - i - 1] << (i * 2);
  }
  lKarr[lIdx++] = kval;

  for (let i = lKSIZE, imax = seq.encodedSeq.length; i < imax; i++) {
    kval = (kval << 2) + seq.encodedSeq[i];
    lKarr[lIdx++] = kval;
  }

  // Compute minimizers

  const lStartStore = wsize - lKSIZE;

  // STORIES-UNTOLD(B11): the store holds at most one minimizer per window, so it is sized to the number of windows.
  // Upstream sized it to twice a random sequence's density (4 / (wsize + 1) of the length); a tandem repeat has a new
  // minimizer every few windows, its writes past the end were dropped while `count` grew, and the seeds after it read
  // undefined positions (NaN diagonals: a wrong alignment, or no seed and the full-matrix fallback).
  const lLength = Math.max(0, lKarr.length - lStartStore);
  const lStoreMinz: TMinzStore = {
    kmer: new Uint16Array(lLength),
    kmerPos: new Int32Array(lLength),
    winPos: new Int32Array(lLength),
    winPosEnd: new Int32Array(lLength),
    count: 0,
  };
  let lCurrMinzPos = NaN;
  let j = 0;
  let lWinStart = -lStartStore - 1;

  // Main loop for minimizer computation.
  // i is indexing kmers. The kmer at index i is considered to be the
  // kmer entering in the sliding window.
  // The sliding window starts at i - window length + kmer length as Kmer(i)
  // is part of the window.
  // Computing the minimizer of the sliding window is finding the kmer with
  // the lowest value in this window.
  // This is done using a deque where kmers entering the window (kmer(i))
  // are appended to the tail and remove all kmers at the tail that are of
  // higher value. This maintains a sorting of kmers in the deque.
  // Hence the deque head contains the kmer with minimal value.

  for (let i = 0; i < lKarr.length; i++) {
    const lKmer = lKarr[i];
    lWinStart++;

    // Remove kmers with lower order than current kmer
    // Note: for sake of symplicity, order is just kmer value.

    while (
      !lQueue.isEmpty &&
      ((j = lQueue.getTail()) < lWinStart || // not in this window anymore
        lKarr[j] >= lKmer) // can't be a minimizer of this and future windows: overthrown by kmer(i)
    ) {
      lQueue.popTail();
    }

    // Append current kmer to tail

    lQueue.pushTail(i);

    if (lWinStart < 0) continue; // compute all kmers in first window first

    // Remove kmers that are not in this window anymore

    while (lQueue.getHead() < lWinStart) {
      lQueue.popHead();
    }

    // Store head kmer in minimizer hash

    let lHead = lQueue.getHead();

    // Is it the same minimizer as in previous loop iteration?
    // If yes, only extend the minimized subarray

    if (lHead === lCurrMinzPos) {
      lStoreMinz.winPosEnd[lStoreMinz.count - 1] = i + lKSIZE;
    } else {
      const lHeadKmer = lKarr[lHead];
      const lIdx = lStoreMinz.count;
      lStoreMinz.count++;

      lStoreMinz.kmer[lIdx] = lHeadKmer;
      lStoreMinz.kmerPos[lIdx] = lHead;
      lStoreMinz.winPos[lIdx] = lWinStart;
      lStoreMinz.winPosEnd[lIdx] = i + lKSIZE;
      lCurrMinzPos = lHead;

      // Add the kmer to the map, either as a new list or as a new
      // item in a previous list (when the same kmer also minimizes
      // another window)

      if (!lMinzMap.has(lHeadKmer)) {
        lMinzMap.set(lHeadKmer, [lIdx]);
      } else {
        let lList = lMinzMap.get(lHeadKmer)!;
        lList.push(lIdx);
      }
    }
  }

  return [lMinzMap, lStoreMinz, lKarr];
}

/** Map Kmer value to a table of all its indices in the sequence */
type TMinzMap = Map<number, number[]>;
type TMinzComp = [TMinzMap, TMinzStore, Uint16Array];

export function noalignPair(
  seqA: TSequence,
  seqB: TSequence,
  pAlignParam: TAlignmentParam,
  pMinzA?: TMinzComp,
  pMinzB?: TMinzComp,
) {
  const lDebugStats: { [k: string]: Partial<{ a: any; b: any; all: any }> } = {};

  const [lMinzA, lMinzAStore, lKmerAArr] = pMinzA ?? extractMinimizers(seqA, KSIZE, WSIZE);
  const [lMinzB, lMinzBStore, lKmerBArr] = pMinzB ?? extractMinimizers(seqB, KSIZE, WSIZE);
  if (DEBUG) {
    Log.add("Extract Minimizers");
    lDebugStats["Nb Minimizers"] = { a: lMinzAStore.count, b: lMinzBStore.count };
    lDebugStats["Dupl. Minimizers"] = { a: lMinzAStore.count - lMinzA.size, b: lMinzBStore.count - lMinzB.size };
  }

  let lNbCommonWindows = 0;
  const lMaxLen = Math.max(seqA.rawSeq.length, seqB.rawSeq.length);
  const lDiagMap = new Map<number, TRange[]>();

  // Loop through all kmers from sequence A to find which ones are matching
  // in sequence B. Doing so, favour the kmers that minimize a common window.
  // Those are better seed candidates for diagonals. They also cover a longer
  // range and thus avoid extra computations.

  for (let i = 0; i < lMinzAStore.count; i++) {
    let kmer = lMinzAStore.kmer[i];

    if (!lMinzB.has(kmer)) continue;
    let listB = lMinzB.get(kmer)!;
    if (listB.length > 4) continue; // poor minimizer. prevent combinatorial explosion

    // Compare minimized string in A with those in listB to take advantage of
    // the ones that share common minimized strings

    let lMinzSubA = seqA.rawSeq.substring(lMinzAStore.winPos[i], lMinzAStore.winPosEnd[i]);

    for (let j = 0; j < listB.length; j++) {
      let lBidx = listB[j];
      let lMinzSubB = seqB.rawSeq.substring(lMinzBStore.winPos[lBidx], lMinzBStore.winPosEnd[lBidx]);
      let lLen = Math.min(lMinzSubA.length, lMinzSubB.length);

      let lHasCommonWindow = false;
      if (lLen == lMinzSubA.length) {
        lHasCommonWindow = lMinzSubB.indexOf(lMinzSubA) === 0;
      } else lHasCommonWindow = lMinzSubA.indexOf(lMinzSubB) === 0;

      // common range to store
      let lDiagId = 0;
      let lRange: TRange;
      if (lHasCommonWindow) {
        lDiagId = lMinzAStore.winPos[i] - lMinzBStore.winPos[lBidx];
        lRange = {
          diagId: lDiagId,
          begin: lMinzAStore.winPos[i],
          end: lMinzAStore.winPos[i] + lLen,
        };
        lNbCommonWindows++;
      } else {
        lDiagId = lMinzAStore.kmerPos[i] - lMinzBStore.kmerPos[lBidx];
        lRange = {
          diagId: lDiagId,
          begin: lMinzAStore.kmerPos[i],
          end: lMinzAStore.kmerPos[i] + KSIZE - 1,
        };
      }

      if (!lDiagMap.has(lDiagId)) {
        lDiagMap.set(lDiagId, [lRange]);
      } else {
        let lDiagList = lDiagMap.get(lDiagId);
        lDiagList!.push(lRange);
      }
    }
  }
  if (DEBUG) {
    Log.add("Filter Common Minimizers");
    lDebugStats["Nb common windows"] = { all: lNbCommonWindows };
  }

  let lDiagList: TRange[] = [];
  // Merge consecutive or overlapping subarrays on the same diagonals
  // Note: the extensions are made on range references
  lDiagMap.forEach((diags) => {
    if (!diags.length) return; // safeguard

    let lCurrentSegment = diags[0];
    lDiagList.push(lCurrentSegment);
    diags.forEach((range) => {
      if (range.begin <= lCurrentSegment.end) {
        lCurrentSegment.end = range.end;
      } else {
        lCurrentSegment = range;
        lDiagList.push(lCurrentSegment);
      }
    });
  });
  if (DEBUG) {
    Log.add("Merge Minimizers");
    lDebugStats["Nb diagonals"] = { all: lDiagList.length };
  }

  // Special case: no diagonal found: fallback to NW alignment
  if (lDiagList.length === 0) {
    if (DEBUG) console.table(lDebugStats);

    return pairwiseAlignment(seqA, seqB, pAlignParam).estrings;
  }

  // Sort ranges in diag list from left to right and from top to bottom
  // STORIES-UNTOLD(B12): by the position in A; `begin + diagId` put a diagonal after a long insertion in B first.
  lDiagList.sort((a, b) => {
    let lDelta = a.begin - b.begin;
    return lDelta === 0 ? a.diagId - b.diagId : lDelta;
  });
  if (DEBUG) Log.add("Sort Minimizers");

  const lLenA = seqA.encodedSeq.length;
  const lLenB = seqB.encodedSeq.length;
  lDiagList = extractLongestPath(lDiagList, lLenA, lLenB);
  // STORIES-UNTOLD(B14): an empty diagonal at the end of both sequences, so the trailing segment is aligned as the
  // leading one is (below); upstream appended the tails, pairing the first letters left in each.
  const lLast = lDiagList.length;
  lDiagList.push({ diagId: lLenA - lLenB, begin: lLenA, end: lLenA });

  if (DEBUG) {
    Log.add("Filter Minimizers - Optimal list");
    let lCoverage = 0;
    lDiagList.forEach((d) => {
      lCoverage += d.end - d.begin;
    });
    lCoverage /= lMaxLen;
    lDebugStats["Strict Coverage"] = { all: lCoverage };

    lDebugStats["Nb filtered diagonals"] = { all: lDiagList.length };
  }

  // Extend consecutive segments on same diagonal. We apply a tolerance for
  // extension at 25% difference between Kmers.
  // This is measured by bitwise operations between 16bits numbers.

  let lPrevSegment: TRange = { diagId: 0, begin: 0, end: 0 }; //lDiagList[0];
  let lExtDiagList = [lPrevSegment];
  let lMissingSegments: TMissingSegment[] = [];

  for (let i = 0; i < lDiagList.length; i++) {
    let lCurrentSegment = lDiagList[i];
    // STORIES-UNTOLD(B12): a diagonal may overlap the previous one by a few letters (it was chained so); trim it.
    lCurrentSegment.begin += Math.max(
      0,
      lPrevSegment.end - lCurrentSegment.begin,
      lPrevSegment.end - lPrevSegment.diagId - lCurrentSegment.begin + lCurrentSegment.diagId,
    );

    // TODO: when difference between end and begin is >> KSIZE (e.g. 400),
    // there are likely unmapped gaps. It might be useful to speed up
    // matches using shorter Kmers/windows. Test the overhead + noise
    // generated.

    // first case: this segment is on the same diagonal as the previous one
    // Try to extend begin position towards previous segment.

    if (lCurrentSegment.diagId === lPrevSegment.diagId && i > 0 && i < lLast) {
      // Compare kmers in forward direction until score drops
      let m = lPrevSegment.end; // end position on seq A
      let n = m - lCurrentSegment.diagId;
      let lWinScore = 0;

      while (m < lCurrentSegment.begin && n < lCurrentSegment.begin - lCurrentSegment.diagId && lWinScore < 2) {
        m += KSIZE;
        n += KSIZE;
        lWinScore = dnaHammingDistance(lKmerAArr[m], lKmerBArr[n]) <= EXTENSION_THRESHOLD ? 0 : lWinScore + 1;
      }

      // To avoid boundaries conditions, remove the tip of the extension
      if (m > lPrevSegment.end) {
        lPrevSegment.end = Math.max(m - EXTENSION_THRESHOLD - KSIZE, lPrevSegment.end);
      }

      // Compare in backward direction until score drops or merge has been
      // achieved

      m = lCurrentSegment.begin - KSIZE;
      n = m - lCurrentSegment.diagId;
      lWinScore = 0;
      // WIP Make stats on full window, not just 1 Kmer
      while (m > lPrevSegment.end && n > lPrevSegment.end - lPrevSegment.diagId && lWinScore < 2) {
        m -= KSIZE;
        n -= KSIZE;
        lWinScore = dnaHammingDistance(lKmerAArr[m], lKmerBArr[n]) <= EXTENSION_THRESHOLD ? 0 : lWinScore + 1;
      }
      //TODO: check this, the Maths look dubious
      if (m < lCurrentSegment.begin - 2 * KSIZE) {
        lCurrentSegment.begin = m + 3 * KSIZE + EXTENSION_THRESHOLD;
      }

      // merge?
      if (lCurrentSegment.begin - lPrevSegment.end <= 3 * KSIZE) {
        lPrevSegment.end = lCurrentSegment.end;
        continue;
      } else {
        lMissingSegments.push({
          beginDiagId: lPrevSegment.diagId,
          endDiagId: lCurrentSegment.diagId,
          begin: lPrevSegment.end,
          end: lCurrentSegment.begin,
        });
        lExtDiagList.push(lCurrentSegment);
        lPrevSegment = lCurrentSegment;
        continue;
      }
    } else {
      // Compare nucleotides in forward direction until they differ
      let m = lPrevSegment.end; // end position on seq A
      let n = m - lPrevSegment.diagId;

      // STORIES-UNTOLD(B14): the empty diagonals at (0, 0) and at the end are not extended: the segments next to them
      // are aligned (below), where extending them paired chance equal letters at the ends (`3M6000D` for `6000D`).
      while (
        i > 0 &&
        seqA.encodedSeq[m] === seqB.encodedSeq[n] &&
        m < lCurrentSegment.begin &&
        n < lCurrentSegment.begin - lCurrentSegment.diagId
      ) {
        m++;
        n++;
      }

      lPrevSegment.end = m; // last before discrepancy. This position is
      // exclusive in the prevSegment, and inclusive
      // in the missing segment if any.

      // Compare in backward direction until sequences differ

      m = lCurrentSegment.begin;
      n = m - lCurrentSegment.diagId;

      // STORIES-UNTOLD(B14): compares the letters before `m`, so the extension reaches the bound when they are equal
      // (it stopped one short, leaving a letter to the DP: `1M6000I` where a leading `6000I` was exact).
      while (
        i < lLast &&
        m > lPrevSegment.end &&
        n > lPrevSegment.end - lPrevSegment.diagId &&
        seqA.encodedSeq[m - 1] === seqB.encodedSeq[n - 1]
      ) {
        m--;
        n--;
      }

      lCurrentSegment.begin = m; // first after discrepancy

      lMissingSegments.push({
        endDiagId: lCurrentSegment.diagId,
        beginDiagId: lPrevSegment.diagId,
        begin: lPrevSegment.end,
        end: lCurrentSegment.begin,
      });
      lExtDiagList.push(lCurrentSegment);
      lPrevSegment = lCurrentSegment;
    }
  }

  if (DEBUG) {
    Log.add("Diagonals extension");
    let lCoverage = 0;
    lExtDiagList.forEach((d) => {
      lCoverage += d.end - d.begin;
    });
    lCoverage /= lMaxLen;
    lDebugStats["Coverage"] = { all: lCoverage };
    lDebugStats["Extended Diagonals"] = { all: lExtDiagList.length };
    console.table(lDebugStats);

    //console.table(lMissingSegments);
  }

  // Make the alignment from the diagonals and fill the gaps
  let lDiag = lExtDiagList[0];
  const lEpathA: number[] = [];
  const lEpathB: number[] = [];
  for (let i = 0; i < lExtDiagList.length; i++) {
    lDiag = lExtDiagList[i];
    lEpathA.push(lDiag.end - lDiag.begin);
    lEpathB.push(lDiag.end - lDiag.begin);

    let lMis = lMissingSegments[i];
    if (lMis) {
      // Special case next diag begins at same col prev diag finishes
      if (lMis.begin === lMis.end) {
        const lDelta = Math.abs(lMis.endDiagId - lMis.beginDiagId);
        lEpathA.push(-lDelta);
        lEpathB.push(lDelta);
        continue;
      }

      // Special case next diag begins at same row prev diag finishes
      if (lMis.begin - lMis.beginDiagId === lMis.end - lMis.endDiagId) {
        const lDelta = lMis.end - lMis.begin;
        lEpathB.push(-lDelta);
        lEpathA.push(lDelta);
        continue;
      }

      let lResult = pairwiseAlignment(
        {
          rawSeq: seqA.rawSeq.substring(lMis.begin, lMis.end),
          type: seqA.type,
          compressedSeq: new Uint8Array(0),
          encodedSeq: seqA.encodedSeq.subarray(lMis.begin, lMis.end),
        },
        {
          rawSeq: seqB.rawSeq.substring(lMis.begin - lMis.beginDiagId, lMis.end - lMis.endDiagId),
          type: seqB.type,
          compressedSeq: new Uint8Array(0),
          encodedSeq: seqB.encodedSeq.subarray(lMis.begin - lMis.beginDiagId, lMis.end - lMis.endDiagId),
        },
        pAlignParam,
        // STORIES-UNTOLD(B14): the leading and trailing segments favor their outer gaps, as a full DP does.
        (i === 0 ? 0 : ALIGNOPT.DISABLE_FAVOR_START_GAP) |
          (i === lExtDiagList.length - 2 ? 0 : ALIGNOPT.DISABLE_FAVOR_END_GAP),
      );

      lEpathA.push(...lResult.estrings[0]);
      lEpathB.push(...lResult.estrings[1]);
    }
  }

  // Append tail if any
  // At this stage both aligned sequences have the same length.
  // For now, just append the missing part. Ultimately, assess, depending
  // on the tail length if a SW alignmnent is deemed necessary.

  sanitizeEpathEnd(lEpathA, seqA.encodedSeq.length);
  sanitizeEpathEnd(lEpathB, seqB.encodedSeq.length);

  let lEstringA = epath2estring(lEpathA, EPATH_2_STRING.NO_REVERSE);
  let lEstringB = epath2estring(lEpathB, EPATH_2_STRING.NO_REVERSE);

  const lEszA = estringLength(lEstringA);
  const lEszB = estringLength(lEstringB);
  if (lEszA > lEszB) lEstringB = estringCat(lEstringB, [lEszB - lEszA]);
  else if (lEszA < lEszB) lEstringA = estringCat(lEstringA, [lEszA - lEszB]);

  if (DEBUG) {
    Log.add("Fill between diagonals");
  }

  // TODO: Finalize alignment here
  return [lEstringA, lEstringB];
}

/**
 * Given an epath, it ensures that the sum of its positive values equals the
 * target length.
 *
 * @param {number[]} epath
 * @param {number} targetLength
 */
function sanitizeEpathEnd(epath: number[], targetLength: number) {
  let lCount = estringCountPositive(epath);
  let lDelta = targetLength - lCount;
  if (lDelta > 0) epath.push(lDelta);
  else if (lDelta < 0) epath[epath.length - 1] += lDelta;
}

/**
 * Returns the number of different letters in DNA alphabet between 2 kmers
 * encoded as numbers, using bitwise operations
 *
 * @param {number} a
 * @param {number} b
 */
export function dnaHammingDistance(a: number, b: number) {
  // DNA alphabet encodes 4 symbols using 2 bits. In the following, we'll
  // count the pairs of bits that are different between a and b.
  let c = a ^ b; // set bit for pair of bits that are different
  c = (c & 0x55555555) | ((c >> 1) & 0x55555555); // intersect with 010101... to retain only the lowest bit in each pair
  return hammingWeight(c);
}

/**
 * Align multiple (long) sequences using the center-star method.
 *
 * @export
 * @param {TSequence[]} pSeq
 * @param {TAlignmentParam} pAlignParam
 * @returns
 */
export function centerStarNoAlign(pSeq: TSequence[], pAlignParam: TAlignmentParam) {
  let lMinz = pSeq.map((seq) => {
    return extractMinimizers(seq, KSIZE, WSIZE);
  });

  // TODO: find center sequence
  let lCenter = 0;
  let lEStrings: number[][] = [];
  let lCenterEString: number[] = [];
  let lCenterEStringRaw: number[][] = [];

  for (let i = 0; i < lMinz.length; i++) {
    if (i === lCenter) continue;

    let lEsAln = noalignPair(pSeq[lCenter], pSeq[i], pAlignParam, lMinz[lCenter], lMinz[i]);

    // Center Edit strings combines the Edit strings from the alignment to the
    // current sequence with its previous edit strings.
    lCenterEString = lCenterEString.length === 0 ? lEsAln[0] : estringMerge(lCenterEString, lEsAln[0]);

    // Keep track of what was this alignment edit string for center sequence
    // This is used to make a diff, at next iteration.
    lCenterEStringRaw[i] = lEsAln[0];

    // Edit string to align this sequence with the center. It will be
    // modified to account for other sequences introduced editions
    lEStrings[i] = lEsAln[1];
  }

  for (let i = 0; i < lMinz.length; i++) {
    if (i === lCenter) {
      lEStrings[i] = lCenterEString;
      continue;
    }
    let lDiff = estringDifference(lCenterEString, lCenterEStringRaw[i]);
    if (!lDiff) {
      console.error(
        "An error occured while computing the center diff. for sequence #" + i,
        lCenterEString,
        lCenterEStringRaw[i],
      );
      console.dir(lEStrings);
      console.dir(lCenterEStringRaw);
      throw new RangeError();
    }
    lEStrings[i] = estringProduct(lDiff, lEStrings[i]);
  }

  return lEStrings;
}

/**
 * Find an optimal order of traversal among a list of diagonals so that they
 * make the longest past covered by diagonals in the alignment matrix between
 * both sequences.
 * @param pDiagList list of diagonals ordered from left to right and top to
 *                  bottom
 * @param lMaxLen   maximum sequence length. Used to compute a minimal gap
 *                  penalty
 * @returns
 */
// STORIES-UNTOLD(B13): exported for the tests.
export function extractLongestPath(pDiagList: TRange[], pLenA: number, pLenB: number) {
  /* Find optimal order of diagonals using Longest Increasing Sequence
       like algorithm. Here, the Word Sequence is replaced by Suite so as to
       avoid confusion with biological sequences.
       In the context of diagonals, a new diagonal can be appended to a list
       iff it is in the SE quadrant from the bottom tip of the last diagonal
       in the list.
       ____________________________________________________________________
       |              \
       |  \            \ B
       |   \ A          \
       |    \                 \
       |                       \ C
       |                        \

       C can extend B but it can't extand A

       A new diagonal can start a new list if:
       - its diag id is > to the one of the current shorter list in the suites
         (the rationale is that any diagonal that would extend a diagonal below
         the smallest one would also be an extension of the smallest one)
       - and its row (defined by begin - diag id) is lower.
         (the rationale is that this new diagonal will be able to catch next
         diagonals that can't be caught by the current minimum)
       - or if this diagonal is, by itself, long enough to be the optimum while
         not being able to be appended to any other list.

       Otherwise, we search for the list of diagonals with maximum size that can
       accept the new diagonal and we duplicate it + append the new diagonal to
       it. The new suite is placed in the list of suites so that the order of
       increasing sizes is maintained.
       When this addition is made, the list of suites is cleaned to remove any
       list that is shorter and that also extends past this new extent (in other
       words a diagonal that extends past, is a diagonal which bottom tip is
       below the lowest point of the new diagonal).
    */

  const lUnitPenalty = 1 / Math.max(pLenA, pLenB);
  const lDividerDiag = pLenA - pLenB;

  /** Collections of list of increasing diagonals */
  let lIncreasingPaths: number[] = [0];
  /** Stores at the index corresponding to each diagonal, the index of its
   * preceding diagonal
   */
  // STORIES-UNTOLD(B8): Int32Array (an Int16Array wrapped past 32,767 diagonals).
  let lTraceBack = new Int32Array(pDiagList.length);
  lTraceBack[0] = -1; // starts a path

  /** The current diag in the list */
  let lDiag = pDiagList[0];
  /** Diagonal length. What it will add to the score of the path it may be appended to */
  let lDiagLen = lDiag.end - lDiag.begin;
  /** Bottom row the diagonal extends to. It's the end posision of the segment in the Sequence B */
  let lDiagBottomRow = lDiag.end - lDiag.diagId;
  /** Top row the diagonal begins at. It's the start posision of the segment in the Sequence B */
  let lDiagTopRow = lDiag.begin - lDiag.diagId;
  let lDiagDistToEnd = Math.min(pLenB - lDiagBottomRow, pLenA - lDiag.end);
  /** Length of each list of increasing diagonals as the sum of diags length */
  let lPathsScores: number[] = [lDiagLen];
  /** Y position of the tip of the path */
  let lPathsBottomRow: number[] = [lDiagBottomRow];

  let lPathsDistToEnd = [lDiagDistToEnd];
  /** Reference to current minimum in list of paths */
  let lMinDiag = pDiagList[0];
  let lMaxDiagId = 0;
  let lMaxDiag = pDiagList[lMaxDiagId];
  let lOptimumScore = lDiagLen;

  let lGCCount = 0;

  for (let i = 1; i < pDiagList.length; i++) {
    lDiag = pDiagList[i];
    lDiagLen = lDiag.end - lDiag.begin;
    lDiagBottomRow = lDiag.end - lDiag.diagId;
    lDiagTopRow = lDiag.begin - lDiag.diagId;
    lDiagDistToEnd = lDiag.diagId < lDividerDiag ? pLenB - lDiagBottomRow : pLenA - lDiag.end;

    // Garbage Collection

    if (lGCCount === 20) {
      lGCCount = 0;
      let lGarbageIdx: number[] = [];
      let k = lIncreasingPaths.length - 1;
      let lPrevBottom = lPathsBottomRow[k];
      let lRemainingDist = pLenA - lDiag.begin;

      // clean paths array by removing all paths:
      // - that extend below higher scoring path
      // - that score too low for ever being in an optimal path

      while (k--) {
        if (lPathsBottomRow[k] > lPrevBottom) {
          lGarbageIdx.push(k);
          continue;
        }

        lPrevBottom = lPathsBottomRow[k];

        if (lPathsScores[k] + Math.min(lPathsDistToEnd[k], lRemainingDist) < lOptimumScore) {
          lGarbageIdx.push(k);
        }
      }

      lGarbageIdx.forEach((k) => {
        lIncreasingPaths.splice(k, 1);
        lPathsScores.splice(k, 1);
        lPathsBottomRow.splice(k, 1);
        lPathsDistToEnd.splice(k, 1);
      });
    }

    if (
      lDiag.diagId >= lMinDiag.diagId && // diagonals below the minimum can't be a candidate for a new minimum
      lDiagBottomRow < lPathsBottomRow[0] // diagonal tip row must be above current minimum
    ) {
      if (lDiagLen + lDiagDistToEnd < lOptimumScore) continue;

      lTraceBack[i] = -1;
      lMinDiag = pDiagList[i];
      lGCCount++;

      if (lDiagLen < lPathsScores[0]) {
        lIncreasingPaths.unshift(i);
        lPathsScores.unshift(lDiagLen);
        lPathsBottomRow.unshift(lDiagBottomRow);
        lPathsDistToEnd.unshift(lDiagDistToEnd);
        continue;
      }

      if (lDiagLen > lOptimumScore) lOptimumScore = lDiagLen;
      lIncreasingPaths[0] = i;
      lPathsScores[0] = lDiagLen;
      lPathsBottomRow[0] = lDiagBottomRow;
      lPathsDistToEnd[0] = lDiagDistToEnd;
      continue;
    }

    // Loop over the list of increasing sequences and see to which one this
    // diagonal could be appended to.
    // STORIES-UNTOLD(B12): the path the diagonal extends best, not the first from the top that accepts it. A jump to
    // another diagonal costs its length up to KSIZE, more than a chance k-mer match (7 letters) adds, so a seed inside
    // a long insertion does not split it in two. The diagonal may overlap the path's tip by fewer letters than its
    // length (equal letters at an insertion's ends), only the rest counted (noalignPair trims it).

    let lBest = -1;
    let lNewSize = lDiagLen; // the diagonal alone
    for (let j = lIncreasingPaths.length - 1; j >= 0 && lPathsScores[j] + lDiagLen > lNewSize; j--) {
      lMaxDiag = pDiagList[lIncreasingPaths[j]];
      const lOverlap = Math.max(0, lMaxDiag.end - lDiag.begin, lPathsBottomRow[j] - lDiagTopRow);
      if (lOverlap >= lDiagLen) continue;
      // This small gap penalty favours diagonals that extend a
      // list with the lowest jump. It should only have an effect
      // when there is a tie between short diagonals
      // STORIES-UNTOLD(B13): the tip's diagonal number; upstream subtracted the path's index in the list.
      const lJump = lDiag.diagId - lMaxDiag.diagId;
      const lGapPenalty = lJump === 0 ? 0 : Math.min(KSIZE, Math.abs(lJump)) + Math.abs(lJump) * lUnitPenalty;
      const lSize = lPathsScores[j] + lDiagLen - lOverlap - lGapPenalty;
      if (lSize > lNewSize) {
        lNewSize = lSize;
        lBest = j;
      }
    }

    if (lBest >= 0 && lNewSize + lDiagDistToEnd >= lOptimumScore) {
      // TODO: binary search to find where it should be inserted now
      let k = lBest;
      while (k > 0 && lPathsScores[k - 1] > lNewSize) k--;
      while (k < lIncreasingPaths.length && lPathsScores[k] < lNewSize) k++;
      if (k === lIncreasingPaths.length || lDiagBottomRow < lPathsBottomRow[lPathsBottomRow.length - 1]) {
        if (k === lIncreasingPaths.length) lOptimumScore = lNewSize;
        lTraceBack[i] = lIncreasingPaths[lBest];
        lPathsDistToEnd.splice(k, 0, lDiagDistToEnd);
        lPathsScores.splice(k, 0, lNewSize);
        lPathsBottomRow.splice(k, 0, lDiagBottomRow);
        lIncreasingPaths.splice(k, 0, i);
        lGCCount++;
      }
    }

    // special case: the diagonal is an optimal path in itself

    if (lDiagLen > lOptimumScore) {
      lTraceBack[i] = -1;
      lOptimumScore = lDiagLen;
      lPathsDistToEnd.push(lDiagDistToEnd);
      lPathsScores.push(lDiagLen);
      lPathsBottomRow.push(lDiagBottomRow);
      lIncreasingPaths.push(i);

      lGCCount++;
    }
  }

  let lPath = [];
  let k = lIncreasingPaths.pop() ?? 0;
  let lSafeguard = k;
  while (k > -1 && lSafeguard > -1) {
    lPath.push(k);
    k = lTraceBack[k];
    lSafeguard--;
  }

  let lLongestPath: TRange[] = [];
  for (let i = lPath.length - 1; i >= 0; i--) {
    lLongestPath.push(pDiagList[lPath[i]]);
  }
  return lLongestPath;
}
