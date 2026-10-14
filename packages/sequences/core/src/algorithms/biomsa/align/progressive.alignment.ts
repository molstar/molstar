/**
 * @file progressive.alignment.ts
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description Controller for the optimal multiple alignment procedure
 * Most of this algorithm is inspired by R.C. Edgar MUSCLE library (a).
 *
 * (a) Edgar, R.C. MUSCLE: a multiple sequence alignment method with reduced
 * time and space complexity. BMC Bioinformatics 5, 113 (2004).
 * https://doi.org/10.1186/1471-2105-5-113
 *
 */

import { DEBUG, type TAlignmentParam } from './params.js';
import { type TSequence, distanceMatrix, sortMSA } from '../sequence/sequence.js';
import { type InternalNode, isLeafNode, makeTree, clustalWeights, setProfiles } from '../sequence/tree.js';
import Log from '../utils/logger.js';
import { pairwiseAlignment, MSASeqAlignment, MSAMSAAlignment } from './align.js';
import { estringProduct } from '../utils/estring.js';
import { mergeProfiles } from '../sequence/profile.js';

// STORIES-UNTOLD(B17): `step` is awaited before each alignment of a node's two children (a pair, a profile and a
// sequence, or two profiles), so the caller can yield to the event loop and throw to stop; `treeAlign` is async for it.
export async function progressiveAlignment(
  seq: TSequence[],
  pParam: TAlignmentParam,
  step?: () => void | Promise<void>,
) {
  // Recursive depth first tree based alignment
  const treeAlign = async (node: InternalNode, tabWeight: number[]) => {
    const nodeA = tree[node.childA],
      nodeB = tree[node.childB];
    let result = {} as { score: number };

    // note: msa property has already been computed when the tree
    // has had subtrees from a previous alignment, copied over.
    // Then, they can be skipped to avoid costly computations.
    if (!isLeafNode(nodeA) && nodeA.msa.length === 0) {
      await treeAlign(nodeA, tabWeight);
    }

    if (!isLeafNode(nodeB) && nodeB.msa.length === 0) {
      await treeAlign(nodeB, tabWeight);
    }

    // Now nodeA and nodeB are either leaf or msa, we can align them
    await step?.();
    if (isLeafNode(nodeA)) {
      //A is a single sequence
      if (isLeafNode(nodeB)) {
        //B is a single sequence
        const lR = pairwiseAlignment(nodeA.seq, nodeB.seq, pParam);
        result.score = lR.score;

        node.numSeq = [nodeA.numSeq[0], nodeB.numSeq[0]];
        node.estring = lR.estrings;
        node.profile = mergeProfiles(nodeA.profile, nodeB.profile, lR.estrings[0], lR.estrings[1], pParam);

        if (DEBUG) Log.add(`seq ${nodeA.numSeq} - seq ${nodeB.numSeq}`);
      } else {
        //B is a MSA
        nodeB.tabWeight = nodeB.numSeq.map((seqidx) => tabWeight[seqidx]);

        // !!! B and A order are inverted

        const lR = MSASeqAlignment(nodeB, nodeA, pParam);
        result.score = lR.score;

        // !!! estrings[0] applies to nodeB

        node.numSeq = [nodeA.numSeq[0], ...nodeB.numSeq];
        node.estring = [lR.estrings[1], ...nodeB.estring.map((es) => estringProduct(lR.estrings[0], es))];
        node.profile = mergeProfiles(nodeA.profile, nodeB.profile, lR.estrings[1], lR.estrings[0], pParam);

        if (DEBUG) Log.add(`seq ${nodeA.numSeq} - MSA ${nodeB.numSeq}`);
      }
    } else if (!isLeafNode(nodeB)) {
      // A & B are both MSA

      nodeB.tabWeight = nodeB.numSeq.map((seqidx) => tabWeight[seqidx]);
      nodeA.tabWeight = nodeA.numSeq.map((seqidx) => tabWeight[seqidx]);

      const lR = MSAMSAAlignment(nodeA, nodeB, pParam);
      result.score = lR.score;

      node.estring = [
        ...nodeA.estring.map((es) => estringProduct(lR.estrings[0], es)),
        ...nodeB.estring.map((es) => estringProduct(lR.estrings[1], es)),
      ];
      node.profile = mergeProfiles(nodeA.profile, nodeB.profile, lR.estrings[0], lR.estrings[1], pParam);

      if (DEBUG) Log.add(`MSA ${nodeA.numSeq} - MSA ${nodeB.numSeq}`);
    }
    node.numSeq = [...nodeA.numSeq, ...nodeB.numSeq];

    return result.score;
  };

  if (DEBUG) Log.add('Start Progressive Alignment');

  // Compute fast pair similarity (see k-mers binary Muscle)
  const lDistMatrixKMers = distanceMatrix(seq);

  if (DEBUG) Log.add('K-mer distance matrix');

  // Compute tree from distance matrix
  const tree = makeTree(lDistMatrixKMers, seq);
  const root = tree[tree.length - 1] as InternalNode;

  if (DEBUG) Log.add('Build Tree');

  // Compute weights (this is used to prevent close sequences to skew
  // the computation toward them, which would result in not disfavouring
  // gaps in their already aligned sequences).
  const weights = clustalWeights(tree);

  if (DEBUG) Log.add('Compute Weights - Start MSA');

  setProfiles(tree, pParam);

  // First alignment following guide tree
  const score1 = await treeAlign(root, weights);
  let lEstrings = sortMSA(root.estring.slice(), root.numSeq);

  if (DEBUG) Log.add('End MSA computation');

  return lEstrings;

  //================== For now, don't refine ====================

  // Refinement.
  // Now that we got a first alignment of all sequences, let see if some
  // sequences have a better likeliness now than what was found in the
  // first place with the fast k-mer based distance.
  // First step: build a distance matrix based on identities, pondered
  // using Kimura method for long branches (effective mutations are >
  // to measured mutations)
  // const lDistMatrixKimura = distanceKimura(msa);

  // if (DEBUG)
  //     Log.add('Distance Kimura');

  // // compute new tree
  // var tree2 = makeTree(lDistMatrixKimura, seq);
  // const root2 = tree2[tree2.length - 1] as InternalNode;

  // if (DEBUG)
  //     Log.add('Build tree 2');

  // // compare trees.
  // // Note: Impure function, it will modify tree2 if differences are found
  // // to take advantage of all the alignments already compared.
  // var treeIdentity = compareTrees(tree, tree2);

  // if (treeIdentity === true) { // No possible improvement
  //     return msa;
  // } else {
  //     //Repeat progressive alignment along this new tree
  //     var tabWeight2 = clustalWeights(tree2);

  //     //$log.debug(tabWeight);
  //     const score2 = treeAlign(root2, tabWeight2);

  //     //$log.debug('scores:', score1, score2);
  //     //var matriceDistancesKimura2 = _utils.distancesKimura(tree2[tree2.length - 1].msa);
  //     if (score2 > score1) {
  //         msa = sortMSA(root2.msa, root2.numSeq);
  //     }
  //     return msa;
  // }
}
