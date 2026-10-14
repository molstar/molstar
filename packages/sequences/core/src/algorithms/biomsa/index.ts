/**
 * @file index
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2020-2021
 */

import { makeSequence, getSequenceType, type TSequence } from './sequence/sequence.js';
import { DEBUG, SEQUENCE_TYPE, getAlignmentParameters, DEFAULT_GAP_CHAR } from './align/params.js';
import { pairwiseAlignment } from './align/align.js';
import { progressiveAlignment } from './align/progressive.alignment.js';
import Log from './utils/logger.js';
import { centerStarNoAlign, noalignPair } from './align/noalign.js';
import { applyEstring } from './utils/estring.js';

type TAlignMethod = 'diag' | 'complete' | 'auto';
type TAlignSequenceType = 'amino' | 'nucleic' | 'auto';
export type TAlignOpt = {
  /**
   * Gap open penalty (generally a negative number)
   * @default -300|-400 set based on the sequence type
   */
  gapopen: number;

  /**
   * Gap extend penalty (generally a negative number)
   * @default -11|-60 set based on the sequence type
   */
  gapextend: number;

  /**
   * Set true to report debug information to the JS console
   * @default false
   */
  debug: boolean;

  /**
   * Alignment method. `auto` is the default value. For sequences length up
   * to 1600 residues, a `complete` Needleman-Wunsch global alignment is
   * performed. Beyond, a diagonal finding algorithm  (`diag`) is run first
   * and global alignment is restricted to regions between diagonals.
   * `complete` and `diag` parameters force the use of one of those
   * strategies.
   * @default 'auto'
   */
  method: TAlignMethod;

  /**
   * Sequence type
   * `auto` performs a detection based on sequences content
   * @default 'auto'
   */
  type: TAlignSequenceType;

  /**
   * Substitution scoring matrix. Residues are sorted lexically by their 1
   * letter code (e.g. "A, C, D, E..." for proteins, not "Ala, Arg, Asn,...")
   * which is a discrepancy with BLOSUM matrices found on NCBI.
   * @default set based on sequence type
   */
  matrix: number[][];

  /**
   * Char that is used as the symbol for gaps in sequences
   * @default '-'
   */
  gapchar: string;
};

type TBioMSAConfig = {
  sequenceType: TAlignSequenceType;
  gapchar: string;
  alignmentMethod: TAlignMethod;
  debug: boolean;
};

const DEFAULT_CONFIG: TBioMSAConfig = {
  sequenceType: 'auto',
  gapchar: DEFAULT_GAP_CHAR,
  alignmentMethod: 'auto',
  debug: false,
};

class BioMSAClass {
  private sequences: TSequence[] = [];
  // STORIES-UNTOLD(B1): typed, as the const enum's member widened to the enum.
  private typeSeq: SEQUENCE_TYPE = SEQUENCE_TYPE.UNSET;
  private config: TBioMSAConfig = {
    ...DEFAULT_CONFIG,
  };

  /**
   * Add sequence to the internal list of sequences of the aligner.
   * All sequences must be of same type.
   * @param {(string|string[])} seq sequence or array of sequences to add
   * @returns void
   * @throws Will throw an exception if sequence added is not of same type as
   *      the other internal sequences.
   */
  private addSequence(seq: string | string[]) {
    // Handle Array parameter overload (recursive)

    if (Array.isArray(seq)) {
      seq.forEach((s) => {
        this.addSequence(s);
      });
      return;
    }

    // Sanity checks

    if (typeof seq !== 'string') {
      throw new TypeError('String type expected for sequences to add.');
    }
    seq = seq.toUpperCase();

    // In sequence type automatic detection mode, some checks are added
    // to infer sequence type and sequence type consistency
    // In user defined mode, we let the user do their own checks and we
    // create the sequences eagerly (unrecognized residues are interpolated
    // randomly)

    if (this.config.sequenceType === 'auto') {
      const type = getSequenceType(seq);

      // First one

      if (this.typeSeq === SEQUENCE_TYPE.UNSET) {
        this.typeSeq = type;

        // Once set, all the others must follow
      } else if (this.typeSeq !== type) {
        throw new Error('All sequences must be of same type.');
      }
    }

    // Add sequence to store

    this.sequences.push(makeSequence(seq, this.typeSeq));
  }

  /**
   * Reset the internal sequence store content
   */
  private reset() {
    this.sequences = [];
    this.typeSeq = SEQUENCE_TYPE.UNSET;
    this.setDefaultConfiguration();
  }

  private setDefaultConfiguration() {
    this.config = {
      ...DEFAULT_CONFIG,
    };
  }

  private setUserConfiguration(opt: Partial<TAlignOpt>) {
    function isValid<K>(value: any, allowed: K[]) {
      return allowed.some((v) => v == value);
    }

    if (opt.method && isValid<TAlignOpt['method']>(opt.method, ['complete', 'diag']))
      this.config.alignmentMethod = opt.method;

    if (opt.type && isValid<TAlignOpt['type']>(opt.type, ['amino', 'nucleic'])) {
      this.config.sequenceType = opt.type;
      this.typeSeq = opt.type === 'amino' ? SEQUENCE_TYPE.PROTEIN : SEQUENCE_TYPE.NUCLEIC;
    }

    if (opt.gapchar?.length === 1) this.config.gapchar = opt.gapchar;

    if (opt.debug !== undefined) this.config.debug = !!opt.debug;
  }

  /**
   * Main alignment function. Takes an array of sequences (strings) as an
   * argument or aligns sequences that were previously added to internal
   * store.
   * This function returns a promise.
   * @param seqArr
   * @returns { Promise<string[]> }
   */
  // STORIES-UNTOLD(B17): `step` (see progressiveAlignment) is awaited between the alignments of a multiple alignment; the
  // executor became an async function whose throws reject (upstream's executor ran the whole alignment in one tick).
  public align(seqArr: string[], opt?: Partial<TAlignOpt>, step?: () => void | Promise<void>) {
    if (DEBUG) Log.start();

    const msa = (async () => {
      // Sequences prerequisites

      if (!Array.isArray(seqArr) || seqArr.some((v) => typeof v !== 'string')) {
        throw 'Array of sequences expected';
      }

      if (seqArr.length < 2) {
        throw 'At least 2 sequences are required.';
      }

      this.reset();
      if (opt) this.setUserConfiguration(opt);

      this.addSequence(seqArr);
      // The singleton's state is replaced by the next call, which may begin while this one waits in `step`.
      const sequences = this.sequences;
      const gapchar = this.config.gapchar;

      if (DEBUG) Log.add('Prepared sequences');

      const lParam = getAlignmentParameters(this.typeSeq, opt);

      if (DEBUG) Log.add('Get sequences type');

      // Set heuristics used:
      // if sequences are long and user does not dislike it, go for
      // fast alignment.

      let doNoAlign = false;
      if (this.config.alignmentMethod === 'auto') doNoAlign = sequences.some((seq) => seq.rawSeq.length > 1600);
      else doNoAlign = this.config.alignmentMethod === 'diag';

      // Start alignment

      let lEStrings: number[][];
      if (sequences.length == 2) {
        if (doNoAlign) {
          lEStrings = noalignPair(sequences[0], sequences[1], lParam);
        } else {
          let lResult = pairwiseAlignment(sequences[0], sequences[1], lParam);
          lEStrings = lResult.estrings;
        }
      } else {
        lEStrings = doNoAlign
          ? centerStarNoAlign(sequences, lParam)
          : await progressiveAlignment(sequences, lParam, step);
      }

      if (DEBUG) Log.summary();

      return applyEstring(sequences, lEStrings, { gapchar });
    })().catch((err) => {
      // Finalize logs collected until this point
      if (DEBUG) Log.summary();
      return Promise.reject(err);
    });

    return msa;
  }
}
const BioMSA = new BioMSAClass();
export default BioMSA;
