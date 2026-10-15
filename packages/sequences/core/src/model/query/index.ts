/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seq-core/model/query: the SeqQL runtime (seqspec:seqql#the-selection, seqspec:seqql#extension): the sequence
// selection and its run kernels, the query context with its limits, the compiler and the base runtimes. The language
// is @molstar/seq-core/query's; the model imports it from this folder only, so code that never queries does not load it.
export {
  type CompileContext,
  type CompiledQuery,
  compileQuery,
  createQueryEnv,
  type QueryEnv,
  type QueryExtension,
  type QueryFn,
  type QueryResult,
  type SymbolRuntime,
} from './compile.js';
export {
  DEFAULT_LIMITS,
  F64,
  QueryAbort,
  QueryContext,
  type QueryLimits,
  type QueryOptions,
  type QuerySlot,
} from './context.js';
export { type BuildGuard, SelectionBuilder } from './kernels.js';
export { type RangesInputOptions, SeqSelection, type ToLociOptions } from './selection.js';
