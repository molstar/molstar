/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export {
  CARRY_OVER_HEADER,
  carryOverInstruction,
  type CompactionSettings,
  checkedCompaction,
  compactionFor,
  DEFAULT_COMPACTION,
  SIZE_NOTE,
} from './settings.js';
export { type CompactionDomain, type CompactionGroup, requestTextOf, toolkitDomain } from './domain.js';
export {
  type CarryOver,
  checkWishes,
  type CompactionPlan,
  type CompactionReport,
  digestOf,
  editTokens,
  type HeadEntry,
  normalizedForQuotes,
  type PlanInput,
  parseCarryOver,
  planCompaction,
  type Turn,
} from './plan.js';
export { type CompactionOptions, runCompaction } from './run.js';
