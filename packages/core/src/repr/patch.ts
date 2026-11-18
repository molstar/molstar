/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/core/repr/patch: revisions, patches and diff (_spec/README.md "Patches"). Revisions (immutable roots, named
// by their hash) are the truth; a patch is derived (`diff`) or written, and applying one gives a new revision with its
// inverse and the net `Changes`. The entry of `repr/patch/` (as `data/int.ts` is of `data/int/`), reached through the
// package's wildcard export; its conformance suite is `@molstar/core/repr/patch/conformance`.

export { applyPatch } from './patch/apply.js';
export { diff } from './patch/diff.js';
export { AnyPatchLanguage, patchFromNode, patchLanguage, patchToNode } from './patch/language.js';
export {
  type ApplyOptions,
  type Changes,
  type InsertOp,
  type MoveOp,
  type Op,
  type Patch,
  PatchIssues,
  type PatchResult,
  type PatchValidator,
  type Place,
  type RemoveOp,
  type ReplaceOp,
  type SetOp,
  type UnsetOp,
  type Where,
} from './patch/patch.js';
