/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `seq.story.*` codes (_spec/README.md "Stories"): what a step's state reports.

import { defineIssues } from '@molstar/core/repr/tree';

// pure: a composition with no story part bundles none of these codes
export const StoryIssues = /* @__PURE__ */ defineIssues('seq.story', {
  'seq.story.step.unknown': {
    severity: 'error',
    description: 'No story step has this id',
    hint: 'Name a step of the `snapshots` container',
  },
  'seq.story.step.stale': {
    severity: 'warning',
    description: 'Some of a step’s ops no longer apply (a target gone, a sibling renamed, a move now refused)',
    hint: 'The step shows the ops that still apply: edit the step, or store it again from what it should show',
  },
  'seq.story.step.from': {
    severity: 'error',
    description: 'A step starts from a step that is not there, or a chain of starts leads back to itself',
    hint: 'Name an existing step as `from`, or drop it: the step then applies over the spec',
  },
  'seq.story.step.chain': {
    severity: 'warning',
    description:
      'A step is more deltas from its keyframe than `story.maxChain`: reaching it costs more applies than the bound',
    hint: 'Store the story again (the encoder keeps every step within the bound), or store this step over the spec',
  },
  'seq.story.step.form': {
    severity: 'error',
    description: 'A step holds neither a patch nor a full state it can read',
    hint: 'Give the step a `patch` (no base) or a `state`',
  },
});
