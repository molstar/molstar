/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's codes (_spec/README.md "Bindings", "SpecModel").

import { defineIssues } from '@molstar/core/repr/tree';

export const ViewerIssues = defineIssues('seq.viewer', {
  'seq.viewer.binding.conflict': {
    severity: 'error',
    description: 'Two bindings of one layer take the same gesture, modifiers and target; the later one is used',
    hint: 'Keep one binding for each gesture, modifiers and target',
  },
  'seq.viewer.binding.action': {
    severity: 'error',
    description:
      'A binding names an action the viewer does not have (no model method, no command defined); it is left out',
    hint: 'Name a method (hover, pan, zoom) or the id of a command the viewer defines',
  },
  'seq.viewer.behavior.missing': {
    severity: 'warning',
    description: 'A gesture needs a behavior the viewer spec does not list; the gesture does nothing',
    hint: 'List the behavior in the viewer spec (the default spec has it)',
  },
  'seq.viewer.binding.key': {
    severity: 'warning',
    description: 'A key gesture in the pointer bindings: keys are bound through the keymap',
    hint: 'Give the command its key in the keymap',
  },
  'seq.viewer.highlight.style': {
    severity: 'warning',
    description: 'A highlight names a style that is not `tint`, `band` or `outline`; the channel is a tint',
    hint: 'tint, band or outline',
  },
  'seq.viewer.selection.dropped': {
    severity: 'info',
    description:
      'A host’s update removed nodes the selection named (a unit, a transcript): those parts are dropped from it, as when a patch removes a node',
  },
  'seq.viewer.preview.dropped': {
    severity: 'warning',
    description: 'A preview no longer applies over the new root (an edit landed during the gesture): it is dropped',
    hint: 'Start the gesture again over what the spec is now',
  },
  'seq.viewer.step.stale': {
    severity: 'warning',
    description: 'The story step shown no longer applies whole over the spec: the ops that still apply are shown',
    hint: 'Edit the step, or leave it',
  },
  'seq.viewer.step.frozen': {
    severity: 'error',
    description:
      'An edit while a full-state story step is shown: the step does not follow the spec, so it would show nothing',
    hint: 'Edit the step itself, or leave it to edit the spec',
  },
  'seq.viewer.step.unknown': {
    severity: 'error',
    description: 'No story step has that id: nothing was entered',
    hint: 'Name a step of the spec (`/snapshots/<id>`)',
  },
  'seq.viewer.behavior.shape': {
    severity: 'warning',
    description: 'A behavior listed under a built-in’s name lacks a method the viewer calls: its gestures do nothing',
    hint: 'Give the provider a name of its own, or a handle with the built-in behavior’s methods',
  },
  'seq.viewer.navigators.none': {
    severity: 'warning',
    description: 'The view declares navigators and the viewer lists none: they are not drawn, the ruler rails stay',
    hint: 'List `NavigatorsPart` of `@molstar/seqx-viewer/navigators` in the viewer spec (`navigators`), as the default spec does',
  },
  'seq.viewer.card.floor': {
    severity: 'warning',
    description:
      'The card’s providers leave out `seq:item` or `seq:place`: a hovered mark may show no card, the hover contract’s floor',
    hint: 'Keep `seq:item` and `seq:place` in the registry (the default spec lists every built-in provider)',
  },
  'seq.viewer.card.provider': {
    severity: 'warning',
    description: 'A card provider threw: its rows are left out of the card, the rest shows',
    hint: 'Fix the provider; the card is composed without it meanwhile',
  },
  'seq.viewer.step.none': {
    severity: 'error',
    description: 'The viewer lists no story part: it enters no story step',
    hint: 'List `StoryPart` of `@molstar/seqx-viewer/story` in the viewer spec (`story`), as the default spec does',
  },
  'seq.viewer.edit.masked': {
    severity: 'info',
    description:
      'The edit sets a key the story step shown sets too: it lands in the spec and shows once the step is left',
    hint: 'Leave the step to see it',
  },
  'seq.viewer.rebase.base': {
    severity: 'error',
    description: 'A stale patch cannot be rebased: the revision it was made against is no longer known',
    hint: 'Make the patch again against the current revision',
  },
  'seq.viewer.rebase.conflict': {
    severity: 'error',
    description: 'A stale patch conflicts with what changed since it was made: nothing was applied',
    hint: 'Review the conflicts, or make the patch again against the current revision',
  },
});
