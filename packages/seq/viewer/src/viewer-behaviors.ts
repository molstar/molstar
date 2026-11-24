/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { State } from '@molstar/core/reactive';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { CladeActions, GroupActions, SelectActions } from '@molstar/seqx-viewer/bindings';
import { ViewerIssues } from './issues.js';
import type { ViewerCommand } from './viewer.js';
import type { BehaviorContext, SeqViewerSpec } from './viewer-spec.js';

// The viewer's behaviors (_spec/README.md "Composition"): made from the spec's list after the models, one of each name,
// their commands and disposers collected; looked up by name with the methods the viewer and its input call; a gesture
// whose behavior is missing or of another shape said once as an issue (`behaviorIssues`).

/** The actions a behavior's bindings may name, by its name: bound only when it is listed. */
const BehaviorActions: Readonly<Record<string, readonly string[]>> = {
  select: SelectActions,
  rows: CladeActions,
  groups: GroupActions,
};

/**
 * The methods the viewer and its input call on a built-in behavior's handle, by the behavior's name: a handle listed
 * under that name without one of them is refused (`seq.viewer.behavior.shape`), never called.
 */
const BehaviorShapes: Readonly<Record<string, readonly string[]>> = {
  hover: ['move', 'leave'],
  navigate: ['drag', 'wheel', 'pinch'],
  select: ['click', 'beginDrag', 'dragTo', 'endDrag', 'cancelDrag'],
  window: ['begin', 'dragTo', 'end', 'cancel'],
  navigator: ['center', 'zoomTo', 'press', 'dragTo', 'release', 'cancel', 'click', 'doubleClick'],
  rows: ['wheel', 'lay', 'railGesture'],
  addressing: ['runLink'],
};

/**
 * The actions of the behaviors listed (the selection's, the rows'), by their names: a behavior's actions are bound
 * only with the behavior that runs them, else their bindings are reported, never dead.
 */
export function listedActions(behaviors: SeqViewerSpec['behaviors']): string[] {
  return behaviors.flatMap((b) => BehaviorActions[b.provider.name] ?? []);
}

/**
 * Makes the behaviors listed, after the models, one of each name: each handle by its provider's name, its commands,
 * its disposer registered. Throws when a name is listed twice (what was made is the caller's to retire).
 */
export function createBehaviors(
  behaviors: SeqViewerSpec['behaviors'],
  context: BehaviorContext,
  handles: Map<string, any>,
  commands: ViewerCommand[],
  register: (dispose: () => void) => void,
): void {
  for (const b of behaviors) {
    if (handles.has(b.provider.name)) throw new Error(`The behavior ${b.provider.name} is listed twice`);
    const handle = b.provider.create(context, b.params);
    handles.set(b.provider.name, handle);
    if (Array.isArray(handle?.commands)) commands.push(...handle.commands);
    if (typeof handle?.dispose === 'function') register(() => handle.dispose());
  }
}

/**
 * The handle of a behavior by its name. One not listed, or one of a built-in's name whose handle lacks a method the
 * viewer calls (`BehaviorShapes`: a host's provider under that name with another shape), is reported once as an issue
 * and undefined is returned: the gesture does nothing.
 */
export function behaviorOf<H>(
  handles: ReadonlyMap<string, any>,
  name: string,
  issues: State<readonly Issue[]>,
): H | undefined {
  const handle = handles.get(name);
  if (handle !== undefined) {
    // a loop, no closure: the input asks on every wheel and pointer event
    const shape = BehaviorShapes[name];
    let lacks: string | undefined;
    for (let i = 0; shape && i < shape.length && lacks === undefined; i++)
      if (typeof handle[shape[i]] !== 'function') lacks = shape[i];
    if (lacks === undefined) return handle as H;
    behaviorIssue(
      issues,
      name,
      'seq.viewer.behavior.shape',
      `The viewer's ${name} behavior has no ${lacks}: its gestures do nothing`,
    );
    return undefined;
  }
  behaviorIssue(
    issues,
    name,
    'seq.viewer.behavior.missing',
    `The viewer has no ${name} behavior: its gestures do nothing`,
  );
  return undefined;
}

/** A behavior's issue, once by its name. */
function behaviorIssue(
  issues: State<readonly Issue[]>,
  name: string,
  code: 'seq.viewer.behavior.missing' | 'seq.viewer.behavior.shape',
  message: string,
): void {
  if (issues.value.some((i) => i.address === name)) return;
  issues.set([
    ...issues.value,
    issueOf(ViewerIssues, code, message, { address: name, hint: ViewerIssues.codes[code].hint }),
  ]);
}
