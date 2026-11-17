/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The commands' change kinds (_spec/README.md "Commands"): `@molstar/ui`'s `CommandDefinition` has no such field, so the
// viewer keeps a table by command id: whether a command writes a spec patch (undoable, one step), changes session
// state, or only navigates. Menus show what Undo will undo from it; agents tell an edit from a look.

/** What a command changes ("Streaming and patches"): a spec patch, session state, or the viewports only. */
export type ChangeKind = 'patch' | 'session' | 'navigate';

/** The model methods a binding may name: continuous, allocation-free, never commands (no history, no outcome). */
export const MethodActions = ['hover', 'pan', 'zoom', 'drag-select'] as const;
export type MethodAction = (typeof MethodActions)[number];

/**
 * The first slice's commands by id and their change kinds: those the default bindings name and the view's moves. The
 * spec's commands (patches) join with the commands that write them.
 */
export const CommandChangeKinds: Readonly<Record<string, ChangeKind>> = {
  'selection.set': 'session',
  'selection.extend': 'session',
  'selection.toggle': 'session',
  'selection.clear': 'session',
  'view.focus': 'navigate',
  'view.fit': 'navigate',
  'view.zoom-in': 'navigate',
  'view.zoom-out': 'navigate',
  'view.pan-left': 'navigate',
  'view.pan-right': 'navigate',
  'view.show': 'session',
};

/** A command's change kind: the table's, then `extra` (an extension's commands); undefined when unknown. */
export function changeKindOf(id: string, extra: Readonly<Record<string, ChangeKind>> = {}): ChangeKind | undefined {
  return Object.hasOwn(extra, id)
    ? extra[id]
    : Object.hasOwn(CommandChangeKinds, id)
      ? CommandChangeKinds[id]
      : undefined;
}

export function isMethodAction(action: string): action is MethodAction {
  return (MethodActions as readonly string[]).includes(action);
}
