/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The commands' change kinds (_spec/README.md "Commands"): `@molstar/ui`'s `CommandDefinition` has no such field, so the
// viewer keeps a table by command id: whether a command writes a spec patch (undoable, one step), changes session
// state, or only navigates. Menus show what Undo will undo from it; agents tell an edit from a look.

/**
 * What a command changes (_spec/README.md "Commands"): a spec patch (undoable, one step), a gesture's preview (ends in
 * one patch or nothing, never a step itself), session state, what is shown without changing the spec (the viewports, a
 * story step), or nothing (`none`: it reads, as copying a link does).
 */
export type ChangeKind = 'patch' | 'preview' | 'session' | 'navigate' | 'none';

/** Every change kind, in the order above. */
export const ChangeKinds: readonly ChangeKind[] = ['patch', 'preview', 'session', 'navigate', 'none'];

/** The model methods a binding may name: continuous, allocation-free, never commands (no history, no outcome). */
export const MethodActions = ['hover', 'pan', 'zoom', 'drag-select'] as const;
export type MethodAction = (typeof MethodActions)[number];

/**
 * The commands by id and their change kinds: those the default bindings name, the view's moves, the actions on
 * addresses and the edits of the spec.
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
  // the actions on addresses: a focus or a reveal moves the view; a selection or a highlight is session state
  'address.focus': 'navigate',
  'address.reveal': 'navigate',
  'address.select': 'session',
  'address.highlight': 'session',
  'address.hover': 'session',
  'address.unhover': 'session',
  // copying a link reads the session and changes nothing
  'address.copy-link': 'none',
  // a bookmark is a ref in the spec: one patch, one undo step
  'refs.add': 'patch',
  // the editing behavior's: each one patch
  'layer.hide': 'patch',
  'layer.show': 'patch',
  'layer.color': 'patch',
  'track.hide': 'patch',
  'track.show': 'patch',
  'track.move': 'patch',
  'refs.remove': 'patch',
  // a story step is navigation: entered and left, never undone
  'step.enter': 'navigate',
  'step.leave': 'navigate',
  // Esc ends a gesture's preview in nothing
  'preview.cancel': 'preview',
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
