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
export const MethodActions = [
  'hover',
  'pan',
  'zoom',
  'drag-select',
  'drag-add',
  'drag-subtract',
  'drag-extend',
] as const;
export type MethodAction = (typeof MethodActions)[number];

/**
 * The commands by id and their change kinds: those the default bindings name, the view's moves, the actions on
 * addresses and the edits of the spec.
 */
export const CommandChangeKinds: Readonly<Record<string, ChangeKind>> = {
  'selection.set': 'session',
  'selection.extend': 'session',
  'selection.toggle': 'session',
  'selection.add': 'session',
  'selection.subtract': 'session',
  'selection.clear': 'session',
  // copying reads: the letters go to the clipboard
  'selection.copy-letters': 'none',
  // the keys panel opens and closes: the reader's, nothing in the models
  'help.keys': 'none',
  // the keyboard cursor and what it selects are session state (every other `cursor.*`: `changeKindOf`); Enter's focus
  // moves the view
  // Home and End scroll the scroller under the pointer and move no cursor: a scroll is a viewport's (a band's top, a
  // section's `top`), so a navigation, as `rows.page-*` (the page keys are the cursor's: a shown cursor moves with them)
  'cursor.scroll-start': 'navigate',
  'cursor.scroll-end': 'navigate',
  'cursor.focus': 'navigate',
  'view.focus': 'navigate',
  'view.fit': 'navigate',
  // the view's history of jumps is session state's kind: it moves the view, never the spec's history
  'view.back': 'navigate',
  'view.forward': 'navigate',
  'view.zoom-in': 'navigate',
  'view.zoom-out': 'navigate',
  'view.pan-left': 'navigate',
  'view.pan-right': 'navigate',
  'view.pan-up': 'navigate',
  'view.pan-down': 'navigate',
  // a window's steps as a slider (its section's viewport)
  'window.pan-back': 'navigate',
  'window.pan-on': 'navigate',
  'window.grow': 'navigate',
  'navigator.center': 'navigate',
  'navigator.zoom-to': 'navigate',
  'window.shrink': 'navigate',
  // the actions on addresses: a focus or a reveal moves the view; a selection or a highlight is session state
  'address.focus': 'navigate',
  'address.reveal': 'navigate',
  'address.select': 'session',
  'address.highlight': 'session',
  'address.hover': 'session',
  'address.unhover': 'session',
  // clearing a highlight channel (or all) is session state, as setting one is
  'highlight.clear': 'session',
  // copying a link reads the session and changes nothing
  'address.copy-link': 'none',
  // a bookmark is a ref in the spec: one patch, one undo step
  'refs.add': 'patch',
  'selection.bookmark': 'patch',
  // the editing behavior's: each one patch
  // a section's direction and pieces are spec params: one patch, one undo step each
  'section.reverse': 'patch',
  'section.collapse': 'patch',
  // a block of a row: the next or the previous one a jump; its direction and the section's `include` a patch each
  'block.next': 'navigate',
  'block.previous': 'navigate',
  'block.reverse': 'patch',
  'section.include': 'patch',
  'block.form': 'patch',
  'layer.hide': 'patch',
  'layer.show': 'patch',
  'layer.color': 'patch',
  // the legend behavior's: a theme edit is one patch; a legend's selection and hover, and the reader's substitution,
  // are session state
  'layer.channel': 'patch',
  'scale.set': 'patch',
  'theme.value-color': 'patch',
  'layer.filter': 'patch',
  'legend.select': 'session',
  'legend.hover': 'session',
  'theme.cvd-safe': 'session',
  'track.hide': 'patch',
  'track.show': 'patch',
  'track.move': 'patch',
  'track.add': 'patch',
  'refs.remove': 'patch',
  'refs.label': 'patch',
  // a story step is navigation: entered and left, never undone
  'step.enter': 'navigate',
  'step.leave': 'navigate',
  // an alignment's rows: folds are session state (a bulk fold a jump), pins, hidden rows, the order and the kept folds
  // spec edits, the rows' zoom a move
  'clades.toggle': 'session',
  'clades.collapseAll': 'session',
  'clades.expandAll': 'session',
  'clades.collapseTo': 'session',
  'clades.reveal': 'session',
  'clades.keep': 'patch',
  // a section's track groups: folds are session state (a toggle a hand's move, a bulk fold a jump), the kept folds a
  // spec edit
  'groups.toggle': 'session',
  'groups.collapseAll': 'session',
  'groups.expandAll': 'session',
  'groups.keep': 'patch',
  'rows.pin': 'patch',
  'rows.hide': 'patch',
  'rows.show': 'patch',
  'rows.order': 'patch',
  'rows.order-next': 'patch',
  'rows.unsort': 'patch',
  'rows.zoom-in': 'navigate',
  'rows.zoom-out': 'navigate',
  'rows.page-down': 'navigate',
  'rows.page-up': 'navigate',
  'rows.tree': 'session',
  // Esc ends a gesture's preview in nothing
  'preview.cancel': 'preview',
  // a track header's menu opens: the reader's, nothing in the models
  'headers.menu': 'none',
  // a bridge's (`bridge-view.ts`, lazy): an align is a jump, select across a selection's, the lock toggle an edit
  'bridge.align': 'navigate',
  'bridge.selectAcross': 'session',
  'bridge.lock': 'patch',
  // the `search` behavior's: a run and a clear are session state, and so is a step (the current hit is its channel; the
  // reveal it runs is a jump of the view's history, as `clades.reveal`'s), a save one patch
  'search.run': 'session',
  'search.next': 'session',
  'search.previous': 'session',
  'search.go': 'session',
  'search.clear': 'session',
  'search.save': 'patch',
};

/** The viewer's own cursor commands that are session state (an extension's `cursor.x` declares its kind). */
const CursorMoves =
  /^cursor\.((prev|next|up|down)(-feature)?|(extend|add)-(prev|next|up|down)|page-(up|down)|(first|last)-row|toggle)$/;

/** A command's change kind: the table's, then `extra` (an extension's commands); undefined when unknown. */
export function changeKindOf(id: string, extra: Readonly<Record<string, ChangeKind>> = {}): ChangeKind | undefined {
  return Object.hasOwn(extra, id)
    ? extra[id]
    : Object.hasOwn(CommandChangeKinds, id)
      ? CommandChangeKinds[id]
      : // the cursor's moves are session state (the ones above, the scrolls and the focus, navigate)
        CursorMoves.test(id)
        ? 'session'
        : undefined;
}

export function isMethodAction(action: string): action is MethodAction {
  return (MethodActions as readonly string[]).includes(action);
}
