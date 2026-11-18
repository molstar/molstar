/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `cursor` behavior (_spec/README.md "Behaviors"): the keyboard cursor (the interaction grammar): the step keys
// (`cursor.keys`) move it a position along its section's axis in display order, or to the section above or below (its
// position carried along the view's link); Shift with a step extends the range the anchor is in (from the cursor's cell
// when there is no anchor); Cmd/Ctrl+Shift with a step starts a new range at the cursor, the others kept; Space toggles
// the cursor's cell; Enter focuses the selection, else the cursor's cell. A cursor moved out of view pans its section
// (`cursor.reveal`). Every move is a command.

import { EmptyElements, rangesLocation } from '@molstar/seqx-core/query';
import { projectOntoView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import {
  type AddressModel,
  axisTotal,
  type Cursor,
  offsetOf,
  positionAt,
  type SelectionModel,
  type SpecModel,
  type ViewModel,
} from '@molstar/seqx-viewer/models';

/** The section a key acts on when there is no cursor yet (the one under the pointer, else the first). */
export interface CursorContext {
  readonly section?: string;
}

export interface CursorBehavior {
  readonly commands: readonly CommandDefinition<CursorContext>[];
}

export function cursorBehavior(
  addresses: AddressModel,
  selection: SelectionModel,
  view: ViewModel,
  spec: SpecModel,
  config: ViewerConfig = DefaultViewerConfig,
): CursorBehavior {
  const keys = config.cursor.keys;
  const sectionOf = (address: string | undefined) => view.resolved.value?.sections.find((s) => s.address === address);

  /** The cursor now, or one placed at the first position the section shows. */
  const current = (ctx: CursorContext): Cursor | undefined => {
    const c = selection.cursor.value;
    if (c && sectionOf(c.section)) return c;
    const section = sectionOf(ctx.section) ?? view.resolved.value?.sections[0];
    const v = section && view.get(section.address);
    if (!section || !v || v.pending) return undefined;
    const p = positionAt(section.axis, Math.ceil(v.start));
    return p === undefined
      ? undefined
      : { section: section.address, system: section.axis.system, position: Math.floor(p) };
  };

  /**
   * A section's viewport pans to show its offset `o` when it is out of view. A follower's does not (it would move its
   * leader): the cursor on a follower keeps the rule it had before followers led, revealing only once detached.
   */
  const reveal = (section: string, o: number) => {
    if (!config.cursor.reveal || !view.canMove(section)) return;
    if (view.isFollowing(section) && config.follow.followerMoves !== 'detach') return;
    const v = view.get(section);
    if (!v || v.pending) return;
    if (o < v.start) view.pan(section, o - v.start);
    else if (o + 1 > v.end) view.pan(section, o + 1 - v.end);
  };

  /** The cursor a step along its axis, in display order, clamped to the axis. */
  const stepped = (c: Cursor, delta: number): Cursor | undefined => {
    const axis = sectionOf(c.section)?.axis;
    const total = axis && axisTotal(axis);
    const o = axis && offsetOf(axis, c.position);
    if (!axis || total === undefined || o === undefined) return undefined;
    const next = Math.min(Math.max(Math.floor(o) + delta, 0), total - 1);
    const p = positionAt(axis, next + 0.5);
    return p === undefined ? undefined : { ...c, position: Math.floor(p) };
  };

  const moveTo = (c: Cursor) => {
    selection.setCursor({ ...c, keyboard: true });
    const axis = sectionOf(c.section)?.axis;
    const o = axis && offsetOf(axis, c.position);
    if (o !== undefined) reveal(c.section, Math.floor(o));
  };

  /** The cursor on the section `delta` above or below, its position carried along the view's link. */
  const across = (c: Cursor, delta: number): Cursor | undefined => {
    const resolved = view.resolved.value;
    if (!resolved) return undefined;
    const i = resolved.sections.findIndex((s) => s.address === c.section);
    const to = resolved.sections[i + delta];
    if (i < 0 || !to) return undefined;
    if (to.axis.system === c.system) return { section: to.address, system: c.system, position: c.position };
    const part = projectOntoView(
      spec.model.value,
      resolved,
      rangesLocation(c.system, [c.position, c.position + 1]),
    ).parts.find((p) => p.system === to.axis.system);
    if (part?.ranges.length) return { section: to.address, system: to.axis.system, position: part.ranges[0] };
    // nothing carries it: the first position the section shows
    const v = view.get(to.address);
    const p = v && !v.pending ? positionAt(to.axis, Math.ceil(v.start)) : undefined;
    return p === undefined ? undefined : { section: to.address, system: to.axis.system, position: Math.floor(p) };
  };

  const placed = (ctx: CursorContext): true | string => (current(ctx) ? true : 'No section to move the cursor in');

  const move = (
    id: string,
    title: string,
    key: string,
    to: (c: Cursor) => Cursor | undefined,
    select: 'move' | 'extend' | 'add' = 'move',
  ): CommandDefinition<CursorContext> => ({
    id,
    title,
    category: 'Cursor',
    shortcut: {
      key,
      ...(select !== 'move' ? { shift: true } : {}),
      ...(select === 'add' ? { mod: true } : {}),
    },
    ...(select === 'move' ? { scope: 'canvas' } : {}),
    when: placed,
    run: (_, ctx): CommandResult => {
      const c = current(ctx)!;
      // a cursor not yet shown is placed first: its first key shows it where it is
      const next = selection.cursor.value ? to(c) : c;
      if (!next) return { message: '' };
      if (select !== 'move') {
        const cell = rangesLocation(c.system, [c.position, c.position + 1]);
        const anchor = { system: c.system, start: c.position, end: c.position + 1 };
        const a = selection.anchor;
        // a new range at the cursor, the others kept; else from the anchor (none on this system: the cursor's cell)
        if (select === 'add') selection.add(cell, anchor);
        else if (!a || a.system !== c.system)
          selection.set({ elements: EmptyElements, location: cell, shown: cell }, anchor);
        selection.extend({ system: next.system, position: next.position }, sectionOf(next.section)?.axis);
      }
      moveTo(next);
      return { message: select === 'move' ? 'Moved' : select === 'add' ? 'Added' : 'Extended' };
    },
  });

  const commands: CommandDefinition<CursorContext>[] = [
    move('cursor.prev', 'Cursor back', keys.prev, (c) => stepped(c, -1)),
    move('cursor.next', 'Cursor on', keys.next, (c) => stepped(c, 1)),
    move('cursor.up', 'Cursor to the section above', keys.up, (c) => across(c, -1)),
    move('cursor.down', 'Cursor to the section below', keys.down, (c) => across(c, 1)),
    move('cursor.extend-prev', 'Extend the selection back', keys.prev, (c) => stepped(c, -1), 'extend'),
    move('cursor.extend-next', 'Extend the selection on', keys.next, (c) => stepped(c, 1), 'extend'),
    move('cursor.add-prev', 'Start a new range back from the cursor', keys.prev, (c) => stepped(c, -1), 'add'),
    move('cursor.add-next', 'Start a new range on from the cursor', keys.next, (c) => stepped(c, 1), 'add'),
    {
      id: 'cursor.toggle',
      title: 'Toggle the cursor’s cell in the selection',
      category: 'Cursor',
      shortcut: { key: ' ' },
      scope: 'canvas',
      when: (ctx) => (selection.cursor.value && current(ctx) ? true : 'No cursor'),
      run: (_, ctx): CommandResult => {
        const c = current(ctx)!;
        selection.toggle({ system: c.system, position: c.position }, sectionOf(c.section)?.axis);
        return { message: 'Toggled' };
      },
    },
    {
      id: 'cursor.focus',
      title: 'Focus the selection, else the cursor',
      category: 'Cursor',
      shortcut: { key: 'Enter' },
      scope: 'canvas',
      when: () => (selection.cursor.value || selection.value.value.shown.parts.length ? true : 'Nothing to focus'),
      run: (): CommandResult => {
        // the selection, else the cursor's cell, focused as an address is: the leading sections move (a cursor on a
        // follower focuses its leader), each to the hull of what it shows of it
        const c = selection.cursor.value;
        const shown = selection.value.value.shown;
        const loc = shown.parts.length ? shown : c ? rangesLocation(c.system, [c.position, c.position + 1]) : undefined;
        if (!loc) return { message: '' };
        const r = addresses.focusLocation(loc);
        return { message: r === 'moved' ? 'Focused' : r === 'pending' ? 'Kept until the view is placed' : '' };
      },
    },
  ];

  return { commands };
}
