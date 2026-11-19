/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `cursor` behavior (_spec/README.md "Behaviors"): the keyboard cursor (the interaction grammar): the step keys
// (`cursor.keys`) move it a position along its section's axis in display order, or a line up or down in a layout of
// lines (a text: the section's layout's step), else to the section above or below (its position carried along the
// view's link); Shift with a step extends the range the anchor is in (from the cursor's cell when there is no anchor),
// up and down by a line; Cmd/Ctrl+Shift with a step starts a new range at the cursor, the others kept; Space toggles
// the cursor's cell; Enter focuses the selection, else the cursor's cell. A cursor moved out of view brings its section
// to it by the layout's least move (`reveal.cursor`). Every move is a command.

import { EmptyElements, rangesLocation } from '@molstar/seqx-core/query';
import { projectOntoView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import {
  type AddressModel,
  axisBreaks,
  axisTotal,
  type Cursor,
  cellOffset,
  elementAt,
  HiddenMessage,
  hiddenLength,
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
    const p = elementAt(section.axis, Math.ceil(v.start));
    return p === undefined ? undefined : { section: section.address, system: section.axis.system, position: p };
  };

  /**
   * A section's viewport pans to show its offset `o` when it is out of view. A follower's does not (it would move its
   * leader): the cursor on a follower keeps the rule it had before followers led, revealing only once detached.
   */
  const reveal = (section: string, o: number) => {
    if (!config.reveal.cursor || !view.canMove(section)) return;
    if (view.isFollowing(section) && config.follow.followerMoves !== 'detach') return;
    const v = view.get(section);
    if (!v || v.pending) return;
    view.reveal(section, o, o + 1);
  };

  /** The cursor a line up or down within its section (a layout of lines), else undefined. */
  const line = (c: Cursor, dir: 'up' | 'down'): Cursor | undefined => {
    const axis = sectionOf(c.section)?.axis;
    const o = axis && cellOffset(axis, c.position);
    if (!axis || o === undefined) return undefined;
    const next = view.step(c.section, o, dir);
    const p = next < 0 ? undefined : elementAt(axis, next);
    return p === undefined ? undefined : { ...c, position: p };
  };

  /**
   * The cursor a step along its axis, in display order (whatever the direction), clamped to the axis; a step from one
   * piece to the next goes over the break between them (_spec/README.md "Behaviors", Cursor).
   */
  const stepped = (c: Cursor, delta: number): Cursor | undefined => {
    const axis = sectionOf(c.section)?.axis;
    // a section of blocks on several systems: its units are no system's positions (stepping across blocks is open)
    if (axis?.pieces === 'blocks' && axis.system !== c.system) return undefined;
    const total = axis && axisTotal(axis);
    const o = axis && cellOffset(axis, c.position);
    if (!axis || total === undefined || o === undefined) return undefined;
    const next = Math.min(Math.max(o + delta, 0), total - 1);
    const p = elementAt(axis, next);
    return p === undefined ? undefined : { ...c, position: p };
  };

  /**
   * What a step from `a` to `b` crossed: the hidden stretch of the break between them, named for the announcement
   * ("crossed a break, 1,240 positions hidden"); empty when it crossed none.
   */
  const crossed = (a: Cursor, b: Cursor): string => {
    const axis = sectionOf(b.section)?.axis;
    if (!axis || a.section !== b.section) return '';
    const oa = cellOffset(axis, a.position);
    const ob = cellOffset(axis, b.position);
    if (oa === undefined || ob === undefined) return '';
    const breaks = axisBreaks(axis);
    let hidden = 0;
    let n = 0;
    for (let i = 0; i + 2 < breaks.length; i += 3) {
      const at = breaks[i];
      if (at > Math.min(oa, ob) && at <= Math.max(oa, ob)) {
        hidden += hiddenLength(breaks[i + 1], breaks[i + 2], axis.length);
        n++;
      }
    }
    if (!n) return '';
    const what = n === 1 ? 'a break' : `${n} breaks`;
    return `, crossed ${what}, ${hidden.toLocaleString('en-US')} positions hidden`;
  };

  const moveTo = (c: Cursor) => {
    selection.setCursor({ ...c, keyboard: true });
    const axis = sectionOf(c.section)?.axis;
    const o = axis && cellOffset(axis, c.position);
    if (o !== undefined) reveal(c.section, o);
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
    const p = v && !v.pending ? elementAt(to.axis, Math.ceil(v.start)) : undefined;
    return p === undefined ? undefined : { section: to.address, system: to.axis.system, position: p };
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
      const verb = select === 'move' ? 'Moved' : select === 'add' ? 'Added' : 'Extended';
      return { message: `${verb}${crossed(c, next)}` };
    },
  });

  const commands: CommandDefinition<CursorContext>[] = [
    move('cursor.prev', 'Cursor back', keys.prev, (c) => stepped(c, -1)),
    move('cursor.next', 'Cursor on', keys.next, (c) => stepped(c, 1)),
    move('cursor.up', 'Cursor a line up, or to the section above', keys.up, (c) => line(c, 'up') ?? across(c, -1)),
    move(
      'cursor.down',
      'Cursor a line down, or to the section below',
      keys.down,
      (c) => line(c, 'down') ?? across(c, 1),
    ),
    move('cursor.extend-up', 'Extend the selection a line up', keys.up, (c) => line(c, 'up'), 'extend'),
    move('cursor.extend-down', 'Extend the selection a line down', keys.down, (c) => line(c, 'down'), 'extend'),
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
        const focused = addresses.movedToBreak ? HiddenMessage : 'Focused';
        return { message: r === 'moved' ? focused : r === 'pending' ? 'Kept until the view is placed' : '' };
      },
    },
  ];

  return { commands };
}
