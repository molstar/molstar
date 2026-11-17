/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The list's rows (commands §5.1): a group to open, a command (or one of its ways in), a prefix's item. Each is made
// once; a `watch` writes what changes while it shows (the highlight, a count, the query's marks).

import type { ReadonlyState } from '@molstar/core/reactive';
import type { CommandManager } from '../../commands/manager.js';
import { matchPositions } from '../../commands/match.js';
import type { CommandMatch } from '../../commands/query.js';
import type { CommandGroup, CommandMenuItem } from '../../commands/types.js';
import { setAttr, setText, type WatchSource, watch } from '../../dom/index.js';
import { Icon } from '../icon/icon.js';
import { classes } from './_generated/command-menu.module.css.js';
import { Keys } from './footer.js';

export type Row =
  | { kind: 'caption'; label: string }
  | { kind: 'group'; group: CommandGroup; index: number }
  /** `index`: its place among the pickable rows. */
  | { kind: 'match'; match: CommandMatch<any>; tail?: string; index: number }
  | { kind: 'item'; item: CommandMenuItem; index: number };
export type Pickable = Exclude<Row, { kind: 'caption' }>;

/** A row's height (`--msui-row`), for the scroll height of the rows not rendered yet. */
export const RowPx = 24;

/** A row's identity across lists: what it shows and where. */
export function rowKey(r: Row): string {
  switch (r.kind) {
    case 'caption':
      return `caption:${r.label}`;
    case 'group':
      return `group:${r.group.id}:${r.index}`;
    case 'item':
      return `item:${r.item.id}:${r.item.label}:${r.index}`;
    case 'match':
      return `match:${r.match.command.id}:${r.match.preset?.id ?? ''}:${r.index}:${r.tail ?? ''}`;
  }
}

/** The same row as before: its view shows the same (the query's marks are its own `watch`'s). */
export function sameRow(a: Row, b: Row): boolean {
  if (a.kind === 'caption' || b.kind === 'caption')
    return a.kind === 'caption' && b.kind === 'caption' && a.label === b.label;
  if (a.kind !== b.kind || a.index !== b.index) return false;
  if (a.kind === 'group') return a.group === (b as typeof a).group;
  if (a.kind === 'item') return a.item === (b as typeof a).item;
  const m = b as typeof a;
  return (
    a.match.command === m.match.command &&
    a.match.preset === m.match.preset &&
    a.match.blocked === m.match.blocked &&
    a.tail === m.tail
  );
}

/** What every option row is given: its id, whether it is the highlighted one, what it does. */
export interface OptionProps {
  id: string;
  /** Told only when it becomes or stops being the highlighted one (a move renders two rows, not all). */
  selected: ReadonlyState<boolean>;
  onPick(): void;
  onHover(): void;
}

/** The highlight as an option's render writes it: `aria-selected`, and into view when it becomes the highlighted one. */
export function writeSelected(el: HTMLElement, selected: boolean, was: boolean): boolean {
  setAttr(el, 'aria-selected', String(selected));
  if (selected && !was) queueMicrotask(() => el.scrollIntoView?.({ block: 'nearest' }));
  return selected;
}

export function GroupRow(props: OptionProps & { manager: CommandManager<any>; group: CommandGroup }) {
  const path = props.manager.groupPath(props.group.parent).map((g) => g.title);
  const el = (
    <div
      id={props.id}
      class={classes.row}
      role="option"
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowSlot}>
        <Icon name={props.group.icon ?? 'list-tree'} />
      </span>
      <span class={classes.rowTitle}>{props.group.title}</span>
      <span class={classes.rowMeta}>{path.join(' › ')}</span>
      <span class={classes.rowEnd}>
        <Icon name="chevron-right" />
      </span>
    </div>
  ) as HTMLElement;
  let was = false;
  watch([props.selected], () => {
    was = writeSelected(el, props.selected.value, was);
  });
  return el;
}

export function ItemRow(props: OptionProps & { item: CommandMenuItem }) {
  const { item } = props;
  const el = (
    <div
      id={props.id}
      class={classes.row}
      role="option"
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowSlot}>{item.icon ? <Icon name={item.icon} /> : null}</span>
      <span class={classes.rowTitle}>{item.label}</span>
      <span class={classes.rowMeta}>{item.detail ?? ''}</span>
      <span />
    </div>
  ) as HTMLElement;
  let was = false;
  watch([props.selected], () => {
    was = writeSelected(el, props.selected.value, was);
  });
  return el;
}

/**
 * A command's row: icon, title with the query's letters marked (and the typed tail after a direct alias), the reason
 * it cannot run or its path or category (highlighted: how much it would do, its count), its key after the keymap.
 */
export function CommandRow(
  props: OptionProps & {
    manager: CommandManager<any>;
    row: Extract<Row, { kind: 'match' }>;
    /** Its group path below the group being browsed, shown instead of the category when it has one. */
    path: readonly string[];
    /** Browsing (not searching): the captions say the categories, so a row without a deeper path shows nothing more. */
    inGroup?: boolean;
    query: ReadonlyState<string>;
  },
) {
  const { manager, query } = props;
  const { match, tail } = props.row;
  const c = match.command;
  const title = match.preset ? `${c.title.replace(/…$/, '')} ${match.preset.title}` : c.title;
  const chord = manager.shortcut(c);
  const meta = match.blocked ?? (props.path.length ? props.path.join(' › ') : props.inGroup ? '' : c.category);
  const titleText = document.createElement('span');
  const metaText = document.createTextNode(meta);
  const el = (
    <div
      id={props.id}
      class={classes.row}
      role="option"
      aria-disabled={match.blocked ? 'true' : undefined}
      onClick={() => props.onPick()}
      onPointerMove={() => props.onHover()}
    >
      <span class={classes.rowSlot}>{c.icon ? <Icon name={c.icon} /> : null}</span>
      <span class={classes.rowTitle}>
        {titleText}
        {tail ? <span class={classes.tail}>{tail}</span> : null}
      </span>
      <span class={classes.rowMeta}>{metaText}</span>
      {chord ? <Keys keys={chord} class={classes.rowEnd} /> : <span />}
    </div>
  ) as HTMLElement;
  let was = false;
  const counted = !!c.count && !match.preset && !match.blocked;
  const sources: WatchSource[] = [props.selected];
  if (counted) sources.push(manager.browsing);
  watch(sources, () => {
    was = writeSelected(el, props.selected.value, was);
    if (counted) {
      const b = manager.browsing.value;
      setText(metaText, (b?.id === c.id ? b.count : undefined) ?? meta);
    }
  });
  // The query's marks in their own watch: the query is a frame state, the highlight a key's.
  let marked: string | undefined;
  const mark = (q: string) => {
    if (q === marked) return;
    marked = q;
    writeMarks(titleText, title, q);
  };
  if (tail) mark('');
  else watch([query], () => mark(query.value));
  return el;
}

/** The text with the query's letters in `mark`s: a few text nodes, written whole. */
function writeMarks(el: HTMLElement, text: string, query: string) {
  const at = new Set(matchPositions(text, query));
  if (!at.size) {
    el.textContent = text;
    return;
  }
  const parts: Node[] = [];
  let run = '';
  let marked = false;
  const flush = () => {
    if (!run) return;
    if (marked) {
      const m = document.createElement('mark');
      m.textContent = run;
      parts.push(m);
    } else parts.push(document.createTextNode(run));
    run = '';
  };
  for (let i = 0; i < text.length; i++) {
    const m = at.has(i);
    if (m !== marked) {
      flush();
      marked = m;
    }
    run += text[i];
  }
  flush();
  el.replaceChildren(...parts);
}
