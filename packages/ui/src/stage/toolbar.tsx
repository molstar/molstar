/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, isState, type ReadonlyState, State } from '@molstar/core/reactive';
import type { IconName } from '../components/icon/icon.js';
import { ScrollArea } from '../components/scroll-area/scroll-area.js';
import {
  activity,
  type Child,
  effect,
  observeResize,
  on,
  requireScope,
  runInScope,
  Show,
  setAttr,
  setStyle,
} from '../dom/index.js';
import { Toolbar, type ToolbarGroup, type ToolbarItem, type ToolbarMenuRow } from '../components/toolbar/toolbar.js';
import { classes } from './_generated/toolbar.module.css.js';
import { StageCompactWidth, type StageModel } from './stage-model.js';

// The stage's toolbar (stage §3.2): `components/toolbar` in the stage's corner. What is the stage's: its menus open
// upward, as tall as the room above (measured once as a menu opens), a narrow row scrolls across, and folded
// (`collapsed`) it is one button whose menu lists the rows.

/** The margin a menu keeps from the stage's top, and its gap to its button (px: `--msui-stage-edge`, `MenuGap`). */
const MenuEdge = 8;
const MenuGap = 8;

/** The stage toolbar's own text (English by default, `DefaultStageToolbarLabels`). */
export interface StageToolbarLabels {
  /** The row's accessible name (`aria-label` takes precedence). */
  tools: string;
  /** The folded button's name and tooltip. */
  folded: string;
  /** The folded toolbar's accessible name (its one button's toolbar). */
  foldedToolbar: string;
}

export const DefaultStageToolbarLabels: StageToolbarLabels = {
  tools: 'Stage tools',
  folded: 'Toolbar',
  foldedToolbar: 'Toolbar',
};

export interface StageToolbarProps {
  /** The toolbar's blocks (`components/toolbar`): buttons, toggles, menus and submenus. */
  groups: readonly ToolbarGroup[] | ReadonlyState<readonly ToolbarGroup[]>;
  /** Before the toolbar, in its row (a mode switch). */
  lead?: Child;
  /** `lead` as rows of the folded toolbar's menu (it is not shown folded). */
  leadRows?: readonly ToolbarMenuRow[] | ReadonlyState<readonly ToolbarMenuRow[]>;
  /** Beside the row, outside its scroll (what opens next to it). */
  aside?: Child;
  /**
   * Folded into one button: `auto` (default) while the stage is compact (`StageModel.compact`), or, without `stage`,
   * while the element it is placed in (the closest `[data-msui-stage]`, else its parent) is narrower than
   * `StageCompactWidth`; or always, or never, or as a state says.
   */
  collapsed?: 'auto' | boolean | ReadonlyState<boolean>;
  /** The stage it sits on (`collapsed: 'auto'`). */
  stage?: StageModel;
  /** The folded button's icon (default `wrench`). */
  foldedIcon?: IconName;
  labels?: Partial<StageToolbarLabels>;
  'aria-label'?: string;
}

/** The folded toolbar's button. */
const FoldedId = 'msui-stage-toolbar';

function valueOf<T>(v: T | ReadonlyState<T>): T {
  return isState(v) ? (v as ReadonlyState<T>).value : (v as T);
}

/**
 * A toolbar item as a row of the folded toolbar's menu: a toggle (checkable) or an action; an item with a menu is a
 * submenu of its rows, led by its own action (unless a row of the menu already is it). Its face, disabled state and
 * menu as they are now (the rows are built as the menu opens). Hidden items have none.
 */
export function toolRow(item: ToolbarItem, group: string, anchor: () => HTMLElement): ToolbarMenuRow | undefined {
  if (item.visible !== undefined && !valueOf(item.visible)) return undefined;
  const face = item.face?.value;
  const label = face?.title ?? item.title;
  const icon = face?.icon ?? item.icon;
  const disabled = item.disabled === undefined ? undefined : valueOf(item.disabled);
  const run = () => item.run?.(anchor());
  const base = {
    id: item.id,
    label,
    icon,
    group,
    ...(disabled ? { disabled } : {}),
  };
  const menu = item.menu === undefined ? undefined : valueOf(item.menu);
  if (!menu?.length) return { ...base, ...(item.pressed ? { pressed: item.pressed } : {}), run };
  const own: ToolbarMenuRow[] =
    item.run && !menu.some((m) => m.label === label)
      ? [
          {
            id: `${item.id}:run`,
            label,
            icon,
            group: '\u0000own',
            ...(item.pressed ? { pressed: item.pressed } : {}),
            run,
          },
        ]
      : [];
  return { ...base, items: [...own, ...menu] };
}

/** The folded toolbar's rows: `lead`'s first (a group of their own), then each block's items, a hairline between. */
export function collapsedRows(
  groups: readonly ToolbarGroup[],
  lead: readonly ToolbarMenuRow[],
  anchor: () => HTMLElement,
): ToolbarMenuRow[] {
  const rows: ToolbarMenuRow[] = lead.map((r) => (r.group === undefined ? { ...r, group: '\u0000lead' } : r));
  for (const g of groups) {
    for (const item of g.items) {
      const row = toolRow(item, g.id, anchor);
      if (row) rows.push(row);
    }
  }
  return rows;
}

/**
 * The stage's toolbar (stage §3.2), for its `corner`: a `Toolbar` whose menus open upward, bounded by the room above it
 * in the stage (measured once as a menu opens: the stage's top to its button, less the menu's gap and the edge); a row
 * wider than the stage scrolls across. Folded, one button (`Toolbar`) whose menu lists the rows (`collapsedRows`);
 * the row stays mounted, hidden, inert and paused (`activity`). On touch a tap opens a button's menu with its own
 * action first (`components/toolbar`).
 */
export function StageToolbar(props: StageToolbarProps) {
  const labels = { ...DefaultStageToolbarLabels, ...props.labels };
  const groups = isState(props.groups)
    ? (props.groups as ReadonlyState<readonly ToolbarGroup[]>)
    : new State(props.groups as readonly ToolbarGroup[]);
  const lead = props.leadRows === undefined ? undefined : props.leadRows;
  const c = props.collapsed ?? 'auto';
  // `auto` without a stage: its container's width, measured once it is in the page (below).
  const narrow = c === 'auto' && !props.stage ? new State(false) : undefined;
  const collapsed: ReadonlyState<boolean> =
    c === 'auto' ? (props.stage?.compact ?? narrow!) : typeof c === 'boolean' ? new State(c) : c;

  // The folded menu's rows, built as it opens (and as the groups change while it is closed).
  const rows = new State<readonly ToolbarMenuRow[]>([], Object.is);
  let foldedButton: HTMLElement | undefined;
  const anchor = () => foldedButton ?? root;
  const refresh = () => rows.set(collapsedRows(groups.value, lead === undefined ? [] : valueOf(lead), anchor));
  const folded: ToolbarItem = { id: FoldedId, icon: props.foldedIcon ?? 'wrench', title: labels.folded, menu: rows };

  const bar = activity(
    derive(collapsed, (f) => !f),
    () => (
      <div class={classes.bar}>
        <ScrollArea axis="x" class={classes.scroll} contentClass={classes.row}>
          {props.lead}
          <Toolbar groups={groups} menuSide="above" elevated aria-label={props['aria-label'] ?? labels.tools} />
        </ScrollArea>
        {props.aside}
      </div>
    ),
  ) as HTMLElement;
  const root = (
    <div class={classes.stageToolbar} data-msui-stage-toolbar="">
      {bar}
      <Show when={collapsed}>
        {() => (
          <div
            class={classes.folded}
            ref={(el: HTMLElement) => {
              queueMicrotask(() => {
                foldedButton = el.querySelector<HTMLElement>('[data-toolbar-button]') ?? undefined;
              });
            }}
          >
            <Toolbar
              groups={[{ id: 'folded', items: [folded] }]}
              menuSide="above"
              elevated
              aria-label={labels.foldedToolbar}
            />
          </div>
        )}
      </Show>
    </div>
  ) as HTMLElement;
  // A menu of this toolbar is open: its rows and its room stay as they were when it opened (rebuilding the folded
  // menu's rows under the pointer would close its submenus; measuring from a row would shrink it).
  const menuOpen = () => !!root.querySelector('[aria-haspopup="menu"][aria-expanded="true"]:not([data-menu-row])');
  effect(collapsed, (f) => setAttr(bar, 'hidden', f));
  effect(groups, () => {
    if (collapsed.value && !menuOpen()) refresh();
  });

  // As a menu is about to open (the pointer comes to a menu button, or the focus or a key or a press does): the room
  // above, written when it changed; folded, the rows as they are now. Not for what happens in an open menu (its rows
  // are the toolbar's descendants too, a submenu's row a menu button).
  let room = -1;
  const opening = (e: Event) => {
    const target = e.target as Element | null;
    if (!target?.closest || target.closest('[role="menu"]')) return;
    const button = target.closest('[aria-haspopup="menu"]') as HTMLElement | null;
    if (!button || button.getAttribute('aria-expanded') === 'true' || !root.contains(button) || menuOpen()) return;
    const stage = root.closest('[data-msui-stage]');
    const top = stage ? stage.getBoundingClientRect().top : 0;
    const next = Math.max(0, Math.round(button.getBoundingClientRect().top - top - MenuGap - MenuEdge));
    if (next !== room) {
      room = next;
      setStyle(root, '--msui-stage-menu-room', `${next}px`);
    }
    if (button.closest(`.${classes.folded}`)) refresh();
  };
  for (const type of ['pointerover', 'pointerdown', 'focusin', 'keydown']) on(root, type, opening, { capture: true });

  // `auto` without a stage (stage §3.2): folded while its container is narrower than `StageCompactWidth`, measured
  // once it is in the page (the container's width does not follow the toolbar's: no folding back and forth).
  if (narrow) {
    const scope = requireScope('StageToolbar');
    queueMicrotask(() => {
      if (scope.disposed) return;
      const container = root.closest('[data-msui-stage]') ?? root.parentElement;
      if (!container) return;
      runInScope(scope, () => {
        observeResize(container, (entries) => {
          const r = entries[entries.length - 1]?.contentRect;
          if (r && r.width > 0) narrow.set(r.width < StageCompactWidth);
        });
      });
    });
  }
  return root;
}
