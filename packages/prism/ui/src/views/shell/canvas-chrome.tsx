/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { activity, For, Show } from '../../dom/control-flow.js';
import { eventTarget, isInside, ownsKey } from '../../dom/events.js';
import {
  deepActiveElement,
  focusables,
  focusedByKeyboard,
  focusFirstItem,
  menuKeys,
  rememberFocus,
  trapTab,
  uniqueId,
} from '../../dom/focus.js';
import {
  effect,
  type MaybeReactive,
  onCleanup,
  requireScope,
  runInScope,
  type Scope,
  toStream,
  type WatchSource,
  watch,
} from '../../dom/scope.js';
import { frame, on } from '../../dom/scoped.js';
import { setAttr, setClass, setStyle, setText } from '../../dom/writes.js';
import { keyedSelector } from '../../reactive/keyed-selector.js';
import type { CanvasMenuItem, CanvasToolbarItem, OpenOverlay, ShellModel } from '../../shell/shell.js';
import type { MessageKind, NotificationsManager, StatusTask } from '../../managers/notifications.js';
import type { Child } from '../../dom/jsx-types.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive, isState, toState } from '../../reactive/state.js';
import { classes } from './_generated/canvas.module.css.js';
import { PinButton } from './dock-view.js';
import { sectionsOf, shownItems } from './rail.js';
import { ViewHost } from './view-host.js';

/**
 * The toolbar: its buttons in sections, each a glass block of its own, the gaps between them the grouping. `part`
 * picks what to show: the rail's items (`all`), split around its tab strip (`lead`: the command menu above it, `rest`
 * below), or the canvas stack's items (`canvas`).
 */

/**
 * A model's stream as a state, for `derive` (shell.ts types an item's `face`, `visible` and `pressed` as streams;
 * most are states): a state as it is, a plain value as one, another stream through `toState` (the boundary,
 * rendering.md §2.2), its subscription ended with the current scope.
 */
function stateOf<T>(source: MaybeReactive<T>, initial: T): ReadonlyState<T> {
  const stream = toStream(source);
  if (!stream) return observableState(source as T).asReadonly;
  if (isState(stream)) return stream as ReadonlyState<T>;
  const state = toState(stream, initial);
  onCleanup(() => state.dispose());
  return state;
}

/** What a toolbar item looks like now: its icon and title. */
export interface ItemFace {
  icon: IconName;
  title: string;
}

/**
 * A toolbar item's face as a state, when it changes (the theme, the mode switch: `face`); undefined for an item that
 * looks the same for good (`staticTitle`, its own icon). Made once a view: one subscription to a stream.
 */
export function faceOf(item: CanvasToolbarItem): ReadonlyState<ItemFace> | undefined {
  return item.face ? stateOf(item.face, { icon: item.icon, title: item.title }) : undefined;
}

/** The title of an item whose face does not change: its own, with its caption. */
export function staticTitle(item: CanvasToolbarItem): string {
  return item.caption ? `${item.title} (${item.caption})` : item.title;
}

/** A toolbar item's icon as shown: its face's, following it (`face`: made by `faceOf`, else here), or its own. */
export function ItemIcon(props: {
  item: CanvasToolbarItem;
  size: number;
  class?: string;
  face?: ReadonlyState<ItemFace>;
}) {
  const { item } = props;
  const face = props.face ?? faceOf(item);
  if (!face) return <Icon name={item.icon} size={props.size} />;
  let name: IconName | undefined;
  const holder = (<span class={props.class} />) as HTMLElement;
  watch([face], () => {
    const next = face.value.icon;
    if (next === name) return;
    name = next;
    holder.replaceChildren(Icon({ name: next, size: props.size }) as Node);
  });
  return holder;
}

export function CanvasToolbar(props: {
  shell: ShellModel;
  part?: 'all' | 'lead' | 'rest' | 'canvas';
  /** Where buttons with a menu open it (the canvas stack); without it, menus are not shown. */
  menus?: ToolMenus;
  /** The rail's left-out buttons (`RailLayout.hidden`): they are not shown, they are under More. */
  hidden?: ReadonlyState<ReadonlySet<string>>;
  /** On the rail: a button with a menu opens it beside the rail on a click (the theme), as More does. */
  railMenu?: (item: CanvasToolbarItem, anchor: HTMLElement) => void;
}) {
  const { shell } = props;
  const part = props.part ?? 'all';
  const canvas = part === 'canvas';
  const items$ = props.hidden
    ? derive([shell.toolbar, props.hidden], (items, hidden) =>
        hidden.size ? items.filter((i) => !hidden.has(`tool:${i.id}`)) : items,
      )
    : shell.toolbar;
  const sections = derive(items$, (items) => {
    const all = sectionsOf(items.filter((i) => (i.place === 'canvas') === canvas));
    return part === 'all' || canvas
      ? all
      : all.filter((s) => (part === 'lead' ? s.group === undefined : s.group !== undefined));
  });
  return (
    <div class={classes.sections}>
      <For each={sections} by={(s) => s.group ?? ''}>
        {(section) => (
          <div
            class={classes.toolbar}
            role="toolbar"
            aria-orientation={props.menus ? 'horizontal' : 'vertical'}
            aria-label={section.group ?? 'commands'}
          >
            {section.items.map((item) =>
              item.visible ? (
                <Show when={item.visible}>
                  {() => <ToolButton shell={shell} item={item} menus={props.menus} railMenu={props.railMenu} />}
                </Show>
              ) : (
                <ToolButton shell={shell} item={item} menus={props.menus} railMenu={props.railMenu} />
              ),
            )}
          </div>
        )}
      </For>
    </div>
  );
}

/** Whether a button shows pressed: its overlay is open, or its own state says so. */
function toolPressed(shell: ShellModel, item: CanvasToolbarItem): ReadonlyState<boolean> {
  const own = stateOf(item.pressed ?? false, false);
  const id = item.overlay;
  return id ? derive([shell.overlay, own], (o, p) => o?.overlay.id === id || p) : own;
}

function ToolButton(props: {
  shell: ShellModel;
  item: CanvasToolbarItem;
  menus?: ToolMenus;
  railMenu?: (item: CanvasToolbarItem, anchor: HTMLElement) => void;
}) {
  const { shell, item, menus } = props;
  const railMenu = !menus && item.menu?.length ? props.railMenu : undefined;
  // A button with a menu shows the menu on hover instead of a tooltip (both would open in the same spot).
  const withMenu = menus && item.menu?.length ? menus : undefined;
  const pressed = item.pressed || item.overlay ? toolPressed(shell, item) : undefined;
  // Its title is its face's when that changes, else fixed.
  const face = faceOf(item);
  const title = face ? undefined : staticTitle(item);
  // The menu button pattern (accessibility.md): ↓ or ↑ opens the menu at its first item (Enter and Space too, on a
  // button that only holds a menu); the pointer opens it on hover. Focus alone does not open it: Tab goes on to the
  // next button. Whether its menu is open: its own key of the open menu.
  const expanded = withMenu?.shown.key(item);
  const button = (
    <button
      type="button"
      class={classes.tool}
      data-captioned={item.caption ? '' : undefined}
      data-tip={withMenu ? undefined : title}
      data-tip-side={item.place === 'canvas' ? 'top' : 'left'}
      aria-label={title}
      aria-haspopup={withMenu || railMenu ? 'menu' : undefined}
      data-overlay-trigger={item.overlay ? '' : undefined}
      data-keeps-card={item.keepsCard ? '' : undefined}
      onClick={(e) => {
        // A button that only holds a menu opens it (from the keyboard, at its first item: a click with no pointer
        // has `detail` 0); otherwise the click runs the button's own action.
        if (withMenu && !item.overlay && !item.run) return withMenu.show(item, e.currentTarget, e.detail === 0);
        // A touch screen has no hover to show the menu: a tap opens it (the button's own action its first row), or
        // closes it.
        if (withMenu && noHover()) {
          if (withMenu.open.value?.item === item) return withMenu.hide();
          return withMenu.show(item, e.currentTarget, e.detail === 0);
        }
        if (railMenu) return railMenu(item, e.currentTarget);
        withMenu?.hide();
        if (item.overlay) shell.toggleOverlay(item.overlay, e.currentTarget);
        else item.run?.(e.currentTarget);
      }}
      onKeyDown={
        withMenu
          ? (e) => {
              if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return;
              e.preventDefault();
              withMenu.show(item, e.currentTarget, true);
            }
          : undefined
      }
      onPointerEnter={withMenu ? (e) => e.pointerType === 'mouse' && withMenu.show(item, e.currentTarget) : undefined}
      onPointerLeave={withMenu ? () => withMenu.hideSoon() : undefined}
      onBlur={withMenu ? () => withMenu.hideSoon() : undefined}
    >
      <ItemIcon item={item} size={item.caption ? 16 : 20} class={classes.face} {...(face ? { face } : {})} />
      {item.caption ? <span class={classes.caption}>{item.caption}</span> : null}
      {withMenu ? <span class={classes.menuMark} aria-hidden="true" /> : null}
    </button>
  ) as HTMLElement;
  // What changes (a face's title, the pressed state, the open menu): one watch; a button with none has no watch. Only a
  // button with a state is a toggle (a button that cycles, like the theme, names its state in its label).
  const sources: WatchSource[] = [];
  if (face) sources.push(face);
  if (pressed) sources.push(pressed);
  if (expanded) sources.push(expanded);
  if (sources.length)
    watch(sources, () => {
      if (face) {
        setAttr(button, 'aria-label', face.value.title);
        if (!withMenu) setAttr(button, 'data-tip', face.value.title);
      }
      if (pressed) setAttr(button, 'aria-pressed', String(pressed.value));
      if (expanded) setAttr(button, 'aria-expanded', String(expanded.value));
    });
  return button;
}

/** No hover (a touch screen): what a pointer over a button would show, a tap shows. */
const noHover = () => typeof matchMedia === 'function' && matchMedia('(hover: none)').matches;

/**
 * A button's menu rows; on a touch screen, where a tap opens the menu instead of running the button, its own action
 * first (the camera's Reset before its variants).
 */
function withOwnAction(item: CanvasToolbarItem, button: Element): readonly CanvasMenuItem[] {
  const rows = item.menu ?? [];
  if (!item.run || !noHover()) return rows;
  return [
    { id: `${item.id}:self`, label: item.title, icon: item.icon, group: '', run: () => item.run?.(button) },
    ...rows,
  ];
}

/** Where a submenu row opens its rows: beside the menu (the canvas stack's), with the row it belongs to. */
export type OpenSubmenu = (item: CanvasMenuItem | undefined, row?: HTMLElement, focus?: boolean) => void;

/**
 * A toolbar button's menu rows (the canvas stack's menus, the rail's): grouped, a toggle's or a choice's check. A row
 * with `items` is a submenu: with `onSubmenu` it opens beside the menu (hover, click, →; any other row's hover closes
 * it), without it its rows show under its name.
 */
export function ToolMenuRows(props: { items: readonly CanvasMenuItem[]; onRun(): void; onSubmenu?: OpenSubmenu }) {
  const { onSubmenu } = props;
  return (
    <>
      {props.items.map((m, i, all) => (
        <>
          {i > 0 && m.group !== all[i - 1].group ? <div class={classes.toolMenuDivider} role="separator" /> : null}
          {m.items && !onSubmenu ? (
            <>
              <div class={classes.toolMenuHeading} role="presentation">
                {m.label}
              </div>
              <Show when={m.items}>{(rows) => <ToolMenuRows items={rows} onRun={props.onRun} />}</Show>
            </>
          ) : m.items && onSubmenu ? (
            <button
              type="button"
              role="menuitem"
              aria-haspopup="menu"
              class={classes.toolMenuRow}
              disabled={!!m.disabled}
              onPointerEnter={(e: PointerEvent) => onSubmenu(m, e.currentTarget as HTMLElement)}
              onClick={(e: MouseEvent) => onSubmenu(m, e.currentTarget as HTMLElement, true)}
              onKeyDown={(e: KeyboardEvent) => {
                if (e.key !== 'ArrowRight') return;
                e.preventDefault();
                onSubmenu(m, e.currentTarget as HTMLElement, true);
              }}
            >
              {m.icon ? <Icon name={m.icon} size={16} /> : <span />}
              <span class={classes.toolMenuLabel}>
                {m.label}
                {m.disabled || m.hint ? <span class={classes.toolMenuHint}>{m.disabled ?? m.hint}</span> : null}
              </span>
              <Icon name="chevron-right" size={14} class={classes.toolMenuChevron} />
            </button>
          ) : (
            <button
              type="button"
              role={m.pressed ? (m.choice ? 'menuitemradio' : 'menuitemcheckbox') : 'menuitem'}
              aria-checked={m.pressed}
              class={classes.toolMenuRow}
              disabled={!!m.disabled}
              onPointerEnter={onSubmenu ? () => onSubmenu(undefined) : undefined}
              onClick={() => {
                props.onRun();
                m.run();
              }}
            >
              {m.icon ? <Icon name={m.icon} size={16} /> : <span />}
              <span class={classes.toolMenuLabel}>
                {m.label}
                {m.disabled || m.hint ? <span class={classes.toolMenuHint}>{m.disabled ?? m.hint}</span> : null}
              </span>
              {m.pressed ? (
                <Show when={m.pressed}>{() => <Icon name="check" size={14} class={classes.toolMenuCheck} />}</Show>
              ) : null}
            </button>
          )}
        </>
      ))}
    </>
  );
}

/**
 * A canvas toolbar item as a row of the collapsed toolbar's menu (shell.md §2): a toggle (pressed) or an action; an
 * item with a menu is a submenu of its own rows, led by its own action (unless a row of the menu is that action
 * already: Expand's first row).
 */
function toolRow(shell: ShellModel, item: CanvasToolbarItem, anchor: () => Element): CanvasMenuItem {
  const toggle = item.pressed || item.overlay ? toolPressed(shell, item) : undefined;
  const act = () => {
    if (item.overlay) shell.toggleOverlay(item.overlay, anchor());
    else item.run?.(anchor());
  };
  const row = { id: item.id, label: item.title, icon: item.icon, group: item.group ?? '' };
  const menu = item.menu ?? [];
  if (!menu.length) return { ...row, ...(toggle ? { pressed: toggle } : {}), run: act };
  const own: CanvasMenuItem[] =
    (item.run || item.overlay) && !menu.some((m) => m.label === item.title)
      ? [
          {
            id: `${item.id}:self`,
            label: item.title,
            icon: item.icon,
            group: '',
            ...(toggle ? { pressed: toggle } : {}),
            run: act,
          },
        ]
      : [];
  return { ...row, run: () => {}, items: observableState<readonly CanvasMenuItem[]>([...own, ...menu]).asReadonly };
}

/**
 * What the collapsed toolbar's menu lists: the app's lead rows (the mode buttons), then the canvas items that show, in
 * the order of the bar's sections, a hairline between sections. Follows the items and their `visible`.
 */
function collapsedRows(
  shell: ShellModel,
  lead: MaybeReactive<readonly CanvasMenuItem[]> | undefined,
  anchor: () => Element,
): ReadonlyState<readonly CanvasMenuItem[]> {
  // An item's row is made once, in a scope of its own under the stack's (its pressed state follows the item); the row
  // of an item that left is dropped with its scope.
  const scope = requireScope('collapsedRows');
  const rows = new Map<CanvasToolbarItem, { row: CanvasMenuItem; scope: Scope }>();
  const rowOf = (i: CanvasToolbarItem) => {
    let made = rows.get(i);
    if (!made) {
      const own = scope.child();
      rows.set(i, (made = { row: runInScope(own, () => toolRow(shell, i, anchor)), scope: own }));
    }
    return made.row;
  };
  const tools = shownItems(
    derive(shell.toolbar, (items) => sectionsOf(items.filter((i) => i.place === 'canvas')).flatMap((s) => s.items)),
  );
  return derive([stateOf(lead ?? [], []), tools], (l, t) => {
    const out = [...l.map((r) => ({ ...r, group: r.group ?? 'lead' })), ...t.map(rowOf)];
    if (rows.size > t.length) {
      const now = new Set(t);
      for (const [item, made] of rows) {
        if (now.has(item)) continue;
        made.scope.dispose();
        rows.delete(item);
      }
    }
    return out;
  });
}

/** What the collapsed toolbar's menu is made of (`ToolMenus.show`): its rows. */
interface StackPopup {
  readonly rows: ReadonlyState<readonly CanvasMenuItem[]>;
}

/** The collapsed toolbar's one button, as a toolbar item (its menu is `StackPopup`'s). */
const CollapsedItem: CanvasToolbarItem = { id: 'toolbar', place: 'canvas', icon: 'wrench', title: 'Toolbar' };

/**
 * The canvas stack at the bottom right of the canvas, beside the rail or the open panel: one frame with the mode buttons
 * (`lead`) and then the toolbar items placed on the canvas (camera, lighting, XR, expand, sequence); `aside` is what
 * opens beside it (the picking level).
 *
 * Compact (`ShellModel.toolbarCollapsed` from `ShellOptions.toolbar`, shell.md §2) it is one button whose menu lists the
 * same tools (`leadMenu`: the rows of what `lead` holds; the items' own menus as submenus). The person has no control
 * to switch; the host's option and, under `auto`, the shell's mode decide. The bar stays mounted while it is folded,
 * hidden, inert and paused (`activity`).
 */
export function CanvasStack(props: {
  shell: ShellModel;
  lead?: Child;
  aside?: Child;
  /** The rows of `lead` for the collapsed toolbar's menu (the mode choices); they follow what `lead` shows. */
  leadMenu?: MaybeReactive<readonly CanvasMenuItem[]>;
  /**
   * The canvas is empty: the stack fades out (`CanvasStackFadeMs`; what is in it should stay that long), and in again
   * when it has something (canvas.module.css, `@starting-style`).
   */
  quiet?: MaybeReactive<boolean>;
}) {
  const { shell } = props;
  const menus = new ToolMenus();
  let stack: HTMLElement | undefined;
  let folded: HTMLElement | undefined;
  const popup: StackPopup = { rows: collapsedRows(shell, props.leadMenu, () => folded ?? stack!) };
  const collapsed = shell.toolbarCollapsed;
  const bar = activity(
    derive(collapsed, (c) => !c),
    () => (
      <div class={classes.bar}>
        <div class={classes.frame}>
          <ScrollArea axis="x" contentClass={classes.stackTools}>
            {props.lead}
            <CanvasToolbar shell={shell} part="canvas" menus={menus} />
          </ScrollArea>
        </div>
        {props.aside}
      </div>
    ),
  ) as HTMLElement;
  const hasTools = derive(popup.rows, (rows) => rows.length > 0);
  const quiet = props.quiet === undefined ? undefined : stateOf(props.quiet, false);
  const root = (
    <div class={classes.stack}>
      {bar}
      <Show when={collapsed}>
        {() => (
          <Show when={hasTools}>
            {() => (
              <div class={`${classes.frame} ${classes.folded}`}>
                <CollapsedButton
                  menus={menus}
                  popup={popup}
                  ref={(el) => {
                    folded = el;
                  }}
                />
              </div>
            )}
          </Show>
        )}
      </Show>
      <ToolMenu menus={menus} />
    </div>
  ) as HTMLElement;
  stack = root;
  // Folded: the bar is hidden and inert (its buttons are neither seen, tabbed to nor read); a menu that went closes.
  watch(quiet ? [collapsed, quiet] : [collapsed], () => {
    setAttr(bar, 'hidden', collapsed.value);
    if (quiet) setAttr(root, 'data-quiet', quiet.value);
  });
  effect(collapsed, () => menus.hide());
  menus.stack = root;
  onCleanup(() => menus.dispose());
  return root;
}

/** The collapsed toolbar's button: opens the menu of the tools as a tool with a menu does (hover, a click, ↓ or ↑). */
function CollapsedButton(props: { menus: ToolMenus; popup: StackPopup; ref(el: HTMLElement): void }) {
  const { menus, popup } = props;
  const expanded = menus.shown.key(CollapsedItem);
  const button = (
    <button
      type="button"
      class={classes.tool}
      data-stack-collapsed=""
      aria-label="Toolbar"
      aria-haspopup="menu"
      ref={props.ref}
      onClick={(e) => {
        // A touch screen has no hover: a tap opens the menu, or closes it. From the keyboard it opens at its first row.
        if (noHover() && menus.open.value?.item === CollapsedItem) return menus.hide();
        menus.show(CollapsedItem, e.currentTarget, e.detail === 0, popup);
      }}
      onKeyDown={(e) => {
        if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return;
        e.preventDefault();
        menus.show(CollapsedItem, e.currentTarget, true, popup);
      }}
      onPointerEnter={(e) => e.pointerType === 'mouse' && menus.show(CollapsedItem, e.currentTarget, false, popup)}
      onPointerLeave={() => menus.hideSoon()}
      onBlur={() => menus.hideSoon()}
    >
      <Icon name="wrench" size={18} />
      <span class={classes.menuMark} aria-hidden="true" />
    </button>
  ) as HTMLElement;
  watch([expanded], () => setAttr(button, 'aria-expanded', String(expanded.value)));
  return button;
}

/**
 * The open tool menu of the canvas stack: which button, at what height. Opens on hover or from the keyboard (its
 * button's ↓, ↑, Enter), closes a moment after the pointer leaves both the button and the menu (so it can travel
 * between them), when the focus leaves it, on a choice, or on Escape (the focus goes back to its button).
 */
/** The margin (px) a menu keeps from the canvas's top: `--msp-canvas-edge`. */
const MenuEdgePx = 8;

export class ToolMenus {
  readonly open = observableState<
    { item: CanvasToolbarItem; right: number; room: number; button: Element; popup?: StackPopup } | undefined
  >(undefined);
  /** Whose menu is open, by key: a button follows its own (a move of the menu wakes two buttons, not all). */
  readonly shown = keyedSelector(derive(this.open.asReadonly, (o) => o?.item));
  stack: HTMLElement | undefined;
  /** The open menu's element (`ToolMenu` sets it). */
  menu: HTMLElement | undefined;
  #timer: ReturnType<typeof setTimeout> | undefined;

  /**
   * Opens `item`'s menu by its button; `focus`: at its first item (opened from the keyboard). `popup`: the collapsed
   * toolbar's menu, whose rows follow its tools in place of the item's `menu`.
   */
  show(item: CanvasToolbarItem, button: Element, focus = false, popup?: StackPopup) {
    this.keep();
    // The menu opens upward from the stack (a row in the canvas' corner), its right edge level with the button's.
    const right = this.stack ? this.stack.getBoundingClientRect().right - button.getBoundingClientRect().right : 0;
    const room = this.roomAbove(this.stack);
    const current = this.open.value;
    if (current?.item !== item || current.right !== right || current.room !== room) {
      this.open.set({ item, right, room, button, ...(popup ? { popup } : {}) });
    }
    if (focus) queueMicrotask(() => this.menu && focusFirstItem(this.menu));
  }

  /**
   * The room a menu has above `el` (px): from the canvas's top (its viewport, else the page's) to `el`'s top, less a
   * margin. A menu is as tall as that and scrolls inside; 0 where there is no layout. Measured once, when a menu opens.
   */
  roomAbove(el: Element | undefined, below?: number): number {
    if (!el) return 0;
    const viewport = el.closest('[data-msp-viewport]');
    const room = (below ?? el.getBoundingClientRect().top) - (viewport?.getBoundingClientRect().top ?? 0) - MenuEdgePx;
    return room > 0 ? Math.round(room) : 0;
  }

  keep() {
    if (this.#timer !== undefined) clearTimeout(this.#timer);
    this.#timer = undefined;
  }

  hideSoon() {
    this.keep();
    this.#timer = setTimeout(() => this.hide(), 160);
  }

  hide() {
    this.keep();
    if (this.open.value) this.open.set(undefined);
  }

  dispose() {
    this.keep();
  }
}

function ToolMenu(props: { menus: ToolMenus }) {
  const { menus } = props;
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented || !menus.open.value) return;
    // Escape in another viewer on the page is not this menu's.
    if (menus.stack && !ownsKey(menus.stack, e)) return;
    menus.hide();
    e.preventDefault();
  };
  // A tap outside the menu and its button closes it (no pointer leaves it on a touch screen).
  const onDown = (e: PointerEvent) => {
    const open = menus.open.value;
    if (!open) return;
    const path = e.composedPath();
    if ((menus.menu && path.includes(menus.menu)) || path.includes(open.button)) return;
    menus.hide();
  };
  on(document, 'keydown', onKey);
  on(document, 'pointerdown', onDown, { capture: true });
  return (
    <Show when={menus.open}>
      {(open) => {
        // Closed while it had the focus (Escape, a choice): the focus goes back to its button, not to the page.
        let hadFocus = false;
        let el: HTMLElement | undefined;
        onCleanup(() => {
          if (menus.menu === el) menus.menu = undefined;
          if (!hadFocus) return;
          const now = deepActiveElement();
          if (!now || el?.contains(now)) (open.button as HTMLElement).focus({ preventScroll: true });
        });
        // The open submenu: its row, and where it sits (its bottom level with the row's, measured once on opening).
        const sub = observableState<
          { item: CanvasMenuItem; row: HTMLElement; bottom: number; room: number } | undefined
        >(undefined);
        let subEl: HTMLElement | undefined;
        const onSubmenu: OpenSubmenu = (item, row, focus = false) => {
          if (!item || !row || !el) {
            if (sub.value) sub.set(undefined);
            return;
          }
          // Rects, not offsets: the rows sit in a scrolling area. Its room: from the canvas's top to its row's bottom.
          const rowBottom = row.getBoundingClientRect().bottom;
          const bottom = Math.round(el.getBoundingClientRect().bottom - rowBottom);
          const room = menus.roomAbove(menus.stack, rowBottom);
          if (sub.value?.item !== item || sub.value.bottom !== bottom || sub.value.room !== room) {
            sub.set({ item, row, bottom, room });
          }
          if (focus) queueMicrotask(() => subEl && focusFirstItem(subEl));
        };
        // ← or Escape in the submenu: back to its row (the menu stays open).
        const subKeys = (e: KeyboardEvent) => {
          if (e.key === 'ArrowLeft' || e.key === 'Escape') {
            e.preventDefault();
            const row = sub.value?.row;
            sub.set(undefined);
            row?.focus({ preventScroll: true });
            return;
          }
          menuKeys(e);
        };
        return (
          <div
            class={classes.toolMenu}
            style={{ right: `${open.right}px`, ...(open.room ? { '--msp-menu-room': `${open.room}px` } : {}) }}
            role="menu"
            aria-label={open.item.title}
            ref={(x) => {
              el = x;
              menus.menu = x;
            }}
            onPointerEnter={() => menus.keep()}
            onPointerLeave={() => menus.hideSoon()}
            onFocusIn={() => {
              hadFocus = true;
              menus.keep();
            }}
            onFocusOut={() => menus.hideSoon()}
          >
            <div class={classes.toolMenuCard} onKeyDown={menuKeys}>
              {open.popup ? (
                <>
                  {/* The rows scroll where the canvas is too short for them (an embedded viewer); the first row stays. */}
                  <Show when={open.popup.rows}>
                    {(rows) => (
                      <ScrollArea class={classes.toolMenuScroll} contentClass={classes.toolMenuScrollContent}>
                        <ToolMenuRows items={rows} onRun={() => menus.hide()} onSubmenu={onSubmenu} />
                      </ScrollArea>
                    )}
                  </Show>
                </>
              ) : (
                <ScrollArea class={classes.toolMenuScroll} contentClass={classes.toolMenuScrollContent}>
                  <ToolMenuRows
                    items={withOwnAction(open.item, open.button)}
                    onRun={() => menus.hide()}
                    onSubmenu={onSubmenu}
                  />
                </ScrollArea>
              )}
            </div>
            <Show when={sub}>
              {(s) => (
                <div
                  class={classes.toolSubmenu}
                  style={{ bottom: `${s.bottom}px`, ...(s.room ? { '--msp-menu-room': `${s.room}px` } : {}) }}
                  role="menu"
                  aria-label={s.item.label}
                  ref={(x) => {
                    subEl = x;
                  }}
                  onKeyDown={subKeys}
                >
                  <div class={classes.toolMenuCard}>
                    <ScrollArea class={classes.toolMenuScroll} contentClass={classes.toolMenuScrollContent}>
                      <Show when={s.item.items!}>
                        {(rows) => <ToolMenuRows items={rows} onRun={() => menus.hide()} />}
                      </Show>
                    </ScrollArea>
                  </div>
                </div>
              )}
            </Show>
          </div>
        );
      }}
    </Show>
  );
}

/**
 * Compact: the end of the bottom tab bar. The command menu, then More: the tabs not in the bar (State tree, Log,
 * Settings) and the rail's actions (selection, animation, canvas settings, save, export), opening just above it.
 */
export function RailMenu(props: { shell: ShellModel }) {
  const { shell } = props;
  const open = observableState(false);
  const right = shell.dock('right');
  // In the bar: the command menu (the rail's leading items), then More with the tabs not in the bar and the actions.
  const lead = derive(shell.toolbar, (all) => all.filter((i) => i.place !== 'canvas' && !i.group));
  const actions = derive(shell.toolbar, (all) => all.filter((i) => i.place !== 'canvas' && i.group));
  const tabs = derive(right.tabs, (all) => (all.some((t) => t.compactBar) ? all.filter((t) => !t.compactBar) : []));
  let list: HTMLElement | undefined;
  const run = (item: CanvasToolbarItem, anchor: Element) => {
    open.set(false);
    if (item.overlay) shell.toggleOverlay(item.overlay, anchor);
    else item.run?.(anchor);
  };
  let more!: HTMLElement;
  const root = (
    <div class={classes.railMenu}>
      <For each={lead} by={(i) => i.id}>
        {(item) => {
          const face = faceOf(item);
          const title = face ? undefined : staticTitle(item);
          const button = (
            <button
              type="button"
              class={classes.railMenuButton}
              aria-label={title}
              data-tip={title}
              onClick={(e) => run(item, e.currentTarget)}
            >
              <ItemIcon item={item} size={20} class={classes.face} {...(face ? { face } : {})} />
            </button>
          ) as HTMLElement;
          if (face)
            watch([face], () => {
              setAttr(button, 'aria-label', face.value.title);
              setAttr(button, 'data-tip', face.value.title);
            });
          return button;
        }}
      </For>
      {/* More and its list: the list opens above the button (its anchor), inside the screen. */}
      <div class={classes.railMenuAnchor}>
        <button
          type="button"
          class={classes.railMenuButton}
          aria-label="More"
          aria-haspopup="menu"
          data-tip="More"
          ref={(el) => (more = el)}
          onClick={(e) => {
            open.set(!open.value);
            // From the keyboard, the menu opens at its first item.
            if (open.value && e.detail === 0) queueMicrotask(() => list && focusFirstItem(list));
          }}
          onKeyDown={(e) => {
            if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return;
            e.preventDefault();
            open.set(true);
            queueMicrotask(() => list && focusFirstItem(list));
          }}
        >
          <Icon name="ellipsis" size={20} />
        </button>
        <Show when={open}>
          {() => {
            const restore = rememberFocus();
            onCleanup(() => restore(list));
            return (
              <div
                class={classes.railMenuList}
                role="menu"
                aria-label="More"
                ref={(el) => {
                  list = el;
                }}
                onKeyDown={(e) => {
                  menuKeys(e);
                  if (e.key !== 'Escape') return;
                  e.preventDefault();
                  e.stopPropagation();
                  open.set(false);
                }}
              >
                <For each={tabs} by={(t) => t.id}>
                  {(tab) => (
                    <button
                      type="button"
                      role="menuitem"
                      class={classes.toolMenuRow}
                      onClick={() => {
                        open.set(false);
                        shell.show(tab.id);
                      }}
                    >
                      <Icon name={tab.icon} size={16} />
                      <span class={classes.toolMenuLabel}>{tab.title}</span>
                      <span />
                    </button>
                  )}
                </For>
                <Show when={derive([tabs, actions], (t, a) => t.length > 0 && a.length > 0)}>
                  {() => <div class={classes.railMenuDivider} />}
                </Show>
                <For each={actions} by={(i) => i.id}>
                  {(item) => {
                    const face = faceOf(item);
                    const visible = item.visible ? stateOf(item.visible, true) : undefined;
                    const label = document.createTextNode(face ? face.value.title : staticTitle(item));
                    const row = (
                      <button
                        type="button"
                        role="menuitem"
                        class={classes.toolMenuRow}
                        onClick={(e) => run(item, e.currentTarget)}
                      >
                        <ItemIcon item={item} size={16} class={classes.face} {...(face ? { face } : {})} />
                        <span class={classes.toolMenuLabel}>{label}</span>
                        <span />
                      </button>
                    ) as HTMLElement;
                    // As on the rail: an item shown only while it applies (Play snapshots with two or more).
                    const sources: WatchSource[] = [];
                    if (face) sources.push(face);
                    if (visible) sources.push(visible);
                    if (sources.length)
                      watch(sources, () => {
                        if (face) setText(label, face.value.title);
                        if (visible) setAttr(row, 'hidden', !visible.value);
                      });
                    return row;
                  }}
                </For>
              </div>
            );
          }}
        </Show>
      </div>
    </div>
  ) as HTMLElement;
  watch([open.asReadonly], () => setAttr(more, 'aria-expanded', String(open.value)));
  const onDown = (e: PointerEvent) => {
    if (open.value && !isInside(e, root)) open.set(false);
  };
  on(document, 'pointerdown', onDown);
  return root;
}

/**
 * A task in overlay mode (rendering an animation) holds the plugin: a scrim over everything and a card with what is
 * running, its progress, and Cancel.
 */
export function BusyOverlay(props: { status: NotificationsManager }) {
  const task = derive(
    props.status.tasks,
    (tasks) => tasks.find((t) => t.overlay),
    (a, b) => a?.id === b?.id,
  );
  return (
    <Show when={task}>
      {(first) => {
        const live = derive(props.status.tasks, (tasks) => tasks.find((t) => t.id === first.id) ?? first);
        // The task's label and progress, in whole percent: the overlay renders when one of them changes.
        const percent = derive(live, (t) => (t.progress === undefined ? undefined : Math.round(t.progress * 100)));
        const label = document.createTextNode('');
        const percentText = document.createTextNode('');
        let track!: HTMLElement;
        let bar!: HTMLElement;
        const card = (
          <div class={classes.busy} role="alertdialog" aria-modal="true" aria-label={first.label}>
            <div class={classes.busyCard}>
              <div class={classes.busyHead}>
                <Icon name="loader-circle" size={18} class={classes.busySpin} />
                <span class={classes.busyLabel}>{label}</span>
                <span class={classes.busyPercent}>{percentText}</span>
              </div>
              <div
                class={classes.busyTrack}
                role="progressbar"
                aria-valuemin={0}
                aria-valuemax={100}
                ref={(el: HTMLElement) => (track = el)}
              >
                <span class={classes.busyBar} ref={(el: HTMLElement) => (bar = el)} />
              </div>
              {first.abort ? (
                <div class={classes.busyActions}>
                  <Button size="sm" onClick={() => first.abort?.()}>
                    Cancel
                  </Button>
                </div>
              ) : null}
            </div>
          </div>
        ) as HTMLElement;
        watch([live, percent], () => {
          const p = percent.value;
          setText(label, live.value.label);
          setText(percentText, p === undefined ? '' : `${p}%`);
          setAttr(track, 'aria-valuenow', p);
          setAttr(bar, 'data-indeterminate', p === undefined);
          setStyle(bar, 'scale', `${(p ?? 30) / 100} 1`);
        });
        return card;
      }}
    </Show>
  );
}

/** Per shell, the place the compact sheet offers a card (`CardSlot`); none: the card is over the canvas. */
const CardSlots = new WeakMap<ShellModel, ObservableState<HTMLElement | undefined>>();
function cardSlot(shell: ShellModel): ObservableState<HTMLElement | undefined> {
  let slot = CardSlots.get(shell);
  // By identity: elements are never shallowly the same.
  if (!slot) CardSlots.set(shell, (slot = observableState<HTMLElement | undefined>(undefined, { compare: Object.is })));
  return slot;
}

/**
 * Where the compact sheet (the right `DockView`, compact) takes the card: in the panel's place, between the panel and
 * the tab bar (shell.md, Compact). Offered while it is mounted; the card moves in and out (`OverlayCard`).
 */
export function CardSlot(props: { shell: ShellModel }) {
  const slot = cardSlot(props.shell);
  const el = (<div class={classes.cardSlot} />) as HTMLElement;
  slot.set(el);
  onCleanup(() => {
    if (slot.value === el) slot.set(undefined);
  });
  return el;
}

/**
 * The open overlay as a card in the right panel's place (the panel steps aside while it is open). One card while any is
 * open: another one opened swaps its title and body, as a tab swaps a panel's content (no closing and opening again).
 * Over the canvas, or in the compact sheet where one offers it (`CardSlot`): one card, moved between the two as the
 * layout changes (never a second one made: its view and state stay, its focus too).
 */
export function OverlayHost(props: { shell: ShellModel }) {
  const { shell } = props;
  const isOpen = derive(shell.overlay, (o) => !!o);
  return <Show when={isOpen}>{() => <OverlayCard shell={shell} />}</Show>;
}

function OverlayCard(props: { shell: ShellModel }) {
  const { shell } = props;
  /** The overlay open now (the card lives while one is). */
  const current = () => shell.overlay.value as OpenOverlay;
  // A card is a dialog (accessibility.md: cards trap and restore focus): opened from the keyboard, the focus goes into
  // it; Tab stays inside while it floats (pinned, it is part of the workspace); closed, the focus goes back to the
  // button of the last card opened.
  let anchor = current().anchor as HTMLElement | undefined;
  const restore = rememberFocus(() => anchor);
  const titleId = uniqueId('msp-card-title');
  const pinnedState = derive(shell.dock('right').state, (s) => s.pinned);
  const title = document.createTextNode('');
  let cancelFocus: (() => void) | undefined;
  const slot = cardSlot(shell);
  /** In the compact sheet (a slot offered), not over the canvas. */
  const docked = () => !!slot.value;
  const card = (
    <section
      class={classes.card}
      role="dialog"
      aria-labelledby={titleId}
      onKeyDown={(e) => {
        if (!inPlace()) trapTab(e, e.currentTarget);
      }}
    >
      <header class={classes.cardHeader}>
        <h2 class={classes.cardTitle} id={titleId}>
          {title}
        </h2>
        <PinButton shell={shell} dock="right" />
        <Button variant="ghost" size="sm" icon={<Icon name="x" />} title="Close" onClick={() => shell.closeOverlay()} />
      </header>
      <Show when={shell.overlay}>
        {(open) => {
          anchor = open.anchor as HTMLElement | undefined;
          // Into the card: its body's first control once the view is there (a frame), else its Close.
          if (focusedByKeyboard(open.anchor)) {
            cancelFocus?.();
            cancelFocus = frame(() => {
              cancelFocus = undefined;
              const body = card.querySelector(`.${classes.cardBody}, .${classes.cardBodyFill}`);
              const all = focusables(card);
              (all.find((el) => body?.contains(el)) ?? all[all.length - 1])?.focus();
            });
          }
          const ctx = { shell, ...(open.focus ? { focus: open.focus } : {}) };
          return open.overlay.fill ? (
            <div class={classes.cardBodyFill}>
              <ViewHost renderer={open.overlay.render} ctx={ctx} />
            </div>
          ) : (
            <ScrollArea class={classes.cardBody} contentClass={classes.cardBodyContent}>
              <ViewHost renderer={open.overlay.render} ctx={ctx} />
            </ScrollArea>
          );
        }}
      </Show>
    </section>
  ) as HTMLElement;
  // Over the canvas (here, where the host put it) or in the sheet's slot: moved, with the focus kept where it was.
  const home = (<div style="display: contents">{card}</div>) as HTMLElement;
  watch([slot.asReadonly], () => {
    const into = slot.value ?? home;
    if (card.parentElement === into) return;
    const focused = deepActiveElement();
    const keep = focused && card.contains(focused) ? (focused as HTMLElement) : undefined;
    into.append(card);
    setClass(card, classes.card, !slot.value);
    setClass(card, classes.cardDocked, !!slot.value);
    keep?.focus({ preventScroll: true });
  });
  onCleanup(() => card.remove());
  // Modal while it floats; its title is the open overlay's (another one swaps it).
  watch([pinnedState, slot.asReadonly, shell.overlay], () => {
    setAttr(card, 'aria-modal', String(!slot.value && !pinnedState.value));
    setText(title, shell.overlay.value?.overlay.title ?? '');
  });
  // Esc, or a pointer down outside the card and its button, closes it; not one on the canvas: the person turns and
  // zooms the view with the card open (setting up an export, trying a canvas setting). Pinned, the card is part of the
  // workspace like the pinned panel: only its ✕ (or another card or tab) closes it; a persistent card (the Story, read
  // beside the scene) is not closed by a click elsewhere either.
  const pinned = () => shell.dock('right').state.value.pinned;
  // Part of the layout, not a dialog over it: pinned, or in the compact sheet (Tab leaves it, as it leaves a panel).
  const inPlace = () => docked() || pinned();
  const doc = document;
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented || pinned() || !ownsKey(card, e)) return;
    shell.closeOverlay();
    e.preventDefault();
  };
  const onPointer = (e: PointerEvent) => {
    const open = shell.overlay.value;
    if (!open) return;
    const target = eventTarget(e) as Element;
    if (pinned() || open.overlay.persistent || isInside(e, card) || isInside(e, open.anchor)) return;
    // Buttons that open another card, and panel tabs, swap the card on their click: closing it here first would show
    // the panel underneath for a frame (a blink), and a tab click would then collapse the panel.
    if (target.closest?.('[data-overlay-trigger], [role=tab]')) return;
    // What acts on the open card or the scene leaves it open: the command menu and its button, the canvas and the bar
    // at its bottom (the status, the selection's Clear, the toolbar).
    if (target.closest?.('[data-keeps-card], [data-msp-canvas]')) return;
    shell.closeOverlay();
  };
  on(doc, 'keydown', onKey);
  on(doc, 'pointerdown', onPointer, { capture: true });
  onCleanup(() => restore(card));
  return home;
}

export const KindIcons: Record<MessageKind, IconName> = {
  info: 'info',
  success: 'circle-check',
  warning: 'triangle-alert',
  error: 'circle-alert',
};

/** The bottom left of the canvas: the latest message, running tasks, and under them `after` (the assistant's turn). */
export function StatusLine(props: { status: NotificationsManager; onOpenLog?: () => void; after?: Child }) {
  const { status } = props;
  // An error opens the Log at the problems; another message, at everything (it may be in no other filter).
  const openLog = (kind: MessageKind) => {
    status.setLogFilter(kind === 'error' ? 'problems' : 'all');
    props.onOpenLog?.();
  };
  const spokenText = document.createTextNode('');
  watch([status.spoken], () => setText(spokenText, status.spoken.value));
  // One polite live region for results (accessibility.md): the latest message only, not the tasks' progress or the
  // assistant's streamed words (its card says when a turn is done).
  return (
    <div class={classes.status}>
      <div class={classes.announce} role="status" aria-live="polite">
        {spokenText}
      </div>
      <Show when={status.current}>
        {(m) =>
          // An error stays until it is dismissed: its ✕ does (the others fade).
          m.kind === 'error' ? (
            <span class={`${classes.pill} ${classes.message}`} data-kind={m.kind}>
              <button
                type="button"
                class={classes.pillMain}
                title="Show the problems in the log"
                onClick={() => openLog(m.kind)}
              >
                <Icon name={KindIcons[m.kind]} size={14} />
                <span class={classes.text}>{m.text}</span>
              </button>
              <button
                type="button"
                class={classes.pillButton}
                aria-label="Dismiss"
                data-tip="Dismiss"
                onClick={() => status.dismiss(m.id)}
              >
                <Icon name="x" size={14} />
              </button>
            </span>
          ) : (
            <button
              type="button"
              class={`${classes.pill} ${classes.message}`}
              data-kind={m.kind}
              title="Show the log"
              onClick={() => openLog(m.kind)}
            >
              <Icon name={KindIcons[m.kind]} size={14} />
              <span class={classes.text}>{m.text}</span>
            </button>
          )
        }
      </Show>
      <For each={derive(status.tasks, (tasks) => tasks.filter((t) => !t.overlay))}>
        {(task) => <TaskPill task={task} />}
      </For>
      {props.after}
    </div>
  );
}

/**
 * One line of the status column (a task's look): an icon, its text, what follows (state, progress), and an ✕; the line
 * itself clickable when `onClick` is given. For lines that are not tasks (the assistant's turn).
 */
export function StatusPill(props: {
  icon: Child;
  text: Child;
  after?: Child;
  label?: string;
  onClick?(): void;
  onDismiss?(): void;
  dismissLabel?: string;
  /**
   * Dismissed on its own after this many ms (a finished turn), a ring around the ✕ showing the time left; the pointer
   * over the line or the focus in it holds the time. Undefined: it stays.
   */
  autoDismiss?: MaybeReactive<number | undefined>;
}) {
  // The time left, held while the pointer is over the line or the focus is in it.
  let timer: ReturnType<typeof setTimeout> | undefined;
  let left = 0;
  let since = 0;
  let held = false;
  const ring = observableState<number | undefined>(undefined);
  const stop = () => {
    if (timer) clearTimeout(timer);
    timer = undefined;
  };
  const run = () => {
    stop();
    if (held || left <= 0) return;
    since = performance.now();
    timer = setTimeout(() => props.onDismiss?.(), left);
  };
  const hold = (on: boolean) => {
    if (on === held) return;
    held = on;
    if (on && timer) left -= performance.now() - since;
    if (on) stop();
    else run();
  };
  if (props.autoDismiss && props.onDismiss) {
    effect(props.autoDismiss, (ms) => {
      stop();
      left = ms ?? 0;
      ring.set(ms);
      run();
    });
    onCleanup(stop);
  }
  const body = (
    <>
      {props.icon}
      <span class={classes.text}>{props.text}</span>
      {props.after}
    </>
  );
  return (
    <span
      class={classes.pill}
      onPointerEnter={() => hold(true)}
      onPointerLeave={() => hold(false)}
      onFocusIn={() => hold(true)}
      onFocusOut={(e) => {
        if (!e.currentTarget.contains(e.relatedTarget as Node | null)) hold(false);
      }}
    >
      {props.onClick ? (
        <button type="button" class={classes.pillMain} aria-label={props.label} onClick={() => props.onClick?.()}>
          {body}
        </button>
      ) : (
        body
      )}
      {props.onDismiss ? (
        <button
          type="button"
          class={classes.pillButton}
          aria-label={props.dismissLabel ?? 'Dismiss'}
          data-tip={props.dismissLabel ?? 'Dismiss'}
          onClick={() => props.onDismiss?.()}
        >
          <Icon name="x" size={14} />
          <Show when={ring.asReadonly}>
            {(ms) => (
              // The time left, drawn down around the ✕ (the line's hover holds it, as the timer).
              <svg
                class={classes.countdown}
                viewBox="0 0 20 20"
                aria-hidden="true"
                style={{ '--msp-countdown': `${ms}ms` }}
              >
                <circle cx="10" cy="10" r="9" pathLength="100" />
              </svg>
            )}
          </Show>
        </button>
      ) : null}
    </span>
  );
}

/** The bottom centre of the canvas: what the person works on now (selection, focus) as chips. */
export function ContextChips(props: { status: NotificationsManager; lead?: Child }) {
  return (
    <div class={classes.context}>
      {props.lead}
      <For each={props.status.chips}>
        {(chip) => (
          <span class={`${classes.pill} ${classes.chip}`}>
            <button type="button" class={classes.text} onClick={() => chip.activate?.()}>
              {chip.label}
            </button>
            {chip.dismiss ? (
              <button type="button" class={classes.pillButton} title="Clear" onClick={() => chip.dismiss?.()}>
                <Icon name="x" size={14} />
              </button>
            ) : null}
          </span>
        )}
      </For>
    </div>
  );
}

function TaskPill(props: { task: StatusTask }) {
  const { task } = props;
  return (
    <span class={classes.pill}>
      <Icon name="loader-circle" size={14} class={classes.spin} />
      <span class={classes.text}>{task.label}</span>
      {task.progress !== undefined ? (
        <span class={classes.progress}>
          <div style={{ width: `${Math.round(task.progress * 100)}%` }} />
        </span>
      ) : null}
      {task.abort ? (
        <button type="button" class={classes.pillButton} title="Cancel" onClick={() => task.abort?.()}>
          <Icon name="x" size={14} />
        </button>
      ) : null}
    </span>
  );
}

/**
 * The compact rail's button (`ShellOptions.compactRail`): the rail's head glyph, in a card of its own at the bottom of
 * the rail's column (centered on its icons, at the viewport's edge); it shows and hides the whole rail column.
 */
export function RailToggle(props: { shell: ShellModel }) {
  const { shell } = props;
  const icon = derive(
    shell.dock('right').tabs,
    (tabs) => tabs.find((t) => t.place === 'head')?.icon ?? ('prism' as const),
  );
  let button!: HTMLElement;
  let face!: HTMLElement;
  const content = (
    <div class={`${classes.stack} ${classes.railToggle}`}>
      <div class={classes.frame}>
        <button
          type="button"
          class={classes.tool}
          data-tip-side="top"
          ref={(el) => (button = el)}
          onClick={() => shell.toggleRail()}
        >
          <span class={classes.face} ref={(el) => (face = el)} />
        </button>
      </div>
    </div>
  );
  let shown: IconName | undefined;
  watch([icon, shell.railOpen], () => {
    const open = shell.railOpen.value;
    const title = open ? 'Hide the panel rail' : 'Show the panel rail';
    setAttr(button, 'data-tip', title);
    setAttr(button, 'aria-label', title);
    setAttr(button, 'aria-pressed', String(open));
    if (icon.value !== shown) {
      shown = icon.value;
      face.replaceChildren(Icon({ name: shown, size: 18 }) as Node);
    }
  });
  return content;
}
