/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, type Observable, of } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { For, Show } from '../../dom/control-flow.js';
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
import { bind, onCleanup } from '../../dom/scope.js';
import type { CanvasMenuItem, CanvasToolbarItem, OpenOverlay, ShellModel } from '../../shell/shell.js';
import type { MessageKind, StatusModel, StatusTask } from '../../shell/status.js';
import type { Child } from '../../dom/jsx-types.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/canvas.module.css.js';
import { PinButton } from './dock-view.js';
import { sectionsOf } from './rail.js';
import { ViewHost } from './view-host.js';

/**
 * The toolbar: its buttons in sections, each a glass block of its own, the gaps between them the grouping. `part`
 * picks what to show: the rail's items (`all`), split around its tab strip (`lead`: the command menu above it, `rest`
 * below), or the canvas stack's items (`canvas`).
 */

/** A toolbar item's title as shown: its face's (it changes: the theme, the mode switch), or its own with its caption. */
export function itemTitle(item: CanvasToolbarItem) {
  return item.face
    ? item.face.pipe(map((f) => f.title))
    : item.caption
      ? `${item.title} (${item.caption})`
      : item.title;
}

/** A toolbar item's icon as shown: its face's, following it, or its own. */
export function ItemIcon(props: { item: CanvasToolbarItem; size: number; class?: string }) {
  const { item } = props;
  return item.face ? (
    <span
      class={props.class}
      ref={(el) => bind(item.face!, (f) => el.replaceChildren(Icon({ name: f.icon, size: props.size })))}
    />
  ) : (
    <Icon name={item.icon} size={props.size} />
  );
}

export function CanvasToolbar(props: {
  shell: ShellModel;
  part?: 'all' | 'lead' | 'rest' | 'canvas';
  /** Where buttons with a menu open it (the canvas stack); without it, menus are not shown. */
  menus?: ToolMenus;
  /** The rail's left-out buttons (`RailLayout.hidden`): they are not shown, they are under More. */
  hidden?: Observable<ReadonlySet<string>>;
  /** On the rail: a button with a menu opens it beside the rail on a click (the theme), as More does. */
  railMenu?: (item: CanvasToolbarItem, anchor: HTMLElement) => void;
}) {
  const { shell } = props;
  const part = props.part ?? 'all';
  const canvas = part === 'canvas';
  const items$ = props.hidden
    ? combineLatest([shell.toolbar, props.hidden]).pipe(
        map(([items, hidden]) => (hidden.size ? items.filter((i) => !hidden.has(`tool:${i.id}`)) : items)),
      )
    : shell.toolbar;
  const sections = items$.pipe(
    map((items) => sectionsOf(items.filter((i) => (i.place === 'canvas') === canvas))),
    map((all) =>
      part === 'all' || canvas
        ? all
        : all.filter((s) => (part === 'lead' ? s.group === undefined : s.group !== undefined)),
    ),
  );
  return (
    <div class={classes.sections}>
      <For each={sections} key={(s) => s.group ?? ''}>
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
  const overlayOpen = item.overlay ? shell.overlay.pipe(map((o) => o?.overlay.id === item.overlay)) : of(false);
  const pressed = combineLatest([overlayOpen, item.pressed ?? of(false)]).pipe(
    map(([a, b]) => a || b),
    distinctUntilChanged(),
  );
  // The menu button pattern (README §3.7): ↓ or ↑ opens the menu at its first item (Enter and Space too, on a button
  // that only holds a menu); the pointer opens it on hover. Focus alone does not open it: Tab goes on to the next button.
  const expanded = withMenu?.open.asReadonly.pipe(
    map((o) => o?.item === item),
    distinctUntilChanged(),
  );
  return (
    <button
      type="button"
      class={classes.tool}
      data-captioned={item.caption ? '' : undefined}
      data-tip={withMenu ? undefined : itemTitle(item)}
      data-tip-side={item.place === 'canvas' ? 'top' : 'left'}
      aria-label={itemTitle(item)}
      // Only a button with a state is a toggle (a button that cycles, like the theme, names its state in its label).
      aria-pressed={item.pressed || item.overlay ? pressed : undefined}
      aria-haspopup={withMenu || railMenu ? 'menu' : undefined}
      aria-expanded={expanded}
      data-overlay-trigger={item.overlay ? '' : undefined}
      data-keeps-card={item.keepsCard ? '' : undefined}
      onClick={(e) => {
        // A button that only holds a menu opens it (from the keyboard, at its first item: a click with no pointer
        // has `detail` 0); otherwise the click runs the button's own action.
        if (withMenu && !item.overlay && !item.run) return withMenu.show(item, e.currentTarget, e.detail === 0);
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
      onPointerEnter={withMenu ? (e) => withMenu.show(item, e.currentTarget) : undefined}
      onPointerLeave={withMenu ? () => withMenu.hideSoon() : undefined}
      onBlur={withMenu ? () => withMenu.hideSoon() : undefined}
    >
      <ItemIcon item={item} size={item.caption ? 16 : 20} class={classes.face} />
      {item.caption ? <span class={classes.caption}>{item.caption}</span> : null}
      {withMenu ? <span class={classes.menuMark} aria-hidden="true" /> : null}
    </button>
  );
}

/** A toolbar button's menu rows (the canvas stack's menus, the rail's): grouped, a toggle's or a choice's check. */
export function ToolMenuRows(props: { items: readonly CanvasMenuItem[]; onRun(): void }) {
  return (
    <>
      {props.items.map((m, i, all) => (
        <>
          {i > 0 && m.group !== all[i - 1].group ? <div class={classes.toolMenuDivider} role="separator" /> : null}
          <button
            type="button"
            role={m.pressed ? (m.choice ? 'menuitemradio' : 'menuitemcheckbox') : 'menuitem'}
            aria-checked={m.pressed}
            class={classes.toolMenuRow}
            disabled={!!m.disabled}
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
        </>
      ))}
    </>
  );
}

/**
 * The canvas stack at the top right of the canvas, below what floats at the top and beside the rail or the open
 * panel: one frame with the mode buttons (`lead`) and then the toolbar items placed on the canvas (camera, lighting,
 * XR, expand, sequence); `aside` is what opens beside it (the picking level).
 */
export function CanvasStack(props: { shell: ShellModel; lead?: Child; aside?: Child }) {
  const menus = new ToolMenus();
  const stack = (
    <div class={classes.stack}>
      <div class={classes.frame}>
        <ScrollArea axis="x" contentClass={classes.stackTools}>
          {props.lead}
          <CanvasToolbar shell={props.shell} part="canvas" menus={menus} />
        </ScrollArea>
      </div>
      {props.aside}
      <ToolMenu menus={menus} />
    </div>
  ) as HTMLElement;
  menus.stack = stack;
  onCleanup(() => menus.dispose());
  return stack;
}

/**
 * The open tool menu of the canvas stack: which button, at what height. Opens on hover or from the keyboard (its
 * button's ↓, ↑, Enter), closes a moment after the pointer leaves both the button and the menu (so it can travel
 * between them), when the focus leaves it, on a choice, or on Escape (the focus goes back to its button).
 */
export class ToolMenus {
  readonly open = observableState<{ item: CanvasToolbarItem; right: number; button: Element } | undefined>(undefined);
  stack: HTMLElement | undefined;
  /** The open menu's element (`ToolMenu` sets it). */
  menu: HTMLElement | undefined;
  #timer: ReturnType<typeof setTimeout> | undefined;

  /** Opens `item`'s menu by its button; `focus`: at its first item (opened from the keyboard). */
  show(item: CanvasToolbarItem, button: Element, focus = false) {
    this.keep();
    // The menu opens upward from the stack (a row in the canvas' corner), its right edge level with the button's.
    const right = this.stack ? this.stack.getBoundingClientRect().right - button.getBoundingClientRect().right : 0;
    const current = this.open.value;
    if (current?.item !== item || current.right !== right) this.open.set({ item, right, button });
    if (focus) queueMicrotask(() => this.menu && focusFirstItem(this.menu));
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
  document.addEventListener('keydown', onKey);
  onCleanup(() => document.removeEventListener('keydown', onKey));
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
        return (
          <div
            class={classes.toolMenu}
            style={{ right: `${open.right}px` }}
            role="menu"
            aria-label={open.item.title}
            ref={(x) => {
              el = x;
              menus.menu = x;
            }}
            onKeyDown={menuKeys}
            onPointerEnter={() => menus.keep()}
            onPointerLeave={() => menus.hideSoon()}
            onFocusIn={() => {
              hadFocus = true;
              menus.keep();
            }}
            onFocusOut={() => menus.hideSoon()}
          >
            <div class={classes.toolMenuCard}>
              <ToolMenuRows items={open.item.menu!} onRun={() => menus.hide()} />
            </div>
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
  const lead = shell.toolbar.pipe(map((all) => all.filter((i) => i.place !== 'canvas' && !i.group)));
  const actions = shell.toolbar.pipe(map((all) => all.filter((i) => i.place !== 'canvas' && i.group)));
  const tabs = right.tabs.pipe(map((all) => (all.some((t) => t.compactBar) ? all.filter((t) => !t.compactBar) : [])));
  let list: HTMLElement | undefined;
  const run = (item: CanvasToolbarItem, anchor: Element) => {
    open.set(false);
    if (item.overlay) shell.toggleOverlay(item.overlay, anchor);
    else item.run?.(anchor);
  };
  const root = (
    <div class={classes.railMenu}>
      <For each={lead} key={(i) => i.id}>
        {(item) => (
          <button
            type="button"
            class={classes.railMenuButton}
            aria-label={itemTitle(item)}
            data-tip={itemTitle(item)}
            onClick={(e) => run(item, e.currentTarget)}
          >
            <ItemIcon item={item} size={20} class={classes.face} />
          </button>
        )}
      </For>
      {/* More and its list: the list opens above the button (its anchor), inside the screen. */}
      <div class={classes.railMenuAnchor}>
        <button
          type="button"
          class={classes.railMenuButton}
          aria-label="More"
          aria-haspopup="menu"
          aria-expanded={open}
          data-tip="More"
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
                <For each={tabs} key={(t) => t.id}>
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
                <Show when={combineLatest([tabs, actions]).pipe(map(([t, a]) => t.length > 0 && a.length > 0))}>
                  {() => <div class={classes.railMenuDivider} />}
                </Show>
                <For each={actions} key={(i) => i.id}>
                  {(item) => (
                    <button
                      type="button"
                      role="menuitem"
                      class={classes.toolMenuRow}
                      // As on the rail: an item shown only while it applies (Play snapshots with two or more).
                      hidden={item.visible ? item.visible.pipe(map((v) => !v)) : false}
                      onClick={(e) => run(item, e.currentTarget)}
                    >
                      <ItemIcon item={item} size={16} class={classes.face} />
                      <span class={classes.toolMenuLabel}>{itemTitle(item)}</span>
                      <span />
                    </button>
                  )}
                </For>
              </div>
            );
          }}
        </Show>
      </div>
    </div>
  ) as HTMLElement;
  const onDown = (e: PointerEvent) => {
    if (open.value && !isInside(e, root)) open.set(false);
  };
  document.addEventListener('pointerdown', onDown);
  onCleanup(() => document.removeEventListener('pointerdown', onDown));
  return root;
}

/**
 * A task in overlay mode (rendering an animation) holds the plugin: a scrim over everything and a card with what is
 * running, its progress, and Cancel.
 */
export function BusyOverlay(props: { status: StatusModel }) {
  const task = props.status.tasks.pipe(
    map((tasks) => tasks.find((t) => t.overlay)),
    distinctUntilChanged((a, b) => a?.id === b?.id),
  );
  return (
    <Show when={task}>
      {(first) => {
        const live = props.status.tasks.pipe(map((tasks) => tasks.find((t) => t.id === first.id) ?? first));
        const percent = live.pipe(map((t) => (t.progress === undefined ? undefined : Math.round(t.progress * 100))));
        return (
          <div class={classes.busy} role="alertdialog" aria-modal="true" aria-label={first.label}>
            <div class={classes.busyCard}>
              <div class={classes.busyHead}>
                <Icon name="loader-circle" size={18} class={classes.busySpin} />
                <span class={classes.busyLabel}>{live.pipe(map((t) => t.label))}</span>
                <span class={classes.busyPercent}>{percent.pipe(map((p) => (p === undefined ? '' : `${p}%`)))}</span>
              </div>
              <div
                class={classes.busyTrack}
                role="progressbar"
                aria-valuemin={0}
                aria-valuemax={100}
                aria-valuenow={percent}
              >
                <span
                  class={classes.busyBar}
                  data-indeterminate={percent.pipe(map((p) => (p === undefined ? '' : undefined)))}
                  style={{ scale: percent.pipe(map((p) => `${(p ?? 30) / 100} 1`)) }}
                />
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
        );
      }}
    </Show>
  );
}

/** The open overlay as a card in the right panel's place (the panel steps aside while it is open). */
export function OverlayHost(props: { shell: ShellModel }) {
  const { shell } = props;
  return <Show when={shell.overlay}>{(open) => <OverlayCard shell={shell} open={open} />}</Show>;
}

function OverlayCard(props: { shell: ShellModel; open: OpenOverlay }) {
  const { shell, open } = props;
  // A card is a dialog (README §3.7: cards trap and restore focus): opened from the keyboard, the focus goes into it;
  // Tab stays inside while it floats (pinned, it is part of the workspace); closed, the focus goes back to its button.
  const restore = rememberFocus(() => open.anchor as HTMLElement | undefined);
  const byKeyboard = focusedByKeyboard(open.anchor);
  const titleId = uniqueId('msp-card-title');
  const pinnedState = shell.dock('right').state.pipe(
    map((s) => s.pinned),
    distinctUntilChanged(),
  );
  const card = (
    <section
      class={classes.card}
      role="dialog"
      aria-modal={pinnedState.pipe(map((p) => !p))}
      aria-labelledby={titleId}
      onKeyDown={(e) => {
        if (!pinned()) trapTab(e, e.currentTarget);
      }}
    >
      <header class={classes.cardHeader}>
        <h2 class={classes.cardTitle} id={titleId}>
          {open.overlay.title}
        </h2>
        <PinButton shell={shell} dock="right" />
        <Button variant="ghost" size="sm" icon={<Icon name="x" />} title="Close" onClick={() => shell.closeOverlay()} />
      </header>
      <ScrollArea class={classes.cardBody} contentClass={classes.cardBodyContent}>
        <ViewHost renderer={open.overlay.render} ctx={{ shell, ...(open.focus ? { focus: open.focus } : {}) }} />
      </ScrollArea>
    </section>
  ) as HTMLElement;
  // Esc, or a pointer down outside the card and its button, closes it; not one on the canvas: the person turns and
  // zooms the view with the card open (setting up an export, trying a canvas setting). Pinned, the card is part of the workspace like the pinned panel: only its ✕ (or another card or tab) closes it.
  const pinned = () => shell.dock('right').state.value.pinned;
  const doc = document;
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented || pinned() || !ownsKey(card, e)) return;
    shell.closeOverlay();
    e.preventDefault();
  };
  const onPointer = (e: PointerEvent) => {
    const target = eventTarget(e) as Element;
    if (pinned() || isInside(e, card) || isInside(e, open.anchor)) return;
    // Buttons that open another card, and panel tabs, swap the card on their click: closing it here first would show
    // the panel underneath for a frame (a blink), and a tab click would then collapse the panel.
    if (target.closest?.('[data-overlay-trigger], [role=tab]')) return;
    // What acts on the open card leaves it open: the command menu and its button.
    if (target.closest?.('[data-keeps-card], [data-msp-canvas]')) return;
    shell.closeOverlay();
  };
  doc.addEventListener('keydown', onKey);
  doc.addEventListener('pointerdown', onPointer, true);
  // Into the card: its body's first control once the view is there (a frame), else its Close.
  let frame = 0;
  if (byKeyboard) {
    frame = requestAnimationFrame(() => {
      frame = 0;
      const body = card.querySelector(`.${classes.cardBody}`);
      const all = focusables(card);
      (all.find((el) => body?.contains(el)) ?? all[all.length - 1])?.focus();
    });
  }
  onCleanup(() => {
    if (frame) cancelAnimationFrame(frame);
    doc.removeEventListener('keydown', onKey);
    doc.removeEventListener('pointerdown', onPointer, true);
    restore(card);
  });
  return card;
}

export const KindIcons: Record<MessageKind, IconName> = {
  info: 'info',
  success: 'circle-check',
  warning: 'triangle-alert',
  error: 'circle-alert',
};

/** The bottom left of the canvas: the latest message, running tasks, and under them `after` (the assistant's turn). */
export function StatusLine(props: { status: StatusModel; onOpenLog?: () => void; after?: Child }) {
  const { status } = props;
  // One polite live region for results (README §3.7): the latest message only, not the tasks' progress or the
  // assistant's streamed words (its card says when a turn is done).
  return (
    <div class={classes.status}>
      <div class={classes.announce} role="status" aria-live="polite">
        {status.spoken}
      </div>
      <Show when={status.current}>
        {(m) => (
          <button
            type="button"
            class={`${classes.pill} ${classes.message}`}
            data-kind={m.kind}
            title="Show the log"
            onClick={() => props.onOpenLog?.()}
          >
            <Icon name={KindIcons[m.kind]} size={14} />
            <span class={classes.text}>{m.text}</span>
          </button>
        )}
      </Show>
      <For each={status.tasks.pipe(map((tasks) => tasks.filter((t) => !t.overlay)))}>
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
}) {
  const body = (
    <>
      {props.icon}
      <span class={classes.text}>{props.text}</span>
      {props.after}
    </>
  );
  return (
    <span class={classes.pill}>
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
        </button>
      ) : null}
    </span>
  );
}

/** The bottom centre of the canvas: what the person works on now (selection, focus) as chips. */
export function ContextChips(props: { status: StatusModel; lead?: Child }) {
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
  const icon = shell.dock('right').tabs.pipe(
    map((tabs) => tabs.find((t) => t.place === 'head')?.icon ?? ('prism' as const)),
    distinctUntilChanged(),
  );
  const title = shell.railOpen.pipe(map((open) => (open ? 'Hide the panel rail' : 'Show the panel rail')));
  return (
    <div class={`${classes.stack} ${classes.railToggle}`}>
      <div class={classes.frame}>
        <button
          type="button"
          class={classes.tool}
          data-tip={title}
          data-tip-side="top"
          aria-label={title}
          aria-pressed={shell.railOpen}
          onClick={() => shell.toggleRail()}
        >
          <span class={classes.face} ref={(el) => bind(icon, (name) => el.replaceChildren(Icon({ name, size: 18 })))} />
        </button>
      </div>
    </div>
  );
}
