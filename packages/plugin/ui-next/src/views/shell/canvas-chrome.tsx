/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, of } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { For, Show } from '../../dom/control-flow.js';
import { onCleanup } from '../../dom/scope.js';
import type { CanvasToolbarItem, OpenOverlay, ShellModel } from '../../shell/shell.js';
import type { MessageKind, StatusModel, StatusTask } from '../../shell/status.js';
import type { Child } from '../../dom/jsx-types.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/canvas.module.css.js';
import { ViewHost } from './view-host.js';

/** A section of the toolbar: its id (none for the leading items: the palette) and its buttons. */
interface ToolSection {
  group: string | undefined;
  items: readonly CanvasToolbarItem[];
}

function sectionsOf(items: readonly CanvasToolbarItem[]): ToolSection[] {
  const sections: ToolSection[] = [];
  for (const item of items) {
    const group = item.group;
    const section = sections.find((s) => s.group === group);
    if (section) (section.items as CanvasToolbarItem[]).push(item);
    else sections.push({ group, items: [item] });
  }
  // Items without a section (the palette) lead.
  return sections.sort((a, b) => (a.group === undefined ? -1 : b.group === undefined ? 1 : 0));
}

/**
 * The toolbar: its buttons in sections, each a glass block of its own, the gaps between them the grouping. `part`
 * picks what to show: the rail's items (`all`), split around its tab strip (`lead`: the palette above it, `rest`
 * below), or the canvas stack's items (`canvas`).
 */
export function CanvasToolbar(props: {
  shell: ShellModel;
  part?: 'all' | 'lead' | 'rest' | 'canvas';
  /** Where buttons with a menu open it (the canvas stack); without it, menus are not shown. */
  menus?: ToolMenus;
}) {
  const { shell } = props;
  const part = props.part ?? 'all';
  const canvas = part === 'canvas';
  const sections = shell.toolbar.pipe(
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
            aria-orientation="vertical"
            aria-label={section.group ?? 'commands'}
          >
            {section.items.map((item) =>
              item.visible ? (
                <Show when={item.visible}>{() => <ToolButton shell={shell} item={item} menus={props.menus} />}</Show>
              ) : (
                <ToolButton shell={shell} item={item} menus={props.menus} />
              ),
            )}
          </div>
        )}
      </For>
    </div>
  );
}

function ToolButton(props: { shell: ShellModel; item: CanvasToolbarItem; menus?: ToolMenus }) {
  const { shell, item, menus } = props;
  // A button with a menu shows the menu on hover instead of a tooltip (both would open in the same spot).
  const withMenu = menus && item.menu?.length ? menus : undefined;
  const overlayOpen = item.overlay ? shell.overlay.pipe(map((o) => o?.overlay.id === item.overlay)) : of(false);
  const pressed = combineLatest([overlayOpen, item.pressed ?? of(false)]).pipe(
    map(([a, b]) => a || b),
    distinctUntilChanged(),
  );
  return (
    <button
      type="button"
      class={classes.tool}
      data-tip={withMenu ? undefined : item.title}
      data-tip-side="left"
      aria-label={item.title}
      aria-pressed={pressed}
      aria-haspopup={withMenu ? 'menu' : undefined}
      data-overlay-trigger={item.overlay ? '' : undefined}
      onClick={(e) => {
        // A button that only holds a menu opens it; otherwise the click runs the button's own action.
        if (withMenu && !item.overlay && !item.run) return withMenu.show(item, e.currentTarget);
        withMenu?.hide();
        if (item.overlay) shell.toggleOverlay(item.overlay, e.currentTarget);
        else item.run?.(e.currentTarget);
      }}
      onPointerEnter={withMenu ? (e) => withMenu.show(item, e.currentTarget) : undefined}
      onPointerLeave={withMenu ? () => withMenu.hideSoon() : undefined}
      onFocus={withMenu ? (e) => withMenu.show(item, e.currentTarget) : undefined}
      onBlur={withMenu ? () => withMenu.hideSoon() : undefined}
    >
      <Icon name={item.icon} size={20} />
      {withMenu ? <span class={classes.menuMark} aria-hidden="true" /> : null}
    </button>
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
        {props.lead}
        <CanvasToolbar shell={props.shell} part="canvas" menus={menus} />
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
 * The open tool menu of the canvas stack: which button, at what height. Opens on hover or focus, closes a moment after
 * the pointer leaves both the button and the menu (so it can travel between them), on a choice, or on Escape.
 */
export class ToolMenus {
  readonly open = observableState<{ item: CanvasToolbarItem; top: number } | undefined>(undefined);
  stack: HTMLElement | undefined;
  #timer: ReturnType<typeof setTimeout> | undefined;

  show(item: CanvasToolbarItem, button: Element) {
    this.keep();
    const top = this.stack ? button.getBoundingClientRect().top - this.stack.getBoundingClientRect().top : 0;
    const current = this.open.value;
    if (current?.item !== item || current.top !== top) this.open.set({ item, top });
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
    menus.hide();
    e.preventDefault();
  };
  document.addEventListener('keydown', onKey);
  onCleanup(() => document.removeEventListener('keydown', onKey));
  return (
    <Show when={menus.open}>
      {(open) => (
        <div
          class={classes.toolMenu}
          style={{ top: `${open.top}px` }}
          role="menu"
          aria-label={open.item.title}
          onPointerEnter={() => menus.keep()}
          onPointerLeave={() => menus.hideSoon()}
          onFocusIn={() => menus.keep()}
          onFocusOut={() => menus.hideSoon()}
        >
          <div class={classes.toolMenuCard}>
            {open.item.menu!.map((m) => (
              <button
                type="button"
                role={m.pressed ? 'menuitemcheckbox' : 'menuitem'}
                aria-checked={m.pressed}
                class={classes.toolMenuRow}
                onClick={() => {
                  menus.hide();
                  m.run();
                }}
              >
                {m.icon ? <Icon name={m.icon} size={16} /> : <span />}
                <span class={classes.toolMenuLabel}>
                  {m.label}
                  {m.hint ? <span class={classes.toolMenuHint}>{m.hint}</span> : null}
                </span>
                {m.pressed ? (
                  <Show when={m.pressed}>{() => <Icon name="check" size={14} class={classes.toolMenuCheck} />}</Show>
                ) : null}
              </button>
            ))}
          </div>
        </div>
      )}
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
  const card = (
    <section class={classes.card} role="dialog" aria-label={open.overlay.title}>
      <header class={classes.cardHeader}>
        <span class={classes.cardTitle}>{open.overlay.title}</span>
        <Button variant="ghost" size="sm" icon={<Icon name="x" />} title="Close" onClick={() => shell.closeOverlay()} />
      </header>
      <div class={classes.cardBody}>
        <ViewHost renderer={open.overlay.render} ctx={{ shell, ...(open.focus ? { focus: open.focus } : {}) }} />
      </div>
    </section>
  ) as HTMLElement;
  // Esc, or a pointer down outside the card and its button, closes it.
  const doc = document;
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented) return;
    shell.closeOverlay();
    e.preventDefault();
  };
  const onPointer = (e: PointerEvent) => {
    const target = e.target as Element;
    if (card.contains(target) || open.anchor?.contains(target)) return;
    // Buttons that open another card, and panel tabs, swap the card on their click: closing it here first would show
    // the panel underneath for a frame (a blink), and a tab click would then collapse the panel.
    if (target.closest?.('[data-overlay-trigger], [role=tab]')) return;
    shell.closeOverlay();
  };
  doc.addEventListener('keydown', onKey);
  doc.addEventListener('pointerdown', onPointer, true);
  onCleanup(() => {
    doc.removeEventListener('keydown', onKey);
    doc.removeEventListener('pointerdown', onPointer, true);
  });
  return card;
}

export const KindIcons: Record<MessageKind, IconName> = {
  info: 'info',
  success: 'circle-check',
  warning: 'triangle-alert',
  error: 'circle-alert',
};

/** The bottom left of the canvas: the latest message and running tasks. */
export function StatusLine(props: { status: StatusModel; onOpenLog?: () => void }) {
  const { status } = props;
  return (
    <div class={classes.status} aria-live="polite">
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
      <For each={status.tasks}>{(task) => <TaskPill task={task} />}</For>
    </div>
  );
}

/** The bottom centre of the canvas: what the person works on now (selection, focus) as chips. */
export function ContextChips(props: { status: StatusModel }) {
  return (
    <div class={classes.context}>
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
