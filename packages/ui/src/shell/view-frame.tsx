/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../components/button/button.js';
import { type ContextMenuItem, openContextMenu } from '../components/context-menu/context-menu.js';
import { Icon, type IconName } from '../components/icon/icon.js';
import { ScrollArea } from '../components/scroll-area/scroll-area.js';
import {
  type Child,
  type ElementProps,
  type MaybeReactive,
  onCleanup,
  requireScope,
  type ViewRenderer,
} from '../dom/index.js';
import { classes } from './_generated/view.module.css.js';
import type { ViewRef } from './layout.js';
import type { ShellModel, ShellViewContext, ShellViewDef, ViewPlacement } from './shell-model.js';

/** The context a view of `ref` is rendered with; `placement` is read where it is now (a view moves). */
export function viewContext(
  shell: ShellModel,
  ref: ViewRef,
  placement: ViewPlacement | (() => ViewPlacement),
): ShellViewContext {
  const i = ref.indexOf('#');
  return {
    shell,
    ref,
    key: i < 0 ? undefined : ref.slice(i + 1),
    subject: shell.subject(ref),
    focus: shell.focus(ref),
    get placement() {
      return typeof placement === 'function' ? placement() : placement;
    },
    close: () => shell.close(ref),
  };
}

/**
 * A view's body (styles §8.4): one `ScrollArea` holding it, or, `fill`, a box it fills and scrolls itself; then the slot
 * its `ViewFooter` moves into, under the scroll (it stays while the content scrolls, out of the scroll's fade and bar).
 * The slot is the body's own, so a footer goes with its view when the view moves (a body is kept across moves).
 */
export function viewBody(def: ShellViewDef, ctx: ShellViewContext): HTMLElement {
  return (
    <div class={classes.bodyWrap} data-view-body="">
      {def.fill ? (
        <div class={classes.fill}>
          <ViewHost renderer={def.render} ctx={ctx} class={classes.fillHost} />
        </div>
      ) : (
        <ScrollArea class={classes.body} contentClass={classes.bodyContent}>
          <ViewHost renderer={def.render} ctx={ctx} class={classes.bodyHost} />
        </ScrollArea>
      )}
      <div class={classes.footerSlot} data-view-footer="" />
    </div>
  ) as HTMLElement;
}

/**
 * Mounts a contributed view (a `ViewRenderer`: this runtime's, a custom element, React's) into an element of its own,
 * in a child scope of the current one: paused with it (a folded section, a covered view), disposed with it.
 */
export function ViewHost<C>(props: { renderer: ViewRenderer<C>; ctx: C; class?: string }) {
  const owner = requireScope('ViewHost');
  const scope = owner.child();
  const host = (<div class={props.class ? `${classes.host} ${props.class}` : classes.host} />) as HTMLElement;
  const dispose = props.renderer.mount(host, props.ctx, scope);
  onCleanup(() => {
    dispose();
    scope.dispose();
  });
  return host;
}

export interface ViewFrameProps {
  readonly shell: ShellModel;
  readonly ref: ViewRef;
  readonly def: ShellViewDef;
  readonly placement: ViewPlacement;
  /** `rail`: the header is a rail button tall (a docked group beside the rail, styles §8); default 32. */
  readonly head?: 'rail' | 'header';
  /** In place of the title: a strip of tabs (a group of several views not at the rail). */
  readonly tabs?: Child;
  /**
   * The layout menu's entries (Float, Dock, Split below, …), asked for when it opens: a right-click or Shift+F10 on the
   * header; with `menuButton`, a ⋯ button too.
   */
  readonly menu?: () => readonly ContextMenuItem[];
  readonly menuButton?: boolean;
  /** Its ✕: closes it (a card, a window, a page) or folds it (a docked view). */
  readonly onClose?: () => void;
  /** The header's element, for a window to drag by. */
  readonly headerRef?: (el: HTMLElement) => void;
  /** The body, made elsewhere (the shell keeps a view's body while it moves); default: made here (`viewBody`). */
  readonly body?: Node;
  readonly class?: string;
  readonly id?: string;
}

/**
 * The frame every view is shown in (styles §8): a header (its title or its tabs, the view's own
 * actions, its menu, ✕) and a body that scrolls (one `ScrollArea`), or that the view fills (`fill`, a view that
 * scrolls its own parts). The same in a docked group, a window, a card, a page and a modal, so a view looks the same
 * wherever a person puts it.
 */
export function ViewFrame(props: ViewFrameProps) {
  const { shell, ref, def } = props;
  const ctx = viewContext(shell, ref, props.placement);
  const titleId = `${props.id ?? `msui-view-${ref.replace(/[^\w-]/g, '_')}`}-title`;
  const heading = def.heading ?? def.title;
  const landmark = props.placement === 'docked' || props.placement === 'sheet';
  let menuButton!: HTMLElement;
  const openMenu = (at?: { x: number; y: number }, from?: Element) => {
    const items = props.menu?.();
    if (!items?.length) return;
    const anchor = from ?? menuButton;
    const r = anchor.getBoundingClientRect();
    openContextMenu(
      anchor.closest<HTMLElement>('.msui-root') ?? document.body,
      at ?? { x: r.left, y: r.bottom },
      items,
    );
  };
  return (
    // A region (a landmark) where nothing around it names it: docked and in the sheet; a window, card, page and modal
    // are named by their host (a dialog, a region), which a second landmark of the same name would only repeat.
    (
      <section
        class={props.class ? `${classes.frame} ${props.class}` : classes.frame}
        role={landmark ? undefined : 'none'}
        aria-labelledby={landmark && !props.tabs ? titleId : undefined}
        aria-label={landmark && props.tabs ? def.title : undefined}
        data-placement={props.placement}
        data-msui-view={ref}
      >
        <header
          class={classes.header}
          data-head={props.head ?? 'header'}
          ref={(el: HTMLElement) => props.headerRef?.(el)}
          onContextMenu={(e: MouseEvent) => {
            // A tab's own menu is the tab's.
            if (!props.menu || (e.target as Element).closest('[role="tab"]')) return;
            e.preventDefault();
            openMenu({ x: e.clientX, y: e.clientY }, e.currentTarget as Element);
          }}
          onKeyDown={(e: KeyboardEvent) => {
            if (!props.menu || !((e.key === 'F10' && e.shiftKey) || e.key === 'ContextMenu')) return;
            if ((e.target as Element).closest('[role="tab"]')) return;
            e.preventDefault();
            openMenu(undefined, e.target as Element);
          }}
        >
          {props.tabs ?? (
            // A heading that is a view (a wordmark) is drawn as it is: not clipped, not faded (its own size).
            <h2
              class={classes.title}
              id={titleId}
              title={def.title}
              data-view={typeof heading === 'function' ? '' : undefined}
              aria-label={typeof heading === 'function' ? def.title : undefined}
            >
              {typeof heading === 'function' ? heading() : heading}
            </h2>
          )}
          {def.actions ? <ViewHost renderer={def.actions} ctx={ctx} class={classes.actions} /> : null}
          {props.menu && props.menuButton ? (
            <span
              class={classes.headerButton}
              ref={(el: HTMLElement) => {
                menuButton = el;
                el.firstElementChild?.setAttribute('aria-haspopup', 'menu');
              }}
            >
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name="ellipsis" />}
                title={`${def.title}: layout`}
                onClick={openMenu}
              />
            </span>
          ) : null}
          {props.onClose ? (
            <Button
              variant="ghost"
              size="sm"
              icon={<Icon name="x" />}
              title={`Close ${def.title}`}
              class={classes.headerButton}
              onClick={() => props.onClose?.()}
            />
          ) : null}
        </header>
        {props.body ?? viewBody(def, ctx)}
      </section>
    ) as HTMLElement
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// What a view is composed of (styles §8): so a view, an extension's too, never sets its own padding,
// gaps or type.

/** A run of a view's content with the view's rhythm: blocks a step apart (`--msui-space-2`), at the view's inset. */
export function ViewStack(props: { children: Child; gap?: 'sections' | 'blocks'; class?: string }) {
  return (
    <div class={props.class ? `${classes.stack} ${props.class}` : classes.stack} data-gap={props.gap ?? 'blocks'}>
      {props.children}
    </div>
  );
}

/** A row of tools for the whole view (a search, filters, Add), above its content: a row high, at the inset. */
export function ViewToolbar(props: { children: Child; label?: string; class?: string }) {
  return (
    <div
      class={props.class ? `${classes.toolbar} ${props.class}` : classes.toolbar}
      role="toolbar"
      aria-label={props.label}
    >
      {props.children}
    </div>
  );
}

/**
 * What an empty view says (styles §8): what goes here and how to get it, one line, an action under it; at the
 * inset, not centered in a void.
 */
export function ViewEmpty(props: { children: Child; action?: Child }) {
  return (
    <div class={classes.empty}>
      <p class={classes.emptyText}>{props.children}</p>
      {props.action ? <div class={classes.emptyAction}>{props.action}</div> : null}
    </div>
  );
}

/** An error in the content (not a toast): a line at the inset, in the danger color, with what to do. */
export function ViewError(props: { children: Child; action?: Child }) {
  return (
    <div class={classes.error} role="alert">
      <Icon name="circle-alert" />
      <span>{props.children}</span>
      {props.action ?? null}
    </div>
  );
}

/**
 * A view's footer of actions (an `ActionFooter`, Apply / Cancel): under the body's scroll, at the view's bottom, while
 * the content scrolls (styles §8.4). Written anywhere in the view, it moves itself into its body's footer slot once in
 * place, and goes with the view's scope. Outside a shell view (no slot) it stays where it is written.
 */
export function ViewFooter(props: { children: Child }) {
  const el = (<div class={classes.footer}>{props.children}</div>) as HTMLElement;
  queueMicrotask(() => {
    const slot = el.closest('[data-view-body]')?.querySelector(':scope > [data-view-footer]');
    if (slot && el.parentElement !== slot) slot.append(el);
  });
  onCleanup(() => el.remove());
  return el;
}

/**
 * A list or tree row (styles §8.5): an icon slot (kept when empty), the label, then a muted description on the same line,
 * what trails (a count, an action). Its label and description fade where they are cut off, and the row carries the
 * shell's side tooltip: while either is cut off, a card beside the view shows both whole (`data-tip-row`, from a pointer
 * over the text or the keyboard's focus on the row). `level` indents it a slot per level (a tree).
 */
export function ViewRow(
  props: {
    icon?: IconName | Child;
    label: MaybeReactive<string>;
    description?: MaybeReactive<string>;
    trailing?: Child;
    level?: number;
    class?: string;
  } & Omit<ElementProps<HTMLDivElement>, 'class' | 'children'>,
) {
  const { icon, label, description, trailing, level, class: cls, ...rest } = props;
  return (
    <div
      {...rest}
      class={cls ? `msui-row ${classes.row} ${cls}` : `msui-row ${classes.row}`}
      data-tip-row=""
      style={level ? `--msui-level: ${level}` : undefined}
    >
      <span class="msui-slot">{typeof icon === 'string' ? <Icon name={icon as IconName} /> : (icon ?? null)}</span>
      <span class={classes.rowText}>
        <span class={classes.rowLabel} data-tip-part="">
          {label}
        </span>
        {description !== undefined ? (
          <span class={classes.rowDescription} data-tip-part="">
            {description}
          </span>
        ) : null}
      </span>
      {trailing !== undefined ? <span class={classes.rowTrailing}>{trailing}</span> : null}
    </div>
  );
}

/** A small heading inside a view (a group of rows without a fold): xs, muted, on the text column. */
export function ViewHeading(props: { children: Child }) {
  return <h3 class={classes.heading}>{props.children}</h3>;
}
