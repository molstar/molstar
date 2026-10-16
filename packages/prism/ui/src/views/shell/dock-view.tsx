/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, of } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { activity, For, Show } from '../../dom/control-flow.js';
import { eventTarget, isInside, ownsKey } from '../../dom/events.js';
import { deepActiveElement, focusFirstItem, menuKeys, rememberFocus, rovingTabs, uniqueId } from '../../dom/focus.js';
import { onCleanup } from '../../dom/scope.js';
import type { Child } from '../../dom/jsx-types.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import type { DockId, DockModel, ShellModel, ShellTab } from '../../shell/shell.js';
import { classes as canvasClasses } from './_generated/canvas.module.css.js';
import { classes } from './_generated/dock.module.css.js';
import { ItemIcon, itemTitle, ToolMenuRows } from './canvas-chrome.js';
import { type RailEntry, RailLayout } from './rail.js';
import { ViewHost } from './view-host.js';

/** A side dock: its tab strip on the outer edge and the active tab's content next to the canvas. */
export function DockView(props: {
  shell: ShellModel;
  dock: DockModel /** After the tabs in the strip (compact: More). */;
  stripEnd?: Child;
}) {
  const { shell, dock } = props;
  const horizontal = dock.id === 'top' || dock.id === 'bottom';
  const hasTabs = dock.tabs.pipe(map((tabs) => tabs.length > 0));
  const compact = shell.mode.pipe(map((m) => m === 'compact' && dock.id === 'right'));
  // Compact: only the tabs that ask for the bar (the rest are under More); elsewhere all of them.
  const compactTabs = combineLatest([dock.tabs, compact]).pipe(
    map(([tabs, c]) => (c && tabs.some((t) => t.compactBar) ? tabs.filter((t) => t.compactBar) : tabs)),
  );
  // A top or bottom dock with one tab needs no strip: its content collapses itself, and a slim bar expands it.
  const showStrip = horizontal
    ? combineLatest([dock.tabs, dock.active]).pipe(map(([tabs, active]) => tabs.length > 1 && !!active))
    : of(true);
  const panelId = uniqueId('msp-panel');
  return (
    <Show when={hasTabs}>
      {() => (
        <div
          class={`${classes.dock} ${classes[dock.id]}`}
          classList={{ [classes.compact]: compact }}
          style={{ '--msp-dock-size': dock.state.pipe(map((s) => `${s.size}px`)) }}
        >
          <Show
            when={dock.active}
            reveal
            fallback={horizontal ? () => <CollapsedBar shell={shell} dock={dock} /> : undefined}
          >
            {(tab) => <DockContent shell={shell} dock={dock} tab={tab} id={panelId} />}
          </Show>
          <Show when={showStrip}>
            {() => (
              <div
                class={classes.strip}
                role="tablist"
                aria-label="Panels"
                aria-orientation={compact.pipe(map((c) => (c || horizontal ? 'horizontal' : 'vertical')))}
                ref={(el) => onCleanup(rovingTabs(el))}
              >
                <For each={compactTabs}>
                  {(tab) => <TabButton shell={shell} dock={dock} tab={tab} panel={panelId} />}
                </For>
                {props.stripEnd}
              </div>
            )}
          </Show>
        </div>
      )}
    </Show>
  );
}

/**
 * A dock floating over the canvas (the right dock by default): a full-height rail at the edge with the tab strip, the
 * rail's toolbar sections below it, and its foot (the command menu) at the bottom; the active tab's content as a glass card
 * beside it. Collapsed, only the rail stays.
 */
export function FloatingDock(props: {
  shell: ShellModel;
  dock: DockModel;
  railLead?: Child;
  /** The rail's toolbar sections; it gets the layout to leave out what does not fit (`CanvasToolbar`'s `hidden`). */
  railItems?: (rail: RailLayout) => Child;
  /** At the bottom of the rail (the command menu). */
  railFoot?: (rail: RailLayout) => Child;
}) {
  const { shell, dock } = props;
  const size = { '--msp-dock-size': dock.state.pipe(map((s) => `${s.size}px`)) };
  // What does not fit the rail's height, even at its smallest size, is under More (rail.ts).
  const rail = new RailLayout(shell, dock);
  // Tabs at the head of the rail (Home), in the strip, and at its foot (Settings, by the command menu).
  const shown = combineLatest([dock.tabs, rail.hidden]).pipe(
    map(([tabs, hidden]) => (hidden.size ? tabs.filter((t) => !hidden.has(`tab:${t.id}`)) : tabs)),
  );
  const head = shown.pipe(map((tabs) => tabs.filter((t) => t.place === 'head')));
  const strip = shown.pipe(map((tabs) => tabs.filter((t) => !t.place || t.place === 'strip')));
  const foot = shown.pipe(map((tabs) => tabs.filter((t) => t.place === 'foot')));
  const panelId = uniqueId('msp-panel');
  if (dock.id === 'top' || dock.id === 'bottom') {
    // A top or bottom card spans the canvas between the edges and what floats on the sides. Collapsed, it is gone: a
    // canvas toolbar button (Sequence) brings it back.
    return (
      <div
        class={`${classes.floating} ${dock.id === 'top' ? classes.floatingTop : classes.floatingBottom}`}
        classList={{ [classes.band]: dock.active.pipe(map((tab) => tab?.frame === 'none')) }}
        style={size}
      >
        <Show when={dock.active} reveal>
          {(tab) => <DockContent shell={shell} dock={dock} tab={tab} id={panelId} />}
        </Show>
      </div>
    );
  }
  const root = (
    <div
      class={`${classes.floating} ${dock.id === 'left' ? classes.floatingLeft : classes.floatingRight}`}
      style={size}
    >
      <div
        class={`${classes.rail} ${canvasClasses.frame}`}
        ref={(el) => {
          onCleanup(rail.observe(el));
        }}
      >
        {props.railLead}
        <Show when={head.pipe(map((tabs) => tabs.length > 0))}>
          {() => (
            <div
              class={classes.strip}
              role="tablist"
              aria-orientation="vertical"
              aria-label="Home"
              ref={(el) => onCleanup(rovingTabs(el))}
            >
              <For each={head}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} panel={panelId} />}</For>
            </div>
          )}
        </Show>
        <Show when={strip.pipe(map((tabs) => tabs.length > 0))}>
          {() => (
            <div
              class={classes.strip}
              role="tablist"
              aria-orientation="vertical"
              aria-label="Panels"
              ref={(el) => onCleanup(rovingTabs(el))}
            >
              <For each={strip}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} panel={panelId} />}</For>
            </div>
          )}
        </Show>
        {props.railItems?.(rail)}
        <RailMore shell={shell} dock={dock} rail={rail} root={() => root} />
        <div class={classes.railFoot}>
          <Show when={foot.pipe(map((tabs) => tabs.length > 0))}>
            {() => (
              <div
                class={classes.strip}
                role="tablist"
                aria-orientation="vertical"
                aria-label="More panels"
                ref={(el) => onCleanup(rovingTabs(el))}
              >
                <For each={foot}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} panel={panelId} />}</For>
              </div>
            )}
          </Show>
          {props.railFoot?.(rail)}
        </div>
      </div>
      <Show when={dock.active} reveal="side">
        {(tab) => <DockContent shell={shell} dock={dock} tab={tab} id={panelId} pinnable coveredBy={shell.overlay} />}
      </Show>
      <RailMoreMenu shell={shell} dock={dock} rail={rail} />
    </div>
  ) as HTMLElement;
  return root;
}

/** Runs a tab's button: selects it, or, with an overlay card in the panel's place, closes that and shows the tab. */
function activateTab(shell: ShellModel, dock: DockModel, tab: ShellTab) {
  if (dock.id === 'right' && shell.overlay.value) {
    shell.closeOverlay();
    shell.show(tab.id);
  } else shell.select(tab.id);
}

/**
 * The rail's More button: where the buttons that do not fit the rail's height go (rail.ts). It is the last button of the
 * rail's sections, above the foot; its menu (`RailMoreMenu`) opens to the left, outside the rail's clipping.
 */
function RailMore(props: { shell: ShellModel; dock: DockModel; rail: RailLayout; root: () => HTMLElement }) {
  const { rail } = props;
  const menu = rail.menu;
  return (
    <Show when={rail.overflow.pipe(map((o) => o.length > 0))}>
      {() => {
        // Gone (the rail grew): the menu is closed, it does not come back open.
        onCleanup(() => menu.set(undefined));
        return (
          <nav class={classes.strip} aria-label="More buttons">
            <button
              type="button"
              class={classes.tab}
              data-tip="More"
              data-tip-side="left"
              aria-label="More"
              aria-haspopup="menu"
              aria-expanded={menu.asReadonly.pipe(map((m) => !!m))}
              onClick={(e) => rail.toggleMenu(e.currentTarget)}
              onKeyDown={(e) => {
                // ↑ or ↓ opens it too (the menu button pattern); it opens at its first item.
                if ((e.key !== 'ArrowDown' && e.key !== 'ArrowUp') || menu.value) return;
                e.preventDefault();
                e.currentTarget.click();
              }}
            >
              <Icon name="ellipsis" size={20} />
            </button>
          </nav>
        );
      }}
    </Show>
  );
}

/** The More menu: a row for each button left out of the rail, in the rail's order; like the canvas menus, a glass card. */
function RailMoreMenu(props: { shell: ShellModel; dock: DockModel; rail: RailLayout }) {
  const { shell, dock, rail } = props;
  const menu = rail.menu;
  // Open only while there is something in it (More's: the rail grew, nothing is left out).
  const open = combineLatest([menu.asReadonly, rail.overflow]).pipe(
    map(([m, o]) => (m && (m.item ? m.item.menu?.length : o.length > 0) ? m : undefined)),
  );
  const run = (entry: RailEntry, anchor: Element) => {
    menu.set(undefined);
    if (entry.tab) return activateTab(shell, dock, entry.tab);
    const item = entry.item!;
    if (item.overlay) shell.toggleOverlay(item.overlay, anchor);
    else item.run?.(anchor);
  };
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || !menu.value || !ownsKey(menu.value.anchor, e)) return;
    const anchor = menu.value.anchor as HTMLElement;
    menu.set(undefined);
    anchor.focus();
    e.preventDefault();
  };
  const onDown = (e: PointerEvent) => {
    const m = menu.value;
    if (!m) return;
    const target = eventTarget(e) as Node;
    if (isInside(e, m.anchor) || (target as Element).closest?.(`.${classes.railMenu}`)) return;
    menu.set(undefined);
  };
  document.addEventListener('keydown', onKey);
  document.addEventListener('pointerdown', onDown);
  onCleanup(() => {
    document.removeEventListener('keydown', onKey);
    document.removeEventListener('pointerdown', onDown);
  });
  return (
    <Show when={open}>
      {(m) => {
        // Closed by a choice or a click elsewhere: a focus left in it goes back to More.
        const restore = rememberFocus(() => m.anchor as HTMLElement);
        let list: HTMLElement | undefined;
        onCleanup(() => restore(list));
        return (
          <div class={classes.railMenu} style={{ bottom: `${m.bottom}px` }}>
            <div
              class={canvasClasses.toolMenuCard}
              role="menu"
              aria-label={m.item ? itemTitle(m.item) : 'More'}
              onKeyDown={menuKeys}
              ref={(el) => {
                list = el;
                queueMicrotask(() => focusFirstItem(el));
              }}
            >
              {m.item ? (
                <ToolMenuRows items={m.item.menu!} onRun={() => menu.set(undefined)} />
              ) : (
                <For each={rail.overflow} key={(e) => e.id}>
                  {(entry) => (
                    <button
                      type="button"
                      role="menuitem"
                      class={canvasClasses.toolMenuRow}
                      onClick={(e) => run(entry, m.anchor)}
                    >
                      {entry.tab ? <Icon name={entry.tab.icon} size={16} /> : <ItemIcon item={entry.item!} size={16} />}
                      <span class={canvasClasses.toolMenuLabel}>
                        {entry.tab ? entry.tab.title : itemTitle(entry.item!)}
                      </span>
                      <span />
                    </button>
                  )}
                </For>
              )}
            </div>
          </div>
        );
      }}
    </Show>
  );
}

/** A collapsed top or bottom dock: its tab names in a slim bar; one click expands. */
function CollapsedBar(props: { shell: ShellModel; dock: DockModel }) {
  const { shell, dock } = props;
  return (
    <div class={classes.collapsedBar}>
      <For each={dock.tabs}>
        {(tab) => (
          <button type="button" class={classes.barTab} onClick={() => shell.show(tab.id)} title={`Show ${tab.title}`}>
            <Icon name={dock.id === 'top' ? 'chevron-down' : 'chevron-up'} size={14} />
            {tab.title}
          </button>
        )}
      </For>
    </div>
  );
}

/** A panel's tab: one of a tab list (`rovingTabs`: one Tab stop, the arrows between tabs); it controls the dock's panel. */
function TabButton(props: { shell: ShellModel; dock: DockModel; tab: ShellTab; panel: string }) {
  const { shell, dock, tab } = props;
  // An overlay card in the right panel's place replaces it: no tab there shows selected, and a click brings the tab back.
  const covered = dock.id === 'right' ? shell.overlay.pipe(map((o) => !!o)) : of(false);
  const selected = combineLatest([dock.state, covered]).pipe(
    map(([s, c]) => !c && !s.collapsed && s.active === tab.id),
    distinctUntilChanged(),
  );
  return (
    <button
      type="button"
      role="tab"
      class={classes.tab}
      data-tip={tab.title}
      data-tip-side="left"
      aria-label={tab.title}
      aria-selected={selected}
      aria-controls={props.panel}
      data-tab={tab.id}
      onClick={() => activateTab(shell, dock, tab)}
    >
      <Icon name={tab.icon} size={20} />
    </button>
  );
}

/**
 * Pins the dock beside the canvas (the canvas makes room) or lets it float over it again; in every panel's and overlay
 * card's header, before the ✕, since they all open in the same place.
 */
export function PinButton(props: { shell: ShellModel; dock: DockId }) {
  const { shell } = props;
  const dock = shell.dock(props.dock);
  return (
    <Button
      variant="ghost"
      size="sm"
      class={classes.pin}
      icon={<Icon name="pin" />}
      title="Pin beside the canvas"
      pressed={dock.state.pipe(
        map((s) => s.pinned),
        distinctUntilChanged(),
      )}
      onClick={() => shell.setPinned(props.dock, !dock.state.value.pinned)}
    />
  );
}

/**
 * A tab's content: a card with its header and the view. It is mounted only while its tab is the dock's active one (`Show`
 * in `DockView`: a closed or replaced panel is gone, not hidden). The one way it stays mounted and hidden is an overlay
 * card taking the panel's place (`coveredBy`, the shell's open overlay): the panel is then paused and inert (`activity`,
 * shell.md "Hidden panels"), and comes back as it was.
 */
function DockContent(props: {
  shell: ShellModel;
  dock: DockModel;
  tab: ShellTab;
  /** The panel's id: its tabs control it (`aria-controls`). */
  id: string;
  pinnable?: boolean;
  coveredBy?: ReadonlyState<unknown>;
}) {
  const { shell, dock, tab } = props;
  const shown = props.coveredBy ? props.coveredBy.pipe(map((o) => !o)) : true;
  let section: HTMLElement | undefined;
  // Whether the focus is in it. A focus lost with its nodes (a view's, removed before this content's cleanup) or with
  // its closing (made inert) still counts; one moved out, to another element or to the page (a click on the canvas),
  // does not.
  let focusIn = false;
  const onFocusOut = (e: FocusEvent) => {
    const to = e.relatedTarget as Node | null;
    if (to) {
      focusIn = !!section?.contains(to);
      return;
    }
    const from = eventTarget(e) as Element;
    queueMicrotask(() => {
      const now = deepActiveElement();
      if (now && section?.contains(now)) return;
      if (from.isConnected && !from.closest('[inert]')) focusIn = false;
    });
  };
  // Closed (✕, its shortcut) with the focus in it: the focus goes to its tab, not to the page (README §3.7). A switch
  // with the focus elsewhere (a script's, the agent's, with the focus on the page) leaves it there.
  onCleanup(() => {
    if (!focusIn) return;
    const now = deepActiveElement();
    if (now && !section?.contains(now)) return;
    const root = (section?.closest('.msp-root') ?? section?.getRootNode()) as ParentNode | undefined;
    const tabs = root ? [...root.querySelectorAll<HTMLElement>('[role="tab"]')] : [];
    tabs.find((t) => t.dataset.tab === tab.id)?.focus({ preventScroll: true });
  });
  return activity(shown, () => (
    <section
      class={tab.frame === 'none' ? `${classes.content} ${classes.frameless}` : classes.content}
      role="tabpanel"
      id={props.id}
      aria-label={tab.title}
      ref={(el) => {
        section = el;
      }}
      onFocusIn={() => {
        focusIn = true;
      }}
      onFocusOut={onFocusOut}
    >
      <Splitter shell={shell} dock={dock} />
      <header class={classes.header}>
        <h2 class={classes.title}>{typeof tab.heading === 'function' ? tab.heading() : (tab.heading ?? tab.title)}</h2>
        {props.pinnable ? <PinButton shell={shell} dock={dock.id} /> : null}
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="x" />}
          title="Close panel"
          onClick={() => shell.toggle(dock.id)}
        />
      </header>
      <ScrollArea class={classes.body} contentClass={classes.bodyContent}>
        {/* Shown at a section: the content again, with that section as its focus. */}
        <Show
          when={shell.tabFocus.pipe(
            map((f) => ({ section: f?.tab === tab.id ? f.section : undefined })),
            distinctUntilChanged((a, b) => a.section === b.section && !b.section),
          )}
        >
          {(f) => (
            <ViewHost
              renderer={tab.render}
              ctx={{ shell, dock: dock.id, ...(f.section ? { focus: f.section } : {}) }}
            />
          )}
        </Show>
      </ScrollArea>
    </section>
  ));
}

/**
 * Drags the dock's inner edge to resize it. Lit while the pointer is on it or dragging (`data-active`, not `:hover`:
 * a captured pointer keeps `:hover` on the handle after release until it moves again); the drag ends on release,
 * cancel or a lost capture, and the light goes unless the pointer was released over the handle.
 */
function Splitter(props: { shell: ShellModel; dock: DockModel }) {
  const { shell, dock } = props;
  const sign = dock.id === 'left' ? 1 : -1;
  let startX = 0;
  let startSize = 0;
  let dragging = false;
  const onMove = (e: PointerEvent) => shell.setSize(dock.id, startSize + sign * (e.clientX - startX));
  return (
    <div
      class={classes.splitter}
      role="separator"
      aria-orientation="vertical"
      // From the keyboard too (WAI-ARIA's window splitter): ← and → resize by 16 px (64 with Shift).
      tabIndex={0}
      aria-label="Panel width"
      aria-valuenow={dock.state.pipe(map((s) => Math.round(s.size)))}
      onKeyDown={(e) => {
        if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight') return;
        e.preventDefault();
        const step = (e.shiftKey ? 64 : 16) * (e.key === 'ArrowRight' ? 1 : -1);
        shell.setSize(dock.id, dock.state.value.size + sign * step);
      }}
      onPointerEnter={(e) => e.currentTarget.setAttribute('data-active', '')}
      onPointerLeave={(e) => {
        if (!dragging) e.currentTarget.removeAttribute('data-active');
      }}
      onPointerDown={(e) => {
        const el = e.currentTarget;
        startX = e.clientX;
        startSize = dock.state.value.size;
        dragging = true;
        el.setPointerCapture(e.pointerId);
        el.setAttribute('data-active', '');
        el.addEventListener('pointermove', onMove);
        const end = (u: PointerEvent) => {
          if (!dragging) return;
          dragging = false;
          el.removeEventListener('pointermove', onMove);
          el.removeEventListener('pointerup', end);
          el.removeEventListener('pointercancel', end);
          el.removeEventListener('lostpointercapture', end);
          if (el.hasPointerCapture(u.pointerId)) el.releasePointerCapture(u.pointerId);
          const r = el.getBoundingClientRect();
          const over = u.clientX >= r.left && u.clientX <= r.right && u.clientY >= r.top && u.clientY <= r.bottom;
          if (!over || u.type !== 'pointerup') el.removeAttribute('data-active');
        };
        el.addEventListener('pointerup', end);
        el.addEventListener('pointercancel', end);
        el.addEventListener('lostpointercapture', end);
      }}
      onDblClick={() => shell.setSize(dock.id, dock.id === 'left' ? 300 : 340)}
    />
  );
}
