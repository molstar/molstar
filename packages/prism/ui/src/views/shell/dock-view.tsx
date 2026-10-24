/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { type ContextMenuItem, openContextMenu } from '../../components/context-menu/context-menu.js';
import { Icon } from '../../components/icon/icon.js';
import { ScrollArea } from '../../components/scroll-area/scroll-area.js';
import { activity, For, Show } from '../../dom/control-flow.js';
import { eventTarget, isInside, ownsKey } from '../../dom/events.js';
import { deepActiveElement, focusFirstItem, menuKeys, rememberFocus, rovingTabs, uniqueId } from '../../dom/focus.js';
import { effect, frameState, onCleanup, watch } from '../../dom/scope.js';
import { setAttr, setClass, setStyle, setText } from '../../dom/writes.js';
import { drag, on } from '../../dom/scoped.js';
import type { Child } from '../../dom/jsx-types.js';
import { type KeyedSelector, keyedSelector } from '../../reactive/keyed-selector.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { sameShallow } from '../../reactive/same-shallow.js';
import { derive } from '../../reactive/state.js';
import type { DockId, DockModel, ShellModel, ShellTab } from '../../shell/shell.js';
import { classes as canvasClasses } from './_generated/canvas.module.css.js';
import { classes } from './_generated/dock.module.css.js';
import { CardSlot, faceOf, ItemIcon, staticTitle, ToolMenuRows } from './canvas-chrome.js';
import { type RailEntry, RailLayout, RailSteps, railStep } from './rail.js';
import { ViewHost } from './view-host.js';

/** What the views of a dock read of its state: one derived value each, shared by them all. */
export interface DockFacets {
  /** Its size, px. */
  readonly size: ReadonlyState<number>;
  /** Its size as the `--msp-dock-size` value. */
  readonly sizePx: ReadonlyState<string>;
  readonly pinned: ReadonlyState<boolean>;
  readonly floating: ReadonlyState<boolean>;
  /**
   * Which tab shows selected (the open one; none while an overlay card takes the right panel's place): a tab's button
   * follows its own key (a switch tells the two tabs it concerns, not every tab).
   */
  readonly selected: KeyedSelector<string | undefined>;
}

const Facets = new WeakMap<DockModel, DockFacets>();

/**
 * The dock's facets (view-models.md §3, the narrowest source): a change of the dock's state (a drag's size, another tab)
 * wakes each facet once, not once a view that reads it; made once a dock.
 */
export function dockFacets(shell: ShellModel, dock: DockModel): DockFacets {
  let f = Facets.get(dock);
  if (f) return f;
  // One gate for what the layout reads (a switch of tabs changes the dock's state, not its layout).
  const layout = derive(
    dock.state,
    (s) => s,
    (a, b) => a.size === b.size && a.pinned === b.pinned && a.floating === b.floating,
  );
  const size = derive(layout, (s) => s.size);
  // An overlay card in the right panel's place replaces it: no tab there shows selected, and a click brings the tab back.
  const open =
    dock.id === 'right'
      ? derive([dock.state, shell.overlay], (s, o) => (o || s.collapsed ? undefined : s.active))
      : derive(dock.state, (s) => (s.collapsed ? undefined : s.active));
  f = {
    size,
    sizePx: derive(size, (px) => `${px}px`),
    pinned: derive(layout, (s) => s.pinned),
    floating: derive(layout, (s) => s.floating),
    selected: keyedSelector(open),
  };
  Facets.set(dock, f);
  return f;
}

/** A side dock: its tab strip on the outer edge and the active tab's content next to the canvas. */
export function DockView(props: {
  shell: ShellModel;
  dock: DockModel /** After the tabs in the strip (compact: More). */;
  stripEnd?: Child;
}) {
  const { shell, dock } = props;
  const horizontal = dock.id === 'top' || dock.id === 'bottom';
  const hasTabs = derive(dock.tabs, (tabs) => tabs.length > 0);
  const compact = derive(shell.mode, (m) => m === 'compact' && dock.id === 'right');
  // A card in the compact sheet, in the panel's place.
  const covered = derive([shell.overlay, compact], (o, c) => (c ? o : undefined));
  // Compact: only the tabs that ask for the bar (the rest are under More); elsewhere all of them.
  const compactTabs = derive([dock.tabs, compact], (tabs, c) =>
    c && tabs.some((t) => t.compactBar) ? tabs.filter((t) => t.compactBar) : tabs,
  );
  // A top or bottom dock with one tab needs no strip: its content collapses itself, and a slim bar expands it.
  const showStrip = horizontal ? derive([dock.tabs, dock.active], (tabs, active) => tabs.length > 1 && !!active) : true;
  const panelId = uniqueId('msp-panel');
  let pinned: HTMLElement | undefined;
  return (
    <Show when={hasTabs}>
      {() => {
        const sizePx = dockFacets(shell, dock).sizePx;
        const root = (
          <div
            class={`${classes.dock} ${classes[dock.id]}`}
            ref={(el) => {
              pinned = el;
            }}
          >
            <Show
              when={dock.active}
              reveal
              // One tab's panel replacing another's: a view transition in the dock (motion.ts `swapIn`).
              swap={() => pinned}
              fallback={horizontal ? () => <CollapsedBar shell={shell} dock={dock} /> : undefined}
            >
              {(tab) => (
                <DockContent
                  shell={shell}
                  dock={dock}
                  tab={tab}
                  id={panelId}
                  // Compact, a card opens in the sheet in the panel's place: the panel stays, paused and inert (out of the
                  // sheet, `data-covered`), and is back as it was when the card closes. Elsewhere a card floats over the
                  // canvas and a docked panel beside it stays live (the floating one is `FloatingDock`'s).
                  {...(dock.id === 'right' ? { coveredBy: covered } : {})}
                />
              )}
            </Show>
            {/* Compact, the place for the card in the sheet (shell.md, Compact): the card moves in while it is offered. */}
            <Show when={compact}>{() => <CardSlot shell={shell} />}</Show>
            <Show when={showStrip}>
              {() => {
                let tablist!: HTMLElement;
                const strip = (
                  <div class={classes.strip}>
                    {/* The tab list holds tabs only; what ends the strip (the compact bar's menu buttons) is beside it. */}
                    <div
                      role="tablist"
                      aria-label="Panels"
                      style="display: contents"
                      ref={(el) => {
                        tablist = el;
                        onCleanup(rovingTabs(el));
                      }}
                    >
                      <For each={compactTabs}>
                        {(tab) => <TabButton shell={shell} dock={dock} tab={tab} panel={panelId} />}
                      </For>
                    </div>
                    {props.stripEnd}
                  </div>
                );
                watch([compact], () =>
                  setAttr(tablist, 'aria-orientation', compact.value || horizontal ? 'horizontal' : 'vertical'),
                );
                return strip;
              }}
            </Show>
          </div>
        ) as HTMLElement;
        // The dock as a function of its size and the compact layout.
        watch([sizePx, compact, covered], () => {
          setClass(root, classes.compact, compact.value);
          setAttr(root, 'data-covered', !!covered.value);
          setStyle(root, '--msp-dock-size', sizePx.value);
        });
        return root;
      }}
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
  const sizePx = dockFacets(shell, dock).sizePx;
  // What does not fit the rail's height, even at its smallest size, is under More (rail.ts).
  const rail = new RailLayout(shell, dock);
  // Tabs at the head of the rail (Home), in the strip, and at its foot (Settings, by the command menu).
  const shown = derive([dock.tabs, rail.hidden], (tabs, hidden) =>
    hidden.size ? tabs.filter((t) => !hidden.has(`tab:${t.id}`)) : tabs,
  );
  // The same tabs in a new array are no change (a layout change that moves no tab renders none of them again).
  const head = derive(shown, (tabs): readonly ShellTab[] => tabs.filter((t) => t.place === 'head'), sameShallow);
  const strip = derive(
    shown,
    (tabs): readonly ShellTab[] => tabs.filter((t) => !t.place || t.place === 'strip'),
    sameShallow,
  );
  const foot = derive(shown, (tabs): readonly ShellTab[] => tabs.filter((t) => t.place === 'foot'), sameShallow);
  const panelId = uniqueId('msp-panel');
  if (dock.id === 'top' || dock.id === 'bottom') {
    // A top or bottom card spans the canvas between the edges and what floats on the sides. Collapsed, it is gone: a
    // canvas toolbar button (Sequence) brings it back.
    // A band while its tab shows and while it fades out (`reveal="fade"` keeps the leaving content): with none active, the
    // last tab's; dropped at once, the card would take a panel's place under the fading band (it moved, then faded).
    let wasBand = false;
    const band = derive(dock.active, (tab) => (wasBand = tab ? tab.frame === 'none' : wasBand));
    const card = (
      <div class={`${classes.floating} ${dock.id === 'top' ? classes.floatingTop : classes.floatingBottom}`}>
        <Show when={dock.active} reveal="fade">
          {(tab) => <DockContent shell={shell} dock={dock} tab={tab} id={panelId} />}
        </Show>
      </div>
    ) as HTMLElement;
    const hidden = derive(dock.active, (tab) => !tab);
    watch([band, sizePx, hidden], () => {
      setClass(card, classes.band, band.value);
      // Hidden, a band fades as a whole (its lines and the canvas-colored backing behind them), where it stands.
      setAttr(card, 'data-hidden', hidden.value);
      setStyle(card, '--msp-dock-size', sizePx.value);
    });
    return card;
  }
  // A rail button's height at the rail's step (rail.ts; the CSS steps the buttons themselves): the panel's header takes
  // it, so the content starts where the rail's first button ends. Before the rail is measured, the normal step.
  const railHead = derive(rail.height.asReadonly, (h) => `${(h > 0 ? railStep(h) : RailSteps[0]).button}px`);
  const root = (
    <div class={`${classes.floating} ${dock.id === 'left' ? classes.floatingLeft : classes.floatingRight}`}>
      <div
        class={`${classes.rail} ${canvasClasses.frame}`}
        ref={(el) => {
          onCleanup(rail.observe(el));
        }}
      >
        {props.railLead}
        <Show when={derive(head, (tabs) => tabs.length > 0)}>
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
        <Show when={derive(strip, (tabs) => tabs.length > 0)}>
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
          <Show when={derive(foot, (tabs) => tabs.length > 0)}>
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
      {/* Opening slides in; one tab's panel replacing another's is a view transition in the dock (motion.ts `swapIn`). */}
      <Show when={dock.active} reveal="side" swap={() => root}>
        {(tab) => <DockContent shell={shell} dock={dock} tab={tab} id={panelId} pinnable coveredBy={shell.overlay} />}
      </Show>
      <RailMoreMenu shell={shell} dock={dock} rail={rail} />
    </div>
  ) as HTMLElement;
  watch([sizePx, railHead], () => {
    setStyle(root, '--msp-dock-size', sizePx.value);
    setStyle(root, '--msp-rail-head', railHead.value);
  });
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
  // Whether its menu is open, a gate: More follows that, not which row opened it.
  const opened = derive(menu.asReadonly, (m) => !!m);
  return (
    <Show when={derive(rail.overflow, (o) => o.length > 0)}>
      {() => {
        // Gone (the rail grew): the menu is closed, it does not come back open.
        onCleanup(() => menu.set(undefined));
        const button = (
          <button
            type="button"
            class={classes.tab}
            data-tip="More"
            data-tip-side="left"
            aria-label="More"
            aria-haspopup="menu"
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
        ) as HTMLElement;
        watch([opened], () => setAttr(button, 'aria-expanded', String(opened.value)));
        return (
          <nav class={classes.strip} aria-label="More buttons">
            {button}
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
  const open = derive([menu.asReadonly, rail.overflow], (m, o) =>
    m && (m.item ? m.item.menu?.length : o.length > 0) ? m : undefined,
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
  on(document, 'keydown', onKey);
  on(document, 'pointerdown', onDown);
  return (
    <Show when={open}>
      {(m) => {
        // Closed by a choice or a click elsewhere: a focus left in it goes back to More.
        const restore = rememberFocus(() => m.anchor as HTMLElement);
        let list: HTMLElement | undefined;
        onCleanup(() => restore(list));
        // The title of the item whose menu it shows (the theme's: its face's).
        const face = m.item ? faceOf(m.item) : undefined;
        const card = (
          <div
            class={canvasClasses.toolMenuCard}
            role="menu"
            aria-label={m.item ? staticTitle(m.item) : 'More'}
            onKeyDown={menuKeys}
            ref={(el) => {
              list = el;
              queueMicrotask(() => focusFirstItem(el));
            }}
          >
            <ScrollArea class={classes.railMenuScroll} contentClass={classes.railMenuRows}>
              {m.item ? (
                <ToolMenuRows items={m.item.menu!} onRun={() => menu.set(undefined)} />
              ) : (
                <For each={rail.overflow} by={(e) => e.id}>
                  {(entry) => <OverflowRow entry={entry} onRun={() => run(entry, m.anchor)} />}
                </For>
              )}
            </ScrollArea>
          </div>
        ) as HTMLElement;
        if (face) watch([face], () => setAttr(card, 'aria-label', face.value.title));
        return (
          <div class={classes.railMenu} style={{ bottom: `${m.bottom}px` }}>
            {card}
          </div>
        );
      }}
    </Show>
  );
}

/** A button left out of the rail, as a row of the More menu: its icon and title (a theme button's follow its face). */
function OverflowRow(props: { entry: RailEntry; onRun(): void }) {
  const { entry } = props;
  const face = entry.item ? faceOf(entry.item) : undefined;
  const label = document.createTextNode(
    entry.tab ? entry.tab.title : face ? face.value.title : staticTitle(entry.item!),
  );
  if (face) watch([face], () => setText(label, face.value.title));
  return (
    <button type="button" role="menuitem" class={canvasClasses.toolMenuRow} onClick={() => props.onRun()}>
      {entry.tab ? (
        <Icon name={entry.tab.icon} size={16} />
      ) : (
        <ItemIcon item={entry.item!} size={16} {...(face ? { face } : {})} />
      )}
      <span class={canvasClasses.toolMenuLabel}>{label}</span>
      <span />
    </button>
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

/**
 * A tab button's menu (its context menu, or Shift+F10 / the menu key): hide this tab, show a hidden one, reset the
 * layout (shell.md §3). The ways back are here and in the app's commands.
 */
export function tabMenuItems(shell: ShellModel, tab: ShellTab): ContextMenuItem[] {
  return [
    { id: 'hide', label: `Hide ${tab.title}`, icon: 'eye-off', run: () => shell.hideTab(tab.id) },
    ...shell.hiddenTabs.value.map(
      (t): ContextMenuItem => ({
        id: `show:${t.id}`,
        label: `Show ${t.title}`,
        icon: 'eye',
        run: () => shell.showTab(t.id),
      }),
    ),
    { id: 'reset', label: 'Reset layout', icon: 'rotate-ccw', run: () => shell.resetLayout() },
  ];
}

/** Opens a tab button's menu at a point (the pointer), or below the button (from the keyboard). */
function openTabMenu(shell: ShellModel, tab: ShellTab, button: HTMLElement, at?: { x: number; y: number }) {
  const root = button.closest<HTMLElement>('.msp-root') ?? document.body;
  const r = button.getBoundingClientRect();
  openContextMenu(root, at ?? { x: r.left, y: r.bottom }, tabMenuItems(shell, tab));
}

/** A panel's tab: one of a tab list (`rovingTabs`: one Tab stop, the arrows between tabs); it controls the dock's panel. */
function TabButton(props: { shell: ShellModel; dock: DockModel; tab: ShellTab; panel: string }) {
  const { shell, dock, tab } = props;
  const selected = dockFacets(shell, dock).selected.key(tab.id);
  const button = (
    <button
      type="button"
      role="tab"
      class={classes.tab}
      data-tip={tab.title}
      data-tip-side="left"
      aria-label={tab.title}
      aria-controls={props.panel}
      data-tab={tab.id}
      onClick={() => activateTab(shell, dock, tab)}
      onContextMenu={(e: MouseEvent) => {
        e.preventDefault();
        openTabMenu(shell, tab, e.currentTarget as HTMLElement, { x: e.clientX, y: e.clientY });
      }}
      onKeyDown={(e: KeyboardEvent) => {
        if (e.key !== 'ContextMenu' && !(e.shiftKey && e.key === 'F10')) return;
        e.preventDefault();
        openTabMenu(shell, tab, e.currentTarget as HTMLElement);
      }}
    >
      <Icon name={tab.icon} size={20} />
    </button>
  ) as HTMLElement;
  // Selected: its own key of the open tab (a switch wakes the two tabs it concerns).
  watch([selected], () => setAttr(button, 'aria-selected', String(selected.value)));
  return button;
}

/**
 * Pins the dock beside the canvas (the canvas makes room) or lets it float over it again; in every panel's and overlay
 * card's header, before the ✕, since they all open in the same place.
 */
export function PinButton(props: { shell: ShellModel; dock: DockId }) {
  const { shell } = props;
  const dock = shell.dock(props.dock);
  // Embedded, nothing is pinned (shell.md §3): no pin to offer.
  // Only on the desktop with the full rail: compact (a bottom sheet), embedded and a compact rail never pin (ShellModel).
  const pinnable = derive([shell.mode, shell.compactRail], (m, c) => !c && m !== 'embedded' && m !== 'compact');
  return (
    <Show when={pinnable}>
      {() => (
        <Button
          variant="ghost"
          size="sm"
          class={classes.pin}
          icon={<Icon name="pin" />}
          title="Pin beside the canvas"
          pressed={dockFacets(shell, dock).pinned}
          onClick={() => shell.setPinned(props.dock, !dock.state.value.pinned)}
        />
      )}
    </Show>
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
  const shown = props.coveredBy ? derive(props.coveredBy, (o) => !o) : true;
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
  // Closed (✕, its shortcut) with the focus in it: the focus goes to its tab, not to the page (accessibility.md). A
  // switch with the focus elsewhere (a script's, the agent's, with the focus on the page) leaves it there.
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
      // Keys scoped to this panel work inside it (commands.md §6).
      data-msp-scope={`panel:${tab.id}`}
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
          when={derive(
            shell.tabFocus,
            (f) => ({ section: f?.tab === tab.id ? f.section : undefined }),
            (a, b) => a.section === b.section && !b.section,
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
  // A drag's size in the frame lane (rendering.md §2.3): pointer events come several a frame, the layout once.
  const dragged = frameState(-1);
  effect(dragged, (size) => {
    if (size >= 0) shell.setSize(dock.id, size);
  });
  let press: { readonly active: boolean } | undefined;
  const size = dockFacets(shell, dock).size;
  const splitter = (
    <div
      class={classes.splitter}
      role="separator"
      aria-orientation="vertical"
      // From the keyboard too (WAI-ARIA's window splitter): ← and → resize by 16 px (64 with Shift).
      tabIndex={0}
      aria-label="Panel width"
      onKeyDown={(e) => {
        if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight') return;
        e.preventDefault();
        const step = (e.shiftKey ? 64 : 16) * (e.key === 'ArrowRight' ? 1 : -1);
        shell.setSize(dock.id, dock.state.value.size + sign * step);
      }}
      onPointerEnter={(e) => e.currentTarget.setAttribute('data-active', '')}
      onPointerLeave={(e) => {
        if (!press?.active) e.currentTarget.removeAttribute('data-active');
      }}
      ref={(el: HTMLElement) => {
        // The drag (dom/scoped.ts `drag`): its listeners attached once, every way it ends ends it.
        press = drag(el, {
          start: (e) => {
            el.setAttribute('data-active', '');
            return { x: e.clientX, size: dock.state.value.size };
          },
          move: (e, s) => dragged.offer(s.size + sign * (e.clientX - s.x)),
          end: (_s, how, u) => {
            if (how === 'up' && u) {
              const r = el.getBoundingClientRect();
              if (u.clientX >= r.left && u.clientX <= r.right && u.clientY >= r.top && u.clientY <= r.bottom) return;
            }
            el.removeAttribute('data-active');
          },
        });
      }}
      onDblClick={() => shell.setSize(dock.id, dock.id === 'left' ? 300 : 340)}
    />
  ) as HTMLElement;
  watch([size], () => setAttr(splitter, 'aria-valuenow', Math.round(size.value)));
  return splitter;
}
