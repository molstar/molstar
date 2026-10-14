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
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import type { DockId, DockModel, ShellModel, ShellTab } from '../../shell/shell.js';
import { classes as canvasClasses } from './_generated/canvas.module.css.js';
import { classes } from './_generated/dock.module.css.js';
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
            {(tab) => <DockContent shell={shell} dock={dock} tab={tab} />}
          </Show>
          <Show when={showStrip}>
            {() => (
              <nav
                class={classes.strip}
                role="tablist"
                aria-orientation={compact.pipe(map((c) => (c || horizontal ? 'horizontal' : 'vertical')))}
              >
                <For each={compactTabs}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} />}</For>
                {props.stripEnd}
              </nav>
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
  railItems?: Child;
  /** At the bottom of the rail (the command menu). */
  railFoot?: Child;
}) {
  const { shell, dock } = props;
  const size = { '--msp-dock-size': dock.state.pipe(map((s) => `${s.size}px`)) };
  // Tabs at the head of the rail (Home), in the strip, and at its foot (Settings, by the command menu).
  const head = dock.tabs.pipe(map((tabs) => tabs.filter((t) => t.place === 'head')));
  const strip = dock.tabs.pipe(map((tabs) => tabs.filter((t) => !t.place || t.place === 'strip')));
  const foot = dock.tabs.pipe(map((tabs) => tabs.filter((t) => t.place === 'foot')));
  if (dock.id === 'top' || dock.id === 'bottom') {
    // A top or bottom card spans the canvas between the edges and what floats on the sides. Collapsed, it is gone: a
    // canvas toolbar button (Sequence) brings it back.
    return (
      <div
        class={`${classes.floating} ${dock.id === 'top' ? classes.floatingTop : classes.floatingBottom}`}
        style={size}
      >
        <Show when={dock.active} reveal>
          {(tab) => <DockContent shell={shell} dock={dock} tab={tab} />}
        </Show>
      </div>
    );
  }
  return (
    <div
      class={`${classes.floating} ${dock.id === 'left' ? classes.floatingLeft : classes.floatingRight}`}
      style={size}
    >
      <div class={`${classes.rail} ${canvasClasses.frame}`}>
        {props.railLead}
        <Show when={head.pipe(map((tabs) => tabs.length > 0))}>
          {() => (
            <nav class={classes.strip} role="tablist" aria-orientation="vertical" aria-label="Home">
              <For each={head}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} />}</For>
            </nav>
          )}
        </Show>
        <Show when={strip.pipe(map((tabs) => tabs.length > 0))}>
          {() => (
            <nav class={classes.strip} role="tablist" aria-orientation="vertical" aria-label="Panels">
              <For each={strip}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} />}</For>
            </nav>
          )}
        </Show>
        {props.railItems}
        <div class={classes.railFoot}>
          <Show when={foot.pipe(map((tabs) => tabs.length > 0))}>
            {() => (
              <nav class={classes.strip} role="tablist" aria-orientation="vertical" aria-label="More panels">
                <For each={foot}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} />}</For>
              </nav>
            )}
          </Show>
          {props.railFoot}
        </div>
      </div>
      <Show when={dock.active} reveal="side">
        {(tab) => <DockContent shell={shell} dock={dock} tab={tab} pinnable />}
      </Show>
    </div>
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

function TabButton(props: { shell: ShellModel; dock: DockModel; tab: ShellTab }) {
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
      onClick={() => {
        if (dock.id === 'right' && shell.overlay.value) {
          shell.closeOverlay();
          shell.show(tab.id);
        } else shell.select(tab.id);
      }}
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

function DockContent(props: { shell: ShellModel; dock: DockModel; tab: ShellTab; pinnable?: boolean }) {
  const { shell, dock, tab } = props;
  return (
    <section
      class={tab.frame === 'none' ? `${classes.content} ${classes.frameless}` : classes.content}
      role="tabpanel"
      aria-label={tab.title}
    >
      <Splitter shell={shell} dock={dock} />
      <header class={classes.header}>
        <span class={classes.title}>
          {typeof tab.heading === 'function' ? tab.heading() : (tab.heading ?? tab.title)}
        </span>
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
  );
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
