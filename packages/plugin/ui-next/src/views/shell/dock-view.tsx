/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map, of } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { For, Show } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import type { DockModel, ShellModel, ShellTab } from '../../shell/shell.js';
import { classes } from './_generated/dock.module.css.js';
import { ViewHost } from './view-host.js';

/** A side dock: its tab strip on the outer edge and the active tab's content next to the canvas. */
export function DockView(props: { shell: ShellModel; dock: DockModel }) {
  const { shell, dock } = props;
  const horizontal = dock.id === 'top' || dock.id === 'bottom';
  const hasTabs = dock.tabs.pipe(map((tabs) => tabs.length > 0));
  const compact = shell.mode.pipe(map((m) => m === 'compact' && dock.id === 'right'));
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
          <Show when={dock.active} fallback={horizontal ? () => <CollapsedBar shell={shell} dock={dock} /> : undefined}>
            {(tab) => <DockContent shell={shell} dock={dock} tab={tab} />}
          </Show>
          <Show when={showStrip}>
            {() => (
              <nav
                class={classes.strip}
                role="tablist"
                aria-orientation={compact.pipe(map((c) => (c || horizontal ? 'horizontal' : 'vertical')))}
              >
                <For each={dock.tabs}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} />}</For>
              </nav>
            )}
          </Show>
        </div>
      )}
    </Show>
  );
}

/**
 * A dock floating over the canvas (the right dock by default): a rail at the edge with the tab strip between any
 * leading items (the palette) and further items (the canvas toolbar's sections), and the active tab's content as a
 * glass card beside it. Collapsed, only the rail stays.
 */
export function FloatingDock(props: { shell: ShellModel; dock: DockModel; railLead?: Child; railItems?: Child }) {
  const { shell, dock } = props;
  const size = { '--msp-dock-size': dock.state.pipe(map((s) => `${s.size}px`)) };
  if (dock.id === 'top' || dock.id === 'bottom') {
    // A top or bottom card spans the canvas between the edges and what floats on the sides. Collapsed, it is gone: a
    // canvas toolbar button (Sequence) brings it back.
    return (
      <div
        class={`${classes.floating} ${dock.id === 'top' ? classes.floatingTop : classes.floatingBottom}`}
        style={size}
      >
        <Show when={dock.active}>{(tab) => <DockContent shell={shell} dock={dock} tab={tab} />}</Show>
      </div>
    );
  }
  return (
    <div
      class={`${classes.floating} ${dock.id === 'left' ? classes.floatingLeft : classes.floatingRight}`}
      style={size}
    >
      <div class={classes.rail}>
        {props.railLead}
        <Show when={dock.tabs.pipe(map((tabs) => tabs.length > 0))}>
          {() => (
            <nav class={classes.strip} role="tablist" aria-orientation="vertical" aria-label="Panels">
              <For each={dock.tabs}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} />}</For>
            </nav>
          )}
        </Show>
        {props.railItems}
      </div>
      <Show when={dock.active}>{(tab) => <DockContent shell={shell} dock={dock} tab={tab} />}</Show>
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
  const selected = dock.state.pipe(
    map((s) => !s.collapsed && s.active === tab.id),
    distinctUntilChanged(),
  );
  return (
    <button
      type="button"
      role="tab"
      class={classes.tab}
      title={tab.title}
      aria-label={tab.title}
      aria-selected={selected}
      onClick={() => shell.select(tab.id)}
    >
      <Icon name={tab.icon} size={20} />
    </button>
  );
}

function DockContent(props: { shell: ShellModel; dock: DockModel; tab: ShellTab }) {
  const { shell, dock, tab } = props;
  return (
    <section class={classes.content} role="tabpanel" aria-label={tab.title}>
      <Splitter shell={shell} dock={dock} />
      <header class={classes.header}>
        <span class={classes.title}>{tab.title}</span>
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="x" />}
          title="Close panel"
          onClick={() => shell.toggle(dock.id)}
        />
      </header>
      <div class={classes.body}>
        <ViewHost renderer={tab.render} ctx={{ shell, dock: dock.id }} />
      </div>
    </section>
  );
}

/** Drags the dock's inner edge to resize it. */
function Splitter(props: { shell: ShellModel; dock: DockModel }) {
  const { shell, dock } = props;
  const sign = dock.id === 'left' ? 1 : -1;
  let startX = 0;
  let startSize = 0;
  const onMove = (e: PointerEvent) => shell.setSize(dock.id, startSize + sign * (e.clientX - startX));
  return (
    <div
      class={classes.splitter}
      role="separator"
      aria-orientation="vertical"
      onPointerDown={(e) => {
        const el = e.currentTarget;
        startX = e.clientX;
        startSize = dock.state.value.size;
        el.setPointerCapture(e.pointerId);
        el.setAttribute('data-dragging', '');
        el.addEventListener('pointermove', onMove);
        el.addEventListener(
          'pointerup',
          () => {
            el.removeEventListener('pointermove', onMove);
            el.removeAttribute('data-dragging');
          },
          { once: true },
        );
      }}
      onDblClick={() => shell.setSize(dock.id, dock.id === 'left' ? 300 : 340)}
    />
  );
}
