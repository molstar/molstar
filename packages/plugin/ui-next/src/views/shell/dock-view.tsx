/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { For, Show } from '../../dom/control-flow.js';
import type { DockModel, ShellModel, ShellTab } from '../../shell/shell.js';
import { classes } from './_generated/dock.module.css.js';
import { ViewHost } from './view-host.js';

/** A side dock: its tab strip on the outer edge and the active tab's content next to the canvas. */
export function DockView(props: { shell: ShellModel; dock: DockModel }) {
  const { shell, dock } = props;
  const hasTabs = dock.tabs.pipe(map((tabs) => tabs.length > 0));
  const compact = shell.mode.pipe(map((m) => m === 'compact' && dock.id === 'right'));
  return (
    <Show when={hasTabs}>
      {() => (
        <div
          class={`${classes.dock} ${classes[dock.id]}`}
          classList={{ [classes.compact]: compact }}
          style={{ '--msp-dock-size': dock.state.pipe(map((s) => `${s.size}px`)) }}
        >
          <Show when={dock.active}>{(tab) => <DockContent shell={shell} dock={dock} tab={tab} />}</Show>
          <nav
            class={classes.strip}
            role="tablist"
            aria-orientation={compact.pipe(map((c) => (c ? 'horizontal' : 'vertical')))}
          >
            <For each={dock.tabs}>{(tab) => <TabButton shell={shell} dock={dock} tab={tab} />}</For>
          </nav>
        </div>
      )}
    </Show>
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
          icon={<Icon name={dock.id === 'left' ? 'panel-left' : 'panel-right'} />}
          title="Collapse panel"
          onClick={() => shell.toggle(dock.id)}
        />
      </header>
      <div class={classes.body}>
        <ViewHost renderer={tab.render} ctx={{ shell }} />
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
