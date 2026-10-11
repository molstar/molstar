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
import { classes } from './_generated/canvas.module.css.js';
import { ViewHost } from './view-host.js';

/** The vertical button strip at the top right of the canvas. */
export function CanvasToolbar(props: { shell: ShellModel }) {
  const { shell } = props;
  return (
    <div class={classes.toolbar} role="toolbar" aria-orientation="vertical" aria-label="Canvas">
      <For each={shell.toolbar}>{(item) => <ToolButton shell={shell} item={item} />}</For>
    </div>
  );
}

function ToolButton(props: { shell: ShellModel; item: CanvasToolbarItem }) {
  const { shell, item } = props;
  const overlayOpen = item.overlay ? shell.overlay.pipe(map((o) => o?.overlay.id === item.overlay)) : of(false);
  const pressed = combineLatest([overlayOpen, item.pressed ?? of(false)]).pipe(
    map(([a, b]) => a || b),
    distinctUntilChanged(),
  );
  return (
    <button
      type="button"
      class={classes.tool}
      title={item.title}
      aria-label={item.title}
      aria-pressed={pressed}
      onClick={(e) => {
        if (item.overlay) shell.toggleOverlay(item.overlay, e.currentTarget);
        else item.run?.(e.currentTarget);
      }}
    >
      <Icon name={item.icon} size={20} />
    </button>
  );
}

/** The open overlay as a card beside the toolbar, level with the button that opened it. */
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
        <ViewHost renderer={open.overlay.render} ctx={{ shell }} />
      </div>
    </section>
  ) as HTMLElement;
  // Level with its button once placed (the card is positioned within the canvas layer).
  queueMicrotask(() => {
    const layer = card.offsetParent;
    if (!open.anchor || !layer) return;
    const top = open.anchor.getBoundingClientRect().top - layer.getBoundingClientRect().top - 4;
    card.style.top = `${Math.max(8, Math.min(top, layer.clientHeight - card.offsetHeight - 8))}px`;
  });
  card.style.top = '12px';
  // Esc, or a pointer down outside the card and its button, closes it.
  const doc = document;
  const onKey = (e: KeyboardEvent) => {
    if (e.key === 'Escape' && !e.defaultPrevented) shell.closeOverlay();
  };
  const onPointer = (e: PointerEvent) => {
    const target = e.target as Node;
    if (!card.contains(target) && !open.anchor?.contains(target)) shell.closeOverlay();
  };
  doc.addEventListener('keydown', onKey);
  doc.addEventListener('pointerdown', onPointer, true);
  onCleanup(() => {
    doc.removeEventListener('keydown', onKey);
    doc.removeEventListener('pointerdown', onPointer, true);
  });
  return card;
}

const KindIcons: Record<MessageKind, IconName> = {
  info: 'info',
  success: 'circle-check',
  warning: 'triangle-alert',
  error: 'circle-alert',
};

/** The line at the bottom of the canvas: the latest message, running tasks, selection and focus chips. */
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
      <For each={status.chips}>
        {(chip) => (
          <span class={classes.pill}>
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
