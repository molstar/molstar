/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ToolMenuRows } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { observableState } from '@molstar/prism-ui/core';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { requireScope, runInScope } from '@molstar/prism-ui/dom/scope';
import { on, timeout } from '@molstar/prism-ui/dom/scoped';
import type { CanvasMenuItem } from '@molstar/prism-ui/shell/shell';
import type { SelectionAction, SelectionModel } from '../../model/selection.js';
import { classes } from './_generated/selection-bar.module.css.js';

/**
 * The selection, bottom right of the canvas while there is one: what is selected, a button that opens what to do with
 * it (Focus, Show in structure tree, Make a component, the extensions' rows, More in the Selection panel: rows as the
 * canvas stack's menus, a title and a line under it), and Clear. Mol*'s selection toolbar actions, as one menu.
 */
export function SelectionBar(props: { selection: SelectionModel }) {
  const { selection } = props;
  const open = observableState(false);
  /** Opened from the keyboard: its first row takes the focus (a hover or a click leaves it on the canvas) */
  let byKeyboard = false;
  let bar!: HTMLElement;
  let button!: HTMLElement;
  const rows = (actions: readonly SelectionAction[]): CanvasMenuItem[] =>
    actions.map((a) => ({
      id: a.id,
      label: a.label,
      icon: a.icon,
      ...(a.hint ? { hint: a.hint } : {}),
      // The viewer's own rows, the extensions', then the way to the panel: a hairline before it
      group: a.id === 'more' ? 'more' : 'actions',
      run: () => a.run(),
    }));
  const close = (refocus?: boolean) => {
    keep();
    open.set(false);
    if (refocus) button.focus();
  };
  // Hover opens it, as the canvas stack's menus open (canvas-chrome.tsx `ToolMenus`): leaving both the button and the
  // menu closes it a moment later (crossing the gap between them keeps it); a click and the keyboard open it too.
  const owner = requireScope('SelectionBar');
  let leaving: (() => void) | undefined;
  const keep = () => {
    leaving?.();
    leaving = undefined;
  };
  const leaveSoon = () => {
    keep();
    leaving = runInScope(owner, () => timeout(() => close(), 160));
  };
  const hoverOpen = (e: PointerEvent) => {
    keep();
    // A touch opens it by its click
    if (e.pointerType !== 'touch' && !open.value) open.set(true);
  };
  return (
    <Show when={selection.summary}>
      {(summary) => (
        <div class={classes.bar} role="toolbar" aria-label="Selection" ref={(e: HTMLElement) => (bar = e)}>
          <button
            type="button"
            class={classes.summary}
            aria-haspopup="menu"
            aria-expanded={open}
            ref={(e: HTMLElement) => (button = e)}
            onPointerEnter={hoverOpen}
            onPointerLeave={(e: PointerEvent) => e.pointerType !== 'touch' && open.value && leaveSoon()}
            onClick={(e: MouseEvent) => {
              keep();
              // `detail` 0: Enter or Space, not a press of the pointer
              byKeyboard = e.detail === 0;
              open.set(!open.value);
            }}
          >
            <span>{summary}</span>
            {/* As the canvas stack's buttons with a menu: a corner mark says there is more on hover */}
            <span class={classes.menuMark} aria-hidden="true" />
          </button>
          <button
            type="button"
            class={classes.action}
            aria-label="Clear selection"
            data-tip="Clear selection"
            data-tip-side="top"
            onClick={() => selection.clear()}
          >
            <Icon name="x" size={16} />
          </button>
          <Show when={open}>
            {() => {
              // Escape, or a press outside the bar and the menu, closes it
              on(
                document,
                'keydown',
                (e: KeyboardEvent) => {
                  if (e.key !== 'Escape' || e.defaultPrevented) return;
                  e.preventDefault();
                  close(true);
                },
                { capture: true },
              );
              on(
                document,
                'pointerdown',
                (e: PointerEvent) => {
                  // The path, not the target: inside a shadow root the document sees the host as the target
                  if (!e.composedPath().includes(bar)) close();
                },
                { capture: true },
              );
              return (
                <div
                  class={classes.menu}
                  role="menu"
                  aria-label="What to do with the selection"
                  ref={(el: HTMLElement) => {
                    if (byKeyboard) queueMicrotask(() => el.querySelector<HTMLElement>('button')?.focus());
                  }}
                  onPointerEnter={keep}
                  onPointerLeave={(e: PointerEvent) => e.pointerType !== 'touch' && leaveSoon()}
                >
                  <div class={classes.menuCard}>
                    <Show when={selection.actions}>
                      {(actions) => <ToolMenuRows items={rows(actions)} onRun={() => close()} />}
                    </Show>
                  </div>
                </div>
              );
            }}
          </Show>
        </div>
      )}
    </Show>
  );
}
