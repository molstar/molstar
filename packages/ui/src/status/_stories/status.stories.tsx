/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Note, page, story } from '../../_stories/support.js';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import { type Child, onCleanup, restartable, Show } from '../../dom/index.js';
import { LogView } from '../../log/log-view.js';
import { state } from '../../state.js';
import { StatusChips, StatusLine as Line, StatusPill as Pill } from '../status-line.js';
import { StatusModel } from '../status-model.js';

const meta: Meta = { title: 'Stage/Status' };
export default meta;

/** A stage's place (a placeholder canvas): the status line sits at its bottom left, the chips at its bottom right. */
function Stage(props: { children: Child; corner?: Child }) {
  return (
    <div
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
      style="position: relative; display: flex; align-items: flex-end; justify-content: space-between; gap: var(--msui-space-4); height: 240px; padding: var(--msui-space-3); border-radius: var(--msui-radius-lg); background: radial-gradient(circle at 60% 40%, var(--msui-surface-raised), var(--msui-surface-sunken) 70%); overflow: hidden;"
    >
      {props.children}
      {props.corner ?? null}
    </div>
  );
}

const row = 'display: flex; flex-wrap: wrap; gap: var(--msui-space-2);';

/** The status line: a message (it fades), an error (it stays, opens the log at the problems), and the log beside. */
export const StatusLine: StoryObj = {
  render: () =>
    story(() => {
      const status = new StatusModel();
      onCleanup(() => status.retire());
      status.post('Opened 1tqn from the PDB', 'success', { detail: '640 ms' });
      return (
        <div style={page}>
          <Stage>
            <Line status={status} onOpenLog={() => status.announce('The log shows the problems')} />
          </Stage>
          <div style={row}>
            <Button onClick={() => status.post('Colored the polymer by chain id', 'info', { detail: '12 ms' })}>
              Info
            </Button>
            <Button onClick={() => status.post('Saved the session', 'success')}>Success</Button>
            <Button onClick={() => status.post('No such chain: B', 'warning')}>Warning</Button>
            <Button variant="danger" onClick={() => status.post('Volume streaming failed: 503', 'error')}>
              Error
            </Button>
          </div>
          <div class="msui-surface" style="display: flex; flex-direction: column; width: 340px; height: 220px;">
            <LogView log={status.log} />
          </div>
          <Note>A message fades after 4 s; an error stays until its ✕. Each is said once by the polite region.</Note>
        </div>
      );
    }),
};

/** A line of its own: dismissed after its time, a ring around its ✕ counting it down, held while hovered. */
export const StatusPill: StoryObj = {
  render: () =>
    story(() => {
      const shown = state(true);
      const error = state(true);
      return (
        <div style={page}>
          <Stage>
            <div style="display: flex; flex-direction: column; align-items: flex-start;">
              <Show when={shown}>
                {() => (
                  <Pill
                    icon={<Icon name="circle-check" />}
                    text="The assistant finished its turn"
                    kind="success"
                    tip="Show the turn"
                    onClick={() => {}}
                    onDismiss={() => shown.set(false)}
                    autoDismiss={8000}
                  />
                )}
              </Show>
              <Show when={error}>
                {() => (
                  <Pill
                    icon={<Icon name="circle-alert" />}
                    text="The structure could not load"
                    kind="error"
                    onDismiss={() => error.set(false)}
                  />
                )}
              </Show>
            </div>
          </Stage>
          <div style={row}>
            <Button
              onClick={() => {
                shown.set(true);
                error.set(true);
              }}
            >
              Show again
            </Button>
          </div>
          <Note>The first goes after 8 s (hover it: the ring and the time hold). The error stays until its ✕.</Note>
        </div>
      );
    }),
};

/** Running tasks: a pill each, made once and written on progress; Cancel while it can be. */
export const TaskPills: StoryObj = {
  render: () =>
    story(() => {
      const status = new StatusModel();
      onCleanup(() => status.retire());
      let progress = 0;
      const tick = restartable(() => {
        progress = progress >= 1 ? 0 : progress + 0.05;
        status.setTasks([
          {
            id: 'download',
            label: 'Downloading 1tqn',
            progress,
            abort: () => status.post('Download cancelled', 'warning'),
          },
          { id: 'compute', label: 'Computing the surface' },
        ]);
        tick.start(400);
      });
      tick.start(0);
      status.putTask({ id: 'extension', label: 'Loading the video export' });
      return (
        <div style={page}>
          <Stage>
            <Line status={status} />
          </Stage>
          <div style={row}>
            <Button onClick={() => status.removeTask('extension')}>Finish the extension’s task</Button>
          </div>
          <Note>
            The download’s pill keeps its spinner and Cancel while its progress moves (one pill per id). A task put by
            one writer stays when the mirrored ones change.
          </Note>
        </div>
      );
    }),
};

/** What the person works on: chips at the stage's bottom right; a click activates one, its ✕ clears it. */
export const Chips: StoryObj = {
  render: () =>
    story(() => {
      const status = new StatusModel();
      onCleanup(() => status.retire());
      const reset = () =>
        status.setChips([
          {
            id: 'selection',
            label: '12 residues selected',
            activate: () => status.announce('Focused the selection'),
            dismiss: () => status.setChips(status.chips.value.filter((c) => c.id !== 'selection')),
          },
          { id: 'focus', label: 'Focus: ligand HEM', activate: () => status.announce('Focused HEM') },
        ]);
      reset();
      return (
        <div style={page}>
          <Stage corner={<StatusChips status={status} />}>
            <Line status={status} />
          </Stage>
          <div style={row}>
            <Button onClick={reset}>Reset the chips</Button>
          </div>
        </div>
      );
    }),
};
