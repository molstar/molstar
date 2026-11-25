/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Note, page, story } from '../../_stories/support.js';
import { Button } from '../../components/button/button.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import { type Child, onCleanup, restartable } from '../../dom/index.js';
import { type LogKind, LogModel, LogProblems } from '../log-model.js';
import { LogView } from '../log-view.js';

const meta: Meta = { title: 'Views/Log' };
export default meta;

/** A log in a panel's place: a surface, a flex column it fills, the tooltips on. */
function Panel(props: { children: Child; height?: number }) {
  return (
    <div
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
      class="msui-surface"
      data-tip-edge="right"
      style={`display: flex; flex-direction: column; width: 340px; height: ${props.height ?? 360}px;`}
    >
      {props.children}
    </div>
  );
}

const start = Date.UTC(2026, 10, 17, 9, 30, 0);

/** A session's log: messages, a load's transforms with their times, a plugin's warnings, debug lines. */
function session(log: LogModel) {
  let t = start;
  const at = () => (t += 1300);
  log.write({ text: 'Opened 1tqn from the PDB', kind: 'success', time: at(), detail: '640 ms' });
  for (const [i, name] of ['Download', 'Parse CIF', 'Trajectory', 'Model', 'Structure', 'Polymer', 'Ligand'].entries())
    log.write({ text: `Created ${name}`, source: 'Transforms', detail: `${3 + i * 7} ms`, time: at() });
  log.write({ text: 'No such chain: B', kind: 'warning', time: at() });
  log.write({ text: 'Colored the polymer by chain id', time: at(), detail: '12 ms' });
  log.write({ text: 'Lost the WebGL context; restoring it', kind: 'warning', source: 'Plugin', time: at() });
  log.write({
    text: 'Volume streaming failed: the server answered 503 Service Unavailable',
    kind: 'error',
    time: at(),
  });
  log.write({ text: 'Frame 3 of 12 rendered in 18 ms', kind: 'debug', source: 'Renderer', time: at() });
  log.write({
    text: 'A long line that does not fit the panel fades at its end, and its row’s side tooltip shows the whole of it',
    time: at(),
  });
  log.flush();
  return log;
}

/** The filters (All, Messages, one per source, Problems; each named in its tooltip) and a text filter (Escape). */
export const Log: StoryObj = {
  render: () =>
    story(() => (
      <div style={page}>
        <Panel>
          <LogView log={session(new LogModel())} sourceIcons={{ Transforms: 'layers', Plugin: 'wrench' }} />
        </Panel>
        <Note>Hover a cut-off row: its side tooltip shows the whole text. Type in the filter; Escape empties it.</Note>
      </div>
    )),
};

/** Thousands of entries: a recycled list of rows, a window of them in the DOM. */
export const Long: StoryObj = {
  render: () =>
    story(() => {
      const log = new LogModel({ size: 5000 });
      const kinds: LogKind[] = ['info', 'info', 'info', 'success', 'debug', 'warning', 'info', 'error'];
      for (let i = 0; i < 5000; i++) {
        log.write({
          text: `Updated cell ${i + 1} of the representation tree`,
          kind: kinds[i % kinds.length],
          source: i % 3 ? 'Transforms' : 'Plugin',
          detail: `${i % 17} ms`,
          time: start + i * 250,
        });
      }
      log.flush();
      return (
        <div style={page}>
          <Panel height={480}>
            <LogView log={log} sourceIcons={{ Transforms: 'layers', Plugin: 'wrench' }} />
          </Panel>
          <Note>5,000 entries: scroll fast, only a screenful of rows is in the DOM.</Note>
        </div>
      );
    }),
};

/** Entries arriving: the list stays at its end, unless scrolled up (then it stays where it is). */
export const Live: StoryObj = {
  render: () =>
    story(() => {
      const log = session(new LogModel());
      let n = 0;
      const running = { on: true };
      const next = restartable(() => {
        n++;
        log.write({
          text: n % 6 ? `Received frame ${n}` : `Frame ${n} was late`,
          kind: n % 6 ? 'info' : 'warning',
          source: 'Stream',
          detail: `${(n * 37) % 90} ms`,
        });
        if (running.on) next.start(700);
      });
      next.start(700);
      return (
        <div style={page}>
          <Panel>
            <LogView log={log} />
          </Panel>
          <div style="display: flex; gap: var(--msui-space-2);">
            <Button
              onClick={() => {
                running.on = !running.on;
                if (running.on) next.start(700);
                else next.cancel();
              }}
            >
              Pause or resume
            </Button>
            <Button variant="ghost" onClick={() => log.clear()}>
              Clear
            </Button>
          </div>
          <Note>An entry every 0.7 s. Scroll up: the list stays; back at the end, it follows again.</Note>
        </div>
      );
    }),
};

/** The problems: warnings and errors only (an error on the status line opens the log here). */
export const Problems: StoryObj = {
  render: () =>
    story(() => {
      const log = session(new LogModel());
      log.setFilter(LogProblems);
      return (
        <div style={page}>
          <Panel>
            <LogView log={log} sourceIcons={{ Transforms: 'layers', Plugin: 'wrench' }} />
          </Panel>
        </div>
      );
    }),
};
