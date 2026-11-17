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
import { onCleanup } from '../../dom/index.js';
import { StatusLine } from '../../status/status-line.js';
import { StatusModel } from '../../status/status-model.js';
import { HistoryModel, type HistoryOrigin } from '../history-model.js';
import { HistoryView } from '../history-view.js';

const meta: Meta = { title: 'Views/History' };
export default meta;

const person: HistoryOrigin = { kind: 'person' };
const agent: HistoryOrigin = { kind: 'agent', turn: 'turn-1' };
const script: HistoryOrigin = { kind: 'script' };

/** A session's history: the person's steps, an agent's turn grouped, a script's group, a change outside it. */
async function session(history: HistoryModel<Record<string, string>>, scene: { value: Record<string, string> }) {
  const set = (key: string, value: string) => async () => {
    scene.value = { ...scene.value, [key]: value };
  };
  await history.run('Load 1tqn', person, set('structure', '1tqn'));
  await history.run('Color by chain', person, set('color', 'chain'));
  history.invalidate('The scene changed outside the history: a plugin call');
  await history.run('Change the background', person, set('background', 'white'));
  const end = history.beginGroup('Make the ligand stand out', agent);
  await history.run('Select HEM', agent, set('selection', 'HEM'));
  await history.run('Ball and stick on HEM', agent, set('ligand', 'ball-and-stick'));
  await history.run('Focus HEM', agent, set('focus', 'HEM'));
  end();
  await history.group('Apply the preset', script, async (ctx) => {
    await history.run('Cartoon', script, set('polymer', 'cartoon'), { within: ctx });
    await history.run('Surface', script, set('surface', 'molecular'), { within: ctx });
  });
  await history.run('Hide the water', person, set('water', 'hidden'));
  await history.undo();
}

/**
 * A session's steps: the position marked, Undo and Redo with their keys, "Undo to here" and "Redo to here" on the rows
 * (on hover or focus), an agent's turn grouped, a step that cannot be undone and why.
 */
export const Steps: StoryObj = {
  render: () =>
    story(() => {
      const scene = { value: {} as Record<string, string> };
      const history = new HistoryModel<Record<string, string>>({
        checkpoint: () => scene.value,
        restore: async (c) => {
          scene.value = c;
        },
      });
      const status = new StatusModel();
      onCleanup(() => {
        status.retire();
        history.dispose();
      });
      void session(history, scene);
      let n = 0;
      return (
        <div style={page}>
          <div
            ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
            class="msui-surface"
            data-tip-edge="right"
            style="display: flex; flex-direction: column; width: 340px; height: 420px;"
          >
            <HistoryView
              history={history}
              onResult={(r) => {
                if (!r.ok) status.post(r.reason, 'warning');
              }}
            />
          </div>
          <div style="display: flex; gap: var(--msui-space-2);">
            <Button
              onClick={() => void history.run(`Change ${++n}`, person, async () => void (scene.value = { n: `${n}` }))}
            >
              A new step
            </Button>
            <Button variant="ghost" onClick={() => history.clear()}>
              Clear
            </Button>
          </div>
          <StatusLine status={status} />
          <Note>
            The first two steps were made before a change outside the history: their warning says why they cannot be
            undone; undo back to them and Undo says so.
          </Note>
        </div>
      );
    }),
};
