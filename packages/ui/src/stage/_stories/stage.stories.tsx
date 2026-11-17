/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/support.js';
import { Button } from '../../components/button/button.js';
import { PlaybackDeck, StepperRow, TimelineRow } from '../../components/playback/playback.js';
import { Stepper } from '../../components/playback/stepper-model.js';
import { Durationless, Timeline } from '../../components/playback/timeline-model.js';
import { onCleanup } from '../../dom/index.js';
import { StatusChips, StatusLine } from '../../status/status-line.js';
import { StatusModel } from '../../status/status-model.js';
import { Stage } from '../stage.js';
import { StageModel, type StageModelOptions } from '../stage-model.js';
import { StageToolbar } from '../toolbar.js';
import { DemoBand, DemoCanvas, demoGroups, Pill, StageBox } from './demo.js';

const meta: Meta = { title: 'Stage/Stage' };
export default meta;

/** A model for the story's lifetime. */
function model(options: StageModelOptions): StageModel {
  const stage = new StageModel(options);
  onCleanup(() => stage.dispose());
  return stage;
}

/** Loaded: the canvas stand-in, the toolbar in the corner; hover a bead: its label shows in the corner above it. */
export const Loaded: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({ empty: false });
      const status = new StatusModel();
      onCleanup(() => status.retire());
      status.post('1CBS loaded', 'success');
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            status={<StatusLine status={status} />}
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};

/** Empty: after 700 ms, the glass wordmark (press it); the canvas is quiet. Load fills it, Clear empties it again. */
export const Empty: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({});
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            lead={
              <div style="display: flex; gap: var(--msui-space-2)">
                <Button size="sm" onClick={() => stage.setEmpty(false)}>
                  Load
                </Button>
                <Button size="sm" onClick={() => stage.setEmpty(true)}>
                  Clear
                </Button>
              </div>
            }
            status={<Pill>{derive(stage.mark, (m) => `mark: ${m}`)}</Pill>}
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};

/** The hover label beside the pointer: one read and one write a frame, flipped near the edges. */
export const HoverAtCursor: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({ empty: false, hoverPlacement: 'cursor' });
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};

/**
 * A band on top (a sequence's shape): the canvas's color behind it, fading out below it (the faded part lets the
 * pointer through to the canvas); the chrome starts under it (`topInset`, shown top left). Hover a residue.
 */
export const BandOnTop: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({ empty: false });
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            top={<DemoBand stage={stage} />}
            lead={
              <div style="display: flex; gap: var(--msui-space-2); align-items: center">
                <Pill>{derive(stage.topInset, (px) => `topInset: ${px}px`)}</Pill>
                <Button size="sm" onClick={() => stage.showBand(!stage.band.value)}>
                  Toggle the band
                </Button>
              </div>
            }
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};

/** Busy: work under way on an empty stage holds the mark back; it comes 700 ms after the work ends. */
export const Busy: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({});
      stage.setBusy(true);
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            lead={
              <Button size="sm" onClick={() => stage.setBusy(!stage.busy.value)}>
                {derive(stage.busy, (b) => (b ? 'End the work' : 'Start work'))}
              </Button>
            }
            status={<Pill>{derive(stage.busy, (b) => (b ? 'Loading…' : 'Idle'))}</Pill>}
          />
        </StageBox>
      );
    }),
};

const Backgrounds: readonly (readonly [string, string])[] = [
  ['White', '#ffffff'],
  ['Light grey', '#e4e7eb'],
  ['Mid grey', '#767676'],
  ['Teal', '#2f6f8f'],
  ['Night', '#0f1115'],
];

/** The ink: pick a background; the chrome over it (toolbar, label, band) takes the theme that reads on it. */
export const Ink: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({ empty: false, background: '#0f1115' });
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            top={<DemoBand stage={stage} rows={1} />}
            lead={
              <div style="display: flex; flex-wrap: wrap; gap: var(--msui-space-1)">
                {Backgrounds.map(([label, color]) => (
                  <Button
                    size="sm"
                    pressed={derive(stage.background, (b) => b === color)}
                    onClick={() => stage.setBackground(color)}
                  >
                    {label}
                  </Button>
                ))}
              </div>
            }
            status={<Pill>{derive(stage.ink, (i) => `${i.color}: ${i.theme} chrome`)}</Pill>}
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};

/** Compact (narrow): the toolbar folds into one button, the status takes a row of its own above. */
export const Compact: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({ empty: false });
      return (
        <StageBox width="360px" height="560px">
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            status={<Pill>1CBS loaded, 1 chain</Pill>}
            chips={<Pill>Residue 42</Pill>}
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};

/** The status and chips: a message that fades, an error that stays, a task with progress and Cancel, a selection chip. */
export const StatusAndChips: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({ empty: false });
      const status = new StatusModel();
      onCleanup(() => status.retire());
      let progress = 0;
      const timer = setInterval(() => {
        progress = (progress + 0.05) % 1.05;
        status.putTask({ id: 'load', label: 'Loading 1CBS', progress, abort: () => status.removeTask('load') });
      }, 300);
      onCleanup(() => clearInterval(timer));
      const select = () =>
        status.setChips([
          {
            id: 'sel',
            label: '12 residues',
            activate: () => status.post('Focused 12 residues'),
            dismiss: () => status.setChips([]),
          },
        ]);
      select();
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            lead={
              <div style="display: flex; flex-wrap: wrap; gap: var(--msui-space-1)">
                <Button size="sm" onClick={() => status.post('Saved the session', 'success')}>
                  Message
                </Button>
                <Button
                  size="sm"
                  onClick={() => status.post('Could not read 1xyz.cif', 'error', { detail: 'line 12' })}
                >
                  Error
                </Button>
                <Button size="sm" onClick={select}>
                  Select
                </Button>
              </div>
            }
            status={<StatusLine status={status} />}
            chips={<StatusChips status={status} />}
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};

/** Playback: a deck top-left, under a band (a tour's steps, a timed animation, a simulation without an end). */
export const Playback: StoryObj = {
  render: () =>
    story(() => {
      const stage = model({ empty: false });
      const steps = ['Overview', 'The pocket', 'The ligand', 'Summary'].map((title) => ({ title }));
      const stepper = new Stepper({ steps, intervalMs: 2000 });
      const timed = new Timeline({ duration: 4 });
      const free = new Durationless({});
      onCleanup(() => {
        stepper.dispose();
        timed.dispose();
        free.dispose();
      });
      return (
        <StageBox>
          <Stage
            model={stage}
            canvas={<DemoCanvas stage={stage} />}
            top={<DemoBand stage={stage} rows={1} />}
            lead={
              <PlaybackDeck>
                <StepperRow model={stepper} label="Tour" />
                <TimelineRow model={timed} />
                <TimelineRow model={free} />
              </PlaybackDeck>
            }
            corner={<StageToolbar groups={demoGroups()} stage={stage} />}
          />
        </StageBox>
      );
    }),
};
