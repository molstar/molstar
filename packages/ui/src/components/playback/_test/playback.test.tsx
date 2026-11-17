/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it } from 'vitest';
import { flush, flushFrame, render } from '../../../dom/index.js';
import { axeViolations } from '../../../testing/axe.js';
import { MuteButton, PlaybackDeck, StepCard, StepperRow, TimelineRow } from '../playback.js';
import { Stepper } from '../stepper-model.js';
import { Durationless, Timeline, type TimelineModel } from '../timeline-model.js';

// The playback rows over their models (playback §3): which rows show, what they say, what a press or a scrub does to the
// model; the deck and the card; axe over them.

const tick = () => new Promise((r) => setTimeout(r, 0));
const steps = (n: number) =>
  Array.from({ length: n }, (_, i) => ({ title: `Step ${i + 1}`, description: `# Step ${i + 1}\n\nText ${i + 1}` }));

let disposers: (() => void)[] = [];
afterEach(() => {
  for (const d of disposers) d();
  disposers = [];
  document.body.innerHTML = '';
});

function mount(content: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  disposers.push(render(content, container));
  const rows = () => [...container.querySelectorAll('[data-row]')];
  const button = (name: string) =>
    container.querySelector<HTMLButtonElement>(`button[data-tip="${name}"], button[aria-label="${name}"]`);
  return { container, rows, button };
}

describe('the playback deck', () => {
  it('has a stepper row from two steps on, a timed row while it has a length, a durationless row always', () => {
    const stepper = new Stepper({ steps: steps(1) });
    const timed = new Timeline({ duration: 0 });
    const m = mount(() => (
      <PlaybackDeck>
        <StepperRow model={stepper} />
        <TimelineRow model={timed} />
      </PlaybackDeck>
    ));
    expect(m.rows()).toHaveLength(0);
    stepper.setSteps(steps(3));
    timed.setLength(2);
    flush();
    expect(m.rows()).toHaveLength(2);
    expect(m.container.querySelector('section')?.getAttribute('aria-label')).toBe('Playback');
    const d = mount(() => <TimelineRow model={new Durationless()} />);
    expect(d.rows()).toHaveLength(1);
  });
});

describe('the stepper row', () => {
  it('shows dots up to seven: each goes to its step, the current one marked', () => {
    const s = new Stepper({ steps: steps(3) });
    const m = mount(() => <StepperRow model={s} />);
    const dots = [...m.container.querySelectorAll<HTMLButtonElement>('[role=group] [role=group] button')];
    expect(dots.map((d) => d.getAttribute('aria-label'))).toEqual(['Step 1', 'Step 2', 'Step 3']);
    expect(dots[0].getAttribute('aria-current')).toBe('step');
    dots[2].click();
    expect(s.index.value).toBe(2);
    expect(dots[2].getAttribute('aria-current')).toBe('step');
    expect(dots[0].hasAttribute('aria-current')).toBe(false);
    m.button('Previous')!.click();
    expect(s.index.value).toBe(1);
    m.button('Next')!.click();
    m.button('Next')!.click();
    expect(s.index.value).toBe(0);
    m.button('Play')!.click();
    expect(s.playing.value).toBe(true);
    expect(m.button('Pause')).not.toBeNull();
    s.dispose();
  });

  it('past seven, "i / N" opens a filterable list; a pick goes there and closes it, Esc closes it', async () => {
    const s = new Stepper({ steps: steps(20) });
    const m = mount(() => <StepperRow model={s} />);
    const position = m.container.querySelector<HTMLButtonElement>('button[aria-expanded]')!;
    expect(position.textContent).toBe('1 / 20');
    expect(position.getAttribute('aria-label')).toBe('Steps: Step 1');
    position.click();
    await tick();
    const filter = m.container.querySelector<HTMLInputElement>('input[aria-label="Filter steps"]')!;
    expect(document.activeElement).toBe(filter);
    filter.value = 'step 1';
    filter.dispatchEvent(new Event('input', { bubbles: true }));
    flush();
    const items = () => [...m.container.querySelectorAll<HTMLButtonElement>('ul button')];
    // Step 1 and Step 10…19.
    expect(items()).toHaveLength(11);
    items()[3].click();
    expect(s.index.value).toBe(11);
    expect(position.textContent).toBe('12 / 20');
    expect(position.getAttribute('aria-expanded')).toBe('false');
    expect(document.activeElement).toBe(position);
    position.click();
    await tick();
    const e = new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, cancelable: true });
    m.container.querySelector('ul button')!.dispatchEvent(e);
    flush();
    expect(e.defaultPrevented).toBe(true);
    expect(m.container.querySelector('ul')).toBeNull();
  });

  it('ends with its actions: a mute', () => {
    const muted = new State(false);
    const m = mount(() => (
      <StepperRow
        model={new Stepper({ steps: steps(2) })}
        actions={() => <MuteButton muted={muted} onClick={() => muted.set(!muted.value)} />}
      />
    ));
    const mute = m.button('Mute')!;
    mute.click();
    flush();
    expect(mute.getAttribute('aria-pressed')).toBe('true');
  });
});

describe('the step card', () => {
  it('shows the title and the description without a heading that repeats it', () => {
    const s = new Stepper({ steps: steps(2) });
    const m = mount(() => <StepCard model={s} />);
    expect(m.container.querySelector('section')?.getAttribute('aria-label')).toBe('Steps');
    expect(m.container.textContent).toContain('Step 1');
    expect(m.container.querySelector('h1')).toBeNull();
    expect(m.container.textContent).toContain('Text 1');
    s.next();
    flush();
    expect(m.container.textContent).toContain('Text 2');
  });

  it('is not there without steps; a single step has no controls', () => {
    const s = new Stepper();
    const m = mount(() => <StepCard model={s} />);
    expect(m.container.querySelector('section')).toBeNull();
    s.setSteps(steps(1));
    flush();
    expect(m.container.querySelector('section')).not.toBeNull();
    expect(m.button('Next')).toBeNull();
  });
});

describe('the timeline row', () => {
  it('timed: scrubs once a frame, says the time, plays', () => {
    const t = new Timeline({ duration: 4, fps: 10 });
    const m = mount(() => <TimelineRow model={t} />);
    const scrubber = m.container.querySelector<HTMLInputElement>('input[type=range]')!;
    expect(scrubber.getAttribute('aria-label')).toBe('Time');
    expect(scrubber.max).toBe('4');
    scrubber.value = '1.2';
    scrubber.dispatchEvent(new Event('input', { bubbles: true }));
    scrubber.value = '1.5';
    scrubber.dispatchEvent(new Event('input', { bubbles: true }));
    // Offered, not sought, until the frame.
    expect(t.time!.value).toBe(0);
    flushFrame();
    flush();
    expect(t.time!.value).toBe(1.5);
    expect(scrubber.getAttribute('aria-valuetext')).toBe('1.5 s of 4.0 s');
    expect(m.container.textContent).toContain('1.5 / 4.0 s');
    // Back to a value it sought before, after the model moved: sought again.
    t.seek(3);
    scrubber.value = '1.5';
    scrubber.dispatchEvent(new Event('input', { bubbles: true }));
    flushFrame();
    expect(t.time!.value).toBe(1.5);
    m.button('Play')!.click();
    expect(t.playing.value).toBe(true);
    t.dispose();
  });

  it('counted: frames as people count, previous and next wrap when periodic', () => {
    const t = new Timeline({ count: 10, unit: 'model' });
    const m = mount(() => <TimelineRow model={t} />);
    const scrubber = m.container.querySelector<HTMLInputElement>('input[type=range]')!;
    expect(scrubber.getAttribute('aria-label')).toBe('Model');
    expect(scrubber.max).toBe('9');
    expect(scrubber.getAttribute('aria-valuetext')).toBe('model 1 of 10');
    m.button('Previous model')!.click();
    flush();
    expect(t.index!.value).toBe(9);
    expect(scrubber.getAttribute('aria-valuetext')).toBe('model 10 of 10');
    m.button('Next model')!.click();
    expect(t.index!.value).toBe(0);
  });

  it('durationless: no scrubber; play, the status said, a step while paused, restart, speed', () => {
    const d = new Durationless({ speeds: [1, 2] });
    const m = mount(() => <TimelineRow model={d} />);
    expect(m.container.querySelector('input[type=range]')).toBeNull();
    const status = m.container.querySelector('[role=progressbar]')!;
    expect(status.getAttribute('aria-valuetext')).toBe('paused at 0 s, step 0');
    m.button('One step')!.click();
    flush();
    expect(status.getAttribute('aria-valuetext')).toBe('paused at 0 s, step 1');
    m.button('Play')!.click();
    flush();
    expect(status.getAttribute('aria-valuetext')).toBe('running, 0 s, step 1');
    expect(m.container.querySelector<HTMLButtonElement>('button[data-tip="One step (pause first)"]')!.disabled).toBe(
      true,
    );
    m.button('Restart')!.click();
    flush();
    expect(status.textContent).toBe('0 s · step 0');
    expect(m.container.querySelector('[aria-label="Speed"]')).not.toBeNull();
    d.setSpeed(2);
    flush();
    expect(m.container.querySelector('[aria-label="Speed"]')?.textContent).toContain('2×');
    d.dispose();
  });

  it('shows only what the model has', () => {
    const playing = new State(true);
    const m: TimelineModel = { label: 'Spin', playing, togglePlay: () => playing.set(!playing.value) };
    const v = mount(() => <TimelineRow model={m} />);
    expect(v.container.querySelector('[role=group]')?.getAttribute('aria-label')).toBe('Spin');
    expect(v.container.querySelector('[role=progressbar]')?.getAttribute('aria-valuetext')).toBe('running');
    expect(v.button('Restart')).toBeNull();
    expect(v.container.querySelector('[aria-label="Speed"]')).toBeNull();
  });
});

describe('playback is accessible (axe)', () => {
  it('a deck of every kind, a card, an open list', async () => {
    const many = new Stepper({ steps: steps(20) });
    const m = mount(() => (
      <div>
        <PlaybackDeck>
          <StepperRow model={new Stepper({ steps: steps(3) })} />
          <StepperRow model={many} label="Scenes" />
          <TimelineRow model={new Timeline({ duration: 3, speeds: [1, 2] })} />
          <TimelineRow model={new Timeline({ count: 12 })} />
          <TimelineRow model={new Durationless({ label: 'Simulation', speeds: [1, 2] })} />
        </PlaybackDeck>
        <StepCard model={new Stepper({ steps: steps(4) })} />
      </div>
    ));
    m.container.querySelector<HTMLButtonElement>('button[aria-expanded]')!.click();
    await tick();
    expect(await axeViolations(m.container)).toBe('');
  });
});
