/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The playback views (playback §3): rows over the DOM-free models, a deck holding them, a step's card. Views render the
// models' state; the ticking is the models' (or the app's).

import { derive, KeyedSelector, type ReadonlyState, State } from '@molstar/core/reactive';
import { type Child, effect, For, frameState, type MaybeReactive, Show } from '../../dom/index.js';
import { state } from '../../state.js';
import { Button } from '../button/button.js';
import { Icon, type IconName } from '../icon/icon.js';
import { Markdown } from '../markdown/markdown.js';
import { ScrollArea } from '../scroll-area/scroll-area.js';
import { Select } from '../select/select.js';
import { TextInput } from '../text-input/text-input.js';
import { classes } from './_generated/playback.module.css.js';
import { type StepperModel, withoutTitle } from './stepper-model.js';
import {
  formatSpeed,
  type TimelineModel,
  timelineKind,
  timelinePosition,
  timelineShown,
  timelineValueText,
} from './timeline-model.js';

/** Up to this many steps, the position is dots; past it, "i / N" opening the list of their titles. */
export const MaxDots = 7;
/** From this many steps, the list has a filter. */
const FilterFrom = 12;

const None: ReadonlyState<undefined> = new State(undefined);

/**
 * The playback controls in one surface, a row each (`StepperRow`, `TimelineRow`, rows of the app's own with
 * `data-row`), a hairline between them. With no row to show, it is not there.
 */
export function PlaybackDeck(props: { children?: Child; label?: string; class?: string }) {
  return (
    <section
      class={props.class ? `${classes.deck} ${props.class}` : classes.deck}
      aria-label={props.label ?? 'Playback'}
    >
      {props.children}
    </section>
  );
}

/**
 * A stepped model's row (from two steps on): where it is (a dot a step, up to `MaxDots`; past it "i / N", which opens
 * a filterable list of the titles), previous, play, next, then `actions` (a mute).
 */
export function StepperRow(props: {
  model: StepperModel;
  /** What the steps are ("Steps", "Scenes"): the row's name. */
  label?: string;
  icon?: IconName;
  /** At the row's end, after next. */
  actions?: () => Child;
}) {
  return (
    <Show when={derive(props.model.count, (c) => c > 1)}>
      {() => (
        <div class={classes.rowWrap} data-row="">
          <StepperBar {...props} lead />
        </div>
      )}
    </Show>
  );
}

/** A step's card: its title, where it is and the controls (from two steps on), then its description. */
export function StepCard(props: {
  model: StepperModel;
  label?: string;
  /** In a surface already (a card's): no surface of its own. */
  bare?: boolean;
  class?: string;
}) {
  const { model } = props;
  // The description without a heading that only repeats the title (shown above it).
  const description = derive([model.description ?? None, model.title], (d, t) => (d ? withoutTitle(d, t ?? '') : d));
  return (
    <Show when={derive(model.count, (c) => c > 0)}>
      {() => (
        <section
          class={[props.bare ? classes.bare : classes.card, props.class ?? ''].join(' ').trim()}
          aria-label={props.label ?? 'Steps'}
        >
          <StepperBar model={model} label={props.label} title />
          <Show when={description}>
            {(source) => (
              <ScrollArea class={classes.description} contentClass={classes.descriptionContent}>
                <Show
                  when={model.plainText ?? false}
                  fallback={() => (
                    <Markdown source={source} onLink={(href) => model.openLink?.(href)} host={model.markdown} />
                  )}
                >
                  {() => <p class={classes.plain}>{source}</p>}
                </Show>
              </ScrollArea>
            )}
          </Show>
        </section>
      )}
    </Show>
  );
}

/** The stepper's bar: a lead icon or the title, the position, the buttons; the list under it when open. */
function StepperBar(props: {
  model: StepperModel;
  label?: string;
  icon?: IconName;
  actions?: () => Child;
  /** The current title at its start (a card), else a lead icon (a row). */
  title?: boolean;
  lead?: boolean;
}) {
  const { model } = props;
  const label = props.label ?? 'Steps';
  const listOpen = state(false);
  let position: HTMLElement | undefined;
  const close = () => {
    listOpen.set(false);
    position?.focus();
  };
  const several = derive(model.count, (c) => c > 1);
  return (
    <>
      <div class={classes.bar} role="group" aria-label={label}>
        {props.lead ? <Lead icon={props.icon ?? 'list'} /> : null}
        {props.title ? <span class={classes.title}>{model.title}</span> : null}
        <Show when={several}>
          {() => (
            <>
              <Show
                when={derive(model.count, (c) => c > MaxDots)}
                fallback={() => <Dots model={model} label={label} />}
              >
                {() => (
                  <button
                    type="button"
                    class={classes.position}
                    aria-expanded={listOpen}
                    aria-label={derive(model.title, (t) => (t ? `${label}: ${t}` : label))}
                    data-tip={derive(model.title, (t) => t ?? `No ${label.toLowerCase()}`)}
                    ref={(el: HTMLElement) => {
                      position = el;
                    }}
                    onClick={() => listOpen.set(!listOpen.value)}
                  >
                    {derive([model.index, model.count], (i, c) => `${i < 0 ? '–' : i + 1} / ${c}`)}
                  </button>
                )}
              </Show>
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name="chevron-left" />}
                title="Previous"
                onClick={() => model.previous()}
              />
              <PlayButton playing={model.playing} title="Play" onClick={() => model.togglePlay()} />
              <Button
                variant="ghost"
                size="sm"
                icon={<Icon name="chevron-right" />}
                title="Next"
                onClick={() => model.next()}
              />
              {props.actions?.()}
            </>
          )}
        </Show>
      </div>
      <Show when={listOpen}>{() => <StepList model={model} label={label} onClose={close} />}</Show>
    </>
  );
}

/** One dot a step, the current one filled (`aria-current`); each goes to its step, its title its name. */
function Dots(props: { model: StepperModel; label: string }) {
  const { model } = props;
  const current = new KeyedSelector(model.index);
  const dots = derive(model.titles, (titles) => titles.map((title, index) => ({ title, index })));
  return (
    <span class={classes.dots} role="group" aria-label={props.label}>
      <For each={dots} by={(d) => d.index}>
        {(d) => (
          <button
            type="button"
            class={classes.dot}
            data-tip={d.title}
            aria-label={d.title}
            aria-current={derive(current.key(d.index), (c) => (c ? 'step' : undefined))}
            onClick={() => model.goTo(d.index)}
          />
        )}
      </For>
    </span>
  );
}

/** Many steps: their titles under the bar, filterable when long, the current one marked; a pick or Esc closes it. */
function StepList(props: { model: StepperModel; label: string; onClose(): void }) {
  const { model } = props;
  const filter = state('');
  const rows = derive([model.titles, filter], (titles, f) => {
    const q = f.trim().toLowerCase();
    const out: { title: string; index: number }[] = [];
    for (let index = 0; index < titles.length; index++)
      if (!q || titles[index].toLowerCase().includes(q)) out.push({ title: titles[index], index });
    return out;
  });
  const current = new KeyedSelector(model.index);
  let first: HTMLElement | undefined;
  queueMicrotask(() => first?.focus());
  return (
    <div
      class={classes.list}
      onKeyDown={(e: KeyboardEvent) => {
        if (e.key !== 'Escape' || e.defaultPrevented) return;
        e.preventDefault();
        e.stopPropagation();
        props.onClose();
      }}
    >
      <Show when={derive(model.count, (c) => c >= FilterFrom)}>
        {() => (
          <span
            ref={(el: HTMLElement) => {
              first = el.querySelector('input') ?? undefined;
            }}
          >
            <TextInput
              search
              size="sm"
              aria-label={`Filter ${props.label.toLowerCase()}`}
              placeholder="Filter"
              value={filter}
              onInput={(v) => filter.set(v)}
            />
          </span>
        )}
      </Show>
      <ScrollArea class={classes.listRows} contentClass={classes.listRowsContent}>
        <ul class={classes.listItems} aria-label={props.label}>
          <For each={rows} by={(r) => r.index}>
            {(r) => (
              <li>
                <button
                  type="button"
                  class={classes.listRow}
                  aria-current={derive(current.key(r.index), (c) => (c ? 'step' : undefined))}
                  ref={(el: HTMLElement) => {
                    if (!first && r.index === Math.max(0, model.index.value)) first = el;
                  }}
                  onClick={() => {
                    model.goTo(r.index);
                    props.onClose();
                  }}
                >
                  <span class={classes.listIndex}>{String(r.index + 1)}</span>
                  <span class={classes.listTitle}>{r.title}</span>
                </button>
              </li>
            )}
          </For>
        </ul>
      </ScrollArea>
    </div>
  );
}

/**
 * A timed, counted or durationless model's row (`timelineKind`). Timed and counted: where it is, a scrubber (said as
 * "1.2 s of 4.0 s", "frame 3 of 10"), previous and next frame (counted), play; shown while it has a length (timed) or
 * more than a frame (counted). Durationless: no scrubber (nothing to seek, no end): play, how long or how far it ran
 * (said as "running, 12 s", "paused at step 340"), one step while paused, restart. Any kind: a speed when it has
 * choices.
 */
export function TimelineRow(props: { model: TimelineModel; icon?: IconName }) {
  const { model } = props;
  const kind = timelineKind(model);
  const shown =
    kind === 'durationless'
      ? true
      : derive(kind === 'timed' ? model.duration! : model.count!, () => timelineShown(model));
  return (
    <Show when={shown}>{() => (kind === 'durationless' ? <Durationless {...props} /> : <Scrubbed {...props} />)}</Show>
  );
}

function Scrubbed(props: { model: TimelineModel; icon?: IconName }) {
  const { model } = props;
  const timed = timelineKind(model) === 'timed';
  const at = (timed ? model.time : model.index) ?? new State(0);
  const length = (timed ? model.duration : model.count)!;
  const unit = model.unit ?? 'frame';
  const Unit = unit.charAt(0).toUpperCase() + unit.slice(1);
  const label = model.label ?? (timed ? 'Timeline' : `${Unit}s`);
  const position = derive([at, length], () => timelinePosition(model));
  const valueText = derive([at, length], () => timelineValueText(model));
  const max = derive(length, (n) => String(timed ? n : Math.max(0, n - 1)));
  const progress = derive([at, length], (a, n) => {
    const end = timed ? n : n - 1;
    return `${end > 0 ? (100 * Math.min(a, end)) / end : 0}%`;
  });
  // Pointer-rate input (dom §3): the scrubber's value is offered at each input and seeks once a frame. Never equal, so
  // a pick of the value it last sought seeks again after the model moved on.
  const seeking = frameState<number | undefined>(undefined, () => false);
  effect(seeking, (v) => {
    if (v !== undefined) model.seek?.(v);
  });
  // One frame back or forth, past the end to the start when periodic.
  const step = (by: number) => {
    const n = length.value;
    if (n < 2) return;
    const i = at.value + by;
    model.seek?.(model.periodic ? ((i % n) + n) % n : Math.min(Math.max(i, 0), n - 1));
  };
  return (
    <div class={classes.bar} role="group" aria-label={label} data-row="">
      <Lead icon={props.icon ?? 'film'} />
      {model.label ? <span class={classes.label}>{model.label}</span> : null}
      <span class={classes.counter}>{position}</span>
      <input
        class={classes.scrubber}
        type="range"
        min="0"
        max={max}
        step={timed ? derive(model.step ?? None, (s) => (s && s > 0 ? String(s) : 'any')) : '1'}
        value={derive(at, String)}
        style={{ '--msui-progress': progress }}
        aria-label={timed ? 'Time' : Unit}
        aria-valuetext={valueText}
        disabled={!model.seek}
        onInput={(e: Event) => seeking.offer(Number((e.currentTarget as HTMLInputElement).value))}
      />
      {timed ? null : (
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="chevron-left" />}
          title={`Previous ${unit}`}
          onClick={() => step(-1)}
        />
      )}
      <PlayButton playing={model.playing} title="Play" onClick={() => model.togglePlay()} />
      {timed ? null : (
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="chevron-right" />}
          title={`Next ${unit}`}
          onClick={() => step(1)}
        />
      )}
      <Speed model={model} />
    </div>
  );
}

function Durationless(props: { model: TimelineModel; icon?: IconName }) {
  const { model } = props;
  const label = model.label ?? 'Playback';
  const sources = [model.playing, model.elapsed ?? None, model.steps ?? None] as const;
  const position = derive(sources, () => timelinePosition(model));
  const valueText = derive(sources, () => timelineValueText(model));
  return (
    <div class={classes.bar} role="group" aria-label={label} data-row="">
      <Lead icon={props.icon ?? 'film'} />
      {model.label ? <span class={classes.label}>{model.label}</span> : null}
      {/* Indeterminate progress (no end): its value said as text, not live (it changes every step). */}
      <span
        class={classes.status}
        role="progressbar"
        aria-label={`${label} status`}
        aria-valuetext={valueText}
        data-playing={model.playing}
      >
        {position}
      </span>
      <PlayButton playing={model.playing} title="Play" onClick={() => model.togglePlay()} />
      {model.stepOnce ? (
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="step-forward" />}
          title={derive(model.playing, (p) => (p ? 'One step (pause first)' : 'One step'))}
          disabled={model.playing}
          onClick={() => model.stepOnce?.()}
        />
      ) : null}
      {model.restart ? (
        <Button
          variant="ghost"
          size="sm"
          icon={<Icon name="rotate-ccw" />}
          title="Restart"
          onClick={() => model.restart?.()}
        />
      ) : null}
      <Speed model={model} />
    </div>
  );
}

/** A speed choice, when the model has choices. */
function Speed(props: { model: TimelineModel }) {
  const { model } = props;
  if (!model.speed || !model.speeds?.length || !model.setSpeed) return null;
  return (
    <span class={classes.speed}>
      <Select
        size="sm"
        aria-label="Speed"
        title="Speed"
        options={model.speeds.map((s) => [String(s), formatSpeed(s)] as const)}
        value={derive(model.speed, String)}
        onChange={(v) => model.setSpeed?.(Number(v))}
      />
    </span>
  );
}

/** A mute button: the speaker, struck through while muted (pressed). */
export function MuteButton(props: { muted: ReadonlyState<boolean>; onClick(): void }) {
  return (
    <Button
      variant="ghost"
      size="sm"
      icon={
        <Show when={props.muted} fallback={() => <Icon name="volume-2" />}>
          {() => <Icon name="volume-x" />}
        </Show>
      }
      title="Mute"
      pressed={props.muted}
      onClick={() => props.onClick()}
    />
  );
}

/** Play and pause in one button: its name says what a press does. */
function PlayButton(props: { playing: ReadonlyState<boolean>; title: MaybeReactive<string>; onClick(): void }) {
  const title = typeof props.title === 'string' ? new State(props.title) : (props.title as ReadonlyState<string>);
  return (
    <Button
      variant="ghost"
      size="sm"
      icon={
        <Show when={props.playing} fallback={() => <Icon name="play" />}>
          {() => <Icon name="pause" />}
        </Show>
      }
      title={derive([props.playing, title], (p, t) => (p ? 'Pause' : t))}
      onClick={() => props.onClick()}
    />
  );
}

/** A row's leading icon, in a button's slot. */
function Lead(props: { icon: IconName }) {
  return (
    <span class={classes.lead} aria-hidden="true">
      <Icon name={props.icon} />
    </span>
  );
}
