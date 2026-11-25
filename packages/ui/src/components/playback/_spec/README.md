# Playback: the specification

`components/playback/`: generic playback controls over small DOM-free models. A model says what it knows; a row shows
only what that allows. Code comments cite its sections (`playback §2`). No app's names: what is stepped or played (a
story's scenes, a structure's models, a spin) is the app's.

| File                | What                                                                                          |
| ------------------- | --------------------------------------------------------------------------------------------- |
| `stepper-model.ts`  | `StepperModel` (interface), `Stepper` (reference), `withoutTitle` (§1, §2)                    |
| `timeline-model.ts` | `TimelineModel`, `timelineKind`, what a row shows and says, `Timeline`, `Durationless` (§1–2) |
| `clock.ts`          | `PlaybackClock`, `TimerClock`: the reference models' ticks (§2)                               |
| `labels.ts`         | `PlaybackLabels`, `DefaultPlaybackLabels`, `formatElapsed`: every word the views say (§3)     |
| `playback.tsx`      | `PlaybackDeck`, `StepperRow`, `StepCard`, `TimelineRow`, `MuteButton` (§3)                    |

## 1. Three kinds

- **Stepped** (`StepperModel`): a count of steps with titles, the current one (-1: none), a description (markdown, or
  plain text), playing (a slideshow); previous, next, go to, toggle play, `openLink`. From none, next is the first and
  previous the last; `goTo` holds to the steps and ignores NaN.
- **Timed** (`TimelineModel` with a `duration` in seconds: `time`, `step`) and **counted** (with a `count` of frames:
  `index`, 0-based): `seek` to a time or a frame, playing, toggle play; `periodic`: past the end to the start.
- **Durationless** (`TimelineModel` with neither): something that runs until stopped (a spin, a live simulation, a
  stream). Play and pause, how long it ran (`elapsed`, s) and/or how far (`steps`), an optional `restart`, `stepOnce`
  (one step while paused). No scrubber, nothing to seek, no end. Why: not every playback has a length, and a scrubber
  over an unknown end lies.
- Any kind may offer a `speed` from a few `speeds` (`setSpeed`). `timelineKind` decides by the fields present, a
  `duration` first; the kind of a model never changes.

## 2. Models

- **Headless** (ui §0.1): models are DOM-free and tested in plain Node; every operation is a method. What a row says is
  a pure function of the model and the labels (`timelinePosition`, `timelineValueText`, `timelineShown`), so a headless
  host and the row agree.
- **Models hold no words**: no name, no unit. A row's name is the view's `label` prop (a string or a state), the same on
  `StepperRow`, `StepCard` and `TimelineRow`; what a frame is called ("model") is the app's `labels`. Why: one way to
  name a row, and every word through `labels`.
- **The reference models** (`Stepper`, `Timeline`, `Durationless`) are for stories, tests and a headless host; an app
  implements the interfaces over what it plays. They tick through an injected `PlaybackClock` (default the host's
  `setInterval`), **only while playing**: paused, they hold no timer. Why: idle is idle, and fake timers drive them.
- A timeline shortened to nothing to play (no length, one frame) stops (`setLength`; a tick finding it so too): no NaN
  time, no tick past it. `seek(NaN)` does nothing.
- A durationless model's `elapsed` is run time: a step takes `stepMs` over the speed it ran at, counted on from the last
  speed change (no drift from summing; half as long a step at 2×).

## 3. Views

- **Views render state; they never tick.** A row binds the model's states; hidden (a paused scope) it holds nothing.
- **Every word is `labels`** (`PlaybackLabels`, `DefaultPlaybackLabels` English; parametrized ones are functions), one
  prop on each view, the screen reader's value text too (`timedValue`, `countedValue`, `running`, `paused`).
- **`PlaybackDeck`** (`section`, named `labels.deck`): a floating surface (raised, hairline, float shadow), a row a
  header high each, a hairline between; not there when no row shows (`:has([data-row])`), so a row of the app's own
  carries `data-row`.
- **`StepperRow`** (from two steps on): a dot a step up to `MaxDots` (7; `aria-current="step"`, each named by its
  title), past it "i / N", a button that opens the titles' list under the row (a filter from 12; a pick goes there and
  closes it; Esc closes it; the focus comes back to the button). Previous, play, next, then `actions` (a `MuteButton`).
- **`StepCard`**: the title, the stepper's controls from two steps on, the description without a first heading that only
  repeats the title (`withoutTitle`), in its own scroll; `bare` inside a surface.
- **`TimelineRow`**: timed and counted: where it is ("0.4 / 2.0 s", "3 / 10"), a range scrubber (`aria-valuetext` "1.2 s
  of 4.0 s", "frame 3 of 10": as people count), previous and next frame (counted), play; shown while it has a length
  (counted: more than one frame). The scrubber offers its value to a frame state that seeks once a frame (dom §3; never
  equal, so a value sought before seeks again after the model moved). Durationless: play, a status (`role="progressbar"`
  with no value: indeterminate, its `aria-valuetext` "running, 12 s", "paused at step 340"; not live, it changes every
  step), one step (disabled while running), restart, speed.
- Play is one button whose name says what a press does (Play, Pause). Buttons are sm ghost, muted over content.
- Not here: a periodic box, or an app's names for what it plays (its `label`, its `labels`).

## Tests

`_test/models.test.ts` (Node, fake timers: stepping, the three kinds, what is said, no timer while paused),
`without-title.test.ts`, `playback.test.tsx` (jsdom: which rows show, presses and scrubs reach the model, axe). Stories:
UI › Playback (Stepper, Timeline, Durationless, StepCard, Deck).
