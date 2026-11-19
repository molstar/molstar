/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '@molstar/seqx-examples';
import { Section } from '@molstar/seqx-examples/p53/story';
import { frameView, sample, storySteps, type TransitionPlan } from '@molstar/seqx-core/story';
import { headlessAnimation, type HeadlessClock } from '@molstar/seqx-viewer/animation';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import type { SeqAnimationPort } from '@molstar/seqx-viewer/ports';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';

// The p53 story's composition, DOM-free: the story and its headless twin call these. The application owns the clock:
// the page's frames in the story, a fake clock in the twin and the export. A stand-in structure camera follows the same
// plans through the animation port (its zoom and the angle its centre turns to), so two views move on one clock.

/** A stand-in for a structure viewer's camera: where it looks, as the plan's view of the protein section says. */
export interface Camera {
  /** Degrees: the protein's centre in view, read round a turn. */
  angle: number;
  /** How far in: the protein's length over the window's. */
  zoom: number;
  /** Frames it was moved. */
  moves: number;
}

export interface StorySession {
  readonly viewer: SeqViewer;
  readonly camera: Camera;
  /** The fake clock, headless. */
  readonly clock?: HeadlessClock;
  /** The step shown, by the story's buttons. */
  step: string | undefined;
  /** The plan a scrub drives (the slider's), until it is released. */
  scrubbing?: TransitionPlan;
  dispose(): void;
}

/** The protein's length (the section's whole window). */
const lengthOf = (v: SeqViewer) => v.view.resolved.value?.sections.find((s) => s.address === Section)?.axis.length ?? 1;

/** Moves the camera to a plan at `t`: what a structure viewer's camera would read of the same plan. */
export function moveCamera(camera: Camera, plan: TransitionPlan, t: number, length: number): void {
  const w = frameView(plan, sample(plan, t))?.sections[Section];
  if (!w) return;
  camera.angle = (((w.start + w.end) / 2 / length) * 360) % 360;
  camera.zoom = length / Math.max(1e-9, w.end - w.start);
  camera.moves++;
}

/** An animation port that runs the viewer's plans and moves the camera with them, on the clock of `inner`. */
function withCamera(inner: SeqAnimationPort, camera: Camera, length: () => number): SeqAnimationPort {
  return (participant) =>
    inner({
      progress: (plan, t) => {
        participant.progress(plan, t);
        moveCamera(camera, plan, t, length());
      },
      show: (plan, frame) => participant.show(plan, frame),
      clear: () => participant.clear(),
      draw: () => participant.draw(),
      drawn: participant.drawn,
      durationOf: (plan) => participant.durationOf(plan),
    });
}

/** A viewer of the story, mounted into `target` when given (the page's clock), else headless on a fake clock. */
export function composeStory(view: ExampleView, options: { readonly target?: HTMLElement } = {}): StorySession {
  const base = DefaultSeqViewerSpec();
  const camera: Camera = { angle: 0, zoom: 1, moves: 0 };
  const headless = options.target ? undefined : headlessAnimation();
  let viewer: SeqViewer | undefined;
  const length = () => (viewer ? lengthOf(viewer) : 1);
  viewer = createSeqViewer(
    {
      ...base,
      ...(options.target ? {} : { ui: undefined }),
      animation: withCamera(headless?.port ?? base.animation!, camera, length),
    },
    options.target,
  );
  const v = viewer;
  return {
    viewer: v,
    camera,
    ...(headless ? { clock: headless.clock } : {}),
    step: undefined,
    dispose: () => v.retire(),
  };
}

/** The story's steps, in order: id and label. */
export function storyRows(s: StorySession): { id: string; label: string }[] {
  const root = s.viewer.spec.root.value;
  return storySteps(root).map((id) => {
    const n = root.children?.find((c) => c.kind === 'snapshots')?.children?.find((c) => c.id === id);
    return { id, label: typeof n?.params?.label === 'string' ? n.params.label : id };
  });
}

/** Plays a step (or leaves the story, `undefined`) on the session's clock; what it says. */
export function play(s: StorySession, id: string | undefined): string {
  s.scrubbing = undefined;
  const { issue } = s.viewer.playStep(id);
  if (issue) return issue.message;
  s.step = id;
  return id === undefined ? 'Left the story' : `Playing ${storyRows(s).find((r) => r.id === id)?.label ?? id}`;
}

/** The step after the one shown (the first when none is), else undefined at the end. */
export function nextStep(s: StorySession): string | undefined {
  const rows = storyRows(s);
  const i = rows.findIndex((r) => r.id === s.step);
  return rows[i + 1]?.id;
}

/** Starts a scrub to a step: its plan is shown at its start, for the slider (the application's clock) to drive. */
export function scrubTo(s: StorySession, id: string): string {
  const { plan, issue } = s.viewer.playStep(id, { run: false });
  if (issue || !plan) return issue?.message ?? 'Nothing to scrub';
  s.scrubbing = plan;
  s.step = id;
  moveCamera(s.camera, plan, 0, lengthOf(s.viewer));
  return `Scrub to ${id}`;
}

/** The slider moved: the scrub's plan at `t`, on the viewer and the camera alike. */
export function scrub(s: StorySession, t: number): void {
  const plan = s.scrubbing;
  if (!plan) return;
  s.viewer.animation.progress(plan, t);
  moveCamera(s.camera, plan, t, lengthOf(s.viewer));
}

/** The slider released: the step's state is shown. */
export function endScrub(s: StorySession): void {
  if (!s.scrubbing) return;
  s.scrubbing = undefined;
  s.viewer.animation.clear();
}

/** The camera as a line of text. */
export function cameraText(c: Camera): string {
  return `camera: ${c.angle.toFixed(1)}°, zoom ${c.zoom.toFixed(2)}×`;
}

/**
 * The story exported headlessly: an SVG per step, its state and view at its end (each step entered as a cut: what the
 * transition ends on). A session of its own, on a fake clock.
 */
export async function exportSteps(
  view: ExampleView,
  options: { readonly width?: number } = {},
): Promise<{ id: string; svg: string }[]> {
  const s = composeStory(view);
  try {
    s.viewer.resize(options.width ?? 720);
    await s.viewer.load(view.build());
    await s.viewer.settle();
    s.viewer.draw();
    const out: { id: string; svg: string }[] = [];
    for (const { id } of storyRows(s)) {
      s.viewer.playStep(id);
      s.clock!.runAll(16);
      // data a step brings (the hotspot's range) is read before its picture
      await s.viewer.settle();
      s.viewer.draw();
      out.push({ id, svg: (await s.viewer.renderSvg({ dpr: 1, idPrefix: `story-${id}-` })).svg });
    }
    return out;
  } finally {
    s.dispose();
  }
}
