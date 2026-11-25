/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import { Domains } from '@molstar/seqx-examples/p53/data';
import { Hotspot, Section, steps } from '@molstar/seqx-examples/p53/story';
import type { RangesLay } from '@molstar/seqx-core/render';
import { frameView, sample, stateOf, stripStory } from '@molstar/seqx-core/story';
import { hashNode, resolveAddress } from '@molstar/core/repr/tree';
import { afterEach, describe, expect, it } from 'vitest';
import { checkMarks } from '../../../examples/_test/draw-checks.js';
import {
  composeStory,
  endScrub,
  exportSteps,
  moveCamera,
  nextStep,
  play,
  scrub,
  scrubTo,
  type StorySession,
  storyRows,
} from '../story-compose.js';

// The p53 story, headless on the animation port's fake clock, composed as the story composes it: each step played to
// its end shows its state and view, drawn against the data (every domain in the window has its marks, every mark is
// hit); a scrub drives the viewer and the stand-in camera from one plan at one time; the export's SVG per step is the
// picture a viewer shows at the step's end.

const view = P53Views.find((v) => v.id === 'story')!;

let session: StorySession | undefined;
afterEach(() => {
  session?.dispose();
  session = undefined;
});

async function open(): Promise<StorySession> {
  const s = composeStory(view);
  session = s;
  s.viewer.resize(720);
  await s.viewer.load(view.build());
  await s.viewer.settle();
  s.viewer.draw();
  return s;
}

/** The domains drawn in the frame: the item keys of their marks. */
function drawnDomains(s: StorySession): Set<string> {
  const out = new Set<string>();
  for (const l of s.viewer.frame!.layout.sections.flatMap((x) => x.tracks.flatMap((t) => t.layers))) {
    if (l.type.name !== 'ranges') continue;
    const lay = l.lay as RangesLay;
    for (let j = 0; j < (lay.hn ?? 0); j++) out.add(String(lay.keys[lay.hitem[j]]));
  }
  return out;
}

describe('the p53 story, headless', () => {
  it('stores four steps as patches over the spec, each giving its state back', () => {
    const spec = view.build();
    const states = steps();
    expect(storyRowsOf(spec)).toEqual(states.map((x) => x.id));
    for (const x of states) {
      const st = stateOf(spec, x.id)!;
      expect(st).toMatchObject({ form: 'patch', applies: true, chain: 0 });
      expect(hashNode(stripStory(st.root))).toBe(hashNode(stripStory(x.root)));
    }
  });

  it('the R175H step marks residue 175: a track of its own over a range of that one residue, drawn there', async () => {
    const st = stateOf(view.build(), 'r175h')!;
    expect(st.applies).toBe(true);
    const track = resolveAddress(st.root, ['views', 'main', 'protein', 'hotspot']);
    expect(track?.params?.label).toBe('R175H');
    const rows = resolveAddress(st.root, ['data', 'hotspot'])?.params?.rows as { start: Int32Array; end: Int32Array };
    expect([rows.start[0], rows.end[0]]).toEqual([Hotspot.start, Hotspot.end]);
    expect(Hotspot.start + 1).toBe(175);
    // drawn: the step played, its mark is hit at residue 175
    const s = await open();
    play(s, 'r175h');
    s.clock!.runAll(16);
    await s.viewer.settle();
    s.viewer.draw();
    expect(drawnDomains(s).has(Hotspot.id)).toBe(true);
    expect(checkMarks(s.viewer).wrong).toEqual([]);
  });

  it('plays each step to its end: its state and view, drawn against the data', async () => {
    const s = await open();
    let id = nextStep(s);
    let played = 0;
    while (id) {
      expect(play(s, id)).toMatch(/^Playing /);
      // the first step is the fitted view the story opens on: nothing to move
      expect(s.viewer.moving).toBe(id !== 'whole');
      s.clock!.runAll(16);
      s.viewer.draw();
      expect(s.viewer.moving).toBe(false);
      const step = steps().find((x) => x.id === id)!;
      const want = (step.view as { sections: Record<string, { start: number; end: number }> }).sections[Section];
      expect(s.viewer.drawnViewports[Section]).toMatchObject(want);
      expect(hashNode(stripStory(s.viewer.spec.shown.value))).toBe(hashNode(stripStory(step.root)));
      // drawn against the data: every domain overlapping the window has its marks, every mark is hit
      const drawn = drawnDomains(s);
      Domains.id.forEach((d, i) => {
        if (Domains.start[i] < want.end && Domains.end[i] > want.start) expect(drawn.has(d)).toBe(true);
      });
      const marks = checkMarks(s.viewer);
      expect(marks.wrong).toEqual([]);
      expect(marks.marks).toBeGreaterThan(0);
      played++;
      id = nextStep(s);
    }
    expect(played).toBe(4);
    // the camera moved with every frame the clock ran
    expect(s.camera.moves).toBeGreaterThan(4);
  });

  it('a scrub drives the viewer and the camera from one plan at one time', async () => {
    const s = await open();
    expect(scrubTo(s, 'tetramer')).toBe('Scrub to tetramer');
    const plan = s.scrubbing!;
    // from the whole protein into the tetramerization domain: the windows overlap, so a pan (a zoom in along)
    expect(plan.channels.find((c) => c.kind === 'view')).toMatchObject({ path: 'pan' });
    expect(plan.channels.map((c) => c.kind)).toContain('leave');
    scrub(s, 0.5);
    s.viewer.draw();
    const at = frameView(plan, sample(plan, 0.5))!.sections[Section];
    expect(s.viewer.drawnViewports[Section]).toMatchObject({ start: at.start, end: at.end });
    const camera = { angle: 0, zoom: 1, moves: 0 };
    moveCamera(camera, plan, 0.5, 393);
    expect(s.camera.angle).toBeCloseTo(camera.angle, 9);
    expect(s.camera.zoom).toBeCloseTo(camera.zoom, 9);
    endScrub(s);
    expect(s.viewer.spec.frame).toBeUndefined();
    expect(storyRows(s).map((r) => r.id)).toEqual(['whole', 'dna-binding', 'r175h', 'tetramer']);
  });

  it('a transition between two non-overlapping windows goes out, across and in', async () => {
    const s = await open();
    play(s, 'r175h');
    s.clock!.runAll(16);
    const { plan } = s.viewer.playStep('tetramer', { run: false });
    expect(plan!.channels.find((c) => c.kind === 'view')).toMatchObject({ path: 'out-across-in' });
    s.viewer.animation.progress(plan!, 0.5);
    s.viewer.draw();
    const mid = s.viewer.drawnViewports[Section];
    // at the middle the window holds both
    expect(mid.start).toBeLessThanOrEqual(160);
    expect(mid.end).toBeGreaterThanOrEqual(363);
  });

  it('the export draws an SVG per step: the picture a viewer shows at each step’s end', {
    tags: ['medium'],
  }, async () => {
    const svgs = await exportSteps(view);
    expect(svgs.map((x) => x.id)).toEqual(['whole', 'dna-binding', 'r175h', 'tetramer']);
    expect(new Set(svgs.map((x) => x.svg)).size).toBe(4);
    // the R175H step's figure draws its mark's label
    expect(svgs.find((x) => x.id === 'r175h')!.svg).toContain('R175H');
    const s = await open();
    for (const { id, svg } of svgs) {
      play(s, id);
      s.clock!.runAll(16);
      await s.viewer.settle();
      s.viewer.draw();
      expect((await s.viewer.renderSvg({ dpr: 1, idPrefix: `story-${id}-` })).svg).toBe(svg);
    }
  });
});

function storyRowsOf(spec: ReturnType<typeof view.build>): string[] {
  return (spec.children?.find((c) => c.kind === 'snapshots')?.children ?? []).map((c) => c.id!);
}
