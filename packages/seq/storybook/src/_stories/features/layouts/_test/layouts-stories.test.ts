/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Layouts story (a tracks section past its cap), headless on a fake clock, composed as the story composes it
// (`compose.ts`: each cap's spec and config, a viewer and its SVG export), light and dark: the box at the cap, the ruler
// fixed while the wheel, the rail and the keys scroll the tracks below it, and the SVG showing every track.

import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { type CapEntry, Caps, capConfig, capSpec, Genes } from '../compose.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open(cap: CapEntry, dark: boolean): Promise<SeqViewer> {
  vi.useFakeTimers({ toFake: ['Date', 'performance'] });
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({ ...base, ui: undefined, config: [...(base.config ?? []), ...capConfig(cap)] });
  viewer = v;
  await v.load(capSpec(cap));
  v.resize(900);
  await v.settle();
  if (dark) v.setTokens(DarkTokens);
  v.draw();
  v.input.mac = false;
  return v;
}

const section = (v: SeqViewer) => v.frame!.layout.sections.find((s) => s.section.address === Genes)!;
const box = { default: 600, host: 300, spec: 200 } as Record<string, number>;

describe('Features/Layouts, a tracks section past its cap, headless', () => {
  for (const cap of Caps)
    for (const dark of [false, true])
      it(`${cap.id} (${dark ? 'dark' : 'light'}): the box at the cap; wheel, rail and keys scroll under the ruler; the SVG every track`, async () => {
        const v = await open(cap, dark);
        expect(v.tokens).toBe(dark ? DarkTokens : LightTokens);
        const s0 = section(v);
        expect(s0.placement.box.height).toBe(box[cap.id]);
        const ruler = s0.tracks[0].box.y;
        const end = s0.placement.box.y + s0.placement.box.height;
        // the wheel near the box's foot (over the lanes band, it scrolls its lanes first: the chain)
        const gene = s0.tracks[s0.tracks.length - 1].box;
        for (let i = 0; i < 40 && !(section(v).placement.top! > 0); i++) {
          expect(v.input.wheel(gene.x + 20, end - 10, { deltaX: 0, deltaY: 40, deltaMode: 0 })).toBe(true);
          // the section's wheel lands at the frame (coalesced as the navigation's)
          await new Promise((r) => setTimeout(r, 20));
          v.draw();
        }
        expect(section(v).placement.top).toBeGreaterThan(0);
        expect(section(v).tracks[0].box.y).toBe(ruler);
        // End with the pointer in the box, then Home
        v.input.move(gene.x + 20, end - 10);
        expect(v.input.key({ key: 'End' })).toBe(true);
        v.draw();
        const s1 = section(v);
        expect(s1.placement.top).toBeCloseTo(s1.placement.contentHeight - s1.placement.box.height, 6);
        // the rail's window dragged back to the top
        const r = s1.sectionRail!;
        const x = r.top + 6;
        const y0 = r.x + ((r.start + r.end) / 2 / r.total) * r.width;
        expect(v.input.down(x, y0)).toBe(true);
        v.input.move(x, y0 - r.width);
        v.input.up();
        v.draw();
        expect(section(v).placement.top).toBe(0);
        // the SVG: every gene, the canvas left where it was
        const svg = (await v.renderSvg()).svg;
        for (const n of [1, 15, 30]) expect(svg).toContain(`>Gene ${n}<`);
        expect(section(v).placement.top).toBe(0);
      });
});
