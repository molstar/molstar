/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Case 9's pixels: case 9 through the tiles provider draws each of a surface's four canvases per view, section, zoom and
// theme (DPR 1) to the hash recorded for this user agent (`case9-pixels.json`, `records` by `navigator.userAgent`). A
// hash is exact and the machine's (Skia, the Chrome version, the OS's rasterizer, the fixture fonts), so each setup
// keeps its own record: a user agent with none fails and says so (never a silent skip), and records its hashes with
// `VITE_SEQ_PIXELS=record` once its render has been checked by eye. Re-recorded on 2026-10-13 (an M4 with 16 GB,
// Chrome 155 headless) after the batch's look settled; SW1's migration record (`parity/baselines/pixels.json`, the old
// and next stacks' hashes on an M5, Chrome 154) stays as the switch's history. Mounted as the parity layer did.

import { describe, expect, it } from 'vitest';
import { commands } from 'vitest/browser';
import { case9 } from '../../../core/_test/format/cases.js';
import { fixtureLoaders, MolstarFixture } from '../../../core/_test/format/load-fixtures.js';
import { mountSection } from '../../src/canvas/index.js';
import { createSeqViewer, DefaultSeqViewerSpec, manualScheduler } from '../../src/index.js';
import recorded from './case9-pixels.json' with { type: 'json' };

type J = any;
const [W, H] = [1000, 420];
const LIGHT =
  '--seq-surface:rgb(255 255 255);--seq-surface-alt:rgb(246 247 249);--seq-ink:rgb(27 31 36);--seq-ink-quiet:rgb(91 99 112);' +
  '--seq-ink-faint:rgb(196 200 207);--seq-grid:rgb(229 231 235);--seq-accent:rgb(181 69 44);--seq-hover:rgb(43 127 255 / 0.14);' +
  '--seq-highlight:rgb(242 180 0 / 0.22);--seq-selection:rgb(43 127 255 / 0.1);--seq-selection-border:rgb(43 127 255 / 0.6);' +
  '--seq-unobserved:rgb(214 217 222);--seq-pending:rgb(226 222 214);--seq-grip:rgb(124 131 142);--seq-grip-hover:rgb(27 31 36);' +
  '--seq-conflict:rgb(192 38 45);--seq-focus-ring:rgb(43 127 255)';
const DARK =
  '--seq-surface:rgb(24 26 30);--seq-surface-alt:rgb(34 37 43);--seq-ink:rgb(230 232 235);--seq-ink-quiet:rgb(156 163 175);' +
  '--seq-ink-faint:rgb(75 80 90);--seq-grid:rgb(48 52 60);--seq-accent:rgb(240 128 96);--seq-hover:rgb(110 160 255 / 0.2);' +
  '--seq-highlight:rgb(250 204 21 / 0.28);--seq-selection:rgb(110 160 255 / 0.14);--seq-selection-border:rgb(110 160 255 / 0.7);' +
  '--seq-unobserved:rgb(60 64 72);--seq-pending:rgb(56 53 48);--seq-grip:rgb(140 147 158);--seq-grip-hover:rgb(240 242 245);' +
  '--seq-conflict:rgb(255 123 114);--seq-focus-ring:rgb(110 160 255)';
const THEMES = { light: LIGHT, dark: DARK };
const FONTS = '--seq-font-ui:Arial;--seq-font-mono:"Courier New"';
const CANVASES = ['main', 'overlay', 'nav', 'navOverlay'] as const;

/** The canvases of a surface as pixels, the vertical overview's and its overlay drawn into the overlay (as the parity layer). */
function pixels(surface: J): Uint8ClampedArray[] {
  const cs = Array.from(surface.element.querySelectorAll('canvas')) as HTMLCanvasElement[];
  const read = (c: HTMLCanvasElement) =>
    c.width * c.height === 0
      ? new Uint8ClampedArray(0)
      : (c.getContext('2d') as CanvasRenderingContext2D).getImageData(0, 0, c.width, c.height).data.slice();
  if (cs.length < 6) return cs.slice(0, 4).map(read);
  const over = document.createElement('canvas');
  over.width = cs[1].width;
  over.height = cs[1].height;
  const g = over.getContext('2d', { willReadFrequently: true }) as CanvasRenderingContext2D;
  const dpr = cs[1].width / Number.parseFloat(cs[1].style.width);
  const left = (c: HTMLCanvasElement) => Math.round(Number.parseFloat(c.style.left) * dpr);
  const top = (c: HTMLCanvasElement) => Math.round(Number.parseFloat(c.style.top) * dpr);
  for (const k of [4, 5]) if (cs[k].width * cs[k].height > 0) g.drawImage(cs[k], left(cs[k]), top(cs[k]));
  g.drawImage(cs[1], 0, 0);
  return [read(cs[0]), read(over), read(cs[2]), read(cs[3])];
}

/** The first 16 hex digits of a SHA-256 of the pixels. */
async function sha(data: Uint8ClampedArray): Promise<string> {
  const d = await crypto.subtle.digest('SHA-256', new Uint8Array(data));
  return Array.from(new Uint8Array(d).slice(0, 8), (b) => b.toString(16).padStart(2, '0')).join('');
}

/** Case 9's every view and section, at the fit and zoomed, light and dark: `view/section@zoom/theme` -> the four hashes. */
async function hashes(): Promise<Record<string, Record<string, string>>> {
  const doc = case9();
  const views = (doc.root.children ?? []).filter((v) => v.kind === 'view' && v.ref).map((v) => v.ref as string);
  const out: Record<string, Record<string, string>> = {};
  for (const view of views)
    for (const theme of Object.keys(THEMES) as (keyof typeof THEMES)[]) {
      const base = DefaultSeqViewerSpec();
      const spec = { ...base, registry: [...(base.registry ?? []), { svs: MolstarFixture }] };
      const manual = manualScheduler(0);
      const ctx = createSeqViewer(spec as never, { loaders: fixtureLoaders(), scheduler: manual });
      const r = await ctx.setDocument(doc);
      if (!r.ok) throw new Error(JSON.stringify(r.issues));
      ctx.views.setCurrent(view);
      await ctx.graph.update();
      const sections: string[] = [...ctx.viewports.sections.value];
      const hosts: HTMLElement[] = [];
      const surfaces = sections.map((s) => {
        const host = document.createElement('div');
        host.style.cssText = `position:absolute;left:0;top:0;width:${W}px;height:${H}px;${FONTS};${THEMES[theme]}`;
        document.body.appendChild(host);
        hosts.push(host);
        return mountSection(ctx, s, host, { dpr: 1, readback: true }) as J;
      });
      const advance = () => {
        for (let k = 0; k < 3; k++) manual.advance(16);
        for (const s of surfaces) s.renderNow();
      };
      // Lazy data arrives on macrotasks: settle on them too.
      for (let k = 0; k < 12; k++) {
        await new Promise((res) => setTimeout(res, 40));
        advance();
      }
      try {
        for (const zoom of [false, true]) {
          for (const s of sections) {
            const vp = ctx.viewports.get(s) as J;
            if (!vp) continue;
            if (zoom) {
              vp.setView(Math.max(vp.live.current.minCell * 4, 14), 0, 'api');
              vp.panBy(300, 0, 'api');
            } else vp.fit();
          }
          advance();
          advance();
          for (const [k, s] of sections.entries()) {
            const px = pixels(surfaces[k]);
            if (px[0].length === 0) continue;
            const key = `${view}/${s}`.replaceAll('.', '_');
            out[`${key}@${zoom ? 'zoomed' : 'fit'}/${theme}`] = Object.fromEntries(
              await Promise.all(CANVASES.map(async (c, i) => [c, await sha(px[i])] as const)),
            );
          }
        }
      } finally {
        for (const s of surfaces) s.dispose();
        ctx.dispose();
        for (const h of hosts) h.remove();
      }
    }
  return out;
}

/** The records: per user agent, when and where its hashes were taken and the hashes (`view/section@zoom/theme`). */
interface Records {
  _about: string;
  records: Record<string, { recorded: string; machine: string; hashes: Record<string, Record<string, string>> }>;
}
/** The record file for `commands.readFile` and `writeFile`, which resolve paths from the package. */
const FILE = '_test/canvas/case9-pixels.json';
const env = (import.meta as { env?: Record<string, string | undefined> }).env ?? {};
const RECORD = env.VITE_SEQ_PIXELS === 'record';
const MACHINE = env.VITE_SEQ_PIXELS_MACHINE ?? '';

describe('case 9 through the tiles provider (the hashes recorded for this user agent)', () => {
  it('draws every canvas of every section to its recorded hash', async () => {
    const ua = navigator.userAgent;
    const records = recorded as Records;
    const got = await hashes();
    if (RECORD) {
      // The file as it is now (another setup may have recorded since this module was bundled), this agent's set.
      const now = JSON.parse(await commands.readFile(FILE)) as Records;
      now.records[ua] = { recorded: new Date().toISOString().slice(0, 10), machine: MACHINE, hashes: got };
      await commands.writeFile(FILE, `${JSON.stringify(now, null, 2)}\n`);
      return;
    }
    const mine = records.records[ua];
    if (!mine)
      throw new Error(
        `No case 9 pixels are recorded for this user agent (${ua}); recorded: ${Object.keys(records.records).join(
          ' | ',
        )}. Check the render, then record with VITE_SEQ_PIXELS=record (and VITE_SEQ_PIXELS_MACHINE="what it runs on").`,
      );
    expect(got).toEqual(mine.hashes);
  }, 120_000);
});
