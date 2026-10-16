/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectCell } from '@molstar/core/state';
import { ColorLists } from '@molstar/core/util/color/lists';
import { story } from '@molstar/prism-ui/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { CategoryThemes } from '../sources/palette.js';
import { RealViewerStory } from '../testing/real-viewer.js';
import type { PrismViewer } from '../viewer.js';

// Candidates for Prism's own data palette (many-distinct-prism, the categorical themes' default: chains, entities, …),
// beside Mol*'s many-distinct and the colorblind-safe Okabe–Ito: swatches in order (chain A takes the first), and a real
// multi-chain structure recolored with the one picked. Made in OKLCH: hues a golden angle apart, lightness alternating,
// chroma kept in gamut, so the first eight (A–H) stay apart. The viewer reaches the network (RCSB).

const hex = (list: readonly string[]) => list.map((h) => Number.parseInt(h, 16));

/** The lists compared, in the order chains take them; the smallest OKLab distance among the first eight in the label. */
const Palettes: Readonly<Record<string, { label: string; colors: readonly number[] }>> = {
  // Mol*'s own order (many-distinct's: Dark2's green, orange, purple, pink, lime, gold, brown), redrawn in OKLCH with
  // even weight, a sky blue where Dark2 has grey, then a second eight filling the hues between, lightness alternating.
  'prism-d': {
    label: 'D · Mol*’s green and orange first, refined (min ΔE 0.103; all 16: 0.102)',
    colors: hex(
      '07a57b e57434 7b74c4 d8588f 78b248 ecb333 915a23 4cbadc cb4644 2a61b1 1a8390 c793df 7f8026 f6a88d 883f73 9bdcb1'.split(
        ' ',
      ),
    ),
  },
  'prism-e': {
    label: 'E · the same, lighter (min ΔE 0.097; all 16: 0.095)',
    colors: hex(
      '46b58f f38c56 8d88d1 e672a1 8fc367 fcc85f a06e3f 6fcceb d9615b 4374bd 4294a0 d8a9ee 909146 ffc0ab 975482 b4efc7'.split(
        ' ',
      ),
    ),
  },
  'prism-a': {
    label: 'A · soft (min ΔE 0.136)',
    colors: hex(
      '55aee8 b96368 60bb83 8273bc d49648 06919f da85b6 758b3a 78a5f0 b8674f 39bda0 976caf c2a140 1b8cb3 e3839a 579153'.split(
        ' ',
      ),
    ),
  },
  'prism-b': {
    label: 'B · vivid (min ΔE 0.146)',
    colors: hex(
      '1fa6ff c34f50 02bb7c 8663c7 d28f09 08899d e96cad 618a02 6e9aff c15625 08b7a1 9e5ab3 ba9c0a 0685b4 f46a86 289142'.split(
        ' ',
      ),
    ),
  },
  'prism-c': {
    label: 'C · the logo’s hues first (min ΔE 0.083)',
    colors: hex(
      '3e8ec1 e7b643 5da37f c06240 ab8be3 a16f04 0eb1ce b35682 83ae53 5875c7 e08256 098e7f c282d0 8e7a07 37aae3 ba5565'.split(
        ' ',
      ),
    ),
  },
  'many-distinct': {
    label: 'Mol*’s many-distinct (min ΔE 0.099; all 16: 0.035)',
    colors: ColorLists['many-distinct'].list as number[],
  },
  'okabe-ito': { label: 'Okabe–Ito, colorblind-safe (min ΔE 0.118)', colors: ColorLists['okabe-ito'].list as number[] },
};
type PaletteId = keyof typeof Palettes;

const css = (c: number) => `#${c.toString(16).padStart(6, '0')}`;

const meta: Meta = { title: 'Viewer/Palettes', parameters: { layout: 'fullscreen' } };
export default meta;

/** Every list's colors in order, on the panel's light and dark surfaces. */
export const Swatches: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; gap: var(--msp-space-5); padding: var(--msp-space-5)">
        {(['light', 'dark'] as const).map((theme) => (
          <div
            data-msp-theme={theme}
            style="display: grid; gap: var(--msp-space-4); padding: var(--msp-space-5); border-radius: var(--msp-radius-lg); background: var(--msp-color-bg); color: var(--msp-color-text)"
          >
            {Object.entries(Palettes).map(([, p]) => (
              <div style="display: grid; gap: var(--msp-space-2)">
                <span style="font-size: var(--msp-font-size-sm); color: var(--msp-color-text-muted)">{p.label}</span>
                <div style="display: flex; gap: 4px; flex-wrap: wrap">
                  {p.colors.map((c, i) => (
                    <div
                      title={`${String.fromCharCode(65 + (i % 26))}: ${css(c)}`}
                      style={`width: 28px; height: 28px; border-radius: 6px; background: ${css(c)}`}
                    />
                  ))}
                </div>
              </div>
            ))}
          </div>
        ))}
      </div>
    )),
};

/** The categorical themes' representations recolored with a list (as the data palette does, `sources/palette.ts`). */
async function recolor(viewer: PrismViewer, colors: readonly number[]): Promise<void> {
  const palette = { name: 'colors', params: { list: { kind: 'set', colors: [...colors] } } };
  const update = viewer.plugin.build();
  let any = false;
  for (const cell of viewer.plugin.state.data.cells.values() as Iterable<StateObjectCell>) {
    const params = cell.transform.params as { colorTheme?: { name: string; params?: Record<string, unknown> } };
    const theme = params?.colorTheme;
    if (!theme || !CategoryThemes.includes(theme.name)) continue;
    update.to(cell.transform.ref).update({ ...params, colorTheme: { ...theme, params: { ...theme.params, palette } } });
    any = true;
  }
  if (any) await update.commit();
}

interface Args {
  palette: PaletteId;
  entry: string;
  theme: 'chain-id' | 'entity-id';
  background: 'light' | 'dark';
}

/**
 * A structure from RCSB colored by chain (or entity) with the picked list: GroEL–GroES (1AON, 21 chains), hemoglobin
 * (4HHB, 4), a nucleosome (1KX5, 10: protein and DNA).
 */
export const OnAStructure: StoryObj<Args> = {
  args: { palette: 'prism-d', entry: '1AON', theme: 'chain-id', background: 'light' },
  argTypes: {
    palette: { control: 'select', options: Object.keys(Palettes) },
    entry: { control: 'select', options: ['1AON', '4HHB', '1KX5'] },
    theme: { control: 'inline-radio', options: ['chain-id', 'entity-id'] },
    background: { control: 'inline-radio', options: ['light', 'dark'] },
  },
  render: (args) =>
    story(() =>
      RealViewerStory({
        setup: async (viewer) => {
          const person = viewer.api.as({ kind: 'person' });
          const response = await fetch(`https://models.rcsb.org/${args.entry}.bcif`);
          if (!response.ok) throw new Error(`${args.entry}: ${response.status}`);
          const file = new File([await response.blob()], `${args.entry.toLowerCase()}.bcif`);
          await person.load.files({ files: [file] as never });
          await person.representation.color({ target: 'polymer', theme: args.theme });
          if (args.background === 'dark') await person.canvas.background({ color: '#101114' });
          await recolor(viewer, Palettes[args.palette].colors);
        },
      }),
    ),
};
