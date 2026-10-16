/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ColorLists } from '@molstar/core/util/color/lists';
import { story } from '@molstar/prism-ui/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import type { DataPalette } from '../sources/palette.js';
import { RealViewerStory } from '../testing/real-viewer.js';

// The data palettes (accessibility.md, viewer-spec.md §2.8): the colors the categorical themes give chains, entities
// and the other categories, in the order they take them (chain A the first). `canvas.palette` picks one (Prism's by
// default); `spec.palettes` gives Prism's and the colorblind-safe lists. The structures come from RCSB (the network).

interface PaletteInfo {
  readonly id: DataPalette;
  readonly title: string;
  readonly about: string;
  readonly colors: readonly number[];
}

const Palettes: readonly PaletteInfo[] = [
  {
    id: 'prism',
    title: 'Prism',
    about:
      'The default. Many-Distinct redrawn in OKLCH: Mol*’s green and orange first, even weight, every pair apart (min ΔE 0.10).',
    colors: ColorLists['many-distinct-prism'].list as number[],
  },
  {
    id: 'classic',
    title: 'Classic',
    about: 'Mol*’s own Many-Distinct (Dark2, Set1, Set2): the themes as they are.',
    colors: ColorLists['many-distinct'].list as number[],
  },
  {
    id: 'colorblind',
    title: 'Colorblind-safe',
    about: 'Okabe–Ito: told apart with any color vision.',
    colors: ColorLists['okabe-ito'].list as number[],
  },
];

const css = (c: number) => `#${c.toString(16).padStart(6, '0')}`;
const letter = (i: number) => String.fromCharCode(65 + (i % 26));

const meta: Meta = { title: 'Viewer/Palettes', parameters: { layout: 'fullscreen' } };
export default meta;

function Swatch(props: { color: number; index: number }) {
  return (
    <div style="display: grid; justify-items: center; gap: 3px" title={`${letter(props.index)}: ${css(props.color)}`}>
      <div style={`width: 30px; height: 30px; border-radius: var(--msp-radius-md); background: ${css(props.color)}`} />
      <span style="font-size: var(--msp-font-size-xs); font-family: var(--msp-font-mono); color: var(--msp-color-text-subtle)">
        {letter(props.index)}
      </span>
    </div>
  );
}

function PaletteRow(props: { palette: PaletteInfo }) {
  const p = props.palette;
  return (
    <section style="display: grid; gap: var(--msp-space-3)">
      <div style="display: flex; align-items: baseline; gap: var(--msp-space-3); flex-wrap: wrap">
        <h3 style="margin: 0; font-size: var(--msp-font-size-md); font-weight: 600">{p.title}</h3>
        <code style="font-size: var(--msp-font-size-xs); font-family: var(--msp-font-mono); color: var(--msp-color-text-subtle)">
          {`canvas.palette: '${p.id}'`}
        </code>
        <span style="font-size: var(--msp-font-size-xs); color: var(--msp-color-text-subtle)">
          {`${p.colors.length} colors`}
        </span>
      </div>
      <p style="margin: 0; max-width: 62ch; font-size: var(--msp-font-size-sm); color: var(--msp-color-text-muted)">
        {p.about}
      </p>
      <div style="display: flex; gap: 6px; flex-wrap: wrap">
        {p.colors.map((c, i) => (
          <Swatch color={c} index={i} />
        ))}
      </div>
    </section>
  );
}

/** The three palettes, each color with the category it goes to (A, B, …), on the light and the dark surface. */
export const Swatches: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 420px), 1fr)); gap: var(--msp-space-4); padding: var(--msp-space-5)">
        {(['light', 'dark'] as const).map((theme) => (
          <div
            data-msp-theme={theme}
            style="display: grid; align-content: start; gap: var(--msp-space-6); padding: var(--msp-space-6); border-radius: var(--msp-radius-lg); background: var(--msp-color-bg); color: var(--msp-color-text)"
          >
            {Palettes.map((p) => (
              <PaletteRow palette={p} />
            ))}
          </div>
        ))}
      </div>
    )),
};

interface Args {
  palette: DataPalette;
  /** A list for `spec.palettes.prism`, comma-separated `#rrggbb` (empty: the default); shown with `prism`. */
  colors: string;
  entry: '1TQN' | '1AON' | '4HHB' | '1KX5';
  theme: 'chain-id' | 'entity-id';
  background: 'light' | 'dark';
}

/**
 * A structure from RCSB in the viewer, colored by chain or entity with the palette the options pick (and the list the
 * spec gives, if any): cytochrome P450 3A4 (1TQN: a chain, its heme and ligand), GroEL–GroES (1AON, 21 chains),
 * hemoglobin (4HHB, 4), a nucleosome (1KX5, 10: protein and DNA).
 */
export const OnAStructure: StoryObj<Args> = {
  name: 'On a structure',
  args: { palette: 'prism', colors: '', entry: '1AON', theme: 'chain-id', background: 'light' },
  argTypes: {
    palette: { control: 'inline-radio', options: Palettes.map((p) => p.id) },
    colors: { control: 'text' },
    entry: { control: 'select', options: ['1TQN', '1AON', '4HHB', '1KX5'] },
    theme: { control: 'inline-radio', options: ['chain-id', 'entity-id'] },
    background: { control: 'inline-radio', options: ['light', 'dark'] },
  },
  render: (args) =>
    story(() => {
      const own = args.colors
        .split(/[\s,]+/)
        .filter((c) => /^#?[0-9a-f]{6}$/i.test(c))
        .map((c) => (c.startsWith('#') ? c : `#${c}`));
      return RealViewerStory({
        spec: (base) => ({
          ...base,
          options: { ...base.options, canvas: { ...base.options?.canvas, palette: args.palette } },
          ...(own.length ? { palettes: { prism: own } } : {}),
        }),
        setup: async (viewer) => {
          const person = viewer.api.as({ kind: 'person' });
          const response = await fetch(`https://models.rcsb.org/${args.entry}.bcif`);
          if (!response.ok) throw new Error(`${args.entry}: ${response.status}`);
          const file = new File([await response.blob()], `${args.entry.toLowerCase()}.bcif`);
          await person.load.files({ files: [file] as never });
          await person.representation.color({ target: 'polymer', theme: args.theme });
          if (args.background === 'dark') await person.canvas.background({ color: '#101114' });
        },
      });
    }),
};

/** Cytochrome P450 3A4 (1TQN) by entity: the protein, its heme and its ligand. */
export const On1TQN: StoryObj<Args> = {
  ...OnAStructure,
  name: 'On 1TQN',
  args: { ...OnAStructure.args, entry: '1TQN', theme: 'entity-id' },
};
