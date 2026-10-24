/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { CommandLanguageExtension } from '@molstar/command-language-extension/prism-viewer';
import { ConsoleTab } from '@molstar/command-language-extension/prism-viewer/parts';
import { story } from '@molstar/prism-ui/testing/index';
import { ExampleFiles, RealViewerStory } from '@molstar/prism-viewer/testing/real-viewer';
import type { PrismViewer } from '@molstar/prism-viewer';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';

// PyMOL and VMD commands in the Prism viewer (the command-language extension's _spec): typed in the command menu after
// `pymol `, `vmd ` or `!`, or in the Console tab, each submission one Undo step through the viewer's operations. 1CBS
// from data/examples; no network.

const meta: Meta = { title: 'Viewer/With extensions/PyMOL and VMD commands', parameters: { layout: 'fullscreen' } };
export default meta;

const withCommands = (base: PrismViewerSpec): PrismViewerSpec => ({
  ...base,
  extensions: [...(base.extensions ?? []), CommandLanguageExtension],
});

const ok = (r: { ok: boolean; message?: string }) => {
  if (!r.ok) throw new Error(r.message);
};

/** The viewer's shell (its shadow root): the story's typing goes there, as a person's would. */
const shell = () => [...document.querySelectorAll('*')].find((e) => e.shadowRoot)?.shadowRoot ?? undefined;

const frame = () => new Promise((r) => requestAnimationFrame(() => r(undefined)));

/** Types a line into a text field and presses Enter (or not), as a person does. */
async function typeInto(field: HTMLInputElement | HTMLTextAreaElement, text: string, enter: boolean) {
  field.focus();
  field.value = text;
  field.dispatchEvent(new Event('input', { bubbles: true, composed: true }));
  await frame();
  if (enter) field.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true, composed: true }));
}

/** The viewer, empty (`loaded`: 1CBS loaded first), then the story's steps. */
function commandsStory(steps?: (v: PrismViewer) => Promise<void>, loaded = false): StoryObj {
  return {
    render: () =>
      story(() => (
        <RealViewerStory spec={withCommands} load={loaded ? [ExampleFiles.cbs] : []} setup={(v) => steps?.(v)} />
      )),
  };
}

/** 1CBS as a file (what a person drops), from data/examples. */
async function cbsFile(): Promise<File> {
  return new File([await (await fetch(ExampleFiles.cbs)).text()], '1cbs.cif');
}

/** Drops files on the console, as a person does. */
async function dropOnConsole(files: readonly File[]) {
  const target = shell()?.querySelector('[aria-label="Console output"]');
  const data = new DataTransfer();
  for (const f of files) data.items.add(f);
  target?.dispatchEvent(new DragEvent('drop', { dataTransfer: data, bubbles: true, composed: true, cancelable: true }));
  await new Promise((r) => setTimeout(r, 1500));
}

/** Opens the command menu (⌘K) and types a line in it. */
async function menuWith(text: string) {
  const root = shell();
  const target = root?.querySelector('canvas') ?? root?.firstElementChild;
  const mac = /Mac|iPhone|iPad/.test(navigator.platform);
  target?.dispatchEvent(
    new KeyboardEvent('keydown', { key: 'k', metaKey: mac, ctrlKey: !mac, bubbles: true, composed: true }),
  );
  await frame();
  const input = shell()?.querySelector<HTMLInputElement>('input[role=combobox]');
  if (input) await typeInto(input, text, false);
}

/** Opens the console's tab (its view is a chunk of its own: the tab is there once the chunk has loaded). */
async function openConsole(v: PrismViewer) {
  for (let i = 0; i < 50 && !(await v.api.view.open({ tab: ConsoleTab })).ok; i++) {
    await new Promise((r) => setTimeout(r, 100));
  }
  await frame();
}

/**
 * The Console as it opens, on an empty viewer: what it is, and example lines of the dialect (a click puts one in the
 * input). The composer at its foot: the input, the Run button, the line that will say what Enter does.
 */
export const ConsoleEmpty: StoryObj = commandsStory(openConsole);

/**
 * The Console tab (the rail's foot, beside the Log), starting empty: 1CBS dropped on it (the line that loads it runs),
 * then four submissions: a PyMOL style (hide everything, cartoon,
 * the ligand as sticks colored by element with green carbons: `util.cbag`), a read (`count_atoms`), a typo the console
 * explains before anything runs, and VMD lines (a VDW rep of the ligand in ColorID 4, yellow). Each submission is one
 * step in History.
 */
export const Console: StoryObj = commandsStory(async (v) => {
  await openConsole(v);
  // 1CBS dropped on the console: the line that loads it runs, and shows in the transcript.
  await dropOnConsole([await cbsFile()]);
  const field = shell()?.querySelector<HTMLTextAreaElement>('textarea');
  if (!field) return;
  for (const line of [
    'hide everything; show cartoon; show sticks, resn REA; util.cbag resn REA',
    'count_atoms resn REA',
    'shw sticks',
    'vmd mol representation VDW; mol selection {resname REA}; mol color ColorID 4; mol addrep top',
  ]) {
    await typeInto(field, line, true);
    await new Promise((r) => setTimeout(r, 300));
  }
});

/**
 * The command menu after `pymol `: the line's Run row says what it will do (Enter runs it as one step), then the
 * completions for the word being typed: here the structure's chains after `chain`.
 */
export const MenuPrefix: StoryObj = commandsStory(async () => {
  await menuWith('pymol show sticks, chain ');
}, true);

/** A line that does not read: the row says why and what was meant, and keeps the menu open to fix it. */
export const MenuError: StoryObj = commandsStory(async () => {
  await menuWith('pymol color slaet, resn REA');
}, true);

/**
 * VMD's way, run as one step on an empty viewer: 1CBS loaded with VMD's look, its Lines rep removed, then a NewCartoon of the protein colored by chain and the
 * ligand as Licorice by Name (VMD's element colors: carbon cyan), each a rep of its own.
 */
export const VmdReps: StoryObj = commandsStory(async (v) => {
  ok(
    await v.api.script!.run({
      dialect: 'vmd',
      files: [await cbsFile()],
      text: [
        // Loaded with VMD's look (one rep: Lines), which the next line removes.
        'mol new 1cbs.cif',
        'mol delrep 0 top',
        'mol representation NewCartoon',
        'mol color Chain',
        'mol selection {protein}',
        'mol addrep top',
        'mol representation Licorice',
        'mol color Name',
        'mol selection {resname REA}',
        'mol addrep top',
      ].join('\n'),
    }),
  );
});
