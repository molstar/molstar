/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { map } from 'rxjs';
import { story } from '../../_stories/story-root.js';
import { formatChord } from '../../commands/keys.js';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { For } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { onCleanup } from '../../dom/scope.js';
import { domView } from '../../dom/view.js';
import { parseColor, themeForColor } from '../../model/canvas-theme.js';
import { observableState } from '../../reactive/observable-state.js';
import type { HoverLabelModel, SnapshotStepperModel, TrajectoryControlsModel } from '../../model/playback.js';
import { ShellModel } from '../../shell/shell.js';
import { StatusModel } from '../../shell/status.js';
import { HoverLabel } from '../canvas/hover-label.js';
import { SnapshotPanel, Timeline } from '../canvas/playback.js';
import { SequenceStrip } from '../sequence/sequence-strip.js';
import { CanvasPlaceholder } from './canvas-placeholder.js';
import { ShellView } from './shell-view.js';

const meta: Meta = { title: 'Shell/Shell', parameters: { layout: 'fullscreen' } };
export default meta;

const pad = 'padding: 12px; display: grid; gap: 12px; align-content: start;';
const muted = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm);';

function Placeholder(props: { title: string; children?: Child }) {
  return (
    <div style={pad}>
      <b>{props.title}</b>
      <p style={muted}>{props.children}</p>
    </div>
  );
}

const Sequence =
  'MALYGTHSHGLFKKLGIPGPTPLPFLGNILSYHKGFCMFDMECHKKYGKVWGFYDGQQPVLAITDPDMIKTVLVKECYSVFTNRRPFGPVGFMKSAISIAEDEEWKRLRSLLSPTFTSGKLKEMVPIIAQYGDVLVRNLRREAETGKPVTLKDVFGAYSMDVITSTSFGVNIDSLNNPQDPFVENTKKLLRFDFLDPFFLSITVFPFLIPILEVLNICVFPREVTNFLRKSVKRMKESRLEDTQKHRVDFLQLMIDSQNSKETESHKALSDLELVAQSIIFIFAGYETTSSVLSFIMYELATHPDVQQKLQEEIDAVLPNKAPPTYDTVLQMEYLDMVVNETLRLFPIAMRLERVCKKDVEINGMFIPKGVVVMIPSYALHRDPKYWTEPEKFLPERFSKKNKDNIDPYIYTPFGSGPRNCIGMRFALMNMKLALIRVLQNFSFKPCKETQIPLKLSLGGLLQPEKPVVLKVESRDGTVSGA';
const ThreeLetter: Record<string, string> = {
  A: 'ALA',
  R: 'ARG',
  N: 'ASN',
  D: 'ASP',
  C: 'CYS',
  Q: 'GLN',
  E: 'GLU',
  G: 'GLY',
  H: 'HIS',
  I: 'ILE',
  L: 'LEU',
  K: 'LYS',
  M: 'MET',
  F: 'PHE',
  P: 'PRO',
  S: 'SER',
  T: 'THR',
  W: 'TRP',
  Y: 'TYR',
  V: 'VAL',
}; // prettier-ignore

function mockSnapshots(status: StatusModel): SnapshotStepperModel {
  const entries = [
    {
      title: 'Overview',
      description: 'Cytochrome **P450 3A4** with its heme cofactor.\n\nSee the [active site](#active-site).',
    },
    {
      title: 'Active site',
      description:
        '### Active site\n- Heme iron coordinated by `Cys442`\n- Substrate channel above the heme\n\nBack to the [overview](#overview).',
    },
    { title: 'Surface', description: undefined },
  ];
  const index = observableState(0);
  const isPlaying = observableState(false);
  const transition = observableState<number | undefined>(undefined);
  const title = observableState<string | undefined>(entries[0].title);
  const description = observableState<string | undefined>(entries[0].description);
  let timer: ReturnType<typeof setInterval> | undefined;
  const go = (i: number) => {
    const n = (i + entries.length) % entries.length;
    index.set(n);
    title.set(entries[n].title);
    description.set(entries[n].description);
    let t = 0;
    const step = setInterval(() => {
      t += 0.1;
      transition.set(t >= 1 ? undefined : t);
      if (t >= 1) clearInterval(step);
    }, 30);
  };
  return {
    count: observableState(entries.length).asReadonly,
    index: index.asReadonly,
    title: title.asReadonly,
    description: description.asReadonly,
    isPlaying: isPlaying.asReadonly,
    transition: transition.asReadonly,
    previous: () => go(index.value - 1),
    next: () => go(index.value + 1),
    togglePlay: () => {
      isPlaying.set(!isPlaying.value);
      clearInterval(timer);
      if (isPlaying.value) timer = setInterval(() => go(index.value + 1), 2500);
    },
    openLink: (href) => {
      const i = entries.findIndex((e) => `#${e.title.toLowerCase().replace(/ /g, '-')}` === href);
      if (i >= 0) go(i);
      else status.post(`Link: ${href}`);
    },
  };
}

function mockTrajectory(): TrajectoryControlsModel {
  const count = 20;
  const index = observableState(0);
  const isPlaying = observableState(false);
  let timer: ReturnType<typeof setInterval> | undefined;
  const seek = (i: number) => index.set(Math.max(0, Math.min(count - 1, i)));
  return {
    count: observableState(count).asReadonly,
    index: index.asReadonly,
    isPlaying: isPlaying.asReadonly,
    first: () => seek(0),
    previous: () => seek(index.value - 1),
    next: () => seek(index.value + 1),
    last: () => seek(count - 1),
    seek,
    togglePlay: () => {
      isPlaying.set(!isPlaying.value);
      clearInterval(timer);
      if (isPlaying.value) timer = setInterval(() => index.set((index.value + 1) % count), 150);
    },
  };
}

/** A mock of everything around the canvas, with the shell's real model. */
function mockShell(options: { loaded: boolean }) {
  const shell = new ShellModel();
  const status = new StatusModel();
  const selectionMode = observableState(false);
  const illumination = observableState(false);
  const snapshots = observableState<readonly string[]>(['Overview', 'Active site']);
  const hoverLines = observableState<readonly string[]>([]);
  const hover: HoverLabelModel = { lines: hoverLines.asReadonly };

  shell.registerTabs([
    {
      id: 'home',
      title: 'Home',
      icon: 'house',
      order: 0,
      render: domView(() => (
        <div style={pad}>
          <b style="font-size: var(--msp-font-size-lg)">Open a structure</b>
          <div style="display: grid; gap: 8px">
            <Button icon={<Icon name="folder-open" />}>Open files…</Button>
            <Button icon={<Icon name="download" />}>Download by ID (PDB, AlphaFold DB, …)</Button>
          </div>
          <b>Examples</b>
          <div style="display: grid; gap: 4px">
            {['1TQN · Cytochrome P450 3A4', '7BV2 · SARS-CoV-2 RdRp', 'AF-Q8W3K0 · AlphaFold model'].map((e) => (
              <Button variant="ghost" onClick={() => status.post(`Loading ${e}…`)}>
                {e}
              </Button>
            ))}
          </div>
          <p style={muted}>
            Press <kbd>{formatChord({ key: 'k', mod: true })}</kbd> for every command.
          </p>
        </div>
      )),
    },
    {
      id: 'structure',
      title: 'Structure',
      icon: 'layers',
      order: 1,
      render: domView(() => (
        <Placeholder title="Scene">
          Structures, components and representations as one tree, with the inspector for the picked node below.
        </Placeholder>
      )),
    },
    {
      id: 'state-tree',
      title: 'State tree',
      icon: 'list-tree',
      order: 2,
      render: domView(() => <Placeholder title="State tree">The full state tree for advanced use.</Placeholder>),
    },
    {
      id: 'settings',
      title: 'Settings',
      icon: 'settings',
      order: 3,
      render: domView(() => (
        <Placeholder title="Settings">Searchable, grouped by intent: appearance, camera, interaction, …</Placeholder>
      )),
    },
  ]);
  if (options.loaded) {
    shell.registerTabs([
      {
        id: 'sequence',
        title: 'Sequence',
        icon: 'sequence',
        dock: 'top',
        render: domView(({ ctx }) => (
          <SequenceStrip
            chain={{ structure: '1TQN', name: 'Chain A', info: 'Cytochrome P450 3A4', sequence: Sequence, start: 1 }}
            onCollapse={() => ctx.shell.toggle(ctx.dock ?? 'top')}
            onHover={(n) =>
              hoverLines.set(n === undefined ? [] : [`${ThreeLetter[Sequence[n - 1]] ?? '?'} ${n}`, 'Chain A · 1TQN'])
            }
            onToggle={(n, on) => status.post(`${on ? 'Selected' : 'Deselected'} residue ${n}`)}
          />
        )),
      },
    ]);
    shell.show('structure');
  }

  shell.registerOverlays([
    {
      id: 'save-state',
      title: 'Save state',
      render: domView(() => {
        const name = observableState(`Snapshot ${snapshots.value.length + 1}`);
        return (
          <div style={pad}>
            <input
              value={name}
              onInput={(e) => name.set(e.currentTarget.value)}
              style="height: 28px; padding: 0 8px; border: 1px solid var(--msp-color-border); border-radius: var(--msp-radius-sm); background: var(--msp-color-bg)"
            />
            <div style="display: flex; gap: 8px; flex-wrap: wrap">
              <Button
                variant="primary"
                icon={<Icon name="plus" />}
                onClick={() => {
                  snapshots.set([...snapshots.value, name.value]);
                  status.post(`Saved snapshot “${name.value}”`, 'success');
                  name.set(`Snapshot ${snapshots.value.length + 1}`);
                }}
              >
                Snapshot
              </Button>
              <Button icon={<Icon name="download" />}>Session</Button>
              <Button icon={<Icon name="download" />}>State</Button>
            </div>
            <div style="display: grid; gap: 2px">
              <For each={snapshots}>
                {(s, i) => (
                  <div style="display: flex; align-items: center; gap: 8px; height: 28px">
                    <span style={muted}>{i.pipe(map((n) => `${n + 1}.`))}</span>
                    <span style="flex: 1">{s}</span>
                    <Button
                      variant="ghost"
                      size="sm"
                      icon={<Icon name="trash-2" />}
                      title="Remove"
                      onClick={() => snapshots.set(snapshots.value.filter((x) => x !== s))}
                    />
                  </div>
                )}
              </For>
            </div>
          </div>
        );
      }),
    },
    {
      id: 'screenshot',
      title: 'Screenshot',
      render: domView(() => <Placeholder title="Preview">Crop, size, copy and download.</Placeholder>),
    },
    {
      id: 'log',
      title: 'Log',
      render: domView(() => (
        <div style={pad}>
          <For each={status.history} fallback={() => <span style={muted}>No messages.</span>}>
            {(m) => (
              <div style="display: flex; gap: 8px">
                <span style={muted}>{m.time.toLocaleTimeString()}</span>
                <span>{m.text}</span>
              </div>
            )}
          </For>
        </div>
      )),
    },
  ]);

  shell.registerToolbar([
    {
      id: 'palette',
      icon: 'command',
      title: `Commands (${formatChord({ key: 'k', mod: true })})`,
      run: () => status.post('The palette comes next.'),
    },
    { id: 'reset', group: 'Canvas', icon: 'rotate-ccw', title: 'Reset camera', run: () => status.post('Camera reset') },
    { id: 'screenshot', group: 'Export', icon: 'camera', title: 'Screenshot', overlay: 'screenshot' },
    { id: 'save', group: 'Export', icon: 'save', title: 'Save state', overlay: 'save-state' },
    {
      id: 'select',
      group: 'Canvas',
      icon: 'mouse-pointer-click',
      title: 'Selection mode',
      pressed: selectionMode.asReadonly,
      run: () => selectionMode.set(!selectionMode.value),
    },
    {
      id: 'illumination',
      group: 'Canvas',
      icon: 'lightbulb',
      title: 'Illumination',
      pressed: illumination.asReadonly,
      run: () => illumination.set(!illumination.value),
    },
    {
      id: 'xr',
      group: 'Canvas',
      icon: 'glasses',
      title: 'Enter XR',
      run: () => status.post('XR is not available in this browser.', 'warning'),
    },
    { id: 'expand', group: 'Canvas', icon: 'maximize', title: 'Expand', run: () => status.post('Expanded') },
  ]);

  if (options.loaded) {
    status.post('Loaded 1TQN (3,624 atoms) in 0.4 s', 'success');
    status.setChips([
      { id: 'selection', label: '12 residues selected', activate: () => status.post('Focused on the selection') },
    ]);
    status.setTasks([{ id: 'surface', label: 'Computing surface', progress: 0.4, abort: () => status.setTasks([]) }]);
  }

  return { shell, status, hover, hoverLines, snapshots: mockSnapshots(status), trajectory: mockTrajectory() };
}

function ShellStory(props: { loaded: boolean; background?: string }) {
  const { shell, status, hover, hoverLines, snapshots, trajectory } = mockShell(props);
  shell.mount();
  status.mount();
  onCleanup(() => {
    shell.unmount();
    status.unmount();
  });
  return (
    <ShellView
      shell={shell}
      status={status}
      canvas={
        <CanvasPlaceholder
          empty={!props.loaded}
          background={props.background}
          onHover={(i) =>
            hoverLines.set(
              i === undefined ? [] : [`${ThreeLetter[Sequence[i * 4]] ?? 'GLY'} ${i * 4 + 1}`, 'Chain A · 1TQN'],
            )
          }
        />
      }
      overlays={
        props.loaded ? (
          <>
            <SnapshotPanel model={snapshots} />
            <HoverLabel model={hover} area={() => document.querySelector('main')} />
          </>
        ) : null
      }
      timeline={props.loaded ? <Timeline model={trajectory} /> : null}
      canvasTheme={themeForColor(parseColor(props.background ?? '#ffffff') ?? 0xffffff)}
      onOpenLog={() => shell.openOverlay('log')}
    />
  );
}

/** Nothing loaded: Home is the open tab. */
export const Empty: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={false} />}</div>),
};

/** A structure loaded: the Structure tab, a status message, a running task and the selection chip. */
export const Loaded: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} />}</div>),
};

/**
 * The chrome over the canvas follows the canvas background automatically (light glass on a light background, dark on
 * a dark one), whatever the UI theme. Pick a background color in the controls.
 */
export const CanvasBackground: StoryObj<{ background: string }> = {
  args: { background: '#14181f' },
  argTypes: { background: { control: 'color' } },
  render: (args) =>
    story(() => <div style="height: 100vh">{<ShellStory loaded={true} background={args.background} />}</div>),
};

/** A narrow host: the panel becomes a bottom sheet with its tabs along the bottom. */
export const Compact: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh; display: grid; place-items: center; background: var(--msp-color-border)">
        <div style="width: 390px; height: 760px; border-radius: 16px; overflow: hidden; box-shadow: var(--msp-shadow-overlay)">
          <ShellStory loaded={true} />
        </div>
      </div>
    )),
};
