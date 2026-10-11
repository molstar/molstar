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
import { Select } from '../../components/select/select.js';
import { TextInput } from '../../components/text-input/text-input.js';
import { For } from '../../dom/control-flow.js';
import type { Child } from '../../dom/jsx-types.js';
import { onCleanup } from '../../dom/scope.js';
import { domView } from '../../dom/view.js';
import { parseColor, themeForColor } from '../../model/canvas-theme.js';
import { observableState } from '../../reactive/observable-state.js';
import type { HoverLabelModel } from '../../model/playback.js';
import { ShellModel } from '../../shell/shell.js';
import { StatusModel } from '../../shell/status.js';
import { HoverLabel } from '../canvas/hover-label.js';
import {
  type CanvasMode,
  type Granularity,
  ModeSwitch,
  type ModeSwitchProps,
  PickLevel,
} from '../canvas/mode-switch.js';
import { PlaybackCard } from '../canvas/playback.js';
import { SelectionBar } from '../canvas/selection-bar.js';
import { SequenceStrip } from '../sequence/sequence-strip.js';
import { CanvasPlaceholder } from './canvas-placeholder.js';
import { LogView } from './log-view.js';
import { HomeTab } from '../home/home-tab.js';
import { SectionList } from './section-list.js';
import { ShellView } from './shell-view.js';
import { StructureTab } from '../structure/structure-tab.js';
import { mockHome } from '../../_stories/mocks/home.js';
import { mockSnapshots, mockTrajectory } from '../../_stories/mocks/playback.js';
import { registerMockAnimation } from '../../_stories/mocks/animation.js';
import { mockSelection } from '../../_stories/mocks/selection.js';
import { registerMockCanvasSettings } from '../../_stories/mocks/canvas-settings.js';
import { registerMockInspectorSections } from '../../_stories/mocks/inspector-sections.js';
import { type MockScenePreset, mockQuickStyles, mockScene } from '../../_stories/mocks/scene.js';

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

/**
 * The Export panel: one rail button, its sections contributed to the `export` slot. Built-in sections (image, session)
 * and extension sections (animation, geometry, models) look alike; each is also a command that opens the panel at it.
 */
function mockExport(shell: ShellModel, status: StatusModel) {
  const rendering = observableState<number | undefined>(undefined);
  const row = 'display: flex; gap: 8px; flex-wrap: wrap';
  shell.registerOverlays([
    {
      id: 'export',
      title: 'Export',
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="export" focus={ctx.focus} />),
    },
  ]);
  shell.registerSections([
    {
      id: 'image',
      slot: 'export',
      title: 'Image',
      icon: 'camera',
      order: 0,
      render: domView(() => (
        <div style={pad}>
          <div style="height: 120px; border-radius: var(--msp-radius-md); background: var(--msp-color-control)" />
          <div style={row}>
            <Button icon={<Icon name="copy" />} onClick={() => status.post('Copied image', 'success')}>
              Copy
            </Button>
            <Button
              variant="primary"
              icon={<Icon name="download" />}
              onClick={() => status.post('Saved image.png', 'success')}
            >
              Download PNG
            </Button>
          </div>
        </div>
      )),
    },
    {
      id: 'session',
      slot: 'export',
      title: 'Session & state',
      icon: 'save',
      order: 1,
      render: domView(() => (
        <div style={row}>
          <Button icon={<Icon name="download" />}>Session (.molx)</Button>
          <Button icon={<Icon name="download" />}>State (.molj)</Button>
        </div>
      )),
    },
    // Contributed by the mp4-export extension.
    {
      id: 'animation',
      slot: 'export',
      title: 'Animation (MP4)',
      icon: 'film',
      keywords: ['mp4', 'video', 'movie'],
      order: 10,
      render: domView(() => (
        <div style={pad}>
          <Select
            aria-label="Animation"
            options={[
              ['spin', 'Spin'],
              ['rock', 'Rock'],
              ['snapshots', 'Snapshots'],
            ]}
          />
          <Button
            variant="primary"
            icon={<Icon name="film" />}
            disabled={rendering.asReadonly.pipe(map((r) => r !== undefined))}
            onClick={() => {
              let p = 0;
              rendering.set(0);
              const t = setInterval(() => {
                p += 0.1;
                rendering.set(p >= 1 ? undefined : p);
                if (p >= 1) {
                  clearInterval(t);
                  status.post('Rendered animation.mp4 (6 s)', 'success');
                }
              }, 150);
            }}
          >
            {rendering.asReadonly.pipe(
              map((r) => (r === undefined ? 'Render MP4' : `Rendering ${Math.round(r * 100)}%`)),
            )}
          </Button>
        </div>
      )),
    },
    // Contributed by the geo-export extension.
    {
      id: 'geometry',
      slot: 'export',
      title: 'Geometry',
      icon: 'box',
      keywords: ['glb', 'stl', 'obj', 'usdz', 'mesh', '3d print', 'ar'],
      order: 11,
      render: domView(() => (
        <div style={row}>
          {['GLB', 'STL', 'OBJ', 'USDZ'].map((f) => (
            <Button
              icon={<Icon name="download" />}
              onClick={() => status.post(`Exported scene.${f.toLowerCase()}`, 'success')}
            >
              {f}
            </Button>
          ))}
        </div>
      )),
    },
    // Contributed by the model-export extension.
    {
      id: 'models',
      slot: 'export',
      title: 'Models',
      icon: 'file',
      keywords: ['mmcif', 'bcif', 'cif', 'coordinates'],
      order: 12,
      render: domView(() => (
        <div style={row}>
          <Button icon={<Icon name="download" />}>mmCIF</Button>
          <Button icon={<Icon name="download" />}>BinaryCIF</Button>
        </div>
      )),
    },
  ]);
}

// Mol*'s granularities (interactivity), labelled and drawn as picking levels.
const Granularities: readonly Granularity[] = [
  { value: 'element', label: 'Atom', icon: 'pick-atom' },
  { value: 'residue', label: 'Residue', icon: 'pick-residue' },
  { value: 'chain', label: 'Chain', icon: 'pick-chain' },
  { value: 'entity', label: 'Entity', icon: 'pick-entity' },
  { value: 'model', label: 'Model', icon: 'pick-model' },
  { value: 'operator', label: 'Operator', icon: 'copy' },
  { value: 'structure', label: 'Structure', icon: 'atom' },
  { value: 'elementInstances', label: 'Atom', icon: 'pick-atom', group: 'All instances' },
  { value: 'residueInstances', label: 'Residue', icon: 'pick-residue', group: 'All instances' },
  { value: 'chainInstances', label: 'Chain', icon: 'pick-chain', group: 'All instances' },
];

/** A mock of everything around the canvas, with the shell's real model. */
function mockShell(options: { loaded: boolean; preset?: MockScenePreset }) {
  const shell = new ShellModel();
  const status = new StatusModel();
  const canvasMode = observableState<CanvasMode>('interact');
  const granularity = observableState('residue');
  const illumination = observableState(false);
  const snapshots = observableState<readonly string[]>(['Overview', 'Active site']);
  const hoverLines = observableState<readonly string[]>([]);
  // The scene; hovering a node highlights it (the hover label stands in for the 3D highlight here).
  const scene = mockScene({
    preset: options.preset,
    onHighlight: (label) => hoverLines.set(label ? [label] : []),
    onAction: (text) => status.post(text),
  });
  // Extensions' inspector sections (volume streaming, assembly symmetry), contributed to the node kind's slot.
  registerMockInspectorSections(shell, scene);
  const hover: HoverLabelModel = { lines: hoverLines.asReadonly };
  const home = mockHome((text, kind) => status.post(text, kind));

  shell.registerTabs([
    {
      id: 'home',
      title: 'Home',
      icon: 'molstar',
      order: 0,
      place: 'head',
      render: domView(() => <HomeTab home={home} onOpenPalette={() => status.post('The palette comes next.')} />),
    },
    {
      id: 'structure',
      title: 'Structure',
      icon: 'layers',
      order: 1,
      render: domView(({ ctx }) => (
        <StructureTab
          shell={ctx.shell}
          scene={scene}
          quickStyles={mockQuickStyles((t) => status.post(t))}
          onAction={(a) =>
            status.post(
              `${a === 'preset' ? 'Apply preset' : a === 'color' ? 'Color by' : 'Add component'}… opens in the palette`,
            )
          }
        />
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
      id: 'log',
      title: 'Log',
      icon: 'scroll-text',
      order: 2.5,
      render: domView(() => <LogView status={status} />),
    },
    {
      id: 'settings',
      title: 'Settings',
      icon: 'settings',
      order: 3,
      place: 'foot',
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
            <TextInput aria-label="Snapshot name" value={name} onInput={(v) => name.set(v)} />
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
  ]);

  shell.registerToolbar([
    {
      id: 'palette',
      icon: 'command',
      title: `Commands (${formatChord({ key: 'k', mod: true })})`,
      run: () => status.post('The palette comes next.'),
    },
    {
      id: 'reset',
      place: 'canvas',
      group: 'canvas',
      icon: 'focus',
      title: 'Reset camera',
      run: () => status.post('Camera reset'),
      menu: [
        {
          id: 'reset-zoom',
          label: 'Reset zoom',
          icon: 'zoom-out',
          hint: 'Fit everything, keep the orientation',
          run: () => status.post('Zoom reset'),
        },
        {
          id: 'orient-axes',
          label: 'Orient axes',
          icon: 'axis-3d',
          hint: 'Principal axes along the screen',
          run: () => status.post('Oriented along the principal axes'),
        },
        {
          id: 'reset-axes',
          label: 'Reset axes',
          icon: 'rotate-ccw',
          hint: 'Back to the default orientation',
          run: () => status.post('Axes reset'),
        },
      ],
    },
    // The scene as a whole: how it plays, how the canvas renders it.
    { id: 'animation', group: 'scene', icon: 'clapperboard', title: 'Animation', overlay: 'animation' },
    {
      id: 'canvas-settings',
      group: 'scene',
      icon: 'aperture',
      title: 'Canvas settings',
      overlay: 'canvas-settings',
    },
    { id: 'save', group: 'export', icon: 'save', title: 'Save state (snapshots)', overlay: 'save-state' },
    { id: 'export', group: 'export', icon: 'download', title: 'Export…', overlay: 'export' },
    // Rendering modes under one button: a click opens the menu (no action of its own); shown pressed while one is on.
    {
      id: 'rendering',
      place: 'canvas',
      group: 'canvas',
      icon: 'sparkles',
      title: 'Illumination & XR',
      pressed: illumination.asReadonly,
      menu: [
        {
          id: 'illumination',
          label: 'Illumination',
          icon: 'lightbulb',
          hint: 'Path-traced global illumination',
          pressed: illumination.asReadonly,
          run: () => illumination.set(!illumination.value),
        },
        {
          id: 'xr',
          label: 'Enter XR',
          icon: 'glasses',
          hint: 'View on a headset',
          run: () => status.post('XR is not available in this browser.', 'warning'),
        },
      ],
    },
    {
      id: 'expand',
      place: 'canvas',
      group: 'canvas',
      icon: 'maximize',
      title: 'Expand',
      run: () => status.post('Expanded'),
    },
    // Panels over the canvas: brought back from the stack when hidden (hidden, they take no room at all).
    {
      id: 'sequence',
      place: 'canvas',
      group: 'view',
      icon: 'sequence',
      title: 'Sequence',
      // Only while the card is hidden: shown, it carries its own collapse button.
      visible: shell.dock('top').active.pipe(map((t) => !t)),
      run: () => shell.toggle('top'),
    },
  ]);

  // One Export button; built-in and extension exports are sections of its panel.
  mockExport(shell, status);
  registerMockCanvasSettings(shell);
  registerMockAnimation(shell, (text) => status.post(text));
  const selection = mockSelection(shell, (text) => status.post(text));

  if (options.loaded) {
    status.post('Loaded 1TQN (3,624 atoms) in 0.4 s', 'success');
    status.setTasks([{ id: 'surface', label: 'Computing surface', progress: 0.4, abort: () => status.setTasks([]) }]);
  }

  // The canvas mode lives on the canvas, not in the rail: what a click does is visible before the first click.
  const modeProps: ModeSwitchProps = {
    mode: canvasMode.asReadonly,
    setMode: (m) => {
      if (m === canvasMode.value) return;
      canvasMode.set(m);
      status.post(m === 'select' ? 'Select: click to select, Shift+click to extend' : 'Interact: click to focus');
    },
    granularity: granularity.asReadonly,
    granularities: Granularities,
    setGranularity: (g) => granularity.set(g),
  };
  const modeSwitch = <ModeSwitch {...modeProps} />;
  const pickLevel = <PickLevel {...modeProps} />;

  return {
    shell,
    status,
    hover,
    hoverLines,
    modeSwitch,
    pickLevel,
    selection,
    snapshots: mockSnapshots((t) => status.post(t)),
    trajectory: mockTrajectory(),
  };
}

function ShellStory(props: {
  loaded: boolean;
  preset?: MockScenePreset;
  background?: string;
  openExportAt?: string;
  /** What the loaded data has to play: a session's snapshots, a trajectory's models (plain 1TQN has neither). */
  playback?: { snapshots?: boolean; trajectory?: boolean };
}) {
  const { shell, status, hover, hoverLines, modeSwitch, pickLevel, selection, snapshots, trajectory } =
    mockShell(props);
  shell.mount();
  status.mount();
  if (props.openExportAt) shell.openOverlay('export', undefined, props.openExportAt);
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
            <HoverLabel model={hover} area={() => document.querySelector('main')} />
          </>
        ) : null
      }
      playback={
        props.loaded && (props.playback?.snapshots || props.playback?.trajectory) ? (
          <PlaybackCard
            snapshots={props.playback.snapshots ? snapshots : undefined}
            trajectory={props.playback.trajectory ? trajectory : undefined}
          />
        ) : null
      }
      canvasColor={props.background ?? '#ffffff'}
      canvasTheme={themeForColor(parseColor(props.background ?? '#ffffff') ?? 0xffffff)}
      context={props.loaded ? <SelectionBar selection={selection} /> : null}
      canvasLead={props.loaded ? modeSwitch : null}
      canvasAside={props.loaded ? pickLevel : null}
      onOpenLog={() => shell.show('log')}
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

/** Several structures: four CYP3A4 structures superposed on 1TQN, collapsed but the first. */
/**
 * A session with snapshots and a multi-model trajectory: playback appears at the top left only because the data has
 * something to play (plain 1TQN has neither).
 */
export const SessionWithPlayback: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory loaded={true} playback={{ snapshots: true, trajectory: true }} />
      </div>
    )),
};

export const ManyStructures: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} preset="many" />}</div>),
};

/** An EMDB entry: the map is a top-level volume next to the fitted model (EMD-2984 and 5A1A). */
export const EmdbMap: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ShellStory loaded={true} preset="em" />}</div>),
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

/** A command opened the Export panel at a section ("Export animation (MP4)…" in the palette). */
export const ExportAtSection: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory loaded={true} openExportAt="animation" />
      </div>
    )),
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
