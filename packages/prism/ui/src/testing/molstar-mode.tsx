/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { combineLatest, map, of } from 'rxjs';
import { ActionFooter } from '../components/action-footer/action-footer.js';
import { MolstarWordmark } from '../components/icon/wordmark.js';
import { ParamsForm } from '../components/params/params-form.js';
import { domView } from '../dom/view.js';
import { ParamsFormModel } from '../forms/params-form.js';
import type { HoverLabelModel } from '../model/playback.js';
import { observableState } from '../reactive/observable-state.js';
import type { ShellModel } from '../shell/shell.js';
import type { StatusModel } from '../shell/status.js';
import { HoverLabel } from '../views/canvas/hover-label.js';
import {
  type CanvasMode,
  type Granularity,
  ModeSwitch,
  type ModeSwitchProps,
  PickLevel,
} from '../views/canvas/mode-switch.js';
import { SnapshotCard, TrajectoryBar } from '../views/canvas/playback.js';
import { SelectionBar } from '../views/canvas/selection-bar.js';
import { ImageExport } from '../views/export/image-export.js';
import { HomeTab } from '../views/home/home-tab.js';
import { type SequenceSelector, SequenceStrip } from '../views/sequence/sequence-strip.js';
import { CanvasPlaceholder } from '../views/shell/canvas-placeholder.js';
import { SectionList } from '../views/shell/section-list.js';
import { RemoteStatesPanel } from '../views/state/remote-states-panel.js';
import { SnapshotsPanel } from '../views/state/snapshots-panel.js';
import { StateTreeTab } from '../views/state-tree/state-tree-tab.js';
import { StructureTab } from '../views/structure/structure-tab.js';
import { registerMockAnimation } from './animation.js';
import { registerMockCanvasSettings } from './canvas-settings.js';
import { mockCommands, mockPrefixes } from './commands.js';
import { mockHome } from './home.js';
import { mockImageExport } from './image-export.js';
import { mockMeasurements, registerMockMeasurements } from './measurements.js';
import { mockSnapshots, mockTrajectory } from './playback.js';
import { type MockScenePreset, mockComponentOptions, mockPresets, mockQuickStyles, mockScene } from './scene.js';
import { mockSelection } from './selection.js';
import type { MockMode } from './shell.js';
import { mockSnapshotsModel } from './snapshots.js';
import { mockStateTree } from './state-tree.js';
import { mockSuperposition, registerMockSuperposition } from './superposition.js';

// The Mol* mode of the mock shell: the plugin's own UI (Home, Structure, Snapshots, State tree, the sequence, the
// canvas chrome, Mol*'s commands) over the mocks, as one mode the frame activates and deactivates (shell.tsx). A story
// mode (@molstar/prism-stories) is another.

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
  const imageExport = mockImageExport((t, k) => status.post(t, k));
  const rendering = observableState<number | undefined>(undefined);
  // Each section: choose in a form, then one action in the footer (Mol*'s action footer).
  const form = 'display: grid; gap: 8px';
  const sessionForm = new ParamsFormModel({
    kind: PD.Select('session', [
      ['session', 'Session (.molx) · with data'],
      ['state', 'State (.molj) · description only'],
    ]),
  });
  const mp4Form = new ParamsFormModel({
    animation: PD.Select('spin', [
      ['spin', 'Camera spin'],
      ['rock', 'Camera rock'],
      ['snapshots', 'State snapshots'],
    ]),
    resolution: PD.Select('full-hd', [
      ['hd', 'HD (1280 × 720)'],
      ['full-hd', 'Full HD (1920 × 1080)'],
      ['ultra-hd', 'Ultra HD (3840 × 2160)'],
    ]),
    transparent: PD.Boolean(false),
  });
  const geometryForm = new ParamsFormModel({
    format: PD.Select('glb', [
      ['glb', 'GLB (glTF binary)'],
      ['stl', 'STL'],
      ['obj', 'OBJ'],
      ['usdz', 'USDZ (AR)'],
    ]),
  });
  const modelsForm = new ParamsFormModel({
    format: PD.Select('cif', [
      ['cif', 'mmCIF'],
      ['bcif', 'BinaryCIF'],
    ]),
  });
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
      render: domView(() => <ImageExport model={imageExport} />),
    },
    {
      id: 'session',
      slot: 'export',
      title: 'Session & state',
      icon: 'save',
      order: 1,
      render: domView(() => (
        <div style={form}>
          <ParamsForm form={sessionForm} />
          <ActionFooter
            icon="download"
            label="Download"
            onApply={() =>
              status.post(
                (sessionForm.values.value as any).kind === 'session' ? 'Saved state.molx' : 'Saved state.molj',
                'success',
              )
            }
          />
        </div>
      )),
    },
    // Contributed by the mp4-export extension.
    {
      id: 'animation',
      slot: 'export',
      // Named for what it makes, not the format: more formats may follow.
      title: 'Animation',
      icon: 'film',
      keywords: ['mp4', 'video', 'movie'],
      order: 10,
      render: domView(() => (
        <div style={form}>
          <ParamsForm form={mp4Form} />
          <ActionFooter
            icon="film"
            label={rendering.asReadonly.pipe(
              map((r) => (r === undefined ? 'Render' : `Rendering ${Math.round(r * 100)}%`)),
            )}
            busy={rendering.asReadonly.pipe(map((r) => r !== undefined))}
            onReset={() => mp4Form.reset()}
            onApply={() => {
              // Rendering holds the plugin: an overlay task, with Cancel.
              let p = 0;
              rendering.set(0);
              const done = (message: string) => {
                clearInterval(t);
                rendering.set(undefined);
                status.setTasks(status.tasks.value.filter((x) => x.id !== 'mp4'));
                status.post(message, 'success');
              };
              const task = () => ({
                id: 'mp4',
                label: 'Rendering animation.mp4',
                progress: p,
                overlay: true,
                abort: () => done('Rendering cancelled'),
              });
              const put = () => status.setTasks([...status.tasks.value.filter((x) => x.id !== 'mp4'), task()]);
              put();
              const t = setInterval(() => {
                p = Math.min(1, p + 0.04);
                rendering.set(p);
                put();
                if (p >= 1) done('Rendered animation.mp4 (6 s)');
              }, 120);
            }}
          />
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
        <div style={form}>
          <ParamsForm form={geometryForm} />
          <ActionFooter
            icon="download"
            label="Download"
            onApply={() => status.post(`Exported scene.${(geometryForm.values.value as any).format}`, 'success')}
          />
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
        <div style={form}>
          <ParamsForm form={modelsForm} />
          <ActionFooter
            icon="download"
            label="Download"
            onApply={() => status.post(`Saved 1tqn.${(modelsForm.values.value as any).format}`, 'success')}
          />
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
  { value: 'structure', label: 'Structure', icon: 'structure' },
  { value: 'elementInstances', label: 'Atom', icon: 'pick-atom', group: 'All instances' },
  { value: 'residueInstances', label: 'Residue', icon: 'pick-residue', group: 'All instances' },
  { value: 'chainInstances', label: 'Chain', icon: 'pick-chain', group: 'All instances' },
];

/** The sequence view's choices for 1TQN (Mol*'s sequence controls): structure, mode, then entity and chain. */
function sequenceSelectors(post: (text: string) => void): SequenceSelector[] {
  const selector = (
    id: string,
    title: string,
    initial: string,
    options: [string, string][],
    extra: Partial<SequenceSelector> = {},
  ): SequenceSelector => {
    const value = observableState(initial);
    return {
      id,
      title,
      options,
      value: value.asReadonly,
      set: (v) => {
        value.set(v);
        post(`Sequence: ${title.toLowerCase()} ${options.find((o) => o[0] === v)?.[1] ?? v}`);
      },
      ...extra,
    };
  };
  const mode = selector('mode', 'Mode', 'single', [
    ['single', 'Chain'],
    ['polymers', 'Polymers'],
    ['all', 'Everything'],
  ]);
  const single = mode.value.pipe(map((m) => m === 'single'));
  return [
    selector('structure', 'Structure', '1tqn', [['1tqn', '1TQN | Crystal structure of human microsomal P450 3A4']], {
      display: () => '1TQN',
      strong: true,
    }),
    mode,
    selector(
      'entity',
      'Entity',
      '1',
      [
        ['1', '1: Cytochrome P450 3A4'],
        ['2', '2: Protoporphyrin IX containing Fe'],
        ['3', '3: Sulfate ion'],
        ['4', '4: Water'],
      ],
      {
        display: (v) =>
          ['Cytochrome P450 3A4', 'Protoporphyrin IX containing Fe', 'Sulfate ion', 'Water'][Number(v) - 1],
        visible: single,
      },
    ),
    selector('chain', 'Chain', 'A', [['A', 'A']], { display: (v) => `Chain ${v}`, visible: single }),
  ];
}

/** Mol*'s own UI as a mode of the mock shell. */
export function mockMolstarMode(
  options: {
    preset?: MockScenePreset;
    /** The spec turns remote states (a state server) on; a number fills the server with that many. */
    remoteStates?: boolean | number;
    /** What the loaded data has to play: a session's snapshots, a trajectory's models (plain 1TQN has neither). */
    playback?: { snapshots?: boolean; trajectory?: boolean };
  } = {},
): MockMode {
  return {
    id: 'molstar',
    title: 'Mol* mode',
    icon: 'molstar',
    activate: ({ shell, status, loaded, background, openCommandMenu }) => {
      const canvasMode = observableState<CanvasMode>('interact');
      const granularity = observableState('residue');
      const illumination = observableState(false);
      const hoverLines = observableState<readonly string[]>([]);
      // The scene; hovering a node highlights it (the hover label stands in for the 3D highlight here).
      const scene = mockScene({
        preset: options.preset,
        onHighlight: (label) => hoverLines.set(label ? [label] : []),
        onAction: (text) => status.post(text),
        onMeasure: () => shell.show('selection', 'measure'),
      });
      // Superposing, with two or more structures: the Selection tab's Superpose section.
      registerMockSuperposition(
        shell,
        mockSuperposition(scene, (t) => status.post(t)),
      );
      // Adding measurements: the Selection tab's Measure section.
      registerMockMeasurements(
        shell,
        mockMeasurements(scene, (t) => status.post(t)),
      );
      const hover: HoverLabelModel = { lines: hoverLines.asReadonly };
      const home = mockHome((text, kind) => status.post(text, kind), { remoteStates: options.remoteStates });
      const snapshotsModel = mockSnapshotsModel((text, kind) => status.post(text, kind));
      const stateTree = mockStateTree((text) => status.post(text));
      const snapshotStepper = mockSnapshots((t) => status.post(t));
      const trajectory = mockTrajectory();
      const commands = mockCommands(shell, () => scene.roots.value.filter((n) => n.kind === 'structure').length);

      shell.registerTabs([
        {
          id: 'home',
          compactBar: true,
          title: 'Home',
          // The Mol* wordmark over Home (the rail's tooltip stays "Home").
          heading: () => <MolstarWordmark live />,
          icon: 'prism',
          order: 0,
          place: 'head',
          render: domView(() => <HomeTab home={home} onOpenCommandMenu={() => openCommandMenu()} />),
        },
        {
          id: 'structure',
          compactBar: true,
          title: 'Structure',
          icon: 'structure',
          order: 1,
          render: domView(({ ctx }) => (
            <StructureTab
              shell={ctx.shell}
              scene={scene}
              presets={mockPresets((t) => status.post(t))}
              options={mockComponentOptions()}
              quickStyles={mockQuickStyles((t) => status.post(t))}
              openCommands={(group) => openCommandMenu(group)}
              onAction={(a) =>
                status.post(
                  `${a === 'preset' ? 'Apply preset' : a === 'styles' ? 'Styles' : 'Add component'}… opens in the command menu`,
                )
              }
            />
          )),
        },
        {
          id: 'snapshots',
          title: 'Snapshots',
          icon: 'images',
          order: 1.5,
          render: domView(() => <SnapshotsPanel model={snapshotsModel} />),
        },
        {
          id: 'state-tree',
          title: 'State tree',
          icon: 'list-tree',
          order: 2,
          render: domView(() => <StateTreeTab model={stateTree} />),
        },
      ]);
      if (loaded) {
        shell.registerTabs([
          {
            id: 'sequence',
            title: 'Sequence',
            icon: 'sequence',
            dock: 'top',
            render: domView(({ ctx }) => (
              <SequenceStrip
                chain={{
                  structure: '1TQN',
                  name: 'Chain A',
                  info: 'Cytochrome P450 3A4',
                  sequence: Sequence,
                  start: 1,
                }}
                selectors={sequenceSelectors((t) => status.post(t))}
                onCollapse={() => ctx.shell.toggle(ctx.dock ?? 'top')}
                onHover={(n) =>
                  hoverLines.set(
                    n === undefined ? [] : [`${ThreeLetter[Sequence[n - 1]] ?? '?'} ${n}`, 'Chain A · 1TQN'],
                  )
                }
                onToggle={(n, on) => status.post(`${on ? 'Selected' : 'Deselected'} residue ${n}`)}
              />
            )),
          },
        ]);
      }

      if (home.remoteStates) {
        const remote = home.remoteStates;
        shell.registerOverlays([
          {
            id: 'remote-states',
            title: 'Remote states',
            render: domView(() => <RemoteStatesPanel model={remote} />),
          },
        ]);
      }

      // The canvas stack's camera and rendering tools act on what is shown: nothing loaded, nothing to show them for.
      const hasData = of(loaded);
      shell.registerToolbar([
        {
          id: 'reset',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
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
          icon: 'monitor-cog',
          title: 'Canvas settings',
          overlay: 'canvas-settings',
        },
        { id: 'export', group: 'export', icon: 'download', title: 'Export…', overlay: 'export' },
        ...(home.remoteStates
          ? [
              {
                id: 'remote-states',
                group: 'export',
                icon: 'cloud-upload' as const,
                title: 'Remote states',
                overlay: 'remote-states',
              },
            ]
          : []),
        // Rendering modes under one button: a click opens the menu (no action of its own); shown pressed while one is on.
        {
          id: 'rendering',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'lightbulb',
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
        // Panels over the canvas: brought back from the stack when hidden (hidden, they take no room at all).
        {
          id: 'sequence',
          place: 'canvas',
          group: 'view',
          icon: 'sequence',
          title: 'Sequence',
          // Only while there is a sequence and its card is hidden (shown, it carries its own collapse button).
          visible: combineLatest([shell.dock('top').tabs, shell.dock('top').active]).pipe(
            map(([tabs, active]) => tabs.length > 0 && !active),
          ),
          run: () => shell.toggle('top'),
        },
      ]);

      // One Export button; built-in and extension exports are sections of its panel.
      mockExport(shell, status);
      registerMockCanvasSettings(shell);
      registerMockAnimation(shell, (text) => status.post(text));
      const selection = mockSelection(shell, (text) => status.post(text));

      if (loaded) {
        status.post('Loaded 1TQN (3,624 atoms) in 0.4 s', 'success');
        // Two tasks at once: they stack.
        const drop = (id: string) => () => status.setTasks(status.tasks.value.filter((t) => t.id !== id));
        status.setTasks([
          { id: 'density', label: 'Downloading density', progress: 0.7, abort: drop('density') },
          { id: 'surface', label: 'Computing surface', progress: 0.4, abort: drop('surface') },
        ]);
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
      const playback = options.playback;

      return {
        start: loaded ? 'structure' : 'home',
        commands,
        prefixes: mockPrefixes({ manager: commands, shell, scene, snapshots: snapshotStepper }),
        canvas: () => (
          <CanvasPlaceholder
            empty={!loaded}
            background={background}
            onHover={(i) =>
              hoverLines.set(
                i === undefined ? [] : [`${ThreeLetter[Sequence[i * 4]] ?? 'GLY'} ${i * 4 + 1}`, 'Chain A · 1TQN'],
              )
            }
          />
        ),
        overlays: () =>
          loaded ? (
            <HoverLabel model={hover} area={() => document.querySelector<HTMLElement>('[data-msp-canvas]')} />
          ) : null,
        playback: () =>
          loaded && (playback?.snapshots || playback?.trajectory) ? (
            <>
              {playback.snapshots ? <SnapshotCard model={snapshotStepper} /> : null}
              {playback.trajectory ? <TrajectoryBar model={trajectory} /> : null}
            </>
          ) : null,
        context: loaded ? () => <SelectionBar selection={selection} /> : undefined,
        canvasLead: () => (loaded ? <ModeSwitch {...modeProps} /> : null),
        canvasAside: () => (loaded ? <PickLevel {...modeProps} /> : null),
        empty: of(!loaded),
        onDropFiles: (files) => home.openFiles(files),
        dropFormats: home.formats,
      };
    },
  };
}
