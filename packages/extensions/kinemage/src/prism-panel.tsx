/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, type ContextMenuItem, domView, For, Icon, openContextMenu, ScrollArea } from '@molstar/prism-ui';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { SceneRow } from '@molstar/prism-viewer/views/structure/scene-row';
import type { KinemageApi, KinemageInfo } from './operations.js';
import {
  animateMenu,
  framesDescription,
  groupDescription,
  type KinemagePanelGroup,
  type KinemagePanelLayout,
  type KinemageItem,
  kinemagePanelLayout,
  kinemageDescription,
  kinemageShown,
  mastersDescription,
  sameLayout,
  sameShown,
  shownKey,
} from './prism-panel-model.js';

// The Kinemage panel (parity.md §6; the v5 right panel's controls, anew), a tab on the rail beside Structure while a
// kinemage is loaded, drawn as the Structure panel is: Next frame as its header's button (in one place, whatever the
// frame), then a tree of the scene tree's own rows (`SceneRow`): each kinemage, its views (a click applies one), its
// animate frames (the groups shown one at a time; the row says which), its masters and its other groups with their
// subgroups, each with an eye.
// Every control calls the area's operations as the person, so a click is a step in the history and Undo takes it back.
// It opens by itself when a kinemage loads (the v5 panel showed with the file), and by `kinemage.show-panel`.

export const KinemageTab = 'kinemage.tab.kinemage';

export function kinemageUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api.kinemage;
  if (!area) return {};
  // Asked again at each change of the state; the layout and the eyes are told on only when they differ (a visibility
  // step changes `shown` alone, a load the layout).
  const list = derive(area.changed, () => area.list());
  const layout = derive(list, kinemagePanelLayout, sameLayout);
  const shown = derive(list, kinemageShown, sameShown);
  const loaded = derive(list, (l) => l.length > 0);

  // Opened by itself when a kinemage loads with none before: the controls come with the file, as v5's panel did. After
  // the tab is on the rail (its `visible` is told the same change, after this).
  loaded.listen?.((has) => {
    if (has) queueMicrotask(() => void ctx.api.view.open({ tab: KinemageTab }));
  });

  return {
    tabs: [
      {
        id: KinemageTab,
        title: 'Kinemage',
        icon: 'kinemage',
        // Beside Structure (1), before Selection (1.2): what is loaded, then what to do with it.
        order: 1.1,
        visible: loaded,
        render: domView(() => (
          <div style="display: flex; flex-direction: column; flex: 1 1 0; height: 100%; min-height: 0; gap: var(--msp-space-2); padding: 0 var(--msp-panel-inset) var(--msp-panel-inset); font-size: var(--msp-font-size-sm)">
            <KinemageActions area={area} layout={layout} />
            <ScrollArea role="tree" aria-label="Kinemages" edge>
              <For
                each={layout}
                by={(k) => k.ref}
                fallback={() => (
                  <p style="margin: 0; padding: var(--msp-space-2); color: var(--msp-color-text-muted)">
                    No kinemage is loaded. Open a .kin file.
                  </p>
                )}
              >
                {(k) => <KinemageTree area={area} kinemage={k} list={list} shown={shown} />}
              </For>
            </ScrollArea>
          </div>
        )),
      },
    ],
  };
}

/** The header's button, as the Structure panel's Preset and Styles: the next animate frame. */
function KinemageActions(props: { area: KinemageApi; layout: ReadonlyState<readonly KinemagePanelLayout[]> }) {
  const { area, layout } = props;
  const frames = derive(layout, animateMenu);
  const step = (f: (typeof frames.value)[number]) => void area.animate({ sequence: f.sequence, kinemage: f.kinemage });
  return (
    <div style="display: flex; gap: var(--msp-space-1); min-width: 0">
      <Button
        size="sm"
        variant="ghost"
        icon={<Icon name="step-forward" />}
        disabled={derive(frames, (f) => f.length === 0)}
        onClick={(e) => {
          // One sequence: a click steps it; several (2animate, more kinemages) ask which.
          if (frames.value.length === 1) return step(frames.value[0]);
          const button = e.currentTarget as HTMLElement;
          const r = button.getBoundingClientRect();
          const items: ContextMenuItem[] = frames.value.map((f) => ({ id: f.id, label: f.label, run: () => step(f) }));
          openContextMenu(
            button.closest<HTMLElement>('.msp-root') ?? document.body,
            { x: r.left, y: r.bottom + 4 },
            items,
          );
        }}
      >
        Next frame
      </Button>
    </div>
  );
}

/** A kinemage as the scene tree draws a structure: its row, then its views, frames, masters and groups a level deeper. */
function KinemageTree(props: {
  area: KinemageApi;
  kinemage: KinemagePanelLayout;
  list: ReadonlyState<readonly KinemageInfo[]>;
  shown: ReadonlyState<ReadonlyMap<string, boolean>>;
}) {
  const { area, kinemage: k, shown } = props;
  const ref = k.ref;
  const info = derive(props.list, (l) => l.find((x) => x.ref === ref));
  const isShown = (item: Pick<KinemageItem, 'kind' | 'key'>) => {
    const key = shownKey(ref, item);
    return derive(shown, (m) => !!m.get(key));
  };
  const toggle = (item: KinemageItem, now: ReadonlyState<boolean>) =>
    void area.visibility({ kind: item.kind, name: item.key, visible: !now.value, kinemage: ref });
  const leaf = (item: KinemageItem, depth: number) => {
    const now = isShown(item);
    return <SceneRow label={item.label} depth={depth} shown={now} onToggle={() => toggle(item, now)} />;
  };
  const group = (g: KinemagePanelGroup, depth: number) => {
    const now = isShown(g);
    return (
      <SceneRow
        label={g.label}
        depth={depth}
        description={derive(info, () => groupDescription(g))}
        shown={now}
        onToggle={() => toggle(g, now)}
        {...(g.subgroups.length ? { children: () => <>{g.subgroups.map((s) => leaf(s, depth + 1))}</> } : {})}
      />
    );
  };
  // An animate sequence: its groups under one row, which says the frame shown (the header's Next frame steps it).
  const frames = (sequence: 'animate' | '2animate') => {
    const groups = k.groups.filter((g) => g.animate === (sequence === 'animate' ? 1 : 2));
    if (!groups.length) return null;
    return (
      <SceneRow
        label={sequence === 'animate' ? 'Frames' : 'Frames (2animate)'}
        depth={1}
        expanded
        description={derive(info, (i) => framesDescription(i, sequence))}
        children={() => <>{groups.map((g) => group(g, 2))}</>}
      />
    );
  };
  return (
    <SceneRow
      label={k.title}
      icon="kinemage"
      depth={0}
      expanded
      description={derive(info, () => kinemageDescription(k))}
      children={() => (
        <>
          {k.views.length ? (
            <SceneRow
              label="Views"
              depth={1}
              expanded
              description={derive(info, () => `${k.views.length} saved`)}
              children={() => (
                <>
                  {k.views.map((view) => (
                    // A click applies it (the camera, the background, as the kinemage says).
                    <SceneRow
                      label={view.name}
                      depth={2}
                      onActivate={() => void area.view({ view: view.key, kinemage: ref })}
                    />
                  ))}
                </>
              )}
            />
          ) : null}
          {frames('animate')}
          {frames('2animate')}
          {k.masters.length ? (
            <SceneRow
              label="Masters"
              depth={1}
              expanded
              description={derive(info, (i) => mastersDescription(i))}
              children={() => <>{k.masters.map((m) => leaf(m, 2))}</>}
            />
          ) : null}
          {k.groups.filter((g) => !g.animate).map((g) => group(g, 1))}
        </>
      )}
    />
  );
}
