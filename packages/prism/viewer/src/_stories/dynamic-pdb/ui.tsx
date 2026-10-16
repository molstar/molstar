/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, For, Icon, Show } from '@molstar/prism-ui';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import type { SceneNode } from '../../model/scene.js';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { SceneRow } from '../../views/structure/scene-row.js';
import { classes } from './_generated/dynamic-pdb.module.css.js';
import { type DpdbEntry, type DpdbMutation, dpdbPage, MetricScales, mutationTarget } from './data.js';
import { type DpdbLoaded, type DpdbLoadedModel, DpdbCard, type DynamicPdbApi } from './extension.js';

// The Dynamic PDB card (as the PAE card, model-archive's `pae/ui.tsx`): the entry (its title, method, crystal, the
// construct and its mutations) and its models, each with its refinement metrics on bars (lower is better) and its eye;
// opened from the toolbar while an entry is loaded, or by `dynamic-pdb.show`. In the Structure panel, a row under each
// of the entry's models opens it. The deposited model's density streams from the volume server by the button on its
// model (`volume.streaming`; the viewer's own Volume streaming row has its controls); the sequence's tracks show and hide
// by their eyes (`dynamic-pdb.track`).

const Quiet = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs)';
const Inset = 'calc(var(--msp-panel-inset) + var(--msp-inset))';

const distinct = <T,>(o: Observable<T>) => o.pipe(distinctUntilChanged());

/** "X-ray crystallography · 1.95 Å · P 32 2 1". */
function facts(e: DpdbEntry): string {
  return [e.method, e.resolution ? `${e.resolution} Å` : '', e.spaceGroup ?? ''].filter(Boolean).join(' · ');
}

function Metric(props: { name: string; value: number }) {
  const scale = MetricScales[props.name] ?? { label: props.name, max: Math.max(1, props.value), digits: 2 };
  const share = Math.min(1, Math.max(0, props.value / scale.max));
  return (
    <div style="display: grid; grid-template-columns: 6.5em 1fr 3.5em; align-items: center; gap: var(--msp-space-2)">
      <span style={Quiet}>{scale.label}</span>
      <div style="height: 4px; border-radius: 2px; background: var(--msp-color-well); overflow: hidden">
        <div style={`height: 100%; width: ${(share * 100).toFixed(1)}%; background: var(--msp-color-accent)`} />
      </div>
      <span style="font-size: var(--msp-font-size-xs); font-variant-numeric: tabular-nums; text-align: right">
        {props.value.toFixed(scale.digits)}
      </span>
    </div>
  );
}

/**
 * A mutation: under the pointer its residues are lit in the 3D view and the sequence (`highlight.set`), a click selects
 * them (`select.set`); on the deposited model, whose residue data maps UniProt's numbering to the author's.
 */
function MutationChip(props: { ctx: ViewerUIContext; loaded: DpdbLoaded; mutation: DpdbMutation }) {
  const { ctx, loaded, mutation } = props;
  const deposited = loaded.models.find((m) => m.model.residueData.length && m.ref);
  const target = deposited ? mutationTarget(deposited.model.residueData, mutation) : '';
  const structure = deposited?.ref;
  const api = ctx.api;
  return (
    <button
      type="button"
      class={classes.chip}
      disabled={!target}
      data-tip={target ? `${mutation.text}: ${target}` : `${mutation.text} is not modeled`}
      onPointerEnter={() => target && void api.highlight.set({ target, structure, sequence: true })}
      onPointerLeave={() => target && void api.highlight.clear()}
      onClick={() => target && void api.select.set({ target, structure })}
    >
      {mutation.text}
    </button>
  );
}

function ModelBlock(props: { ctx: ViewerUIContext; area: DynamicPdbApi; m: DpdbLoadedModel }) {
  const { ctx, area, m } = props;
  const ref = m.ref;
  const shown = distinct(area.changed.pipe(map(() => !!ref && area.isShown(ref))));
  const density = distinct(area.changed.pipe(map(() => !!ref && area.density(ref).on)));
  const details = [
    m.model.modelType,
    m.model.atoms ? `${m.model.atoms.toLocaleString('en-US')} atoms` : '',
    m.model.altlocFraction ? `${Math.round(m.model.altlocFraction * 100)}% with alternates` : '',
  ].filter(Boolean);
  const metrics = Object.keys(MetricScales).filter((k) => typeof m.model.metrics[k] === 'number');
  return (
    <section style="display: grid; gap: var(--msp-space-2)">
      <div style="display: flex; align-items: center; gap: var(--msp-space-2); min-width: 0">
        <div style="display: grid; flex: 1; min-width: 0">
          <span style="font-weight: 600; font-size: var(--msp-font-size-sm)">{m.model.title}</span>
          <span style={Quiet}>{details.join(' · ')}</span>
        </div>
        {ref && area.hasDensity(ref) ? (
          <Button
            size="sm"
            variant="ghost"
            icon={<Icon name="scan" />}
            pressed={density}
            title={density.pipe(map((on) => (on ? 'Stop the density map' : 'Stream the density map (volume server)')))}
            onClick={() => void ctx.api.run('volume.streaming', { structure: ref, on: !area.density(ref).on })}
          />
        ) : null}
        {ref ? (
          <Button
            size="sm"
            variant="ghost"
            icon={<Show when={shown.pipe(map((v) => (v ? 'eye' : 'eye-off')))}>{(n) => <Icon name={n} />}</Show>}
            title={shown.pipe(map((v) => (v ? 'Hide the model' : 'Show the model')))}
            onClick={() => void ctx.api.run('structure.visibility', { structure: ref, visible: !area.isShown(ref) })}
          />
        ) : (
          <span style={Quiet}>No file</span>
        )}
      </div>
      {metrics.map((k) => (
        <Metric name={k} value={m.model.metrics[k]} />
      ))}
    </section>
  );
}

/** The tracks the sequence shows for the entry, each with its eye (`dynamic-pdb.track`, as the person). */
function TracksView(props: { ctx: ViewerUIContext; area: DynamicPdbApi }) {
  const { ctx, area } = props;
  return (
    <Show when={area.tracks.pipe(map((t) => (t.length ? t : undefined)))}>
      {(tracks) => (
        <section style="display: grid; gap: var(--msp-space-1)">
          <span style="font-weight: 600; font-size: var(--msp-font-size-sm)">Sequence tracks</span>
          {tracks.map((t) => {
            const shown = distinct(area.hiddenTracks.pipe(map((h) => !h.has(t.id))));
            return (
              <div style="display: flex; align-items: center; gap: var(--msp-space-2); min-width: 0">
                <span style={`width: 10px; height: 3px; border-radius: 2px; flex: none; background: ${t.color}`} />
                <span style="flex: 1; min-width: 0; font-size: var(--msp-font-size-xs)">{t.label}</span>
                <Button
                  size="xs"
                  variant="ghost"
                  icon={
                    <Show when={shown.pipe(map((v) => (v ? 'eye' : 'eye-off')))}>
                      {(n) => <Icon name={n} size={14} />}
                    </Show>
                  }
                  title={shown.pipe(map((v) => (v ? `Hide ${t.label}` : `Show ${t.label}`)))}
                  onClick={() =>
                    void ctx.api.run('dynamic-pdb.track', { track: t.id, shown: area.hiddenTracks.value.has(t.id) })
                  }
                />
              </div>
            );
          })}
        </section>
      )}
    </Show>
  );
}

function EntryView(props: { ctx: ViewerUIContext; area: DynamicPdbApi }) {
  const { ctx, area } = props;
  return (
    <Show
      when={area.loaded}
      fallback={() => (
        <div style="display: grid; gap: var(--msp-space-3); justify-items: start">
          <p style={`margin: 0; ${Quiet}`}>No Dynamic PDB entry is loaded.</p>
          <Button size="sm" icon={<Icon name="download" />} onClick={() => void ctx.api.run('dynamic-pdb.load', {})}>
            Load dpdb_gijqn2no
          </Button>
        </div>
      )}
    >
      {(l) => {
        const e = l.entry;
        const entity = e.entities[0];
        return (
          <div style="display: grid; gap: var(--msp-space-4); min-width: 0">
            <div style="display: grid; gap: var(--msp-space-1)">
              <span style="font-size: var(--msp-font-size-sm); line-height: 1.35; text-wrap: balance">{e.title}</span>
              <span style={Quiet}>{facts(e)}</span>
              {entity ? (
                <span style={Quiet}>
                  {[entity.description, entity.construct, entity.organism, entity.uniprot].filter(Boolean).join(' · ')}
                </span>
              ) : null}
              {entity?.mutations.length ? (
                <div style="display: flex; flex-wrap: wrap; gap: var(--msp-space-2); padding-top: var(--msp-space-1)">
                  {entity.mutations.map((mu) => (
                    <MutationChip ctx={ctx} loaded={l} mutation={mu} />
                  ))}
                </div>
              ) : null}
            </div>
            <For each={l.models}>{(m) => <ModelBlock ctx={ctx} area={area} m={m} />}</For>
            <TracksView ctx={ctx} area={area} />
            <a
              href={dpdbPage(e.id)}
              target="_blank"
              rel="noopener"
              style="display: inline-flex; align-items: center; gap: 6px; font-size: var(--msp-font-size-xs); color: var(--msp-color-accent)"
            >
              <Icon name="link" size={12} />
              {`${e.id} on dynamicpdb.com`}
            </a>
          </div>
        );
      }}
    </Show>
  );
}

/** The models' titles by state ref, for the rows. */
const titleOf = (area: DynamicPdbApi, ref: string) =>
  area.loaded.pipe(map((l) => l?.models.find((m) => m.ref === ref)?.model.title));

export function dynamicPdbUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api['dynamic-pdb'];
  if (!area) return {};
  return {
    overlays: [
      {
        id: DpdbCard,
        title: 'Dynamic PDB',
        render: domView(() => (
          <div style={`display: grid; min-width: 0; padding: 0 ${Inset} ${Inset}`}>
            <EntryView ctx={ctx} area={area} />
          </div>
        )),
      },
    ],
    rows: [
      {
        id: 'dynamic-pdb.row.entry',
        when: (node: SceneNode) => node.kind === 'structure' && area.isOfEntry(node.id),
        render: domView(({ ctx: view }) => {
          const ref = (view.subject as SceneNode).id;
          return (
            <SceneRow
              label="Dynamic PDB"
              description={titleOf(area, ref)}
              onOpen={() => void ctx.api.view.open({ card: DpdbCard })}
            />
          );
        }),
      },
    ],
    toolbar: [
      {
        id: 'dynamic-pdb.toolbar.entry',
        icon: 'layers',
        title: 'Dynamic PDB',
        group: 'scene',
        overlay: DpdbCard,
        visible: distinct(area.loaded.pipe(map((l) => !!l))),
      },
    ],
  };
}
