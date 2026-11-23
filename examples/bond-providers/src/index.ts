/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  BallAndStickParams,
  BallAndStickRepresentationProvider,
} from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { syncCameras } from '@molstar/prism-viewer';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import {
  type Bonds,
  download,
  type Removals,
  type Source,
  showServed,
  showTested,
  statesOrders,
  type Tested,
  updateRemovals,
} from './load.js';
import { BondProvidersViewerSpec } from './spec.js';
import './index.html';

// The bond providers example page (the bond-order benchmark's plan §7): two viewers side by side, cameras linked. Left,
// a file as served with its built-in bonds; right, the same file with the removals asked and the bonds of the provider
// under test. The URL follows the state, so a view is a link (the benchmark report's links fill it).

const Perception = 'bond-order-perception';

/** The right side's providers, perception first; the example's own are proofs of concept. */
const Providers: Record<string, { label: string; params?: Record<string, unknown> }> = {
  [Perception]: { label: 'Bond-order perception' },
  default: { label: 'Default (built-in)' },
  none: { label: 'None (format only)' },
  'naive-distance': { label: 'Naive distance', params: { cutoff: 1.9 } },
  'naive-orders': { label: 'Naive orders' },
  random: { label: 'Random', params: { probability: 0.1, maxDistance: 4, seed: 42 } },
  // naive orders for the ligands, the built-in bonds for the rest and between units
  'per-unit': {
    label: 'By unit (ligands: naive orders)',
    params: {
      rules: [
        { match: { name: 'entity-type', params: ['non-polymer'] }, provider: { name: 'naive-orders', params: {} } },
      ],
      inter: { provider: { name: 'default', params: {} } },
    },
  },
};

/** The page's state, as the URL holds it (defaults left out). */
interface PageState {
  /** the PDB / CCD switch: `auto` guesses from the code */
  kind: 'auto' | Source['kind'];
  id: string;
  coords: Source['coords'];
  /** a copy to focus: `label_asym_id`, `auth_seq_id` */
  asym?: string;
  seq?: number;
  provider: string;
  mode: 'auto' | 'forceCompute';
  removeHydrogens: boolean;
  removeCharges: boolean;
}

const Defaults: PageState = {
  kind: 'auto',
  id: '1HSG',
  coords: 'ideal',
  provider: Perception,
  mode: 'forceCompute',
  removeHydrogens: false,
  removeCharges: true,
};

/** 4 characters starting with a digit: a PDB ID; anything else a CCD code. */
function guessKind(id: string): Source['kind'] {
  return /^\d[A-Z0-9]{3}$/i.test(id) ? 'pdb' : 'ccd';
}

function resolve(state: PageState): Source['kind'] {
  return state.kind === 'auto' ? guessKind(state.id) : state.kind;
}

function stateOf(search: string): PageState {
  const q = new URLSearchParams(search);
  const pdb = q.get('pdb')?.trim();
  const ccd = q.get('ccd')?.trim();
  const id = (pdb || ccd || Defaults.id).toUpperCase();
  const given: Source['kind'] | undefined = pdb ? 'pdb' : ccd ? 'ccd' : void 0;
  const seq = q.get('seq');
  const provider = q.get('bonds') ?? Defaults.provider;
  return {
    // the switch shows `auto` when the guess agrees with the link
    kind: given && given !== guessKind(id) ? given : 'auto',
    id,
    coords: q.get('coords') === 'model' ? 'model' : 'ideal',
    asym: q.get('asym') || void 0,
    seq: seq !== null && seq !== '' && Number.isFinite(Number(seq)) ? Number(seq) : void 0,
    provider: provider in Providers ? provider : Defaults.provider,
    mode: q.get('mode') === 'auto' ? 'auto' : 'forceCompute',
    removeHydrogens: q.get('removeH') === '1',
    removeCharges: q.get('removeCharges') !== '0',
  };
}

function searchOf(state: PageState): string {
  const kind = resolve(state);
  const q = new URLSearchParams({ [kind]: state.id });
  if (kind === 'ccd' && state.coords !== Defaults.coords) q.set('coords', state.coords);
  if (state.asym) q.set('asym', state.asym);
  if (state.seq !== void 0) q.set('seq', String(state.seq));
  if (state.provider !== Defaults.provider) q.set('bonds', state.provider);
  if (state.provider === Perception && state.mode !== Defaults.mode) q.set('mode', state.mode);
  if (state.removeHydrogens) q.set('removeH', '1');
  if (!state.removeCharges) q.set('removeCharges', '0');
  return `?${q}`;
}

function element<T extends HTMLElement>(selector: string): T {
  const e = document.querySelector<T>(selector);
  if (!e) throw new Error(`The page has no ${selector}.`);
  return e;
}

// Kekulé orders rather than dashed aromatic bonds in every ball-and-stick (presets, the focus, one added by hand): no
// plugin setting gives a representation's params a default, so the definition and the provider's defaults are patched
// before the viewers start.
BallAndStickParams.aromaticBonds = PD.Boolean(false, { description: BallAndStickParams.aromaticBonds.description });
(BallAndStickRepresentationProvider.defaultValues as { aromaticBonds: boolean }).aromaticBonds = false;

async function main() {
  const form = element<HTMLFormElement>('#load');
  const code = element<HTMLInputElement>('#code');
  const kind = element<HTMLSelectElement>('#kind');
  const coords = element<HTMLSelectElement>('#coords');
  const provider = element<HTMLSelectElement>('#provider');
  const mode = element<HTMLSelectElement>('#mode');
  const autoMode = mode.querySelector<HTMLOptionElement>('option[value="auto"]')!;
  const removeH = element<HTMLInputElement>('#removeH');
  const removeCharges = element<HTMLInputElement>('#removeCharges');
  const status = element<HTMLElement>('#status');

  provider.innerHTML = Object.entries(Providers)
    .map(([name, p]) => `<option value="${name}">${p.label}</option>`)
    .join('');

  let state = stateOf(location.search);
  const showState = () => {
    code.value = state.id;
    kind.value = state.kind;
    coords.value = state.coords;
    coords.hidden = resolve(state) !== 'ccd';
    provider.value = state.provider;
    mode.value = state.mode;
    mode.disabled = state.provider !== Perception;
    removeH.checked = state.removeHydrogens;
    removeCharges.checked = state.removeCharges;
  };
  const writeUrl = () => history.replaceState(null, '', `${location.pathname}${searchOf(state)}`);
  const fail = (e: unknown) => {
    console.error(e);
    status.textContent = `Failed: ${e instanceof Error ? e.message : String(e)}`;
  };
  showState();

  // two viewers on a page: each its own persist key (Prism README decision 29)
  const make = (id: string) => {
    const spec = BondProvidersViewerSpec();
    return createPrismViewer({
      target: element(`#${id}`),
      spec: { ...spec, persist: { ...spec.persist, key: `bond-providers-${id}` } },
    });
  };
  const [served, tested] = await Promise.all([make('served'), make('tested')]);
  (window as unknown as { viewers: PrismViewer[] }).viewers = [served, tested];
  syncCameras([served, tested], served);

  const bondsOf = (s: PageState): Bonds =>
    s.provider === Perception
      ? { provider: Perception, params: { mode: s.mode } }
      : { provider: s.provider, params: Providers[s.provider].params ?? {} };
  const removalsOf = (s: PageState): Removals => ({
    removeHydrogens: s.removeHydrogens,
    removeCharges: s.removeCharges,
  });

  let current: Tested | undefined;

  const focus = async (viewer: PrismViewer) => {
    if (!state.asym || state.seq === void 0) return;
    const r = await viewer.api.focus.set({ target: { label_asym_id: state.asym, auth_seq_id: state.seq } });
    if (!r.ok) throw new Error(r.message);
  };

  // `auto` keeps the orders a model states (plan §7): disabled, and `forceCompute` taken, when the model has them
  const checkMode = async () => {
    const stated = !!current && (await statesOrders(current));
    autoMode.disabled = stated;
    mode.title = stated ? 'auto keeps the orders the file states: nothing to perceive here' : '';
    if (stated && state.mode === 'auto') {
      state = { ...state, mode: 'forceCompute' };
      showState();
      writeUrl();
      await applyBonds();
      return 'The file states its orders: forceCompute instead of auto.';
    }
    return '';
  };

  const describe = () => {
    const p = Providers[state.provider].label + (state.provider === Perception ? ` (${state.mode})` : '');
    const removed = [state.removeHydrogens && 'H', state.removeCharges && 'charges'].filter(Boolean).join(' and ');
    return `Right: ${p}${removed ? `, ${removed} removed` : ''}. Cameras linked.`;
  };

  const load = async () => {
    const source: Source = { kind: resolve(state), id: state.id, coords: state.coords };
    writeUrl();
    status.textContent = `Loading ${source.id}…`;
    current = void 0;
    try {
      const data = await download(source);
      const [, t] = await Promise.all([
        showServed(served, source, data),
        showTested(tested, source, data, removalsOf(state), bondsOf(state)),
      ]);
      current = t;
      await Promise.all([focus(served), focus(tested)]);
      const note = await checkMode();
      status.textContent = `${note} ${describe()}`.trim();
    } catch (e) {
      fail(e);
    }
  };

  // the provider and the mode: the right side's bonds switched in place
  const applyBonds = async () => {
    if (!current) return;
    const bonds = bondsOf(state);
    const r = await tested.api.structure.bonds({ provider: bonds.provider, params: bonds.params });
    status.textContent = r.ok ? describe() : `${r.message}${r.hint ? ` ${r.hint}` : ''}`;
  };

  // the removals: the right side parsed again, the camera kept
  const applyRemovals = async () => {
    if (!current) return;
    status.textContent = 'Parsing again…';
    try {
      await updateRemovals(tested, current, removalsOf(state));
      await focus(tested);
      const note = await checkMode();
      status.textContent = `${note} ${describe()}`.trim();
    } catch (e) {
      fail(e);
    }
  };

  const change = (next: Partial<PageState>) => {
    state = { ...state, ...next };
    showState();
    writeUrl();
  };

  code.addEventListener('input', () => {
    coords.hidden = resolve({ ...state, kind: kind.value as PageState['kind'], id: code.value.trim() }) !== 'ccd';
  });
  kind.addEventListener('change', () => {
    coords.hidden = resolve({ ...state, kind: kind.value as PageState['kind'], id: code.value.trim() }) !== 'ccd';
  });
  form.addEventListener('submit', (e) => {
    e.preventDefault();
    const id = code.value.trim().toUpperCase();
    if (!id) return;
    // a new code drops the focus of the previous one
    const sameEntry = id === state.id;
    change({
      id,
      kind: kind.value as PageState['kind'],
      coords: coords.value as Source['coords'],
      asym: sameEntry ? state.asym : void 0,
      seq: sameEntry ? state.seq : void 0,
    });
    void load();
  });
  provider.addEventListener('change', () => {
    change({ provider: provider.value });
    void applyBonds();
  });
  mode.addEventListener('change', () => {
    change({ mode: mode.value as PageState['mode'] });
    void applyBonds();
  });
  removeH.addEventListener('change', () => {
    change({ removeHydrogens: removeH.checked });
    void applyRemovals();
  });
  removeCharges.addEventListener('change', () => {
    change({ removeCharges: removeCharges.checked });
    void applyRemovals();
  });

  await load();
}

main().catch((e) => {
  console.error(e);
  const p = document.createElement('p');
  p.style.cssText = 'padding: 24px';
  p.textContent = `The viewer could not start: ${e instanceof Error ? e.message : String(e)}`;
  document.body.appendChild(p);
});
