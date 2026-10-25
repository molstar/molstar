/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { StructureElement, StructureProperties as P } from '@molstar/model/model/structure';
import type { ApiContext } from '@molstar/prism-viewer/api/area';
import { messageOf } from '@molstar/prism-viewer/api/run';
import type { ExecResult, ScriptHost } from '../language/script.js';
import type { Command, Sel } from '../language/ir.js';
import { selectionOn } from '../language/selection.js';
import { VmdLinesPresetId } from '../plugin/vmd-preset.js';
import {
  atomCount,
  atomProperty,
  compactExpression,
  eachAtom,
  geometryOf,
  labelText,
  positionOf,
  ReadableProperties,
  shapeOf,
} from './atoms.js';
import { cartoon, display } from './display.js';
import { HostFailure, Ops, type Where } from './ops.js';
import { colorTheme, hex } from './reps.js';
import { envOf, type StateHolder } from './state.js';
import { addRep, modifyRep, removeRep, showRep } from './vmd.js';

// The IR on Prism's API (_spec §5.1): each command through the viewer's operations, called with `script.run`'s
// context, so every change joins its one step (Undo takes the submission back; a failure rolls it back). Selections are
// evaluated per structure, and handed to operations as `{ molql }` targets with the structure's ref.

/** Labels made at most per command: more is a selection to narrow. */
const MaxLabels = 100;
/** Distances `distance … cutoff` adds at most (each a measurement in the scene), and the pairs it checks. */
const MaxContacts = 50;
const MaxContactChecks = 1_000_000;
/** Values `$sel get` and reads print at most. */
const MaxValues = 50;

export class PrismScriptHost implements ScriptHost {
  readonly #ops: Ops;

  constructor(
    ctx: ApiContext,
    readonly holder: StateHolder,
    readonly files: ReadonlyMap<string, File>,
    /** The texts of the scripts among the files (what `@file.pml` and `source file.tcl` run). */
    readonly scripts: ReadonlyMap<string, string> = new Map(),
  ) {
    this.#ops = new Ops(ctx, holder);
  }

  env() {
    return { ...envOf(this.#ops.ctx, this.holder.value, [...this.files.keys()]), scripts: this.scripts };
  }

  /** Names every loaded structure in the state (an object name and a molecule id), so they last once used. */
  keepObjects(): void {
    const objects = this.#ops.objects();
    const s = this.holder.value;
    const known = { ...s.objects };
    let changed = false;
    for (const o of objects) {
      if (!known[o.ref]) {
        known[o.ref] = { name: o.name, molecule: o.molecule };
        changed = true;
      }
    }
    if (changed) {
      this.holder.update({
        objects: known,
        nextMolecule: Math.max(s.nextMolecule, ...Object.values(known).map((x) => x.molecule + 1)),
      });
    }
  }

  async execute(cmd: Command): Promise<ExecResult> {
    try {
      const r = await this.#execute(cmd);
      return typeof r === 'string' ? { ok: true, did: r } : r;
    } catch (e) {
      if (e instanceof HostFailure) return { ok: false, message: e.message, ...(e.hint ? { hint: e.hint } : {}) };
      return { ok: false, message: messageOf(e) };
    }
  }

  async #execute(cmd: Command): Promise<string | ExecResult> {
    const ops = this.#ops;
    switch (cmd.kind) {
      case 'load':
        return this.#load(cmd);
      case 'object.remove': {
        const refs = cmd.objects === 'all' ? ops.objects().map((o) => o.ref) : cmd.objects.map((o) => ops.ref(o));
        for (const ref of refs) await ops.run('structure.remove', { structure: ref });
        const objects = { ...this.holder.value.objects };
        for (const ref of refs) delete objects[ref];
        this.holder.update({
          objects,
          ...(this.holder.value.top && refs.includes(this.holder.value.top) ? { top: undefined } : {}),
        });
        return `Removed ${refs.length} structure${refs.length === 1 ? '' : 's'}.`;
      }
      case 'object.visibility': {
        const refs = cmd.objects === 'all' ? ops.objects().map((o) => o.ref) : cmd.objects.map((o) => ops.ref(o));
        for (const ref of refs) await ops.run('structure.visibility', { visible: cmd.visible, structure: ref });
        return `${cmd.visible ? 'Showed' : 'Hid'} ${refs.length} structure${refs.length === 1 ? '' : 's'}.`;
      }
      case 'object.top':
        this.holder.update({ top: ops.ref(cmd.object) });
        return `${cmd.object} is the top molecule.`;
      case 'name.set':
        this.holder.update({ names: { ...this.holder.value.names, [cmd.name]: cmd.sel } });
        return `Named '${cmd.sel.text}' ${cmd.name}.`;
      case 'name.remove': {
        const names =
          cmd.names === 'all'
            ? {}
            : Object.fromEntries(Object.entries(this.holder.value.names).filter(([k]) => !cmd.names.includes(k)));
        this.holder.update({ names });
        return cmd.names === 'all' ? 'Forgot the named selections.' : `Forgot ${cmd.names.join(', ')}.`;
      }
      case 'select': {
        const where = cmd.mode === 'remove' ? ops.where(cmd.sel) : ops.some(cmd.sel);
        if (cmd.mode === 'set') await ops.run('select.clear', {});
        for (const w of where) {
          await ops.run('select.set', {
            target: target(cmd.sel, w),
            mode: cmd.mode === 'set' ? 'add' : cmd.mode,
            structure: w.ref,
          });
        }
        return `Selected ${count(where)} atoms${where.length > 1 ? ` in ${where.length} structures` : ''}.`;
      }
      case 'select.clear':
        await ops.run('select.clear', {});
        return 'Cleared the selection.';
      case 'display':
        return display(ops, cmd);
      case 'cartoon':
        return cartoon(ops, cmd);
      case 'rep.add':
        return addRep(ops, cmd);
      case 'rep.remove':
        return removeRep(ops, cmd);
      case 'rep.modify':
        return modifyRep(ops, cmd);
      case 'rep.visibility':
        return showRep(ops, cmd);
      case 'rep.defaults':
        this.holder.update({ repDefaults: { ...this.holder.value.repDefaults, ...cmd.patch } });
        // Says which command uses it: the style alone changes nothing on screen (VMD's way), which reads as a failure.
        return `The next rep (\`mol addrep\`) takes ${Object.entries(cmd.patch)
          .map(([k, v]) => `${k} ${typeof v === 'object' ? ('rgb' in v ? hex(v.rgb) : v.theme) : v}`)
          .join(', ')}.`;
      case 'color':
        return this.#color(cmd);
      case 'require':
        ops.some(cmd.sel);
        return { ok: true };
      case 'transparency': {
        const where = ops.some(cmd.sel);
        for (const w of where)
          await ops.run('representation.transparency', {
            target: target(cmd.sel, w),
            alpha: cmd.alpha,
            structure: w.ref,
          });
        return `Made ${count(where)} atoms ${Math.round(cmd.alpha * 100)}% transparent.`;
      }
      case 'size': {
        const where = ops.some(cmd.sel);
        const factor = Math.min(10, Math.max(0.05, cmd.factor));
        for (const w of where)
          await ops.run('representation.size', { target: target(cmd.sel, w), factor, structure: w.ref });
        const warnings = [
          ...(cmd.reps ? ['Mol* sizes the atoms in every representation, not only one kind.'] : []),
          ...(factor !== cmd.factor
            ? [`Sizes go from 0.05 to 10 times the default: ${cmd.factor} became ${factor}.`]
            : []),
        ];
        return { ok: true, did: `Sized ${count(where)} atoms ×${factor}.`, ...(warnings.length ? { warnings } : {}) };
      }
      case 'hide-contents':
        return this.#hideContents(cmd.sel);
      case 'hydrogens':
        await ops.run('look.set', { hydrogens: cmd.visible ? 'all' : 'hide-all' });
        return cmd.visible ? 'Showed hydrogens.' : 'Hid hydrogens.';
      case 'camera':
        return this.#camera(cmd);
      case 'camera.zoom': {
        const factor = Math.min(10, Math.max(0.1, cmd.factor));
        await ops.run('camera.zoom', { factor });
        return `Zoomed by ${factor}.`;
      }
      case 'camera.turn':
      case 'camera.move':
        return this.#turnOrMove(cmd);
      case 'measure':
        return this.#measure(cmd);
      case 'measure.clear':
        return this.#clearMeasurements(cmd.measure);
      case 'label':
        return this.#label(cmd);
      case 'label.clear':
        return this.#clearMeasurements('label', cmd.sel);
      case 'superpose':
        return this.#superpose(cmd);
      case 'canvas':
        return this.#canvas(cmd);
      case 'trajectory':
        if (cmd.play === 'stop') await ops.run('trajectory.stop', {});
        else if (cmd.play) await ops.run('trajectory.play', {});
        else if (cmd.frame === -1) {
          await ops.run('trajectory.frame', { index: 0 });
          await ops.run('trajectory.frame', { by: -1 });
        } else await ops.run('trajectory.frame', { index: cmd.frame ?? 0 });
        return cmd.play
          ? cmd.play === 'stop'
            ? 'Stopped the trajectory.'
            : 'Playing the trajectory.'
          : `Showed frame ${cmd.frame === -1 ? 'last' : cmd.frame}.`;
      case 'export.image': {
        const r = await ops.run('export.image', {
          ...(cmd.name ? { name: cmd.name } : {}),
          ...(cmd.width && cmd.height
            ? { size: { name: 'custom', params: { width: cmd.width, height: cmd.height } } }
            : {}),
          ...(cmd.transparent ? { transparent: true } : {}),
        });
        return r.did;
      }
      case 'history':
        throw new HostFailure(`${cmd.action} runs alone, not inside a script.`, `Type ${cmd.action} on its own.`);
      case 'read':
        return { ok: true, output: this.#read(cmd.read) };
      case 'print':
        return { ok: true, output: cmd.text };
    }
  }

  async #load(cmd: Extract<Command, { kind: 'load' }>): Promise<string | ExecResult> {
    const ops = this.#ops;
    const before = new Set(ops.structures().map((s) => s.cell.transform.ref));
    const vmd = cmd.look === 'vmd-lines';
    let r;
    if ('entry' in cmd.from) {
      r = await ops.run(
        'load.entry',
        { source: cmd.from.entry.db, id: cmd.from.entry.id, ...(vmd ? { preset: VmdLinesPresetId } : {}) },
        { copy: cmd.copy },
      );
    } else if ('url' in cmd.from) {
      r = await ops.run(
        'load.url',
        { url: cmd.from.url, ...(cmd.format ? { format: cmd.format } : {}) },
        { copy: cmd.copy },
      );
    } else {
      const name = cmd.from.file;
      const file =
        this.files.get(name) ??
        [...this.files.entries()].find(([k]) => k.split('/').pop() === name.split('/').pop())?.[1];
      if (!file) throw new HostFailure(`There is no file '${name}' with this script.`);
      r = await ops.run('load.data', { file, ...(cmd.format ? { format: cmd.format } : {}) }, { copy: cmd.copy });
    }
    const made = ops.structures().filter((s) => !before.has(s.cell.transform.ref));
    if (!made.length) return r.did;
    // The compiler's name, unless another structure has it (`fetch 1cbs` when 1cbs came from a file: the entry is not
    // the same, so it loaded): then a free one, said.
    const taken = new Set(
      ops
        .objects()
        .filter((o) => !made.some((s) => s.cell.transform.ref === o.ref))
        .map((o) => o.name.toLowerCase()),
    );
    const free = (wish: string) => {
      let n = wish;
      for (let k = 2; taken.has(n.toLowerCase()); k++) n = `${wish}_${k}`;
      taken.add(n.toLowerCase());
      return n;
    };
    const objects = { ...this.holder.value.objects };
    const names: string[] = [];
    // Molecule ids are never reused: one another structure has (the compiler expected the loaded one to stay) is
    // replaced by the next free id.
    const usedIds = new Set(
      ops
        .objects()
        .filter((o) => !made.some((s) => s.cell.transform.ref === o.ref))
        .map((o) => o.molecule),
    );
    let nextId = Math.max(this.holder.value.nextMolecule, cmd.molecule, ...[...usedIds].map((m) => m + 1));
    made.forEach((s, i) => {
      const name = free(i === 0 ? cmd.name : `${cmd.name}_${i + 1}`);
      names.push(name);
      const wanted = cmd.molecule + i;
      const molecule = usedIds.has(wanted) ? nextId++ : wanted;
      usedIds.add(molecule);
      nextId = Math.max(nextId, molecule + 1);
      objects[s.cell.transform.ref] = { name, molecule };
    });
    this.holder.update({
      objects,
      nextMolecule: Math.max(this.holder.value.nextMolecule, nextId),
      ...(vmd ? { top: made[made.length - 1].cell.transform.ref } : {}),
    });
    if (vmd && !('entry' in cmd.from)) {
      for (const s of made)
        await ops.run('structure.preset', { preset: VmdLinesPresetId, structure: s.cell.transform.ref });
    }
    const did = `${r.did.replace(/\.$/, '')} as ${names[0]}${vmd ? ` (molecule ${cmd.molecule})` : ''}.`;
    if (names[0] === cmd.name) return did;
    return { ok: true, did, warnings: [`Another structure is named ${cmd.name}: this one is ${names[0]}.`] };
  }

  async #color(cmd: Extract<Command, { kind: 'color' }>): Promise<string | ExecResult> {
    const ops = this.#ops;
    // One of a group (a color per element): selecting nothing is no failure; the group's `require` said whether the
    // selection selects anything.
    const where = cmd.optional ? ops.where(cmd.sel) : ops.some(cmd.sel);
    if (!where.length) return { ok: true };
    if ('rgb' in cmd.color) {
      for (const w of where) {
        await ops.run('representation.color', {
          target: target(cmd.sel, w),
          color: hex(cmd.color.rgb),
          structure: w.ref,
        });
      }
      const n = count(where);
      return `Colored ${n} atom${n === 1 ? '' : 's'} ${hex(cmd.color.rgb)}.`;
    }
    // A theme goes on the representations of the components the selection covers (as `representation.color` does),
    // under the colors set per atom (the color layer). On a whole structure those colors are cleared first (PyMOL's
    // `util.cbc` after `color red` shows the chains); on part of one, the theme still covers whole components: said.
    const theme = colorTheme(cmd.color);
    const nodes: string[] = [];
    let partial = 0;
    for (const w of where) {
      if (StructureElement.Loci.isWholeStructure(w.loci)) {
        await ops.run('representation.reset', { layer: 'color', structure: w.ref });
      }
      for (const c of w.structure.components) {
        const data = c.cell.obj?.data;
        if (!data) continue;
        const own = StructureElement.Loci.remap(w.loci, data);
        if (StructureElement.Loci.isEmpty(own)) continue;
        if (atomCount(own) < data.elementCount) partial++;
        for (const r of c.representations) nodes.push(r.cell.transform.ref);
      }
    }
    if (!nodes.length) throw new HostFailure(`Nothing shows '${cmd.sel.text}'.`, 'Show it first: show cartoon, …');
    for (const node of nodes) await ops.run('representation.update', { node, params: { colorTheme: theme } });
    const did = `Colored ${nodes.length} representation${nodes.length === 1 ? '' : 's'} by ${theme.name}.`;
    if (!partial) return did;
    return {
      ok: true,
      did,
      warnings: [
        `Mol* colors whole components by a theme: ${partial} ${partial === 1 ? 'component holds' : 'components hold'} more than '${cmd.sel.text}', and colors set with color stay on top.`,
      ],
    };
  }

  /** PyMOL `remove`: whole chains leave every component (`contents.visibility`); other atoms are not removed. */
  async #hideContents(sel: Sel): Promise<string> {
    const ops = this.#ops;
    let chains = 0;
    for (const w of ops.some(sel)) {
      const whole = StructureElement.Loci.extendToWholeChains(w.loci);
      if (atomCount(whole) !== atomCount(w.loci)) {
        throw new HostFailure(
          `'${sel.text}' is not whole chains: Mol* leaves out whole entities or chains only.`,
          `Hide the atoms instead: hide everything, ${sel.text}`,
        );
      }
      const pairs = new Set<string>();
      StructureElement.Loci.forEachLocation(w.loci, (loc) => {
        pairs.add(`${P.chain.label_entity_id(loc)}\u0000${P.chain.auth_asym_id(loc)}`);
      });
      for (const pair of pairs) {
        const [entity, chain] = pair.split('\u0000');
        await ops.run('contents.visibility', { entity, chain, visible: false, structure: w.ref });
        chains++;
      }
    }
    return `Left ${chains} chain${chains === 1 ? '' : 's'} out of the scene.`;
  }

  async #camera(cmd: Extract<Command, { kind: 'camera' }>): Promise<string> {
    const ops = this.#ops;
    const animate = cmd.animate;
    if (cmd.action === 'reset') {
      await ops.run('camera.reset', { animate });
      return 'Reset the camera.';
    }
    if (!cmd.sel) {
      if (cmd.action === 'orient') {
        for (const o of ops.objects()) await ops.run('camera.orient', { structure: o.ref, animate });
        return 'Oriented the view along the principal axes.';
      }
      await ops.run('camera.reset', { animate });
      return 'Framed everything.';
    }
    const where = ops.some(cmd.sel);
    if (cmd.action === 'orient') await ops.run('camera.orient', { structure: where[0].ref, animate: false });
    if (cmd.action === 'center') return this.#center(where, animate);
    if (where.length === 1 && cmd.buffer === undefined) {
      await ops.run('camera.frame', { target: target(cmd.sel, where[0]), structure: where[0].ref, animate });
    } else {
      // Several structures, or a margin: the camera's own framing of the atoms (view state, as camera.frame is).
      ops.plugin.managers.camera.focusLoci(
        where.map((w) => w.loci),
        { extraRadius: cmd.buffer ?? 0, durationMs: animate ? undefined : 0 },
      );
    }
    return `Framed ${count(where)} atoms.`;
  }

  /** The view moved to look at the atoms' center, at the same distance (PyMOL `center`). */
  async #center(where: readonly Where[], animate: boolean): Promise<string> {
    const camera = this.#ops.plugin.canvas3d?.camera;
    if (!camera) throw new HostFailure('There is no canvas.');
    const centers = where.map((w) => shapeOf(w.loci).center);
    const center = Vec3();
    for (const c of centers) Vec3.add(center, center, c);
    Vec3.scale(center, center, 1 / centers.length);
    const s = camera.getSnapshot();
    const offset = Vec3.sub(Vec3(), s.position, s.target);
    const position = Vec3.add(Vec3(), center, offset);
    await this.#ops.run('camera.set', { position: [...position], target: [...center], up: [...s.up], animate });
    return 'Centered the view.';
  }

  /** PyMOL `turn`/`move`, VMD `rotate`/`translate`: about or along the screen's axes. */
  async #turnOrMove(cmd: Extract<Command, { kind: 'camera.turn' | 'camera.move' }>): Promise<string> {
    const camera = this.#ops.plugin.canvas3d?.camera;
    if (!camera) throw new HostFailure('There is no canvas.');
    const s = camera.getSnapshot();
    const forward = Vec3.normalize(Vec3(), Vec3.sub(Vec3(), s.target, s.position));
    const up = Vec3.normalize(Vec3(), s.up);
    const right = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), forward, up));
    const axis = cmd.axis === 'x' ? right : cmd.axis === 'y' ? up : Vec3.negate(Vec3(), forward);
    let position = Vec3.clone(s.position);
    let target_ = Vec3.clone(s.target);
    let newUp = Vec3.clone(s.up);
    if (cmd.kind === 'camera.turn') {
      // Turning the scene by an angle is the camera turning the other way about the target.
      const m = Mat4.fromRotation(Mat4(), (-cmd.angle * Math.PI) / 180, axis);
      const offset = Vec3.transformMat4(Vec3(), Vec3.sub(Vec3(), s.position, s.target), m);
      position = Vec3.add(Vec3(), s.target, offset);
      newUp = Vec3.transformMat4(Vec3(), s.up, m);
    } else {
      const d = Vec3.scale(Vec3(), axis, -cmd.distance);
      position = Vec3.add(Vec3(), s.position, d);
      target_ = Vec3.add(Vec3(), s.target, d);
    }
    await this.#ops.run('camera.set', {
      position: [...position],
      target: [...target_],
      up: [...newUp],
      animate: false,
    });
    return cmd.kind === 'camera.turn'
      ? `Turned the view ${cmd.angle}° about ${cmd.axis}.`
      : `Moved the view ${cmd.distance} Å along ${cmd.axis}.`;
  }

  /** One atom per selection, all in one structure. */
  #atoms(sels: readonly Sel[]): { ref: string; atoms: Where[] } {
    const atoms = sels.map((sel) => {
      const where = this.#ops.some(sel);
      const n = count(where);
      if (n !== 1 || where.length !== 1) {
        throw new HostFailure(
          `'${sel.text}' selects ${n} atoms${where.length > 1 ? ` in ${where.length} structures` : ''}: a measurement takes one atom each.`,
          'Name one atom: /1cbs//A/45/CA, or chain A and resi 45 and name CA.',
        );
      }
      return where[0];
    });
    const ref = atoms[0].ref;
    if (atoms.some((a) => a.ref !== ref))
      throw new HostFailure('The atoms are in different structures: Mol* measures within one.');
    return { ref, atoms };
  }

  async #measure(cmd: Extract<Command, { kind: 'measure' }>): Promise<string> {
    if (cmd.cutoff !== undefined) return this.#contacts(cmd.atoms, cmd.cutoff);
    const { ref, atoms } = this.#atoms(cmd.atoms);
    const keys = ['a', 'b', 'c', 'd'];
    const values: Record<string, unknown> = { structure: ref };
    atoms.forEach((w, i) => {
      values[keys[i]] = { molql: compactExpression(w.loci) };
    });
    const r = await this.#ops.run(`measure.${cmd.measure}`, values);
    return r.did;
  }

  /** `distance name, s1, s2, cutoff`: a distance for every pair of the two sets within the cutoff, in one structure. */
  async #contacts(sels: readonly Sel[], cutoff: number): Promise<string> {
    const [one, two] = sels.map((sel) => this.#ops.some(sel));
    const pairs: { ref: string; a: StructureElement.Loci; b: StructureElement.Loci }[] = [];
    for (const w1 of one) {
      for (const w2 of two) {
        if (w1.ref !== w2.ref) continue;
        if (atomCount(w1.loci) * atomCount(w2.loci) > MaxContactChecks) {
          throw new HostFailure(
            `'${sels[0].text}' and '${sels[1].text}' are too many atoms to pair up.`,
            'Narrow one: the ligand and the residues around it.',
          );
        }
        const xs = eachAtom(w1.loci, Infinity).map((x) => ({ loci: x.loci, at: positionOf(x.loci)! }));
        const ys = eachAtom(w2.loci, Infinity).map((y) => ({ loci: y.loci, at: positionOf(y.loci)! }));
        for (const x of xs) {
          for (const y of ys) {
            const d = Vec3.distance(x.at, y.at);
            if (d > 0 && d <= cutoff) pairs.push({ ref: w1.ref, a: x.loci, b: y.loci });
          }
        }
      }
    }
    if (!pairs.length) return `No atoms of '${sels[0].text}' and '${sels[1].text}' within ${cutoff} Å.`;
    if (pairs.length > MaxContacts) {
      throw new HostFailure(
        `${pairs.length} pairs within ${cutoff} Å: Mol* measures at most ${MaxContacts} at once.`,
        'Narrow a selection, or lower the cutoff.',
      );
    }
    for (const p of pairs) {
      await this.#ops.run('measure.distance', {
        structure: p.ref,
        a: { molql: compactExpression(p.a) },
        b: { molql: compactExpression(p.b) },
      });
    }
    return `Added ${pairs.length} distance${pairs.length === 1 ? '' : 's'} within ${cutoff} Å.`;
  }

  async #label(cmd: Extract<Command, { kind: 'label' }>): Promise<string> {
    const where = this.#ops.some(cmd.sel);
    const total = count(where);
    if (total > MaxLabels) {
      throw new HostFailure(
        `'${cmd.sel.text}' has ${total} atoms: Mol* labels at most ${MaxLabels} at once.`,
        'Narrow it: name CA and resi 40-60.',
      );
    }
    for (const w of where) {
      for (const atom of eachAtom(w.loci, MaxLabels)) {
        await this.#ops.run('measure.label', {
          a: { molql: compactExpression(atom.loci) },
          text: labelText(atom.location, cmd.text),
          structure: w.ref,
        });
      }
    }
    return `Labelled ${total} atom${total === 1 ? '' : 's'}.`;
  }

  async #clearMeasurements(kind: 'distance' | 'angle' | 'dihedral' | 'label', sel?: Sel): Promise<string> {
    const list = await this.#ops.run('measure.list', {});
    const all = (list.value as readonly { ref: string; kind: string }[] | undefined) ?? [];
    let doomed = all.filter((m) => m.kind === kind);
    if (sel && sel.text !== 'all') {
      // The labels on atoms of the selection: by their atoms' places.
      const where = this.#ops.where(sel);
      doomed = doomed.filter((m) => labelledIn(this.#ops.plugin, m.ref, where));
    }
    for (const m of doomed) await this.#ops.run('measure.remove', { measurement: m.ref });
    const what = kind === 'label' ? 'label' : `${kind} measurement`;
    return `Removed ${doomed.length} ${what}${doomed.length === 1 ? '' : 's'}.`;
  }

  async #superpose(cmd: Extract<Command, { kind: 'superpose' }>): Promise<string> {
    const ops = this.#ops;
    const one = (sel: Sel, what: string) => {
      const where = ops.some(sel);
      if (where.length !== 1) {
        throw new HostFailure(
          `The ${what} '${sel.text}' is in ${where.length} structures: give one structure.`,
          `Name it: /1cbs and ${sel.text}.`,
        );
      }
      return where[0];
    };
    const target_ = one(cmd.target, 'target');
    const mobile = one(cmd.mobile, 'mobile');
    if (target_.ref === mobile.ref) throw new HostFailure('The mobile and the target are in the same structure.');
    if (cmd.pairs && cmd.pairs.length > 1) this.#checkPairs(cmd.pairs);
    const targets = [target_, mobile].map((w) => ({ target: { molql: compactExpression(w.loci) }, structure: w.ref }));
    const id =
      cmd.method === 'sequence'
        ? 'superpose.chains'
        : cmd.method === 'structure'
          ? 'superpose.tm-align'
          : 'superpose.atoms';
    const r = await ops.run(id, { targets, ...(cmd.method === 'sequence' ? { alignSequences: true } : {}) });
    return r.did;
  }

  /**
   * `pair_fit`'s pairs: Mol* pairs the atoms of each side in the order the structures hold them, so the pairs given
   * must be in that order on both sides, with as many atoms on each side of a pair.
   */
  #checkPairs(pairs: readonly { readonly mobile: Sel; readonly target: Sel }[]): void {
    const first = (l: StructureElement.Loci) => {
      const loc = StructureElement.Loci.getFirstLocation(l);
      return loc ? P.atom.sourceIndex(loc) : -1;
    };
    const sides = pairs.map((p, i) => {
      const m = this.#ops.some(p.mobile);
      const t = this.#ops.some(p.target);
      const mn = count(m);
      const tn = count(t);
      if (mn !== tn)
        throw new HostFailure(
          `Pair ${i + 1} has ${mn} mobile and ${tn} target atoms.`,
          'Give as many atoms on each side.',
        );
      return { m: first(m[0].loci), t: first(t[0].loci) };
    });
    const order = (key: 'm' | 't') =>
      sides
        .map((_, i) => i)
        .sort((a, b) => sides[a][key] - sides[b][key])
        .join();
    if (order('m') !== order('t')) {
      throw new HostFailure(
        'The pairs are not in the same order in both structures: Mol* pairs atoms in the order the structures hold them.',
        'Give the pairs in that order (by residue number on both sides).',
      );
    }
  }

  async #canvas(cmd: Extract<Command, { kind: 'canvas' }>): Promise<string> {
    const ops = this.#ops;
    const did: string[] = [];
    if (cmd.background !== undefined)
      did.push((await ops.run('canvas.background', { color: hex(cmd.background) })).did);
    if (cmd.projection !== undefined || cmd.fov !== undefined) {
      did.push(
        (
          await ops.run('canvas.camera', {
            ...(cmd.projection ? { mode: cmd.projection } : {}),
            ...(cmd.fov !== undefined ? { fov: Math.min(130, Math.max(10, cmd.fov)) } : {}),
          })
        ).did,
      );
    }
    if (cmd.fog !== undefined)
      did.push((await ops.run('canvas.lighting', { fog: cmd.fog ? { name: 'on', params: {} } : { name: 'off' } })).did);
    return did.join(' ');
  }

  #read(read: Extract<Command, { kind: 'read' }>['read']): string {
    const ops = this.#ops;
    switch (read.kind) {
      case 'count': {
        const where = ops.where(read.sel);
        const n = count(where);
        const per = where.length > 1 ? ` (${where.map((w) => `${w.object}: ${atomCount(w.loci)}`).join(', ')})` : '';
        return `${n}${per}`;
      }
      case 'names': {
        const objects = ops.objects().map((o) => o.name);
        const names = Object.keys(this.holder.value.names);
        return [...objects, ...names].join(' ') || '(none)';
      }
      case 'objects':
        return (
          ops
            .objects()
            .map((o) => `${o.molecule} ${o.name}`)
            .join('\n') || '(none)'
        );
      case 'view': {
        const camera = ops.plugin.canvas3d?.camera;
        if (!camera) return 'There is no canvas.';
        const s = camera.getSnapshot();
        const v = (x: Vec3) => [...x].map((n) => n.toFixed(2)).join(', ');
        return `position ${v(s.position)}; target ${v(s.target)}; up ${v(s.up)}; ${s.mode}, fov ${((s.fov * 180) / Math.PI).toFixed(0)}°`;
      }
      case 'values': {
        for (const p of read.properties) {
          if (!ReadableProperties.includes(p)) {
            throw new HostFailure(
              `The property '${p}' is not read here.`,
              `Properties: ${ReadableProperties.join(', ')}.`,
            );
          }
        }
        const lines: string[] = [];
        let total = 0;
        for (const w of ops.where(read.sel)) {
          total += atomCount(w.loci);
          for (const atom of eachAtom(w.loci, MaxValues - lines.length)) {
            const vs = read.properties.map((p) => atomProperty(atom.location, p));
            lines.push(read.properties.length > 1 ? `{${vs.join(' ')}}` : vs[0]);
          }
        }
        return `${lines.join(' ')}${total > lines.length ? ` … (${total - lines.length} more)` : ''}`;
      }
      case 'geometry': {
        if (read.measure === 'center' || read.measure === 'rgyr' || read.measure === 'minmax') {
          const where = ops.some(read.atoms[0]);
          const shape = shapeOf(where[0].loci);
          const v = (x: Vec3) => [...x].map((n) => n.toFixed(3)).join(' ');
          return read.measure === 'center'
            ? v(shape.center)
            : read.measure === 'rgyr'
              ? shape.rgyr.toFixed(3)
              : `{${v(shape.min)}} {${v(shape.max)}}`;
        }
        const { atoms } = this.#atoms(read.atoms);
        const points = atoms.map((w) => positionOf(w.loci)!);
        return geometryOf(read.measure, points).toFixed(3);
      }
    }
  }
}

/** A selection's target on one structure: its own MolQL there (no second parse). */
function target(sel: Sel, w: Where) {
  return { molql: selectionOn(sel, w.object) };
}

const count = (where: readonly Where[]) => where.reduce((n, w) => n + atomCount(w.loci), 0);

/** Whether a label's atom is in the selection (its cell's data, as Prism's measurement source reads it). */
function labelledIn(plugin: Ops['plugin'], ref: string, where: readonly Where[]): boolean {
  const data = plugin.state.data.cells.get(ref)?.obj?.data as
    | { sourceData?: { infos?: { loci: unknown }[] } }
    | undefined;
  const loci = data?.sourceData?.infos?.[0]?.loci;
  if (!StructureElement.Loci.is(loci)) return false;
  return where.some((w) => StructureElement.Loci.areIntersecting(StructureElement.Loci.remap(loci, w.data), w.loci));
}
