/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Expression } from '@molstar/query-language/language/expression';
import type { Span } from './diagnostics.js';

// The command IR (_spec §4): what to do in Mol* terms, whoever typed it. Both dialects compile to it; the host runs it.
// JSON-safe (selections are MolQL), so a run can be recorded and replayed without a parser. One person's command may
// be several of these (`as` is a hide and a show), never the reverse.

export type Dialect = 'pymol' | 'vmd';

export const Dialects: readonly Dialect[] = ['pymol', 'vmd'];

export const DialectNames: Readonly<Record<Dialect, string>> = { pymol: 'PyMOL', vmd: 'VMD' };

/** One part of a selection: MolQL evaluated on one structure (_spec §4.3). */
export interface SelPart {
  /** A named object's own part, or `others`: every structure the selection names no part for. */
  readonly scope: { readonly object: string } | 'others';
  readonly molql: Expression;
}

/** A selection, compiled: its text (for messages) and its parts per structure. */
export interface Sel {
  readonly text: string;
  readonly dialect: Dialect;
  readonly parts: readonly SelPart[];
}

/**
 * What a representation looks like, in the people's words; the host maps each to a Mol* representation type and its
 * params (_spec §5.2).
 */
export type RepType =
  | 'cartoon'
  | 'putty'
  | 'tube'
  | 'backbone'
  | 'sticks'
  | 'ball-and-stick'
  | 'lines'
  | 'spheres'
  | 'surface'
  | 'mesh'
  | 'gaussian-surface'
  | 'points'
  | 'ellipsoids';

export const RepTypes: readonly RepType[] = [
  'cartoon',
  'putty',
  'tube',
  'backbone',
  'sticks',
  'ball-and-stick',
  'lines',
  'spheres',
  'surface',
  'mesh',
  'gaussian-surface',
  'points',
  'ellipsoids',
];

/**
 * A color: RGB as `0xRRGGBB`, or a Mol* color theme by name (`chain-id`, `sequence-id`, …). Its params are the
 * language's, which the host turns into the theme's: `element-symbol` takes `carbon` (RGB: PyMOL's `util.cbag`) and
 * `palette: 'vmd'` (VMD's Name colors); a theme's own `palette` (`rainbow`, `blue_white_red`) is PyMOL's `spectrum`.
 */
export type ColorSpec =
  | { readonly rgb: number }
  | { readonly theme: string; readonly params?: Readonly<Record<string, unknown>> };

/** A VMD representation: style, coloring, selection, material (_spec §5.3). */
export interface VmdRep {
  readonly style: RepType;
  readonly color: ColorSpec;
  /** The selection text, VMD's; compiled per molecule when the rep is made. */
  readonly selection: string;
  readonly material: 'opaque' | 'transparent';
}

export type LoadSource =
  | { readonly entry: { readonly db: 'pdb' | 'afdb'; readonly id: string } }
  | { readonly url: string }
  /** A file given with the submission, by name. */
  | { readonly file: string };

/** A label's text: literals and atom properties, joined. */
export type LabelPart =
  | { readonly text: string }
  | { readonly property: 'resn' | 'resi' | 'chain' | 'name' | 'elem' | 'b' | 'index' };

export type Read =
  | { readonly kind: 'count'; readonly sel: Sel }
  | { readonly kind: 'names' }
  | { readonly kind: 'objects' }
  | { readonly kind: 'view' }
  | { readonly kind: 'values'; readonly sel: Sel; readonly properties: readonly string[] }
  | {
      readonly kind: 'geometry';
      readonly measure: 'distance' | 'angle' | 'dihedral' | 'center' | 'rgyr' | 'minmax';
      readonly atoms: readonly Sel[];
    };

export type CommandBody =
  | {
      readonly kind: 'load';
      readonly from: LoadSource;
      /** The object's name (`fetch 1cbs` → `1cbs`), given or made from the source. */
      readonly name: string;
      /** Its VMD molecule id: the next one (ids are never reused, as in VMD), in either dialect. */
      readonly molecule: number;
      readonly format?: string;
      /** VMD's dialect loads with VMD's look (one Lines rep), PyMOL's with Mol*'s preset (_spec §5.3). */
      readonly look: 'default' | 'vmd-lines';
      /** Another copy when it is loaded already (VMD's `mol new`, a new object name); else the loaded one stays. */
      readonly copy: boolean;
    }
  | { readonly kind: 'object.remove'; readonly objects: readonly string[] | 'all' }
  | { readonly kind: 'object.visibility'; readonly objects: readonly string[] | 'all'; readonly visible: boolean }
  /** VMD's top molecule. */
  | { readonly kind: 'object.top'; readonly object: string }
  /** A named selection (PyMOL `select name, …`, VMD `atomselect macro`). */
  | { readonly kind: 'name.set'; readonly name: string; readonly sel: Sel }
  | { readonly kind: 'name.remove'; readonly names: readonly string[] | 'all' }
  /** The interactive selection. */
  | { readonly kind: 'select'; readonly mode: 'set' | 'add' | 'remove'; readonly sel: Sel }
  | { readonly kind: 'select.clear' }
  /** PyMOL's per-atom show flags (_spec §5.2); `only` is `as`. */
  | {
      readonly kind: 'display';
      readonly mode: 'show' | 'hide' | 'only';
      readonly reps: readonly RepType[] | 'all';
      readonly sel: Sel;
    }
  /** PyMOL `cartoon putty|tube|automatic`: what cartoon-flagged atoms in the selection show as. */
  | { readonly kind: 'cartoon'; readonly style: 'automatic' | 'putty' | 'tube'; readonly sel: Sel }
  /** VMD's reps of a molecule (_spec §5.3). */
  | { readonly kind: 'rep.add'; readonly object: string; readonly rep: VmdRep }
  | {
      readonly kind: 'rep.modify';
      readonly object: string;
      readonly index: number;
      readonly patch: Partial<VmdRep>;
    }
  | { readonly kind: 'rep.remove'; readonly object: string; readonly index: number }
  | { readonly kind: 'rep.visibility'; readonly object: string; readonly index: number; readonly visible: boolean }
  /** `mol representation|color|selection|material`: what the next `mol addrep` takes. */
  | { readonly kind: 'rep.defaults'; readonly patch: Partial<VmdRep> }
  /** `optional`: one of several (a color per element): selecting nothing is no failure. */
  | { readonly kind: 'color'; readonly sel: Sel; readonly color: ColorSpec; readonly optional?: true }
  /** The selection selects atoms, else the command fails (before a group of optional commands). */
  | { readonly kind: 'require'; readonly sel: Sel }
  | { readonly kind: 'transparency'; readonly sel: Sel; readonly alpha: number; readonly reps?: readonly RepType[] }
  | { readonly kind: 'size'; readonly sel: Sel; readonly factor: number; readonly reps?: readonly RepType[] }
  /** PyMOL `remove`: whole entities or chains left out of every component (`contents.visibility`). */
  | { readonly kind: 'hide-contents'; readonly sel: Sel }
  | { readonly kind: 'hydrogens'; readonly visible: boolean }
  | {
      readonly kind: 'camera';
      readonly action: 'frame' | 'orient' | 'center' | 'reset';
      readonly sel?: Sel;
      readonly animate: boolean;
      /** Ångströms added around what is framed (PyMOL's `zoom` buffer). */
      readonly buffer?: number;
    }
  | { readonly kind: 'camera.zoom'; readonly factor: number }
  | { readonly kind: 'camera.turn'; readonly axis: 'x' | 'y' | 'z'; readonly angle: number }
  | { readonly kind: 'camera.move'; readonly axis: 'x' | 'y' | 'z'; readonly distance: number }
  | {
      readonly kind: 'measure';
      readonly measure: 'distance' | 'angle' | 'dihedral';
      /** One atom each; with a `cutoff` (a distance), sets: every pair within it (PyMOL's mode 0). */
      readonly atoms: readonly Sel[];
      readonly cutoff?: number;
    }
  | { readonly kind: 'label'; readonly sel: Sel; readonly text: readonly LabelPart[] }
  | { readonly kind: 'label.clear'; readonly sel?: Sel }
  /** VMD `label delete Bonds|Angles|Dihedrals`: the measurements of a kind. */
  | { readonly kind: 'measure.clear'; readonly measure: 'distance' | 'angle' | 'dihedral' }
  | {
      readonly kind: 'superpose';
      readonly method: 'sequence' | 'structure' | 'atoms';
      readonly mobile: Sel;
      readonly target: Sel;
      /** `atoms`: the pairs as given (PyMOL's `pair_fit`), which the structures' order must keep. */
      readonly pairs?: readonly { readonly mobile: Sel; readonly target: Sel }[];
    }
  | {
      readonly kind: 'canvas';
      /** `0xRRGGBB`. */
      readonly background?: number;
      readonly projection?: 'perspective' | 'orthographic';
      /** Degrees. */
      readonly fov?: number;
      readonly fog?: boolean;
    }
  /** A frame (0-based; -1 is the last), or playing. */
  | { readonly kind: 'trajectory'; readonly frame?: number; readonly play?: 'forward' | 'stop' }
  | {
      readonly kind: 'export.image';
      readonly name?: string;
      readonly width?: number;
      readonly height?: number;
      readonly transparent?: boolean;
    }
  | { readonly kind: 'history'; readonly action: 'undo' | 'redo' }
  | { readonly kind: 'read'; readonly read: Read }
  | { readonly kind: 'print'; readonly text: string };

/**
 * A command and where it came from: its place in the submission and its text, for messages and the record; a command
 * of an included script also its file and line there (its place is the include's).
 */
export type Command = CommandBody & {
  readonly at: Span;
  readonly source: string;
  readonly file?: { readonly name: string; readonly line: number };
};
