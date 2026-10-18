/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { QueryLanguageEntry } from '../languages.js';
import { JmolLanguage } from './jmol.js';
import { MolScriptLanguage } from './molscript.js';
import { PyMolLanguage } from './pymol.js';
import { VmdLanguage } from './vmd.js';

/** Every query language, PyMOL first (the default for text): what `DefaultPrismViewerSpec` lists. */
export const AllQueryLanguages: readonly QueryLanguageEntry[] = [
  PyMolLanguage,
  VmdLanguage,
  JmolLanguage,
  MolScriptLanguage,
];
