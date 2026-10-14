#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `svs-print-schema`: the reference, generated from the parameter tables and the kernel's and MVS's registries, as
// Markdown on stdout. Runs with Node's type stripping.

import { createRegistries, KernelRegistry, MvsRegistry, printReference, register } from '../index.js';

const reg = createRegistries();
register(reg, [KernelRegistry, MvsRegistry]);
process.stdout.write(printReference(reg));
