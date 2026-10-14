/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The platform-free entry, the host's contract: `materialize`, its limit, the supported set and the description's types,
// and `extentLoci` for an event's extent. `mount`, the only importer of Nightingale, is the
// `@molstar/seq-nightingale/mount` entry, so this one loads in Node.

export { type AxisBlock, type AxisMeta, type AxisRuns, type BlockRun, extentLoci } from "./axis.js";
export { type ElementDesc, NightingaleSupported } from "./element.js";
export { materialize, NIGHTINGALE_MAX_POSITIONS } from "./materialize.js";
