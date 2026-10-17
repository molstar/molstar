/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The guard harness for the viewer's tiers (packages/prism/AGENTS.md, Tests): prism-ui's (Seq*'s `quickCase` and heap
// scaling, and `retainedScaling` for steps that allocate by design). Test code only, by path as the Mol* adapter's
// guards reach Seq*'s.

export * from '../../../ui/src/_test/perf.js';
