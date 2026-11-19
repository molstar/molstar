#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * seqx-render: a SeqViewSpec file to an SVG or a PNG (README.md). Run: node cli/seqx-render/lib/seqx-render.js <args>
 */

import { run } from './cli.js';

process.exitCode = await run(process.argv.slice(2));
