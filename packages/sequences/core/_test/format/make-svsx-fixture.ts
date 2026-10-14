/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// By hand: `node lib/seqstar/core/_test/format/make-svsx-fixture.ts` rewrites fixtures/kinase.svsx from svsx-fixture.ts.
import { writeFileSync } from "node:fs";
import { bundleUris, writeSvsx } from "@molstar/seq-core/builder/svsx";
import { absoluteDoc, fixtureFetch } from "./svsx-fixture.js";

const { doc, files } = await bundleUris(absoluteDoc(), fixtureFetch().fetch);
writeFileSync(new URL("fixtures/kinase.svsx", import.meta.url), await writeSvsx(doc, files));
