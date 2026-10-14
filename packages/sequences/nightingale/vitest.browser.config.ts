/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { playwright } from "@vitest/browser-playwright";
import { defineConfig } from "vitest/config";
import { browserTests } from "../vitest.shared.ts";

export default defineConfig(browserTests(playwright));
