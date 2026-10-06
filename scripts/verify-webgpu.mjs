/** Native WebGPU verification in an installed Chrome browser. Run with Bun. */
import { build } from 'esbuild';
import { chromium } from '@playwright/test';
import { createServer } from 'node:http';
import { readFile, mkdir } from 'node:fs/promises';

async function verifyInteractions(page, scale) {
    await page.evaluate(() => window.webgpuPrepareInteraction());
    try {
        const before = await page.evaluate(() => window.webgpuInteraction.snapshot());
        await page.mouse.move(180, 130);
        await page.mouse.down();
        await page.mouse.move(230, 165, { steps: 12 });
        await page.mouse.up();
        await page.waitForFunction(async previous => {
            const state = await window.webgpuInteraction.snapshot();
            return state.position.some((v, i) => Math.abs(v - previous.position[i]) > 0.1) && state.frameHash !== previous.frameHash;
        }, before);
        const rotated = await page.evaluate(() => window.webgpuInteraction.snapshot());
        await page.mouse.move(180, 130);
        await page.mouse.down({ button: 'right' });
        await page.mouse.move(200, 140, { steps: 8 });
        await page.mouse.up({ button: 'right' });
        await page.waitForFunction(async previous => {
            const state = await window.webgpuInteraction.snapshot();
            return state.target.some((v, i) => Math.abs(v - previous.target[i]) > 0.05) && state.frameHash !== previous.frameHash;
        }, rotated);
        const panned = await page.evaluate(() => window.webgpuInteraction.snapshot());
        if (Math.abs(panned.distance - rotated.distance) > 0.001) throw new Error('Panning must preserve the camera distance to its target.');
        await page.mouse.move(180, 130);
        await page.mouse.wheel(0, -240);
        await page.waitForFunction(async distance => (await window.webgpuInteraction.snapshot()).distance < distance - 0.1, panned.distance);
        const zoomed = await page.evaluate(() => window.webgpuInteraction.snapshot());
        if (zoomed.frameHash === panned.frameHash) throw new Error('Mouse wheel zoom must change the molecular frame.');
        await page.mouse.move(800, 600);
        await page.waitForFunction(async () => !(await window.webgpuInteraction.snapshot()).hover);
        const unmarked = await page.evaluate(() => window.webgpuInteraction.snapshot());
        const target = await page.evaluate(() => window.webgpuInteraction.target());
        await page.mouse.move(target.x, target.y);
        await page.waitForFunction(async () => (await window.webgpuInteraction.snapshot()).hover);
        const hovered = await page.evaluate(() => window.webgpuInteraction.snapshot());
        if (hovered.frameHash === unmarked.frameHash) throw new Error('Atom hover must visibly highlight the molecular frame.');
        await page.evaluate(() => window.webgpuInteraction.setSelectionMode(true));
        await page.mouse.click(target.x, target.y);
        await page.waitForFunction(async () => {
            const state = await window.webgpuInteraction.snapshot();
            return state.clicks > 0 && state.selectedAtoms === 1;
        });
        await page.mouse.move(800, 600);
        await page.waitForFunction(async () => !(await window.webgpuInteraction.snapshot()).hover);
        const selected = await page.evaluate(() => window.webgpuInteraction.snapshot());
        if (selected.frameHash === unmarked.frameHash) throw new Error('Atom selection must visibly mark the molecular frame after the cursor leaves.');
        await page.mouse.click(target.x, target.y);
        await page.waitForFunction(async () => (await window.webgpuInteraction.snapshot()).selectedAtoms === 0);
        await page.mouse.click(target.x, target.y);
        await page.waitForFunction(async () => (await window.webgpuInteraction.snapshot()).selectedAtoms === 1);
        const emptyTarget = await page.evaluate(() => window.webgpuInteraction.emptyTarget());
        await page.mouse.click(emptyTarget.x, emptyTarget.y);
        await page.waitForFunction(async () => (await window.webgpuInteraction.snapshot()).selectedAtoms === 0);
        await page.mouse.move(800, 600);
        await page.waitForFunction(async () => !(await window.webgpuInteraction.snapshot()).hover);
        await page.waitForFunction(async hash => (await window.webgpuInteraction.snapshot()).frameHash === hash, unmarked.frameHash);
        const final = await page.evaluate(() => window.webgpuInteraction.snapshot());
        if (final.errors.length) throw new Error(final.errors.join('\n'));
        await page.screenshot({ path: `tmp/webgpu/interaction-${scale}x.png` });
        return `trusted browser rotation/panning, wheel zoom, asynchronous atom hover, single-atom selection/toggle, empty-background deselection and visible marking at ${scale}x display scaling`;
    } finally {
        await page.evaluate(() => window.webgpuInteraction.dispose());
    }
}

async function main() {
await mkdir('tmp/webgpu', { recursive: true });
await build({ entryPoints: ['src/tests/browser/webgpu.ts'], outfile: 'tmp/webgpu/verify.js', bundle: true, platform: 'browser', sourcemap: 'inline' });
const script = await readFile('tmp/webgpu/verify.js');
const protein = await readFile('examples/1crn.cif');
const server = createServer((request, response) => {
    if (request.url === '/fixtures/1crn.cif') {
        response.setHeader('Content-Type', 'text/plain'); response.end(protein);
    } else if (request.url === '/verify.js') {
        response.setHeader('Content-Type', 'text/javascript'); response.end(script);
    } else {
        response.setHeader('Content-Type', 'text/html');
        response.end('<!doctype html><html><body style="background:#222"><canvas width="256" height="256"></canvas><script src="/verify.js"></script></body></html>');
    }
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
let browser;
try {
    browser = await chromium.launch({ channel: 'chrome', headless: true, args: ['--enable-unsafe-webgpu'] });
    const page = await browser.newPage();
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.goto(`http://127.0.0.1:${server.address().port}`);
    const result = await page.evaluate(() => window.webgpuVerification);
    result.results.push(await verifyInteractions(page, 1));
    const highDensityPage = await browser.newPage({ deviceScaleFactor: 2 });
    highDensityPage.on('pageerror', error => errors.push(error.message));
    try {
        await highDensityPage.goto(`http://127.0.0.1:${server.address().port}/?interaction-only=1`);
        result.results.push(await verifyInteractions(highDensityPage, 2));
    } finally { await highDensityPage.close(); }
    if (errors.length) throw new Error(errors.join('\n'));
    await page.screenshot({ path: 'tmp/webgpu/verification.png', fullPage: true });
    console.log(JSON.stringify(result, null, 2));
} finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
}

}

main().catch(error => { console.error(error); process.exitCode = 1; });
