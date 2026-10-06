import assert from 'node:assert/strict';
import { chmod, mkdtemp, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import spawn from 'cross-spawn';

test('tool commands preserve paths and arguments through platform wrappers', async () => {
  const dir = await mkdtemp(path.join(tmpdir(), 'molstar process '));
  try {
    const script = path.join(dir, 'fixture.mjs');
    await writeFile(script, '#!/usr/bin/env node\nconsole.log(JSON.stringify(process.argv.slice(2)));\n');
    await chmod(script, 0o755);
    let command = script;
    if (process.platform === 'win32') {
      command = path.join(dir, 'fixture.cmd');
      await writeFile(command, `@echo off\r\n"${process.execPath}" "%~dp0fixture.mjs" %*\r\n`);
    }
    const args = ['a file & more.txt', 'literal $value', '(parentheses)'];
    const synchronous = spawn.sync(command, args, { cwd: dir, encoding: 'utf8' });
    assert.ifError(synchronous.error);
    assert.equal(synchronous.status, 0, synchronous.stderr);
    assert.deepEqual(JSON.parse(synchronous.stdout), args);
    const asynchronous = await new Promise((resolve, reject) => {
      const child = spawn(command, args, { cwd: dir, stdio: ['ignore', 'pipe', 'pipe'] });
      let stdout = '',
        stderr = '';
      child.stdout.setEncoding('utf8').on('data', (chunk) => (stdout += chunk));
      child.stderr.setEncoding('utf8').on('data', (chunk) => (stderr += chunk));
      child.once('error', reject);
      child.once('close', (code) => (code === 0 ? resolve(stdout) : reject(new Error(stderr))));
    });
    assert.deepEqual(JSON.parse(asynchronous), args);
  } finally {
    await rm(dir, { recursive: true, force: true });
  }
});
