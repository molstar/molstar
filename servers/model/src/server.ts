#!/usr/bin/env node
/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import type * as http from 'http';
import compression from 'compression';
import express from 'express';
import { ConsoleLogger } from '@molstar/core/util/console-logger';
import { PerformanceMonitor } from '@molstar/core/util/performance-monitor';
import { configureServer, ModelServerConfig as ServerConfig } from '@molstar/model-server/config';
import { drainWebApi, initWebApi } from '@molstar/model-server/server/api-web';
import { VERSION } from '@molstar/model-server/version';

let shuttingDown = false;

/**
 * Responses not sent completely. A response is still being sent after its job is done (what the client has not taken
 * yet is buffered), and Node counts a connection idle once `end` was called: neither an empty queue nor `close`
 * (`closeIdleConnections`) waits for it. 'finish' is emitted once all of it is handed to the OS.
 */
const openResponses = new Set<http.ServerResponse>();
let onResponsesSent: (() => void) | undefined;

function trackResponses(app: express.Express) {
  app.use((_req, res, next) => {
    openResponses.add(res);
    const done = () => {
      openResponses.delete(res);
      if (openResponses.size === 0) onResponsesSent?.();
    };
    res.once('finish', done);
    res.once('close', done);
    next();
  });
}

function responsesSent() {
  if (openResponses.size === 0) return Promise.resolve();
  return new Promise<void>((resolve) => (onResponsesSent = resolve));
}

/**
 * Stops accepting requests, lets the queued and running ones finish and their responses be sent (for at most
 * `requestTimeoutMs`), then exits. A second signal exits right away.
 */
function shutdown(reason: string) {
  if (shuttingDown) {
    ConsoleLogger.log('Server', `${reason} again, exiting now.`);
    process.exit(1);
  }
  shuttingDown = true;
  ConsoleLogger.log('Server', `Shutting down (${reason}), finishing the requests in flight.`);

  setTimeout(() => {
    ConsoleLogger.log('Server', 'Requests still running after the grace period, exiting.');
    process.exit(0);
  }, ServerConfig.requestTimeoutMs);

  // server keeps listening: new queries are answered 503 and `health-check` reports unhealthy meanwhile
  drainWebApi()
    .then(responsesSent)
    .then(() => {
      ConsoleLogger.log('Server', 'Shut down.');
      process.exit(0);
    });
}

function setupShutdown() {
  if (ServerConfig.shutdownTimeoutVarianceMinutes > ServerConfig.shutdownTimeoutMinutes) {
    ConsoleLogger.log('Server', 'Shutdown timeout variance is greater than the timer itself, ignoring.');
  } else {
    let tVar = 0;
    if (ServerConfig.shutdownTimeoutVarianceMinutes > 0) {
      tVar = 2 * (Math.random() - 0.5) * ServerConfig.shutdownTimeoutVarianceMinutes;
    }
    const tMs = (ServerConfig.shutdownTimeoutMinutes + tVar) * 60 * 1000;

    console.log(`----------------------------------------------------------------------------`);
    console.log(`  The server will shut down in ${PerformanceMonitor.format(tMs)} to prevent slow performance.`);
    console.log(`  Please make sure a daemon is running that will automatically restart it.`);
    console.log(`----------------------------------------------------------------------------`);
    console.log();

    setTimeout(() => shutdown('timeout'), tMs).unref();
  }
}

configureServer();

function startServer() {
  const app = express();
  trackResponses(app);
  app.use(
    compression({
      level: 6,
      memLevel: 9,
      chunkSize: 16 * 16384,
      filter: (req, res) => {
        const ct = res.getHeader('Content-Type');
        if (typeof ct === 'string' && ct.indexOf('tar+gzip') > 0) return false;
        return true;
      },
    }),
  );

  initWebApi(app);

  app.get('/robots.txt', function (req, res) {
    res.type('text/plain').send(ServerConfig.robots);
  });

  const port = process.env.port || ServerConfig.defaultPort;
  app.listen(port).setTimeout(ServerConfig.requestTimeoutMs);

  console.log(`Mol* ModelServer ${VERSION}`);
  console.log(``);
  console.log(`The server is running on port ${port}.`);
  console.log(``);
}

startServer();

process.on('SIGTERM', () => shutdown('SIGTERM'));
process.on('SIGINT', () => shutdown('SIGINT'));

if (ServerConfig.shutdownTimeoutMinutes > 0) {
  setupShutdown();
}
