/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * setImmediate polyfill adapted from https://github.com/YuzuJS/setImmediate
 * Copyright (c) 2012 Barnesandnoble.com, llc, Donavon West, and Domenic Denicola
 * MIT license.
 */

function createActions() {
    type Callback = (...args: any[]) => void;
    type Task = { callback: Callback, args: any[] }

    const tasksByHandle: { [handle: number]: Task } = { };
    const doc = typeof document !== 'undefined' ? document : void 0;

    let currentlyRunningATask = false;
    let nextHandle = 1; // Spec says greater than zero
    let registerImmediate: ((handle: number) => void);

    function setImmediate(callback: Callback, ...args: any[]) {
      // Callback can either be a function or a string
      if (typeof callback !== 'function') {
        callback = new Function('' + callback) as Callback;
      }
      // Store and register the task
      const task = { callback: callback, args: args };
      tasksByHandle[nextHandle] = task;
      registerImmediate(nextHandle);
      return nextHandle++;
    }

    function clearImmediate(handle: number) {
        delete tasksByHandle[handle];
    }

    function run(task: Task) {
        const callback = task.callback;
        const args = task.args;
        switch (args.length) {
        case 0:
            callback();
            break;
        case 1:
            callback(args[0]);
            break;
        case 2:
            callback(args[0], args[1]);
            break;
        case 3:
            callback(args[0], args[1], args[2]);
            break;
        default:
            callback.apply(undefined, args);
            break;
        }
    }

    function runIfPresent(handle: number) {
        // From the spec: 'Wait until any invocations of this algorithm started before this one have completed.'
        // So if we're currently running a task, we'll need to delay this invocation.
        if (currentlyRunningATask) {
            // Delay by doing a setTimeout. setImmediate was tried instead, but in Firefox 7 it generated a
            // 'too much recursion' error.
            setTimeout(runIfPresent, 0, handle);
        } else {
            const task = tasksByHandle[handle];
            if (task) {
                currentlyRunningATask = true;
                try {
                    run(task);
                } finally {
                    clearImmediate(handle);
                    currentlyRunningATask = false;
                }
            }
        }
    }

    function installNextTickImplementation() {
        registerImmediate = function(handle) {
            process.nextTick(function () { runIfPresent(handle); });
        };
    }

    function canUsePostMessage() {
        // The test against `importScripts` prevents this implementation from being installed inside a web worker,
        // where `global.postMessage` means something completely different and can't be used for this purpose.
        const global = window as any;
        if (global && global.postMessage && !global.importScripts) {
            let postMessageIsAsynchronous = true;
            const oldOnMessage = global.onmessage;
            global.onmessage = function() {
                postMessageIsAsynchronous = false;
            };
            global.postMessage('', '*');
            global.onmessage = oldOnMessage;
            return postMessageIsAsynchronous;
        }
    }

    function installPostMessageImplementation() {
        // Installs an event handler on `global` for the `message` event: see
        // * https://developer.mozilla.org/en/DOM/window.postMessage
        // * http://www.whatwg.org/specs/web-apps/current-work/multipage/comms.html#crossDocumentMessages

        const messagePrefix = 'setImmediate$' + Math.random() + '$';
        const onGlobalMessage = function(event: any) {
            if (event.source === global &&
                typeof event.data === 'string' &&
                event.data.indexOf(messagePrefix) === 0) {
                runIfPresent(+event.data.slice(messagePrefix.length));
            }
        };

        if (window.addEventListener) {
            window.addEventListener('message', onGlobalMessage, false);
        } else {
            (window as any).attachEvent('onmessage', onGlobalMessage);
        }

        registerImmediate = function(handle) {
            window.postMessage(messagePrefix + handle, '*');
        };
    }

    function installMessageChannelImplementation() {
        const channel = new MessageChannel();
        channel.port1.onmessage = function(event) {
            const handle = event.data;
            runIfPresent(handle);
        };

        registerImmediate = function(handle) {
            channel.port2.postMessage(handle);
        };
    }

    function installReadyStateChangeImplementation() {
        const html = doc!.documentElement;
        registerImmediate = function(handle) {
            // Create a <script> element; its readystatechange event will be fired asynchronously once it is inserted
            // into the document. Do so, thus queuing up the task. Remember to clean up once it's been called.
            let script = doc!.createElement('script') as any;
            script.onreadystatechange = function () {
                runIfPresent(handle);
                script.onreadystatechange = null;
                html.removeChild(script);
                script = null;
            };
            html.appendChild(script);
        };
    }

    function installSetTimeoutImplementation() {
        registerImmediate = function(handle) {
            setTimeout(runIfPresent, 0, handle);
        };
    }

    // If supported, we should attach to the prototype of global, since that is where setTimeout et al. live.
    //const attachTo = Object.getPrototypeOf && Object.getPrototypeOf(global);
    //attachTo = attachTo && attachTo.setTimeout ? attachTo : global;

    // Don't get fooled by e.g. browserify environments.
    if (typeof process !== 'undefined' && {}.toString.call(process) === '[object process]') {
        // For Node.js before 0.9
        installNextTickImplementation();
    } else if (canUsePostMessage()) {
        // For non-IE10 modern browsers
        installPostMessageImplementation();
    } else if (typeof MessageChannel !== 'undefined') {
        // For web workers, where supported
        installMessageChannelImplementation();
    } else if (doc && 'onreadystatechange' in doc.createElement('script')) {
        // For IE 6–8
        installReadyStateChangeImplementation();
    } else {
        // For older browsers
        installSetTimeoutImplementation();
    }

    return {
        setImmediate,
        clearImmediate
    };
}

const actions = (function () {
    if (typeof setImmediate !== 'undefined') {
        return { setImmediate, clearImmediate };
    }
    return createActions();
}());

function resolveImmediate(res: () => void) {
    actions.setImmediate(res);
}

export default {
    immediate: actions.setImmediate,
    clearImmediate: actions.clearImmediate,

    immediatePromise() { return new Promise<void>(resolveImmediate); }
};
