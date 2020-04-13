/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * setImmediate polyfill adapted from https://github.com/YuzuJS/setImmediate
 * Copyright (c) 2012 Barnesandnoble.com, llc, Donavon West, and Domenic Denicola
 * MIT license.
 */

declare const WorkerGlobalScope: any;
function createImmediateActions() {
    const global: any = (function () {
        const _window = typeof window !== 'undefined' && window;
        const _self = typeof self !== 'undefined' && typeof WorkerGlobalScope !== 'undefined' && self instanceof WorkerGlobalScope && self;
        const _global = typeof global !== 'undefined' && global;
        return _window || _global || _self;
    })();

    type Callback = (...args: any[]) => void;
    type Task = { callback: Callback, args: any[] }

    const tasksByHandle: { [handle: number]: Task } = { };
    const doc = typeof document !== 'undefined' ? document : void 0;

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
        const task = tasksByHandle[handle];
        clearImmediate(handle);
        run(task);
    }

    function installNextTickImplementation() {
        registerImmediate = function(handle) {
            process.nextTick(function () { runIfPresent(handle); });
        };
    }

    function canUsePostMessage() {
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
        const html = doc!.documentElement!;
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

const immediateActions = (function () {
    if (typeof setImmediate !== 'undefined') {
        if (typeof window !== 'undefined') {
            return {
                setImmediate: (handler: any, ...args: any[]) => (window as any).setImmediate(handler, ...args as any) as number,
                clearImmediate: (handle: any) => (window as any).clearImmediate(handle)
            };
        } else {
            return { setImmediate, clearImmediate };
        }
    }
    return createImmediateActions();
}());

function resolveImmediate(res: () => void) {
    immediateActions.setImmediate(res);
}

const Scheduler = {
    setImmediate: immediateActions.setImmediate,
    clearImmediate: immediateActions.clearImmediate,
    immediatePromise() { return new Promise<void>(resolveImmediate); },

    delay<T>(timeout: number, value: T | undefined = void 0): Promise<T> { return new Promise(r => setTimeout(r, timeout, value)); }
};

export { Scheduler };
