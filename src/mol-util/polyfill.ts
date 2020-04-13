/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

if (typeof window !== 'undefined') {
    (function () {
        // Console-polyfill. MIT license.
        // https://github.com/paulmillr/console-polyfill
        // Make it safe to do console.log() always.

        (window as any).console = window.console || {};
        let con = window.console;
        let prop, method;
        let empty = {};
        let dummy = function () {};
        let properties = 'memory'.split(',');
        let methods = (
            'assert,clear,count,debug,dir,dirxml,error,exception,group,' +
        'groupCollapsed,groupEnd,info,log,markTimeline,profile,profiles,profileEnd,' +
        'show,table,time,timeEnd,timeline,timelineEnd,timeStamp,trace,warn'
        ).split(',');

        while ((prop = properties.pop())) if (!(con as any)[prop]) (con as any)[prop] = empty;
        while ((method = methods.pop())) if (!(con as any)[method]) (con as any)[method] = dummy;
    })();
}

if (typeof window.HTMLCanvasElement !== 'undefined' && !window.HTMLCanvasElement.prototype.toBlob) {
    // http://code.google.com/p/chromium/issues/detail?id=67587#57
    Object.defineProperty(window.HTMLCanvasElement.prototype, 'toBlob', {

        value: function (callback: any, type: any, quality: any) {
            let bin = window.atob(this.toDataURL(type, quality).split(',')[ 1 ]);
            let len = bin.length;
            let len32 = len >> 2;
            let a8 = new Uint8Array(len);
            let a32 = new Uint32Array(a8.buffer, 0, len32);

            let j = 0;
            for (let i = 0; i < len32; i++) {
                a32[i] = (
                    bin.charCodeAt(j++) |
                    bin.charCodeAt(j++) << 8 |
                    bin.charCodeAt(j++) << 16 |
                    bin.charCodeAt(j++) << 24
                );
            }

            let tailLength = len & 3;

            while (tailLength--) {
                a8[ j ] = bin.charCodeAt(j++);
            }

            callback(new window.Blob([ a8 ], { 'type': type || 'image/png' }));
        }

    });
}

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/cbrt
Math.cbrt = Math.cbrt || function (x) {
    let y = Math.pow(Math.abs(x), 1 / 3);
    return x < 0 ? -y : y;
};

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/sign
if (!Math.sign) {
    Math.sign = function (x) {
        // If x is NaN, the result is NaN.
        // If x is -0, the result is -0.
        // If x is +0, the result is +0.
        // If x is negative and not -0, the result is -1.
        // If x is positive and not +0, the result is +1.
        x = +x; // convert to a number
        if (x === 0 || isNaN(x)) {
            return Number(x);
        }
        return x > 0 ? 1 : -1;
    };
}

if (!Number.isInteger) {
    // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Number/isInteger
    Number.isInteger = function isInteger (nVal) {
        return typeof nVal === 'number' && isFinite(nVal) && nVal > -9007199254740992 && nVal < 9007199254740992 && Math.floor(nVal) === nVal;
    };
}

if (!Number.isNaN) {
    // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Number/isNaN
    Number.isNaN = function isNaN (value) {
        return value !== value; // eslint-disable-line no-self-compare
    };
}

if (!Object.assign) {
    // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Object/assign
    Object.defineProperty(Object, 'assign', {

        enumerable: false,
        configurable: true,
        writable: true,

        value: function (target: any/* , firstSource */) {
            'use strict';
            if (target === undefined || target === null) { throw new TypeError('Cannot convert first argument to object'); }

            let to = Object(target);

            let hasPendingException = false;
            let pendingException;

            for (let i = 1; i < arguments.length; i++) {
                let nextSource = arguments[i];
                if (nextSource === undefined || nextSource === null) { continue; }

                let keysArray = Object.keys(Object(nextSource));
                for (let nextIndex = 0, len = keysArray.length; nextIndex < len; nextIndex++) {
                    let nextKey = keysArray[nextIndex];
                    try {
                        let desc = Object.getOwnPropertyDescriptor(nextSource, nextKey);
                        if (desc !== undefined && desc.enumerable) { to[nextKey] = nextSource[nextKey]; }
                    } catch (e) {
                        if (!hasPendingException) {
                            hasPendingException = true;
                            pendingException = e;
                        }
                    }
                }

                if (hasPendingException) { throw pendingException; }
            }

            return to;
        }

    });
}

if (!String.prototype.startsWith) {
    /* ! https://mths.be/startswith v0.2.0 by @mathias */

    (function () {
        'use strict'; // needed to support `apply`/`call` with `undefined`/`null`
        let defineProperty = (function () {
        // IE 8 only supports `Object.defineProperty` on DOM elements
            let result;
            try {
                let object = {};
                let $defineProperty = Object.defineProperty;
                result = $defineProperty(object, object as any, object) && $defineProperty;
            } catch (error) {} // eslint-disable-line no-empty
            return result;
        }());
        let toString = {}.toString;
        let startsWith = function (this: any, search: any) {
            if (this === null) {
                throw TypeError();
            }
            let string = String(this);
            if (search && toString.call(search) === '[object RegExp]') {
                throw TypeError();
            }
            let stringLength = string.length;
            let searchString = String(search);
            let searchLength = searchString.length;
            let position = arguments.length > 1 ? arguments[1] : undefined;
            // `ToInteger`
            let pos = position ? Number(position) : 0;
            if (Number.isNaN(pos)) {
                pos = 0;
            }
            let start = Math.min(Math.max(pos, 0), stringLength);
            // Avoid the `indexOf` call if no match is possible
            if (searchLength + start > stringLength) {
                return false;
            }
            let index = -1;
            while (++index < searchLength) {
                if (string.charCodeAt(start + index) !== searchString.charCodeAt(index)) {
                    return false;
                }
            }
            return true;
        };
        if (defineProperty) {
            defineProperty(String.prototype, 'startsWith', {
                'value': startsWith,
                'configurable': true,
                'writable': true
            });
        } else {
        // eslint-disable-next-line no-extend-native
            String.prototype.startsWith = startsWith;
        }
    }());
}

if (!String.prototype.endsWith) {
    // eslint-disable-next-line no-extend-native
    String.prototype.endsWith = function (searchString, position) {
        let subjectString = this.toString();
        if (typeof position !== 'number' || !isFinite(position) || Math.floor(position) !== position || position > subjectString.length) {
            position = subjectString.length;
        }
        position -= searchString.length;
        let lastIndex = subjectString.indexOf(searchString, position);
        return lastIndex !== -1 && lastIndex === position;
    };
}

if (!String.prototype.repeat) {
    // eslint-disable-next-line no-extend-native
    String.prototype.repeat = function (count) {
        'use strict';
        if (this === null) {
            throw new TypeError('can\'t convert ' + this + ' to object');
        }
        let str = '' + this;
        count = +count;
        if (Number.isNaN(count)) {
            count = 0;
        }
        if (count < 0) {
            throw new RangeError('repeat count must be non-negative');
        }
        if (count === Infinity) {
            throw new RangeError('repeat count must be less than infinity');
        }
        count = Math.floor(count);
        if (str.length === 0 || count === 0) {
            return '';
        }
        // Ensuring count is a 31-bit integer allows us to heavily optimize the
        // main part. But anyway, most current (August 2014) browsers can't handle
        // strings 1 << 28 chars or longer, so:
        if (str.length * count >= 1 << 28) {
            throw new RangeError('repeat count must not overflow maximum string size');
        }
        let rpt = '';
        for (;;) {
            if ((count & 1) === 1) {
                rpt += str;
            }
            count >>>= 1;
            if (count === 0) {
                break;
            }
            str += str;
        }
        // Could we try:
        // return Array(count + 1).join(this);
        return rpt;
    };
}

if (!String.prototype.includes) {
    // eslint-disable-next-line no-extend-native
    String.prototype.includes = function (search, start) {
        'use strict';
        if (typeof start !== 'number') {
            start = 0;
        }

        if (start + search.length > this.length) {
            return false;
        } else {
            return this.indexOf(search, start) !== -1;
        }
    };
}

if (!Array.prototype.includes) {
    // eslint-disable-next-line no-extend-native
    Array.prototype.includes = function (searchElement /* , fromIndex */) {
        'use strict';
        if (this == null) {
            throw new TypeError('Array.prototype.includes called on null or undefined');
        }

        let O = Object(this);
        let len = parseInt(O.length, 10) || 0;
        if (len === 0) {
            return false;
        }
        let n = parseInt(arguments[1], 10) || 0;
        let k;
        if (n >= 0) {
            k = n;
        } else {
            k = len + n;
            if (k < 0) { k = 0; }
        }
        let currentElement;
        while (k < len) {
            currentElement = O[k];
            if (searchElement === currentElement ||
            (Number.isNaN(searchElement) && Number.isNaN(currentElement))
            ) {
                return true;
            }
            k++;
        }
        return false;
    };
}

// Production steps of ECMA-262, Edition 6, 22.1.2.1
// Reference: https://people.mozilla.org/~jorendorff/es6-draft.html#sec-array.from
if (!Array.from) {
    Array.from = (function () {
        let toStr = Object.prototype.toString;
        let isCallable = function (fn: any) {
            return typeof fn === 'function' || toStr.call(fn) === '[object Function]';
        };
        let toInteger = function (value: any) {
            let number = Number(value);
            if (isNaN(number)) { return 0; }
            if (number === 0 || !isFinite(number)) { return number; }
            return (number > 0 ? 1 : -1) * Math.floor(Math.abs(number));
        };
        let maxSafeInteger = Math.pow(2, 53) - 1;
        let toLength = function (value: any) {
            let len = toInteger(value);
            return Math.min(Math.max(len, 0), maxSafeInteger);
        };

        // The length property of the from method is 1.
        return function from (this: any, arrayLike: any/* , mapFn, thisArg */) {
        // 1. Let C be the this value.
            let C = this;

            // 2. Let items be ToObject(arrayLike).
            let items = Object(arrayLike);

            // 3. ReturnIfAbrupt(items).
            if (arrayLike == null) {
                throw new TypeError('Array.from requires an array-like object - not null or undefined');
            }

            // 4. If mapfn is undefined, then let mapping be false.
            let mapFn = arguments.length > 1 ? arguments[1] : void undefined;
            let T;
            if (typeof mapFn !== 'undefined') {
                // 5. else
                // 5. a If IsCallable(mapfn) is false, throw a TypeError exception.
                if (!isCallable(mapFn)) {
                    throw new TypeError('Array.from: when provided, the second argument must be a function');
                }

                // 5. b. If thisArg was supplied, let T be thisArg; else let T be undefined.
                if (arguments.length > 2) {
                    T = arguments[2];
                }
            }

            // 10. Let lenValue be Get(items, "length").
            // 11. Let len be ToLength(lenValue).
            let len = toLength(items.length);

            // 13. If IsConstructor(C) is true, then
            // 13. a. Let A be the result of calling the [[Construct]] internal method of C with an argument list containing the single item len.
            // 14. a. Else, Let A be ArrayCreate(len).
            let A = isCallable(C) ? Object(new C(len)) : new Array(len);

            // 16. Let k be 0.
            let k = 0;
            // 17. Repeat, while k < len… (also steps a - h)
            let kValue;
            while (k < len) {
                kValue = items[k];
                if (mapFn) {
                    A[k] = typeof T === 'undefined' ? mapFn(kValue, k) : mapFn.call(T, kValue, k);
                } else {
                    A[k] = kValue;
                }
                k += 1;
            }
            // 18. Let putStatus be Put(A, "length", len, true).
            A.length = len;
            // 20. Return A.
            return A;
        };
    }());
}

if (typeof window !== 'undefined') {
    (function () {
        // http://paulirish.com/2011/requestanimationframe-for-smart-animating/
        // http://my.opera.com/emoller/blog/2011/12/20/requestanimationframe-for-smart-er-animating

        // requestAnimationFrame polyfill by Erik Möller. fixes from Paul Irish and Tino Zijdel

        // MIT license

        let lastTime = 0;
        let vendors = [ 'ms', 'moz', 'webkit', 'o' ];

        for (let x = 0; x < vendors.length && !window.requestAnimationFrame; ++x) {
            window.requestAnimationFrame = (
                (window as any)[ vendors[ x ] + 'RequestAnimationFrame' ]
            );

            window.cancelAnimationFrame = (
                (window as any)[ vendors[ x ] + 'CancelAnimationFrame' ] ||
                (window as any)[ vendors[ x ] + 'CancelRequestAnimationFrame' ]
            );
        }

        if (!window.requestAnimationFrame) {
            window.requestAnimationFrame = function (callback/* , element */) {
                let currTime = new Date().getTime();
                let timeToCall = Math.max(0, 16 - (currTime - lastTime));

                let id = window.setTimeout(function () {
                    let time = currTime + timeToCall;
                    callback(time);
                }, timeToCall);

                lastTime = currTime + timeToCall;

                return id;
            };
        }

        if (!window.cancelAnimationFrame) {
            window.cancelAnimationFrame = function (id) {
                clearTimeout(id);
            };
        }
    }());
}

if (Function.prototype.name === undefined && Object.defineProperty !== undefined) {
    // Missing in IE9-11.
    // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Function/name

    // eslint-disable-next-line no-extend-native
    Object.defineProperty(Function.prototype, 'name', {

        get: function () {
            return this.toString().match(/^\s*function\s*(\S*)\s*\(/)[ 1 ];
        }

    });
}

if (typeof window !== 'undefined') {
    if (window.performance === undefined) {
        /* global self */
        (window as any).performance = {};
    }

    if (window.performance.now === undefined) {
        (function () {
            let start = Date.now();

            window.performance.now = function () {
                return Date.now() - start;
            };
        })();
    }
}

if (Object.defineProperty !== undefined) {
    // Missing in IE < 13
    // MIT license
    // Copyright (c) 2016 Financial Times
    // https://github.com/Financial-Times/polyfill-service
    if (Number.MAX_SAFE_INTEGER === undefined) {
        Object.defineProperty(Number, 'MAX_SAFE_INTEGER', {
            enumerable: false,
            configurable: false,
            writable: false,
            value: Math.pow(2, 53) - 1
        });
    }
    if (Number.MIN_SAFE_INTEGER === undefined) {
        Object.defineProperty(Number, 'MIN_SAFE_INTEGER', {
            enumerable: false,
            configurable: false,
            writable: false,
            value: -(Math.pow(2, 53) - 1)
        });
    }
}

if (!Object.entries) {
    Object.entries = function(obj: any){
        let ownProps = Object.keys( obj ),
            i = ownProps.length,
            resArray = new Array(i); // preallocate the Array
        while (i--)
            resArray[i] = [ownProps[i], obj[ownProps[i]]];

        return resArray;
    };
}

// from https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/find
// https://tc39.github.io/ecma262/#sec-array.prototype.find
if (!Array.prototype.find) {
    Object.defineProperty(Array.prototype, 'find', {
        value: function(predicate: any) {
            // 1. Let O be ? ToObject(this value).
            if (this == null) {
                throw TypeError('"this" is null or not defined');
            }

            let o = Object(this);

            // 2. Let len be ? ToLength(? Get(O, "length")).
            let len = o.length >>> 0;

            // 3. If IsCallable(predicate) is false, throw a TypeError exception.
            if (typeof predicate !== 'function') {
                throw TypeError('predicate must be a function');
            }

            // 4. If thisArg was supplied, let T be thisArg; else let T be undefined.
            let thisArg = arguments[1];

            // 5. Let k be 0.
            let k = 0;

            // 6. Repeat, while k < len
            while (k < len) {
                // a. Let Pk be ! ToString(k).
                // b. Let kValue be ? Get(O, Pk).
                // c. Let testResult be ToBoolean(? Call(predicate, T, « kValue, k, O »)).
                // d. If testResult is true, return kValue.
                let kValue = o[k];
                if (predicate.call(thisArg, kValue, k, o)) {
                    return kValue;
                }
                // e. Increase k by 1.
                k++;
            }

            // 7. Return undefined.
            return undefined;
        },
        configurable: true,
        writable: true
    });
}

// from https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/fill
if (!Array.prototype.fill) {
    Object.defineProperty(Array.prototype, 'fill', {
        value: function(value: any) {

            // Steps 1-2.
            if (this == null) {
                throw new TypeError('this is null or not defined');
            }

            let O = Object(this);

            // Steps 3-5.
            let len = O.length >>> 0;

            // Steps 6-7.
            let start = arguments[1];
            let relativeStart = start >> 0;

            // Step 8.
            let k = relativeStart < 0 ?
                Math.max(len + relativeStart, 0) :
                Math.min(relativeStart, len);

            // Steps 9-10.
            let end = arguments[2];
            let relativeEnd = end === undefined ?
                len : end >> 0;

            // Step 11.
            let finalValue = relativeEnd < 0 ?
                Math.max(len + relativeEnd, 0) :
                Math.min(relativeEnd, len);

            // Step 12.
            while (k < finalValue) {
                O[k] = value;
                k++;
            }

            // Step 13.
            return O;
        }
    });
}

if (!Array.prototype.copyWithin) {
    Object.defineProperty(Array.prototype, 'copyWithin', {
        value: function(target: any, start: any/* , end*/) {
            // Steps 1-2.
            if (this == null) {
                throw new TypeError('this is null or not defined');
            }

            let O = Object(this);

            // Steps 3-5.
            let len = O.length >>> 0;

            // Steps 6-8.
            let relativeTarget = target >> 0;

            let to = relativeTarget < 0 ?
                Math.max(len + relativeTarget, 0) :
                Math.min(relativeTarget, len);

            // Steps 9-11.
            let relativeStart = start >> 0;

            let from = relativeStart < 0 ?
                Math.max(len + relativeStart, 0) :
                Math.min(relativeStart, len);

            // Steps 12-14.
            let end = arguments[2];
            let relativeEnd = end === undefined ? len : end >> 0;

            let final = relativeEnd < 0 ?
                Math.max(len + relativeEnd, 0) :
                Math.min(relativeEnd, len);

            // Step 15.
            let count = Math.min(final - from, len - to);

            // Steps 16-17.
            let direction = 1;

            if (from < to && to < (from + count)) {
                direction = -1;
                from += count - 1;
                to += count - 1;
            }

            // Step 18.
            while (count > 0) {
                if (from in O) {
                    O[to] = O[from];
                } else {
                    delete O[to];
                }

                from += direction;
                to += direction;
                count--;
            }

            // Step 19.
            return O;
        },
        configurable: true,
        writable: true
    });
}

if (!Int8Array.prototype.copyWithin) {
    Object.defineProperty(Int8Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}

if (!Int16Array.prototype.copyWithin) {
    Object.defineProperty(Int16Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}

if (!Int32Array.prototype.copyWithin) {
    Object.defineProperty(Int32Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}

if (!Uint8Array.prototype.copyWithin) {
    Object.defineProperty(Uint8Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}

if (!Uint16Array.prototype.copyWithin) {
    Object.defineProperty(Uint16Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}

if (!Uint32Array.prototype.copyWithin) {
    Object.defineProperty(Uint32Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}

if (!Float32Array.prototype.copyWithin) {
    Object.defineProperty(Float32Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}

if (!Float64Array.prototype.copyWithin) {
    Object.defineProperty(Float64Array.prototype, 'copyWithin', { value: Array.prototype.copyWithin });
}