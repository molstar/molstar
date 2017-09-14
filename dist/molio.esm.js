/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */
/**
 * Efficient integer and float parsers.
 *
 * For the purposes of parsing numbers from the mmCIF data representations,
 * up to 4 times faster than JS parseInt/parseFloat.
 */

function parseInt(str, start, end) {
    var ret = 0, neg = 1;
    if (str.charCodeAt(start) === 45 /* - */) {
        neg = -1;
        start++;
    }
    for (; start < end; start++) {
        var c = str.charCodeAt(start) - 48;
        if (c > 9 || c < 0)
            { return (neg * ret) | 0; }
        else
            { ret = (10 * ret + c) | 0; }
    }
    return neg * ret;
}
function parseScientific(main, str, start, end) {
    // handle + in '1e+1' separately.
    if (str.charCodeAt(start) === 43 /* + */)
        { start++; }
    return main * Math.pow(10.0, parseInt(str, start, end));
}

function parseFloat(str, start, end) {
    var neg = 1.0, ret = 0.0, point = 0.0, div = 1.0;
    if (str.charCodeAt(start) === 45) {
        neg = -1.0;
        ++start;
    }
    while (start < end) {
        var c = str.charCodeAt(start) - 48;
        if (c >= 0 && c < 10) {
            ret = ret * 10 + c;
            ++start;
        }
        else if (c === -2) {
            ++start;
            while (start < end) {
                c = str.charCodeAt(start) - 48;
                if (c >= 0 && c < 10) {
                    point = 10.0 * point + c;
                    div = 10.0 * div;
                    ++start;
                }
                else if (c === 53 || c === 21) {
                    return parseScientific(neg * (ret + point / div), str, start + 1, end);
                }
                else {
                    return neg * (ret + point / div);
                }
            }
            return neg * (ret + point / div);
        }
        else if (c === 53 || c === 21) {
            return parseScientific(neg * ret, str, start + 1, end);
        }
        else
            { break; }
    }
    return neg * ret;
}

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */
/**
 * Eat everything until a newline occurs.
 */
function eatLine(state) {
    while (state.position < state.length) {
        switch (state.data.charCodeAt(state.position)) {
            case 10:// \n
                state.currentTokenEnd = state.position;
                ++state.position;
                ++state.currentLineNumber;
                return;
            case 13:// \r
                state.currentTokenEnd = state.position;
                ++state.position;
                ++state.currentLineNumber;
                if (state.data.charCodeAt(state.position) === 10) {
                    ++state.position;
                }
                return;
            default:
                ++state.position;
        }
    }
    state.currentTokenEnd = state.position;
}
/**
 * Eat everything until a whitespace/newline occurs.
 */
function eatValue(state) {
    while (state.position < state.length) {
        switch (state.data.charCodeAt(state.position)) {
            case 9: // \t
            case 10: // \n
            case 13: // \r
            case 32:// ' '
                state.currentTokenEnd = state.position;
                return;
            default:
                ++state.position;
                break;
        }
    }
    state.currentTokenEnd = state.position;
}
/**
 * Skips all the whitespace - space, tab, newline, CR
 * Handles incrementing line count.
 */
function skipWhitespace(state) {
    var prev = 10;
    while (state.position < state.length) {
        var c = state.data.charCodeAt(state.position);
        switch (c) {
            case 9: // '\t'
            case 32:// ' '
                prev = c;
                ++state.position;
                break;
            case 10:// \n
                // handle \r\n
                if (prev !== 13) {
                    ++state.currentLineNumber;
                }
                prev = c;
                ++state.position;
                break;
            case 13:// \r
                prev = c;
                ++state.position;
                ++state.currentLineNumber;
                break;
            default:
                return prev;
        }
    }
    return prev;
}

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
var Tokens;
(function (Tokens) {
    function resize(tokens) {
        // scale the size using golden ratio, because why not.
        var newBuffer = new Int32Array((1.61 * tokens.indices.length) | 0);
        newBuffer.set(tokens.indices);
        tokens.indices = newBuffer;
        tokens.indicesLenMinus2 = (newBuffer.length - 2) | 0;
    }
    function add(tokens, start, end) {
        if (tokens.count > tokens.indicesLenMinus2) {
            resize(tokens);
        }
        tokens.indices[tokens.count++] = start;
        tokens.indices[tokens.count++] = end;
    }
    Tokens.add = add;
    function addUnchecked(tokens, start, end) {
        tokens.indices[tokens.count++] = start;
        tokens.indices[tokens.count++] = end;
    }
    Tokens.addUnchecked = addUnchecked;
    function create(size) {
        return {
            indicesLenMinus2: (size - 2) | 0,
            count: 0,
            indices: new Int32Array(size)
        };
    }
    Tokens.create = create;
})(Tokens || (Tokens = {}));

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */
/**
 * Represents a column that is not present.
 */
var _UndefinedColumn = function _UndefinedColumn() {
    this.isDefined = false;
};
_UndefinedColumn.prototype.getString = function getString (row) { return null; };

_UndefinedColumn.prototype.getInteger = function getInteger (row) { return 0; };
_UndefinedColumn.prototype.getFloat = function getFloat (row) { return 0.0; };
_UndefinedColumn.prototype.getValuePresence = function getValuePresence (row) { return 1 /* NotSpecified */; };
_UndefinedColumn.prototype.areValuesEqual = function areValuesEqual (rowA, rowB) { return true; };
_UndefinedColumn.prototype.stringEquals = function stringEquals (row, value) { return value === null; };
var UndefinedColumn = new _UndefinedColumn();

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */
var ShortStringPool;
(function (ShortStringPool) {
    function create() { return Object.create(null); }
    ShortStringPool.create = create;
    function get(pool, str) {
        if (str.length > 6)
            { return str; }
        var value = pool[str];
        if (value !== void 0)
            { return value; }
        pool[str] = str;
        return str;
    }
    ShortStringPool.get = get;
})(ShortStringPool || (ShortStringPool = {}));

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
/**
 * Represents a single column.
 */
var TextColumn = function TextColumn(table, data, name, index) {
    this.data = data;
    this.name = name;
    this.index = index;
    this.stringPool = ShortStringPool.create();
    this.isDefined = true;
    this.indices = table.indices;
    this.columnCount = table.columnCount;
};
/**
 * Returns the string value at given row.
 */
TextColumn.prototype.getString = function getString (row) {
    var i = (row * this.columnCount + this.index) * 2;
    return ShortStringPool.get(this.stringPool, this.data.substring(this.indices[i], this.indices[i + 1]));
};
/**
 * Returns the integer value at given row.
 */
TextColumn.prototype.getInteger = function getInteger (row) {
    var i = (row * this.columnCount + this.index) * 2;
    return parseInt(this.data, this.indices[i], this.indices[i + 1]);
};
/**
 * Returns the float value at given row.
 */
TextColumn.prototype.getFloat = function getFloat (row) {
    var i = (row * this.columnCount + this.index) * 2;
    return parseFloat(this.data, this.indices[i], this.indices[i + 1]);
};
/**
 * Returns true if the token has the specified string value.
 */
TextColumn.prototype.stringEquals = function stringEquals (row, value) {
        var this$1 = this;

    var aIndex = (row * this.columnCount + this.index) * 2, s = this.indices[aIndex], len = value.length;
    if (len !== this.indices[aIndex + 1] - s)
        { return false; }
    for (var i = 0; i < len; i++) {
        if (this$1.data.charCodeAt(i + s) !== value.charCodeAt(i))
            { return false; }
    }
    return true;
};
/**
 * Determines if values at the given rows are equal.
 */
TextColumn.prototype.areValuesEqual = function areValuesEqual (rowA, rowB) {
        var this$1 = this;

    var aIndex = (rowA * this.columnCount + this.index) * 2;
    var bIndex = (rowB * this.columnCount + this.index) * 2;
    var aS = this.indices[aIndex];
    var bS = this.indices[bIndex];
    var len = this.indices[aIndex + 1] - aS;
    if (len !== this.indices[bIndex + 1] - bS)
        { return false; }
    for (var i = 0; i < len; i++) {
        if (this$1.data.charCodeAt(i + aS) !== this$1.data.charCodeAt(i + bS)) {
            return false;
        }
    }
    return true;
};
TextColumn.prototype.getValuePresence = function getValuePresence (row) {
    var index = 2 * (row * this.columnCount + this.index);
    if (this.indices[index] === this.indices[index + 1]) {
        return 1 /* NotSpecified */;
    }
    return 0 /* Present */;
};
var CifColumn = (function (TextColumn) {
    function CifColumn () {
        TextColumn.apply(this, arguments);
    }

    if ( TextColumn ) CifColumn.__proto__ = TextColumn;
    CifColumn.prototype = Object.create( TextColumn && TextColumn.prototype );
    CifColumn.prototype.constructor = CifColumn;

    CifColumn.prototype.getString = function getString (row) {
        var ret = TextColumn.prototype.getString.call(this, row);
        if (ret === '.' || ret === '?')
            { return null; }
        return ret;
    };
    /**
     * Returns true if the value is not defined (. or ? token).
     */
    CifColumn.prototype.getValuePresence = function getValuePresence (row) {
        var index = 2 * (row * this.columnCount + this.index);
        var s = this.indices[index];
        if (this.indices[index + 1] - s !== 1)
            { return 0 /* Present */; }
        var v = this.data.charCodeAt(s);
        if (v === 46 /* . */)
            { return 1 /* NotSpecified */; }
        if (v === 63 /* ? */)
            { return 2 /* Unknown */; }
        return 0 /* Present */;
    };

    return CifColumn;
}(TextColumn));

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
/**
 * Represents a category backed by a string.
 */
var TextCategory = function TextCategory(data, name, columns, tokens) {
    this.name = name;
    this.indices = tokens.indices;
    this.data = data;
    this.columnCount = columns.length;
    this.rowCount = (tokens.count / 2 / columns.length) | 0;
    this.initColumns(columns);
};

var prototypeAccessors = { columnNames: {} };

prototypeAccessors.columnNames.get = function () {
    return this.columnNameList;
};
/**
 * Get a column object that makes accessing data easier.
 */
TextCategory.prototype.getColumn = function getColumn (name) {
    var i = this.columnIndices.get(name);
    if (i !== void 0)
        { return new TextColumn(this, this.data, name, i); }
    return UndefinedColumn;
};
TextCategory.prototype.initColumns = function initColumns (columns) {
        var this$1 = this;

    this.columnIndices = new Map();
    this.columnNameList = [];
    for (var i = 0; i < columns.length; i++) {
        this$1.columnIndices.set(columns[i], i);
        this$1.columnNameList.push(columns[i]);
    }
};

Object.defineProperties( TextCategory.prototype, prototypeAccessors );
var CifCategory = (function (TextCategory) {
    function CifCategory () {
        TextCategory.apply(this, arguments);
    }

    if ( TextCategory ) CifCategory.__proto__ = TextCategory;
    CifCategory.prototype = Object.create( TextCategory && TextCategory.prototype );
    CifCategory.prototype.constructor = CifCategory;

    CifCategory.prototype.getColumn = function getColumn (name) {
        var i = this.columnIndices.get(name);
        if (i !== void 0)
            { return new CifColumn(this, this.data, name, i); }
        return UndefinedColumn;
    };
    CifCategory.prototype.initColumns = function initColumns (columns) {
        var this$1 = this;

        this.columnIndices = new Map();
        this.columnNameList = [];
        for (var i = 0; i < columns.length; i++) {
            var colName = columns[i].substr(this$1.name.length + 1);
            this$1.columnIndices.set(colName, i);
            this$1.columnNameList.push(colName);
        }
    };

    return CifCategory;
}(TextCategory));

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */
var ParserResult;
(function (ParserResult) {
    function error(message, line) {
        if ( line === void 0 ) line = -1;

        return new ParserError(message, line);
    }
    ParserResult.error = error;
    function success(result, warnings) {
        if ( warnings === void 0 ) warnings = [];

        return new ParserSuccess(result, warnings);
    }
    ParserResult.success = success;
})(ParserResult || (ParserResult = {}));
var ParserError = function ParserError(message, line) {
    this.message = message;
    this.line = line;
    this.isError = true;
};
ParserError.prototype.toString = function toString () {
    if (this.line >= 0) {
        return ("[Line " + (this.line) + "] " + (this.message));
    }
    return this.message;
};
var ParserSuccess = function ParserSuccess(result, warnings) {
    this.result = result;
    this.warnings = warnings;
    this.isError = false;
};

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
var GroFile = function GroFile(data) {
    this.blocks = [];
    this.data = data;
};
var GroBlock = function GroBlock(data) {
    this.data = data;
    this.categoryMap = new Map();
    this.categoryList = [];
};

GroBlock.prototype.getCategory = function getCategory (name) {
    return this.categoryMap.get(name);
};
/**
 * Adds a category.
 */
GroBlock.prototype.addCategory = function addCategory (category) {
    this.categoryList[this.categoryList.length] = category;
    this.categoryMap.set(category.name, category);
};
function createTokenizer(data) {
    return {
        data: data,
        position: 0,
        length: data.length,
        currentLineNumber: 1,
        currentTokenStart: 0,
        currentTokenEnd: 0,
        numberOfAtoms: 0,
        hasVelocities: false,
        numberOfDecimalPlaces: 3
    };
}
/**
 * title string (free format string, optional time in ps after 't=')
 */
function handleTitleString(state, tokens) {
    eatLine(state);
    // console.log('title', state.data.substring(state.currentTokenStart, state.currentTokenEnd))
    var start = state.currentTokenStart;
    var end = state.currentTokenEnd;
    var valueStart = state.currentTokenStart;
    var valueEnd = start;
    while (valueEnd < end && !isTime(state.data, valueEnd))
        { ++valueEnd; }
    if (isTime(state.data, valueEnd)) {
        var timeStart = valueEnd + 2;
        while (valueEnd > start && isSpaceOrComma(state.data, valueEnd - 1))
            { --valueEnd; }
        Tokens.add(tokens, valueStart, valueEnd); // title
        while (timeStart < end && state.data.charCodeAt(timeStart) === 32)
            { ++timeStart; }
        while (valueEnd > timeStart && state.data.charCodeAt(valueEnd - 1) === 32)
            { --valueEnd; }
        Tokens.add(tokens, timeStart, end); // time
    }
    else {
        Tokens.add(tokens, valueStart, valueEnd); // title
        Tokens.add(tokens, valueEnd, valueEnd); // empty token for time
    }
}
function isSpaceOrComma(data, position) {
    var c = data.charCodeAt(position);
    return c === 32 || c === 44;
}
function isTime(data, position) {
    // T/t
    var c = data.charCodeAt(position);
    if (c !== 84 && c !== 116)
        { return false; }
    // =
    if (data.charCodeAt(position + 1) !== 61)
        { return false; }
    return true;
}
// function isDot(state: TokenizerState): boolean {
//     // .
//     if (state.data.charCodeAt(state.currentTokenStart) !== 46) return false;
//     return true;
// }
// function numberOfDecimalPlaces (state: TokenizerState) {
//     // var ndec = firstLines[ 2 ].length - firstLines[ 2 ].lastIndexOf('.') - 1
//     const start = state.currentTokenStart
//     const end = state.currentTokenEnd
//     for (let i = end; start < i; --i) {
//         // .
//         if (state.data.charCodeAt(i) === 46) return end - start - i
//     }
//     throw new Error('Could not determine number of decimal places')
// }
/**
 * number of atoms (free format integer)
 */
function handleNumberOfAtoms(state, tokens) {
    skipWhitespace(state);
    state.currentTokenStart = state.position;
    eatValue(state);
    state.numberOfAtoms = parseInt(state.data, state.currentTokenStart, state.currentTokenEnd);
    Tokens.add(tokens, state.currentTokenStart, state.currentTokenEnd);
    eatLine(state);
}
// function checkForVelocities (state: GroState) {
// }
/**
 * This format is fixed, ie. all columns are in a fixed position.
 * Optionally (for now only yet with trjconv) you can write gro files
 * with any number of decimal places, the format will then be n+5
 * positions with n decimal places (n+1 for velocities) in stead
 * of 8 with 3 (with 4 for velocities). Upon reading, the precision
 * will be inferred from the distance between the decimal points
 * (which will be n+5). Columns contain the following information
 * (from left to right):
 *     residue number (5 positions, integer)
 *     residue name (5 characters)
 *     atom name (5 characters)
 *     atom number (5 positions, integer)
 *     position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
 *     velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
 */
function handleAtoms(state, block) {
    console.log('MOINMOIN');
    var name = 'atoms';
    var columns = ['residueNumber', 'residueName', 'atomName', 'atomNumber', 'x', 'y', 'z'];
    if (state.hasVelocities) {
        columns.push('vx', 'vy', 'vz');
    }
    var fieldSizes = [5, 5, 5, 5, 8, 8, 8, 8, 8, 8];
    var columnCount = columns.length;
    var tokens = Tokens.create(state.numberOfAtoms * 2 * columnCount);
    var start;
    var end;
    var valueStart;
    var valueEnd = state.position;
    for (var i = 0; i < state.numberOfAtoms; ++i) {
        state.currentTokenStart = state.position;
        end = state.currentTokenStart;
        for (var j = 0; j < columnCount; ++j) {
            start = end;
            end = start + fieldSizes[j];
            // trim
            valueStart = start;
            valueEnd = end;
            while (valueStart < valueEnd && state.data.charCodeAt(valueStart) === 32)
                { ++valueStart; }
            while (valueEnd > valueStart && state.data.charCodeAt(valueEnd - 1) === 32)
                { --valueEnd; }
            Tokens.addUnchecked(tokens, valueStart, valueEnd);
        }
        state.position = valueEnd;
        eatLine(state);
    }
    block.addCategory(new TextCategory(state.data, name, columns, tokens));
}
/**
 * box vectors (free format, space separated reals), values:
 * v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
 * the last 6 values may be omitted (they will be set to zero).
 * Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
 */
function handleBoxVectors(state, tokens) {
    // just read the first three values, ignore any remaining
    for (var i = 0; i < 3; ++i) {
        skipWhitespace(state);
        state.currentTokenStart = state.position;
        eatValue(state);
        Tokens.add(tokens, state.currentTokenStart, state.currentTokenEnd);
    }
}
/**
 * Creates an error result.
 */
// function error(line: number, message: string) {
//     return ParserResult.error<GroFile>(message, line);
// }
/**
 * Creates a data result.
 */
function result(data) {
    return ParserResult.success(data);
}
function parseInternal(data) {
    var state = createTokenizer(data);
    var file = new GroFile(data);
    var block = new GroBlock(data);
    file.blocks.push(block);
    var headerColumns = ['title', 'timeInPs', 'numberOfAtoms', 'boxX', 'boxY', 'boxZ'];
    var headerTokens = Tokens.create(2 * headerColumns.length);
    var header = new TextCategory(state.data, 'header', headerColumns, headerTokens);
    block.addCategory(header);
    handleTitleString(state, headerTokens);
    handleNumberOfAtoms(state, headerTokens);
    handleAtoms(state, block);
    handleBoxVectors(state, headerTokens);
    return result(file);
}
function parse(data) {
    return parseInternal(data);
}

/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export { parse as groReader };
//# sourceMappingURL=molio.esm.js.map
