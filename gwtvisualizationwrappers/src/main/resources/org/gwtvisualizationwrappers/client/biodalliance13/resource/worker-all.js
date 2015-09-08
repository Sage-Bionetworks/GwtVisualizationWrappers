(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);throw new Error("Cannot find module '"+o+"'")}var f=n[o]={exports:{}};t[o][0].call(f.exports,function(e){var n=t[o][1][e];return s(n?n:e)},f,f.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bam.js: indexed binary alignments
//

"use strict";

if (typeof(require) !== 'undefined') {
    var spans = require('./spans');
    var Range = spans.Range;
    var union = spans.union;
    var intersection = spans.intersection;

    var bin = require('./bin');
    var readInt = bin.readInt;
    var readShort = bin.readShort;
    var readByte = bin.readByte;
    var readInt64 = bin.readInt64;
    var readFloat = bin.readFloat;

    var lh3utils = require('./lh3utils');
    var readVob = lh3utils.readVob;
    var unbgzf = lh3utils.unbgzf;
    var reg2bins = lh3utils.reg2bins;
    var Chunk = lh3utils.Chunk;
}


var BAM_MAGIC = 0x14d4142;
var BAI_MAGIC = 0x1494142;

var BamFlags = {
    MULTIPLE_SEGMENTS:       0x1,
    ALL_SEGMENTS_ALIGN:      0x2,
    SEGMENT_UNMAPPED:        0x4,
    NEXT_SEGMENT_UNMAPPED:   0x8,
    REVERSE_COMPLEMENT:      0x10,
    NEXT_REVERSE_COMPLEMENT: 0x20,
    FIRST_SEGMENT:           0x40,
    LAST_SEGMENT:            0x80,
    SECONDARY_ALIGNMENT:     0x100,
    QC_FAIL:                 0x200,
    DUPLICATE:               0x400,
    SUPPLEMENTARY:           0x800
};

function BamFile() {
}


// Calculate the length (in bytes) of the BAI ref starting at offset.
// Returns {nbin, length, minBlockIndex}
function _getBaiRefLength(uncba, offset) {
    var p = offset;
    var nbin = readInt(uncba, p); p += 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(uncba, p);
        var nchnk = readInt(uncba, p+4);
        p += 8 + (nchnk * 16);
    }
    var nintv = readInt(uncba, p); p += 4;

    var minBlockIndex = 1000000000;
    var q = p;
    for (var i = 0; i < nintv; ++i) {
        var v = readVob(uncba, q); q += 8;
        if (v) {
            var bi = v.block;
            if (v.offset > 0)
                bi += 65536;

            if (bi < minBlockIndex)
                minBlockIndex = bi;
            break;
        }
    }
    p += (nintv * 8);

    return {
        minBlockIndex: minBlockIndex,
        nbin: nbin,
        length: p - offset
    };
}


function makeBam(data, bai, indexChunks, callback, attempted) {
    // Do an initial probe on the BAM file to catch any mixed-content errors.
    data.slice(0, 10).fetch(function(header) {
        if (header) {
            return makeBam2(data, bai, indexChunks, callback, attempted);
        } else {
            return callback(null, "Couldn't access BAM.");
        }
    }, {timeout: 5000});
}

function makeBam2(data, bai, indexChunks, callback, attempted) {
    var bam = new BamFile();
    bam.data = data;
    bam.bai = bai;
    bam.indexChunks = indexChunks;

    var minBlockIndex = bam.indexChunks ? bam.indexChunks.minBlockIndex : 1000000000;

    // Fills out bam.chrToIndex and bam.indexToChr based on the first few bytes of the BAM.
    function parseBamHeader(r) {
        if (!r) {
            return callback(null, "Couldn't access BAM");
        }

        var unc = unbgzf(r, r.byteLength);
        var uncba = new Uint8Array(unc);

        var magic = readInt(uncba, 0);
        if (magic != BAM_MAGIC) {
            return callback(null, "Not a BAM file, magic=0x" + magic.toString(16));
        }
        var headLen = readInt(uncba, 4);
        var header = '';
        for (var i = 0; i < headLen; ++i) {
            header += String.fromCharCode(uncba[i + 8]);
        }

        var nRef = readInt(uncba, headLen + 8);
        var p = headLen + 12;

        bam.chrToIndex = {};
        bam.indexToChr = [];
        for (var i = 0; i < nRef; ++i) {
            var lName = readInt(uncba, p);
            var name = '';
            for (var j = 0; j < lName-1; ++j) {
                name += String.fromCharCode(uncba[p + 4 + j]);
            }
            var lRef = readInt(uncba, p + lName + 4);
            bam.chrToIndex[name] = i;
            if (name.indexOf('chr') == 0) {
                bam.chrToIndex[name.substring(3)] = i;
            } else {
                bam.chrToIndex['chr' + name] = i;
            }
            bam.indexToChr.push(name);

            p = p + 8 + lName;
        }

        if (bam.indices) {
            return callback(bam);
        }
    }

    function parseBai(header) {
        if (!header) {
            return "Couldn't access BAI";
        }

        var uncba = new Uint8Array(header);
        var baiMagic = readInt(uncba, 0);
        if (baiMagic != BAI_MAGIC) {
            return callback(null, 'Not a BAI file, magic=0x' + baiMagic.toString(16));
        }

        var nref = readInt(uncba, 4);

        bam.indices = [];

        var p = 8;
        for (var ref = 0; ref < nref; ++ref) {
            var blockStart = p;
            var o = _getBaiRefLength(uncba, blockStart);
            p += o.length;

            minBlockIndex = Math.min(o.minBlockIndex, minBlockIndex);

            var nbin = o.nbin;

            if (nbin > 0) {
                bam.indices[ref] = new Uint8Array(header, blockStart, p - blockStart);
            }
        }

        return true;
    }

    if (!bam.indexChunks) {
        bam.bai.fetch(function(header) {   // Do we really need to fetch the whole thing? :-(
            var result = parseBai(header);
            if (result !== true) {
                if (bam.bai.url && typeof(attempted) === "undefined") {
                    // Already attempted x.bam.bai not there so now trying x.bai
                    bam.bai.url = bam.data.url.replace(new RegExp('.bam$'), '.bai');
                    
                     // True lets us know we are making a second attempt
                    makeBam2(data, bam.bai, indexChunks, callback, true);
                }
                else {
                    // We've attempted x.bam.bai & x.bai and nothing worked
                    callback(null, result);
                }
            } else {
              bam.data.slice(0, minBlockIndex).fetch(parseBamHeader);
            }
        });   // Timeout on first request to catch Chrome mixed-content error.
    } else {
        var chunks = bam.indexChunks.chunks;
        bam.indices = []
        for (var i = 0; i < chunks.length; i++) {
           bam.indices[i] = null;  // To be filled out lazily as needed
        }
        bam.data.slice(0, minBlockIndex).fetch(parseBamHeader);
    }
}



BamFile.prototype.blocksForRange = function(refId, min, max) {
    var index = this.indices[refId];
    if (!index) {
        return [];
    }

    var intBinsL = reg2bins(min, max);
    var intBins = [];
    for (var i = 0; i < intBinsL.length; ++i) {
        intBins[intBinsL[i]] = true;
    }
    var leafChunks = [], otherChunks = [];

    var nbin = readInt(index, 0);
    var p = 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(index, p);
        var nchnk = readInt(index, p+4);
//        dlog('bin=' + bin + '; nchnk=' + nchnk);
        p += 8;
        if (intBins[bin]) {
            for (var c = 0; c < nchnk; ++c) {
                var cs = readVob(index, p);
                var ce = readVob(index, p + 8);
                (bin < 4681 ? otherChunks : leafChunks).push(new Chunk(cs, ce));
                p += 16;
            }
        } else {
            p +=  (nchnk * 16);
        }
    }
    // console.log('leafChunks = ' + miniJSONify(leafChunks));
    // console.log('otherChunks = ' + miniJSONify(otherChunks));

    var nintv = readInt(index, p);
    var lowest = null;
    var minLin = Math.min(min>>14, nintv - 1), maxLin = Math.min(max>>14, nintv - 1);
    for (var i = minLin; i <= maxLin; ++i) {
        var lb =  readVob(index, p + 4 + (i * 8));
        if (!lb) {
            continue;
        }
        if (!lowest || lb.block < lowest.block || lb.offset < lowest.offset) {
            lowest = lb;
        }
    }
    // console.log('Lowest LB = ' + lowest);
    
    var prunedOtherChunks = [];
    if (lowest != null) {
        for (var i = 0; i < otherChunks.length; ++i) {
            var chnk = otherChunks[i];
            if (chnk.maxv.block >= lowest.block && chnk.maxv.offset >= lowest.offset) {
                prunedOtherChunks.push(chnk);
            }
        }
    }
    // console.log('prunedOtherChunks = ' + miniJSONify(prunedOtherChunks));
    otherChunks = prunedOtherChunks;

    var intChunks = [];
    for (var i = 0; i < otherChunks.length; ++i) {
        intChunks.push(otherChunks[i]);
    }
    for (var i = 0; i < leafChunks.length; ++i) {
        intChunks.push(leafChunks[i]);
    }

    intChunks.sort(function(c0, c1) {
        var dif = c0.minv.block - c1.minv.block;
        if (dif != 0) {
            return dif;
        } else {
            return c0.minv.offset - c1.minv.offset;
        }
    });
    var mergedChunks = [];
    if (intChunks.length > 0) {
        var cur = intChunks[0];
        for (var i = 1; i < intChunks.length; ++i) {
            var nc = intChunks[i];
            if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) { // no point splitting mid-block
                cur = new Chunk(cur.minv, nc.maxv);
            } else {
                mergedChunks.push(cur);
                cur = nc;
            }
        }
        mergedChunks.push(cur);
    }
    // console.log('mergedChunks = ' + miniJSONify(mergedChunks));

    return mergedChunks;
}

BamFile.prototype.fetch = function(chr, min, max, callback, opts) {
    var thisB = this;
    opts = opts || {};

    var chrId = this.chrToIndex[chr];
    var chunks;
    if (chrId === undefined) {
        chunks = [];
    } else {
        // Fetch this portion of the BAI if it hasn't been loaded yet.
        if (this.indices[chrId] === null && this.indexChunks.chunks[chrId]) {
            var start_stop = this.indexChunks.chunks[chrId];
            return this.bai.slice(start_stop[0], start_stop[1]).fetch(function(data) {
                var buffer = new Uint8Array(data);
                this.indices[chrId] = buffer;
                return this.fetch(chr, min, max, callback, opts);
            }.bind(this));
        }

        chunks = this.blocksForRange(chrId, min, max);
        if (!chunks) {
            callback(null, 'Error in index fetch');
        }
    }
    
    var records = [];
    var index = 0;
    var data;

    function tramp() {
        if (index >= chunks.length) {
            return callback(records);
        } else if (!data) {
            var c = chunks[index];
            var fetchMin = c.minv.block;
            var fetchMax = c.maxv.block + (1<<16); // *sigh*
            // console.log('fetching ' + fetchMin + ':' + fetchMax);
            thisB.data.slice(fetchMin, fetchMax - fetchMin).fetch(function(r) {
                data = unbgzf(r, c.maxv.block - c.minv.block + 1);
                return tramp();
            });
        } else {
            var ba = new Uint8Array(data);
            var finished = thisB.readBamRecords(ba, chunks[index].minv.offset, records, min, max, chrId, opts);
            data = null;
            ++index;
            if (finished)
                return callback(records);
            else
                return tramp();
        }
    }
    tramp();
}

var SEQRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
var CIGAR_DECODER = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];

function BamRecord() {
}

BamFile.prototype.readBamRecords = function(ba, offset, sink, min, max, chrId, opts) {
    while (true) {
        var blockSize = readInt(ba, offset);
        var blockEnd = offset + blockSize + 4;
        if (blockEnd >= ba.length) {
            return false;
        }

        var record = new BamRecord();

        var refID = readInt(ba, offset + 4);
        var pos = readInt(ba, offset + 8);
        
        var bmn = readInt(ba, offset + 12);
        var bin = (bmn & 0xffff0000) >> 16;
        var mq = (bmn & 0xff00) >> 8;
        var nl = bmn & 0xff;

        var flag_nc = readInt(ba, offset + 16);
        var flag = (flag_nc & 0xffff0000) >> 16;
        var nc = flag_nc & 0xffff;
    
        var lseq = readInt(ba, offset + 20);
        
        var nextRef  = readInt(ba, offset + 24);
        var nextPos = readInt(ba, offset + 28);
        
        var tlen = readInt(ba, offset + 32);
    
        record.segment = this.indexToChr[refID];
        record.flag = flag;
        record.pos = pos;
        record.mq = mq;
        if (opts.light)
            record.seqLength = lseq;

        if (!opts.light) {
            if (nextRef >= 0) {
                record.nextSegment = this.indexToChr[nextRef];
                record.nextPos = nextPos;
            }

            var readName = '';
            for (var j = 0; j < nl-1; ++j) {
                readName += String.fromCharCode(ba[offset + 36 + j]);
            }
            record.readName = readName;
        
            var p = offset + 36 + nl;

            var cigar = '';
            for (var c = 0; c < nc; ++c) {
                var cigop = readInt(ba, p);
                cigar = cigar + (cigop>>4) + CIGAR_DECODER[cigop & 0xf];
                p += 4;
            }
            record.cigar = cigar;
        
            var seq = '';
            var seqBytes = (lseq + 1) >> 1;
            for (var j = 0; j < seqBytes; ++j) {
                var sb = ba[p + j];
                seq += SEQRET_DECODER[(sb & 0xf0) >> 4];
                if (seq.length < lseq)
                    seq += SEQRET_DECODER[(sb & 0x0f)];
            }
            p += seqBytes;
            record.seq = seq;

            var qseq = '';
            for (var j = 0; j < lseq; ++j) {
                qseq += String.fromCharCode(ba[p + j] + 33);
            }
            p += lseq;
            record.quals = qseq;

            while (p < blockEnd) {
                var tag = String.fromCharCode(ba[p], ba[p + 1]);
                var type = String.fromCharCode(ba[p + 2]);
                var value;

                if (type == 'A') {
                    value = String.fromCharCode(ba[p + 3]);
                    p += 4;
                } else if (type == 'i' || type == 'I') {
                    value = readInt(ba, p + 3);
                    p += 7;
                } else if (type == 'c' || type == 'C') {
                    value = ba[p + 3];
                    p += 4;
                } else if (type == 's' || type == 'S') {
                    value = readShort(ba, p + 3);
                    p += 5;
                } else if (type == 'f') {
                    value = readFloat(ba, p + 3);
                    p += 7;
                } else if (type == 'Z' || type == 'H') {
                    p += 3;
                    value = '';
                    for (;;) {
                        var cc = ba[p++];
                        if (cc == 0) {
                            break;
                        } else {
                            value += String.fromCharCode(cc);
                        }
                    }
                } else if (type == 'B') {
                    var atype = String.fromCharCode(ba[p + 3]);
                    var alen = readInt(ba, p + 4);
                    var elen;
                    var reader;
                    if (atype == 'i' || atype == 'I' || atype == 'f') {
                        elen = 4;
                        if (atype == 'f')
                            reader = readFloat;
                        else
                            reader = readInt;
                    } else if (atype == 's' || atype == 'S') {
                        elen = 2;
                        reader = readShort;
                    } else if (atype == 'c' || atype == 'C') {
                        elen = 1;
                        reader = readByte;
                    } else {
                        throw 'Unknown array type ' + atype;
                    }

                    p += 8;
                    value = [];
                    for (var i = 0; i < alen; ++i) {
                        value.push(reader(ba, p));
                        p += elen;
                    }
                } else {
                    throw 'Unknown type '+ type;
                }
                record[tag] = value;
            }
        }

        if (!min || record.pos <= max && record.pos + lseq >= min) {
            if (chrId === undefined || refID == chrId) {
                sink.push(record);
            }
        }
        if (record.pos > max) {
            return true;
        }
        offset = blockEnd;
    }

    // Exits via top of loop.
};

if (typeof(module) !== 'undefined') {
    module.exports = {
        makeBam: makeBam,
        BAM_MAGIC: BAM_MAGIC,
        BAI_MAGIC: BAI_MAGIC,
        BamFlags: BamFlags
    };
}

},{"./bin":3,"./lh3utils":8,"./spans":10}],2:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// bigwig.js: indexed binary WIG (and BED) files
//

"use strict";


if (typeof(require) !== 'undefined') {
    var spans = require('./spans');
    var Range = spans.Range;
    var union = spans.union;
    var intersection = spans.intersection;

    var das = require('./das');
    var DASFeature = das.DASFeature;
    var DASGroup = das.DASGroup;

    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;

    var bin = require('./bin');
    var readInt = bin.readInt;

    var jszlib = require('jszlib');
    var jszlib_inflate_buffer = jszlib.inflateBuffer;
    var arrayCopy = jszlib.arrayCopy;
}

var BIG_WIG_MAGIC = 0x888FFC26;
var BIG_WIG_MAGIC_BE = 0x26FC8F88;
var BIG_BED_MAGIC = 0x8789F2EB;
var BIG_BED_MAGIC_BE = 0xEBF28987;


var BIG_WIG_TYPE_GRAPH = 1;
var BIG_WIG_TYPE_VSTEP = 2;
var BIG_WIG_TYPE_FSTEP = 3;
  
var M1 = 256;
var M2 = 256*256;
var M3 = 256*256*256;
var M4 = 256*256*256*256;

var BED_COLOR_REGEXP = new RegExp("^[0-9]+,[0-9]+,[0-9]+");

function bwg_readOffset(ba, o) {
    var offset = ba[o] + ba[o+1]*M1 + ba[o+2]*M2 + ba[o+3]*M3 + ba[o+4]*M4;
    return offset;
}

function BigWig() {
}

BigWig.prototype.readChromTree = function(callback) {
    var thisB = this;
    this.chromsToIDs = {};
    this.idsToChroms = {};
    this.maxID = 0;

    var udo = this.unzoomedDataOffset;
    var eb = (udo - this.chromTreeOffset) & 3;
    udo = udo + 4 - eb;

    this.data.slice(this.chromTreeOffset, udo - this.chromTreeOffset).fetch(function(bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        var bptReadNode = function(offset) {
            var nodeType = ba[offset];
            var cnt = sa[(offset/2) + 1];
            offset += 4;
            for (var n = 0; n < cnt; ++n) {
                if (nodeType == 0) {
                    offset += keySize;
                    var childOffset = bwg_readOffset(ba, offset);
                    offset += 8;
                    childOffset -= thisB.chromTreeOffset;
                    bptReadNode(childOffset);
                } else {
                    var key = '';
                    for (var ki = 0; ki < keySize; ++ki) {
                        var charCode = ba[offset++];
                        if (charCode != 0) {
                            key += String.fromCharCode(charCode);
                        }
                    }
                    var chromId = (ba[offset+3]<<24) | (ba[offset+2]<<16) | (ba[offset+1]<<8) | (ba[offset+0]);
                    var chromSize = (ba[offset + 7]<<24) | (ba[offset+6]<<16) | (ba[offset+5]<<8) | (ba[offset+4]);
                    offset += 8;

                    thisB.chromsToIDs[key] = chromId;
                    if (key.indexOf('chr') == 0) {
                        thisB.chromsToIDs[key.substr(3)] = chromId;
                    }
                    thisB.idsToChroms[chromId] = key;
                    thisB.maxID = Math.max(thisB.maxID, chromId);
                }
            }
        };
        bptReadNode(rootNodeOffset);

        callback(thisB);
    });
}

function BigWigView(bwg, cirTreeOffset, cirTreeLength, isSummary) {
    this.bwg = bwg;
    this.cirTreeOffset = cirTreeOffset;
    this.cirTreeLength = cirTreeLength;
    this.isSummary = isSummary;
}



BigWigView.prototype.readWigData = function(chrName, min, max, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.
        return callback([]);
    } else {
        this.readWigDataById(chr, min, max, callback);
    }
}

BigWigView.prototype.readWigDataById = function(chr, min, max, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function(result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.readWigDataById(chr, min, max, callback);
        });
        return;
    }

    var blocksToFetch = [];
    var outstanding = 0;

    var beforeBWG = Date.now();

    var filter = function(chromId, fmin, fmax, toks) {
        return ((chr < 0 || chromId == chr) && fmin <= max && fmax >= min);
    }

    var cirFobRecur = function(offset, level) {
        if (thisB.bwg.instrument)
            console.log('level=' + level + '; offset=' + offset + '; time=' + (Date.now()|0));

        outstanding += offset.length;

        if (offset.length == 1 && offset[0] - thisB.cirTreeOffset == 48 && thisB.cachedCirRoot) {
            cirFobRecur2(thisB.cachedCirRoot, 0, level);
            --outstanding;
            if (outstanding == 0) {
                thisB.fetchFeatures(filter, blocksToFetch, callback);
            }
            return;
        }

        var maxCirBlockSpan = 4 +  (thisB.cirBlockSize * 32);   // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }
        
        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    }

    var cirFobStartFetch = function(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function(resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);

                    if (offset[i] - thisB.cirTreeOffset == 48 && offset[i] - fr.min() == 0)
                        thisB.cachedCirRoot = resultBuffer;

                    --outstanding;
                    if (outstanding == 0) {
                        thisB.fetchFeatures(filter, blocksToFetch, callback);
                    }
                }
            }
        });
    }

    var cirFobRecur2 = function(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset/2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                var blockSize = bwg_readOffset(ba, offset+24);
                if (((chr < 0 || startChrom < chr) || (startChrom == chr && startBase <= max)) &&
                    ((chr < 0 || endChrom   > chr) || (endChrom == chr && endBase >= min)))
                {
                    blocksToFetch.push({offset: blockOffset, size: blockSize});
                }
                offset += 32;
            }
        } else {
            var recurOffsets = [];
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                if ((chr < 0 || startChrom < chr || (startChrom == chr && startBase <= max)) &&
                    (chr < 0 || endChrom   > chr || (endChrom == chr && endBase >= min)))
                {
                    recurOffsets.push(blockOffset);
                }
                offset += 24;
            }
            if (recurOffsets.length > 0) {
                cirFobRecur(recurOffsets, level + 1);
            }
        }
    };

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
}


BigWigView.prototype.fetchFeatures = function(filter, blocksToFetch, callback) {
    var thisB = this;

    blocksToFetch.sort(function(b0, b1) {
        return (b0.offset|0) - (b1.offset|0);
    });

    if (blocksToFetch.length == 0) {
        callback([]);
    } else {
        var features = [];
        var createFeature = function(chr, fmin, fmax, opts) {
            if (!opts) {
                opts = {};
            }
        
            var f = new DASFeature();
            f._chromId = chr;
            f.segment = thisB.bwg.idsToChroms[chr];
            f.min = fmin;
            f.max = fmax;
            f.type = thisB.bwg.type;
            
            for (var k in opts) {
                f[k] = opts[k];
            }
            
            features.push(f);
        };

        var tramp = function() {
            if (blocksToFetch.length == 0) {
                var afterBWG = Date.now();
                // dlog('BWG fetch took ' + (afterBWG - beforeBWG) + 'ms');
                callback(features);
                return;  // just in case...
            } else {
                var block = blocksToFetch[0];
                if (block.data) {
                    thisB.parseFeatures(block.data, createFeature, filter);
                    blocksToFetch.splice(0, 1);
                    tramp();
                } else {
                    var fetchStart = block.offset;
                    var fetchSize = block.size;
                    var bi = 1;
                    while (bi < blocksToFetch.length && blocksToFetch[bi].offset == (fetchStart + fetchSize)) {
                        fetchSize += blocksToFetch[bi].size;
                        ++bi;
                    }

                    thisB.bwg.data.slice(fetchStart, fetchSize).fetch(function(result) {
                        var offset = 0;
                        var bi = 0;
                        while (offset < fetchSize) {
                            var fb = blocksToFetch[bi];
                        
                            var data;
                            if (thisB.bwg.uncompressBufSize > 0) {
                                data = jszlib_inflate_buffer(result, offset + 2, fb.size - 2);
                            } else {
                                var tmp = new Uint8Array(fb.size);    // FIXME is this really the best we can do?
                                arrayCopy(new Uint8Array(result, offset, fb.size), 0, tmp, 0, fb.size);
                                data = tmp.buffer;
                            }
                            fb.data = data;
                            
                            offset += fb.size;
                            ++bi;
                        }
                        tramp();
                    });
                }
            }
        }
        tramp();
    }
}

BigWigView.prototype.parseFeatures = function(data, createFeature, filter) {
    var ba = new Uint8Array(data);

    if (this.isSummary) {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var itemCount = data.byteLength/32;
        for (var i = 0; i < itemCount; ++i) {
            var chromId =   la[(i*8)];
            var start =     la[(i*8)+1];
            var end =       la[(i*8)+2];
            var validCnt =  la[(i*8)+3];
            var minVal    = fa[(i*8)+4];
            var maxVal    = fa[(i*8)+5];
            var sumData   = fa[(i*8)+6];
            var sumSqData = fa[(i*8)+7];
            
            if (filter(chromId, start + 1, end)) {
                var summaryOpts = {type: 'bigwig', score: sumData/validCnt, maxScore: maxVal};
                if (this.bwg.type == 'bigbed') {
                    summaryOpts.type = 'density';
                }
                createFeature(chromId, start + 1, end, summaryOpts);
            }
        }
    } else if (this.bwg.type == 'bigwig') {
        var sa = new Int16Array(data);
        var la = new Int32Array(data);
        var fa = new Float32Array(data);

        var chromId = la[0];
        var blockStart = la[1];
        var blockEnd = la[2];
        var itemStep = la[3];
        var itemSpan = la[4];
        var blockType = ba[20];
        var itemCount = sa[11];
        
        if (blockType == BIG_WIG_TYPE_FSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var score = fa[i + 6];
                var fmin = blockStart + (i*itemStep) + 1, fmax = blockStart + (i*itemStep) + itemSpan;
                if (filter(chromId, fmin, fmax))
                    createFeature(chromId, fmin, fmax, {score: score});
            }
        } else if (blockType == BIG_WIG_TYPE_VSTEP) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[(i*2) + 6] + 1;
                var end = start + itemSpan - 1;
                var score = fa[(i*2) + 7];
                if (filter(chromId, start, end))
                    createFeature(chromId, start, end, {score: score});
            }
        } else if (blockType == BIG_WIG_TYPE_GRAPH) {
            for (var i = 0; i < itemCount; ++i) {
                var start = la[(i*3) + 6] + 1;
                var end   = la[(i*3) + 7];
                var score = fa[(i*3) + 8];
                if (start > end) {
                    start = end;
                }
                if (filter(chromId, start, end))
                    createFeature(chromId, start, end, {score: score});
            }
        } else {
            console.log('Currently not handling bwgType=' + blockType);
        }
    } else if (this.bwg.type == 'bigbed') {
        var offset = 0;
        var dfc = this.bwg.definedFieldCount;
        var schema = this.bwg.schema;

        while (offset < ba.length) {
            var chromId = (ba[offset+3]<<24) | (ba[offset+2]<<16) | (ba[offset+1]<<8) | (ba[offset+0]);
            var start = (ba[offset+7]<<24) | (ba[offset+6]<<16) | (ba[offset+5]<<8) | (ba[offset+4]);
            var end = (ba[offset+11]<<24) | (ba[offset+10]<<16) | (ba[offset+9]<<8) | (ba[offset+8]);
            offset += 12;
            var rest = '';
            while (true) {
                var ch = ba[offset++];
                if (ch != 0) {
                    rest += String.fromCharCode(ch);
                } else {
                    break;
                }
            }

            var featureOpts = {};
            
            var bedColumns;
            if (rest.length > 0) {
                bedColumns = rest.split('\t');
            } else {
                bedColumns = [];
            }
            if (bedColumns.length > 0 && dfc > 3) {
                featureOpts.label = bedColumns[0];
            }
            if (bedColumns.length > 1 && dfc > 4) {
                var score = parseInt(bedColumns[1]);
                if (!isNaN(score))
                    featureOpts.score = score;
            }
            if (bedColumns.length > 2 && dfc > 5) {
                featureOpts.orientation = bedColumns[2];
            }
            if (bedColumns.length > 5 && dfc > 8) {
                var color = bedColumns[5];
                if (BED_COLOR_REGEXP.test(color)) {
                    featureOpts.itemRgb = 'rgb(' + color + ')';
                }
            }

            if (bedColumns.length > dfc-3 && schema) {
                for (var col = dfc - 3; col < bedColumns.length; ++col) {
                    featureOpts[schema.fields[col+3].name] = bedColumns[col];
                }
            }

            if (filter(chromId, start + 1, end, bedColumns)) {
                if (dfc < 12) {
                    createFeature(chromId, start + 1, end, featureOpts);
                } else {
                    var thickStart = bedColumns[3]|0;
                    var thickEnd   = bedColumns[4]|0;
                    var blockCount = bedColumns[6]|0;
                    var blockSizes = bedColumns[7].split(',');
                    var blockStarts = bedColumns[8].split(',');

                    if (featureOpts.exonFrames) {
                        var exonFrames = featureOpts.exonFrames.split(',');
                        featureOpts.exonFrames = undefined;
                    }
                    
                    featureOpts.type = 'transcript'
                    var grp = new DASGroup();
                    for (var k in featureOpts) {
                        grp[k] = featureOpts[k];
                    }
                    grp.id = bedColumns[0];
                    grp.segment = this.bwg.idsToChroms[chromId];
                    grp.min = start + 1;
                    grp.max = end;
                    grp.notes = [];
                    featureOpts.groups = [grp];

                    // Moving towards using bigGenePred model, but will
                    // still support old Dalliance-style BED12+gene-name for the
                    // foreseeable future.
                    if (bedColumns.length > 9) {
                        var geneId = featureOpts.geneName || bedColumns[9];
                        var geneName = geneId;
                        if (bedColumns.length > 10) {
                            geneName = bedColumns[10];
                        }
                        if (featureOpts.geneName2)
                            geneName = featureOpts.geneName2;

                        var gg = shallowCopy(grp);
                        gg.id = geneId;
                        gg.label = geneName;
                        gg.type = 'gene';
                        featureOpts.groups.push(gg);
                    }

                    var spanList = [];
                    for (var b = 0; b < blockCount; ++b) {
                        var bmin = (blockStarts[b]|0) + start;
                        var bmax = bmin + (blockSizes[b]|0);
                        var span = new Range(bmin, bmax);
                        spanList.push(span);
                    }
                    var spans = union(spanList);
                    
                    var tsList = spans.ranges();
                    for (var s = 0; s < tsList.length; ++s) {
                        var ts = tsList[s];
                        createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                    }

                    if (thickEnd > thickStart) {
                        var codingRegion = (featureOpts.orientation == '+') ?
                            new Range(thickStart, thickEnd + 3) :
                            new Range(thickStart - 3, thickEnd);
                            // +/- 3 to account for stop codon

                        var tl = intersection(spans, codingRegion);
                        if (tl) {
                            featureOpts.type = 'translation';
                            var tlList = tl.ranges();
                            var readingFrame = 0;

                            var tlOffset = 0;
                            while (tlList[0].min() > tsList[tlOffset].max())
                                tlOffset++;

                            for (var s = 0; s < tlList.length; ++s) {
                                // Record reading frame for every exon
                                var index = s;
                                if (featureOpts.orientation == '-')
                                    index = tlList.length - s - 1;
                                var ts = tlList[index];
                                featureOpts.readframe = readingFrame;
                                if (exonFrames) {
                                    var brf = parseInt(exonFrames[index + tlOffset]);
                                    if (typeof(brf) === 'number' && brf >= 0 && brf <= 2) {
                                        featureOpts.readframe = brf;
                                        featureOpts.readframeExplicit = true;
                                    }
                                }
                                var length = ts.max() - ts.min();
                                readingFrame = (readingFrame + length) % 3;
                                createFeature(chromId, ts.min() + 1, ts.max(), featureOpts);
                            }
                        }
                    }
                }
            }
        }
    } else {
        throw Error("Don't know what to do with " + this.bwg.type);
    }
}

//
// nasty cut/paste, should roll back in!
//

BigWigView.prototype.getFirstAdjacent = function(chrName, pos, dir, callback) {
    var chr = this.bwg.chromsToIDs[chrName];
    if (chr === undefined) {
        // Not an error because some .bwgs won't have data for all chromosomes.
        return callback([]);
    } else {
        this.getFirstAdjacentById(chr, pos, dir, callback);
    }
}

BigWigView.prototype.getFirstAdjacentById = function(chr, pos, dir, callback) {
    var thisB = this;
    if (!this.cirHeader) {
        this.bwg.data.slice(this.cirTreeOffset, 48).fetch(function(result) {
            thisB.cirHeader = result;
            var la = new Int32Array(thisB.cirHeader);
            thisB.cirBlockSize = la[1];
            thisB.getFirstAdjacentById(chr, pos, dir, callback);
        });
        return;
    }

    var blockToFetch = null;
    var bestBlockChr = -1;
    var bestBlockOffset = -1;

    var outstanding = 0;

    var beforeBWG = Date.now();

    var cirFobRecur = function(offset, level) {
        outstanding += offset.length;

        var maxCirBlockSpan = 4 +  (thisB.cirBlockSize * 32);   // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], offset[i] + maxCirBlockSpan);
            spans = spans ? union(spans, blockSpan) : blockSpan;
        }
        
        var fetchRanges = spans.ranges();
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            cirFobStartFetch(offset, fr, level);
        }
    }

    var cirFobStartFetch = function(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        thisB.bwg.data.slice(fr.min(), fr.max() - fr.min()).fetch(function(resultBuffer) {
            for (var i = 0; i < offset.length; ++i) {
                if (fr.contains(offset[i])) {
                    cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);
                    --outstanding;
                    if (outstanding == 0) {
                        if (!blockToFetch) {
                            if (dir > 0 && (chr != 0 || pos > 0)) {
                                return thisB.getFirstAdjacentById(0, 0, dir, callback);
                            } else if (dir < 0 && (chr != thisB.bwg.maxID || pos < 1000000000)) {
                                return thisB.getFirstAdjacentById(thisB.bwg.maxID, 1000000000, dir, callback);
                            }
                            return callback([]);
                        }

                        thisB.fetchFeatures(function(chrx, fmin, fmax, toks) {
                            return (dir < 0 && (chrx < chr || fmax < pos)) || (dir > 0 && (chrx > chr || fmin > pos));
                        }, [blockToFetch], function(features) {
                            var bestFeature = null;
                            var bestChr = -1;
                            var bestPos = -1;
                            for (var fi = 0; fi < features.length; ++fi) {
                                var f = features[fi];
                                var chrx = f._chromId, fmin = f.min, fmax = f.max;
                                if (bestFeature == null || ((dir < 0) && (chrx > bestChr || fmax > bestPos)) || ((dir > 0) && (chrx < bestChr || fmin < bestPos))) {
                                    bestFeature = f;
                                    bestPos = (dir < 0) ? fmax : fmin;
                                    bestChr = chrx;
                                }
                            }

                            if (bestFeature != null) 
                                return callback([bestFeature]);
                            else
                                return callback([]);
                        });
                    }
                }
            }
        });
    }

    var cirFobRecur2 = function(cirBlockData, offset, level) {
        var ba = new Uint8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset/2 + 1];
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = bwg_readOffset(ba, offset+16);
                var blockSize = bwg_readOffset(ba, offset+24);
                if ((dir < 0 && ((startChrom < chr || (startChrom == chr && startBase <= pos)))) ||
                    (dir > 0 && ((endChrom > chr || (endChrom == chr && endBase >= pos)))))
                {
                    // console.log('Got an interesting block: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                    if (/_random/.exec(thisB.bwg.idsToChroms[startChrom])) {
                        // dlog('skipping random: ' + thisB.bwg.idsToChroms[startChrom]);
                    } else if (blockToFetch == null || ((dir < 0) && (endChrom > bestBlockChr || (endChrom == bestBlockChr && endBase > bestBlockOffset)) ||
                                                 (dir > 0) && (startChrom < bestBlockChr || (startChrom == bestBlockChr && startBase < bestBlockOffset))))
                    {
                        //                        dlog('best is: startBase=' + startChrom + ':' + startBase + '; endBase=' + endChrom + ':' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                        blockToFetch = {offset: blockOffset, size: blockSize};
                        bestBlockOffset = (dir < 0) ? endBase : startBase;
                        bestBlockChr = (dir < 0) ? endChrom : startChrom;
                    }
                }
                offset += 32;
            }
        } else {
            var bestRecur = -1;
            var bestPos = -1;
            var bestChr = -1;
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = (la[lo + 4]<<32) | (la[lo + 5]);
                if ((dir < 0 && ((startChrom < chr || (startChrom == chr && startBase <= pos)) &&
                                 (endChrom   >= chr))) ||
                     (dir > 0 && ((endChrom > chr || (endChrom == chr && endBase >= pos)) &&
                                  (startChrom <= chr))))
                {
                    if (bestRecur < 0 || endBase > bestPos) {
                        bestRecur = blockOffset;
                        bestPos = (dir < 0) ? endBase : startBase;
                        bestChr = (dir < 0) ? endChrom : startChrom;
                    }
                }
                offset += 24;
            }
            if (bestRecur >= 0) {
                cirFobRecur([bestRecur], level + 1);
            }
        }
    };
    

    cirFobRecur([thisB.cirTreeOffset + 48], 1);
}

BigWig.prototype.readWigData = function(chrName, min, max, callback) {
    this.getUnzoomedView().readWigData(chrName, min, max, callback);
}

BigWig.prototype.getUnzoomedView = function() {
    if (!this.unzoomedView) {
        var cirLen = 4000;
        var nzl = this.zoomLevels[0];
        if (nzl) {
            cirLen = this.zoomLevels[0].dataOffset - this.unzoomedIndexOffset;
        }
        this.unzoomedView = new BigWigView(this, this.unzoomedIndexOffset, cirLen, false);
    }
    return this.unzoomedView;
}

BigWig.prototype.getZoomedView = function(z) {
    var zh = this.zoomLevels[z];
    if (!zh.view) {
        zh.view = new BigWigView(this, zh.indexOffset, /* this.zoomLevels[z + 1].dataOffset - zh.indexOffset */ 4000, true);
    }
    return zh.view;
}

function makeBwg(data, callback, name) {
    var bwg = new BigWig();
    bwg.data = data;
    bwg.name = name;
    bwg.data.slice(0, 512).salted().fetch(function(result) {
        if (!result) {
            return callback(null, "Couldn't fetch file");
        }

        var header = result;
        var ba = new Uint8Array(header);
        var sa = new Int16Array(header);
        var la = new Int32Array(header);
        var magic = ba[0] + (M1 * ba[1]) + (M2 * ba[2]) + (M3 * ba[3]);
        if (magic == BIG_WIG_MAGIC) {
            bwg.type = 'bigwig';
        } else if (magic == BIG_BED_MAGIC) {
            bwg.type = 'bigbed';
        } else if (magic == BIG_WIG_MAGIC_BE || magic == BIG_BED_MAGIC_BE) {
            return callback(null, "Currently don't support big-endian BBI files");
            
        } else {
            return callback(null, "Not a supported format, magic=0x" + magic.toString(16));
            
        }

        bwg.version = sa[2];             // 4
        bwg.numZoomLevels = sa[3];       // 6
        bwg.chromTreeOffset = bwg_readOffset(ba, 8);
        bwg.unzoomedDataOffset = bwg_readOffset(ba, 16);
        bwg.unzoomedIndexOffset = bwg_readOffset(ba, 24);
        bwg.fieldCount = sa[16];         // 32
        bwg.definedFieldCount = sa[17];  // 34
        bwg.asOffset = bwg_readOffset(ba, 36);
        bwg.totalSummaryOffset = bwg_readOffset(ba, 44);
        bwg.uncompressBufSize = la[13];  // 52
        bwg.extHeaderOffset = bwg_readOffset(ba, 56);

        bwg.zoomLevels = [];
        for (var zl = 0; zl < bwg.numZoomLevels; ++zl) {
            var zlReduction = la[zl*6 + 16]
            var zlData = bwg_readOffset(ba, zl*24 + 72);
            var zlIndex = bwg_readOffset(ba, zl*24 + 80);
            bwg.zoomLevels.push({reduction: zlReduction, dataOffset: zlData, indexOffset: zlIndex});
        }

        bwg.readChromTree(function() {
            bwg.getAutoSQL(function(as) {
                bwg.schema = as;
                return callback(bwg);
            });
        });
    }, {timeout: 5000});    // Potential timeout on first request to catch mixed-content errors on
                            // Chromium.
}


BigWig.prototype._tsFetch = function(zoom, chr, min, max, callback) {
    var bwg = this;
    if (zoom >= this.zoomLevels.length - 1) {
        if (!this.topLevelReductionCache) {
            this.getZoomedView(this.zoomLevels.length - 1).readWigDataById(-1, 0, 300000000, function(feats) {
                bwg.topLevelReductionCache = feats;
                return bwg._tsFetch(zoom, chr, min, max, callback);
            });
        } else {
            var f = [];
            var c = this.topLevelReductionCache;
            for (var fi = 0; fi < c.length; ++fi) {
                if (c[fi]._chromId == chr) {
                    f.push(c[fi]);
                }
            }
            return callback(f);
        }
    } else {
        var view;
        if (zoom < 0) {
            view = this.getUnzoomedView();
        } else {
            view = this.getZoomedView(zoom);
        }
        return view.readWigDataById(chr, min, max, callback);
    }
}

BigWig.prototype.thresholdSearch = function(chrName, referencePoint, dir, threshold, callback) {
    dir = (dir<0) ? -1 : 1;
    var bwg = this;
    var initialChr = this.chromsToIDs[chrName];
    var candidates = [{chrOrd: 0, chr: initialChr, zoom: bwg.zoomLevels.length - 4, min: 0, max: 300000000, fromRef: true}]
    for (var i = 1; i <= this.maxID + 1; ++i) {
        var chrId = (initialChr + (dir*i)) % (this.maxID + 1);
        if (chrId < 0) 
            chrId += (this.maxID + 1);
        candidates.push({chrOrd: i, chr: chrId, zoom: bwg.zoomLevels.length - 1, min: 0, max: 300000000})
    }
       
    function fbThresholdSearchRecur() {
    	if (candidates.length == 0) {
    	    return callback(null);
    	}
    	candidates.sort(function(c1, c2) {
    	    var d = c1.zoom - c2.zoom;
    	    if (d != 0)
    		    return d;

            d = c1.chrOrd - c2.chrOrd;
            if (d != 0)
                return d;
    	    else
    		    return c1.min - c2.min * dir;
    	});

	    var candidate = candidates.splice(0, 1)[0];
        bwg._tsFetch(candidate.zoom, candidate.chr, candidate.min, candidate.max, function(feats) {
            var rp = dir > 0 ? 0 : 300000000;
            if (candidate.fromRef)
                rp = referencePoint;
            
            for (var fi = 0; fi < feats.length; ++fi) {
    	        var f = feats[fi];
                var score;
                if (f.maxScore != undefined)
                    score = f.maxScore;
                else
                    score = f.score;

                if (dir > 0) {
    	            if (score > threshold) {
        		        if (candidate.zoom < 0) {
        		            if (f.min > rp)
                                return callback(f);
        		        } else if (f.max > rp) {
        		            candidates.push({chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef});
        		        }
                    }
                } else {
                    if (score > threshold) {
            		    if (candidate.zoom < 0) {
                	        if (f.max < rp)
                			    return callback(f);
                        } else if (f.min < rp) {
                            candidates.push({chr: candidate.chr, chrOrd: candidate.chrOrd, zoom: candidate.zoom - 2, min: f.min, max: f.max, fromRef: candidate.fromRef});
                        }
    	            }
                }
    	    }
            fbThresholdSearchRecur();
        });
    }
    
    fbThresholdSearchRecur();
}

BigWig.prototype.getAutoSQL = function(callback) {
    var thisB = this;
    if (!this.asOffset)
        return callback(null);


    this.data.slice(this.asOffset, 2048).fetch(function(result) {
        var ba = new Uint8Array(result);
        var s = '';
        for (var i = 0; i < ba.length; ++i) {
            if (ba[i] == 0)
                break;
            s += String.fromCharCode(ba[i]);
        }
        
        /* 
         * Quick'n'dirty attempt to parse autoSql format.
         * See: http://www.linuxjournal.com/files/linuxjournal.com/linuxjournal/articles/059/5949/5949l2.html
         */

        var header_re = /(\w+)\s+(\w+)\s+("([^"]+)")?\s+\(\s*/;
        var field_re = /([\w\[\]]+)\s+(\w+)\s*;\s*("([^"]+)")?\s*/g;

        var headerMatch = header_re.exec(s);
        if (headerMatch) {
            var as = {
                declType: headerMatch[1],
                name: headerMatch[2],
                comment: headerMatch[4],

                fields: []
            };

            s = s.substring(headerMatch[0]);
            for (var m = field_re.exec(s); m != null; m = field_re.exec(s)) {
                as.fields.push({type: m[1],
                             name: m[2],
                             comment: m[4]});
            }

            return callback(as);
        }
    });
}

BigWig.prototype.getExtraIndices = function(callback) {
    var thisB = this;
    if (this.version < 4 || this.extHeaderOffset == 0 || this.type != 'bigbed') {
        return callback(null);
    } else {
        this.data.slice(this.extHeaderOffset, 64).fetch(function(result) {
            if (!result) {
                return callback(null, "Couldn't fetch extension header");
            }

            var ba = new Uint8Array(result);
            var sa = new Int16Array(result);
            var la = new Int32Array(result);
            
            var extHeaderSize = sa[0];
            var extraIndexCount = sa[1];
            var extraIndexListOffset = bwg_readOffset(ba, 4);

            if (extraIndexCount == 0) {
                return callback(null);
            }

            // FIXME 20byte records only make sense for single-field indices.
            // Right now, these seem to be the only things around, but the format
            // is actually more general.
            thisB.data.slice(extraIndexListOffset, extraIndexCount * 20).fetch(function(eil) {
                if (!eil) {
                    return callback(null, "Couldn't fetch index info");
                }

                var ba = new Uint8Array(eil);
                var sa = new Int16Array(eil);
                var la = new Int32Array(eil);

                var indices = [];
                for (var ii = 0; ii < extraIndexCount; ++ii) {
                    var eiType = sa[ii*10];
                    var eiFieldCount = sa[ii*10 + 1];
                    var eiOffset = bwg_readOffset(ba, ii*20 + 4);
                    var eiField = sa[ii*10 + 8]
                    var index = new BBIExtraIndex(thisB, eiType, eiFieldCount, eiOffset, eiField);
                    indices.push(index);
                }
                callback(indices);
            });
        });
    }
}

function BBIExtraIndex(bbi, type, fieldCount, offset, field) {
    this.bbi = bbi;
    this.type = type;
    this.fieldCount = fieldCount;
    this.offset = offset;
    this.field = field;
}

BBIExtraIndex.prototype.lookup = function(name, callback) {
    var thisB = this;

    this.bbi.data.slice(this.offset, 32).fetch(function(bpt) {
        var ba = new Uint8Array(bpt);
        var sa = new Int16Array(bpt);
        var la = new Int32Array(bpt);
        var bptMagic = la[0];
        var blockSize = la[1];
        var keySize = la[2];
        var valSize = la[3];
        var itemCount = bwg_readOffset(ba, 16);
        var rootNodeOffset = 32;

        function bptReadNode(nodeOffset) {
            thisB.bbi.data.slice(nodeOffset, 4 + (blockSize * (keySize + valSize))).fetch(function(node) {
                var ba = new Uint8Array(node);
                var sa = new Uint16Array(node);
                var la = new Uint32Array(node);

                var nodeType = ba[0];
                var cnt = sa[1];

                var offset = 4;
                if (nodeType == 0) {
                    var lastChildOffset = null;
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }

                        var childOffset = bwg_readOffset(ba, offset);
                        offset += 8;
                        
                        if (name.localeCompare(key) < 0 && lastChildOffset) {
                            bptReadNode(lastChildOffset);
                            return;
                        }
                        lastChildOffset = childOffset;
                    }
                    bptReadNode(lastChildOffset);
                } else {
                    for (var n = 0; n < cnt; ++n) {
                        var key = '';
                        for (var ki = 0; ki < keySize; ++ki) {
                            var charCode = ba[offset++];
                            if (charCode != 0) {
                                key += String.fromCharCode(charCode);
                            }
                        }
                        
                        // Specific for EI case.
                        if (key == name) {
                            var start = bwg_readOffset(ba, offset);
                            var length = readInt(ba, offset + 8);

                            return thisB.bbi.getUnzoomedView().fetchFeatures(
                                function(chr, min, max, toks) {
                                    if (toks && toks.length > thisB.field - 3)
                                        return toks[thisB.field - 3] == name;
                                }, 
                                [{offset: start, size: length}], 
                                callback);
                        }
                        offset += valSize;
                    }
                    return callback([]);
                }
            });
        }

        bptReadNode(thisB.offset + rootNodeOffset);
    });
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        makeBwg: makeBwg,
        BIG_BED_MAGIC: BIG_BED_MAGIC,
        BIG_WIG_MAGIC: BIG_WIG_MAGIC
    }
}

},{"./bin":3,"./das":5,"./spans":10,"./utils":11,"jszlib":24}],3:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// bin.js general binary data support
//

"use strict";

if (typeof(require) !== 'undefined') {
    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;

    var sha1 = require('./sha1');
    var b64_sha1 = sha1.b64_sha1;

    var Promise = require('es6-promise').Promise;
}

function BlobFetchable(b) {
    this.blob = b;
}

BlobFetchable.prototype.slice = function(start, length) {
    var b;

    if (this.blob.slice) {
        if (length) {
            b = this.blob.slice(start, start + length);
        } else {
            b = this.blob.slice(start);
        }
    } else {
        if (length) {
            b = this.blob.webkitSlice(start, start + length);
        } else {
            b = this.blob.webkitSlice(start);
        }
    }
    return new BlobFetchable(b);
}

BlobFetchable.prototype.salted = function() {return this;}

if (typeof(FileReader) !== 'undefined') {
    // console.log('defining async BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function(callback) {
        var reader = new FileReader();
        reader.onloadend = function(ev) {
            callback(bstringToBuffer(reader.result));
        };
        reader.readAsBinaryString(this.blob);
    }

} else {
    // if (console && console.log)
    //    console.log('defining sync BlobFetchable.fetch');

    BlobFetchable.prototype.fetch = function(callback) {
        var reader = new FileReaderSync();
        try {
            var res = reader.readAsArrayBuffer(this.blob);
            callback(res);
        } catch (e) {
            callback(null, e);
        }
    }
}

function URLFetchable(url, start, end, opts) {
    if (!opts) {
        if (typeof start === 'object') {
            opts = start;
            start = undefined;
        } else {
            opts = {};
        }
    }

    this.url = url;
    this.start = start || 0;
    if (end) {
        this.end = end;
    }
    this.opts = opts;
}

URLFetchable.prototype.slice = function(s, l) {
    if (s < 0) {
        throw 'Bad slice ' + s;
    }

    var ns = this.start, ne = this.end;
    if (ns && s) {
        ns = ns + s;
    } else {
        ns = s || ns;
    }
    if (l && ns) {
        ne = ns + l - 1;
    } else {
        ne = ne || l - 1;
    }
    return new URLFetchable(this.url, ns, ne, this.opts);
}

var seed=0;
var isSafari = navigator.userAgent.indexOf('Safari') >= 0 && navigator.userAgent.indexOf('Chrome') < 0 ;

URLFetchable.prototype.fetchAsText = function(callback) {
    try {
        var req = new XMLHttpRequest();
        var length;
        var url = this.url;
        if ((isSafari || this.opts.salt) && url.indexOf('?') < 0) {
            url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++seed));
        }
        req.open('GET', url, true);

        if (this.end) {
            if (this.end - this.start > 100000000) {
                throw 'Monster fetch!';
            }
            req.setRequestHeader('Range', 'bytes=' + this.start + '-' + this.end);
            length = this.end - this.start + 1;
        }

        req.onreadystatechange = function() {
            if (req.readyState == 4) {
                if (req.status == 200 || req.status == 206) {
                    return callback(req.responseText);
                } else {
                    return callback(null);
                }
            }
        };
        if (this.opts.credentials) {
            req.withCredentials = true;
        }
        req.send('');
    } catch (e) {
        return callback(null);
    }
}

URLFetchable.prototype.salted = function() {
    var o = shallowCopy(this.opts);
    o.salt = true;
    return new URLFetchable(this.url, this.start, this.end, o);
}

URLFetchable.prototype.getURL = function() {
    if (this.opts.resolver) {
        return this.opts.resolver(this.url).then(function (urlOrObj) {
            if (typeof urlOrObj === 'string') {
                return urlOrObj;
            } else {
                return urlOrObj.url;
            }
        });
    } else {
        return Promise.resolve(this.url);
    }
}

URLFetchable.prototype.fetch = function(callback, opts) {
    var thisB = this;
 
    opts = opts || {};
    var attempt = opts.attempt || 1;
    var truncatedLength = opts.truncatedLength;
    if (attempt > 3) {
        return callback(null);
    }

    this.getURL().then(function(url) {
        try {
            var timeout;
            if (opts.timeout && !thisB.opts.credentials) {
                timeout = setTimeout(
                    function() {
                        console.log('timing out ' + url);
                        req.abort();
                        return callback(null, 'Timeout');
                    },
                    opts.timeout
                );
            }
            
            var req = new XMLHttpRequest();
            var length;
            if ((isSafari || thisB.opts.salt) && url.indexOf('?') < 0) {
                url = url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++seed));
            }
            req.open('GET', url, true);
            req.overrideMimeType('text/plain; charset=x-user-defined');
            if (thisB.end) {
                if (thisB.end - thisB.start > 100000000) {
                    throw 'Monster fetch!';
                }
                req.setRequestHeader('Range', 'bytes=' + thisB.start + '-' + thisB.end);
                length = thisB.end - thisB.start + 1;
            }
            req.responseType = 'arraybuffer';
            req.onreadystatechange = function() {
                if (req.readyState == 4) {
                    if (timeout)
                        clearTimeout(timeout);
                    if (req.status == 200 || req.status == 206) {
                        if (req.response) {
                            var bl = req.response.byteLength;
                            if (length && length != bl && (!truncatedLength || bl != truncatedLength)) {
                                return thisB.fetch(callback, {attempt: attempt + 1, truncatedLength: bl});
                            } else {
                                return callback(req.response);
                            }
                        } else if (req.mozResponseArrayBuffer) {
                            return callback(req.mozResponseArrayBuffer);
                        } else {
                            var r = req.responseText;
                            if (length && length != r.length && (!truncatedLength || r.length != truncatedLength)) {
                                return thisB.fetch(callback, {attempt: attempt + 1, truncatedLength: r.length});
                            } else {
                                return callback(bstringToBuffer(req.responseText));
                            }
                        }
                    } else {
                        return thisB.fetch(callback, {attempt: attempt + 1});
                    }
                }
            };
            if (thisB.opts.credentials) {
                req.withCredentials = true;
            }
            req.send('');
        } catch (e) {
            return callback(null);
        }
    }).catch(function(err) {
        console.log(err);
        return callback(null, err);
    });
}
                       
function bstringToBuffer(result) {
    if (!result) {
        return null;
    }

    var ba = new Uint8Array(result.length);
    for (var i = 0; i < ba.length; ++i) {
        ba[i] = result.charCodeAt(i);
    }
    return ba.buffer;
}

// Read from Uint8Array

(function(global) {
    var convertBuffer = new ArrayBuffer(8);
    var ba = new Uint8Array(convertBuffer);
    var fa = new Float32Array(convertBuffer);


    global.readFloat = function(buf, offset) {
        ba[0] = buf[offset];
        ba[1] = buf[offset+1];
        ba[2] = buf[offset+2];
        ba[3] = buf[offset+3];
        return fa[0];
    };
 }(this));

function readInt64(ba, offset) {
    return (ba[offset + 7] << 24) | (ba[offset + 6] << 16) | (ba[offset + 5] << 8) | (ba[offset + 4]);
}

function readInt(ba, offset) {
    return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
}

function readShort(ba, offset) {
    return (ba[offset + 1] << 8) | (ba[offset]);
}

function readByte(ba, offset) {
    return ba[offset];
}

function readIntBE(ba, offset) {
    return (ba[offset] << 24) | (ba[offset + 1] << 16) | (ba[offset + 2] << 8) | (ba[offset + 3]);
}

// Exports if we are being used as a module

if (typeof(module) !== 'undefined') {
    module.exports = {
        BlobFetchable: BlobFetchable,
        URLFetchable: URLFetchable,

        readInt: readInt,
        readIntBE: readIntBE,
        readInt64: readInt64,
        readShort: readShort,
        readByte: readByte,
        readFloat: this.readFloat
    }
}

},{"./sha1":9,"./utils":11,"es6-promise":12}],4:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// color.js
//

"use strict";

function DColour(red, green, blue, name) {
    this.red = red|0;
    this.green = green|0;
    this.blue = blue|0;
    if (name) {
        this.name = name;
    }
}

DColour.prototype.toSvgString = function() {
    if (!this.name) {
        this.name = "rgb(" + this.red + "," + this.green + "," + this.blue + ")";
    }

    return this.name;
}

function hex2(x) {
    var y = '00' + x.toString(16);
    return y.substring(y.length - 2);
}

DColour.prototype.toHexString = function() {
    return '#' + hex2(this.red) + hex2(this.green) + hex2(this.blue);
}

var palette = {
    red: new DColour(255, 0, 0, 'red'),
    green: new DColour(0, 255, 0, 'green'),
    blue: new DColour(0, 0, 255, 'blue'),
    yellow: new DColour(255, 255, 0, 'yellow'),
    white: new DColour(255, 255, 255, 'white'),
    black: new DColour(0, 0, 0, 'black'),
    gray: new DColour(180, 180, 180, 'gray'),
    grey: new DColour(180, 180, 180, 'grey'),
    lightskyblue: new DColour(135, 206, 250, 'lightskyblue'),
    lightsalmon: new DColour(255, 160, 122, 'lightsalmon'),
    hotpink: new DColour(255, 105, 180, 'hotpink')
};

var COLOR_RE = new RegExp('^#([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})$');
var CSS_COLOR_RE = /rgb\(([0-9]+),([0-9]+),([0-9]+)\)/

function dasColourForName(name) {
    var c = palette[name];
    if (!c) {
        var match = COLOR_RE.exec(name);
        if (match) {
            c = new DColour(('0x' + match[1])|0, ('0x' + match[2])|0, ('0x' + match[3])|0, name);
            palette[name] = c;
        } else {
    	    match = CSS_COLOR_RE.exec(name);
    	    if (match) {
        		c = new DColour(match[1]|0, match[2]|0, match[3]|0, name);
        		palette[name] = c;
	       } else {
		      console.log("couldn't handle color: " + name);
		      c = palette.black;
		      palette[name] = c;
	       }
        }
    }
    return c;
}

function makeColourSteps(steps, stops, colours) {
    var dcolours = [];
    for (var ci = 0; ci < colours.length; ++ci) {
        dcolours.push(dasColourForName(colours[ci]));
    }

    var grad = [];
  STEP_LOOP:
    for (var si = 0; si < steps; ++si) {
        var rs = (1.0 * si) / (steps-1);
        var score = stops[0] + (stops[stops.length -1] - stops[0]) * rs;
        for (var i = 0; i < stops.length - 1; ++i) {
            if (score >= stops[i] && score <= stops[i+1]) {
                var frac = (score - stops[i]) / (stops[i+1] - stops[i]);
                var ca = dcolours[i];
                var cb = dcolours[i+1];

                var fill = new DColour(
                    ((ca.red * (1.0 - frac)) + (cb.red * frac))|0,
                    ((ca.green * (1.0 - frac)) + (cb.green * frac))|0,
                    ((ca.blue * (1.0 - frac)) + (cb.blue * frac))|0
                ).toSvgString();
                grad.push(fill);

                continue STEP_LOOP;
            }
        }
        throw 'Bad step';
    }

    return grad;
}

function makeGradient(steps, color1, color2, color3) {
    if (color3) {
        return makeColourSteps(steps, [0, 0.5, 1], [color1, color2, color3]);
    } else {
        return makeColourSteps(steps, [0, 1], [color1, color2]);
    }
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        makeColourSteps: makeColourSteps,
        makeGradient: makeGradient,
        dasColourForName: dasColourForName
    };
}

},{}],5:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// das.js: queries and low-level data model.
//

"use strict";

if (typeof(require) !== 'undefined') {
    var utils = require('./utils');
    var shallowCopy = utils.shallowCopy;
    var pusho = utils.pusho;

    var color = require('./color');
    var makeColourSteps = color.makeColourSteps;
}

var dasLibErrorHandler = function(errMsg) {
    alert(errMsg);
}
var dasLibRequestQueue = new Array();

function DASSegment(name, start, end, description) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.description = description;
}
DASSegment.prototype.toString = function() {
    return this.name + ':' + this.start + '..' + this.end;
};
DASSegment.prototype.isBounded = function() {
    return this.start && this.end;
}
DASSegment.prototype.toDASQuery = function() {
    var q = 'segment=' + this.name;
    if (this.start && this.end) {
        q += (':' + this.start + ',' + this.end);
    }
    return q;
}


function DASSource(a1, a2) {
    var options;
    if (typeof a1 == 'string') {
        this.uri = a1;
        options = a2 || {};
    } else {
        options = a1 || {};
    }
    for (var k in options) {
        this[k] = options[k];
    }

    if (!this.coords) {
        this.coords = [];
    }
    if (!this.props) {
        this.props = {};
    }

    this.dasBaseURI = this.uri;
    if (this.dasBaseURI && this.dasBaseURI.substr(this.uri.length - 1) != '/') {
        this.dasBaseURI = this.dasBaseURI + '/';
    }
}

DASSource.prototype.getURI = function(uri) {
    if (this.resolver) {
        return this.resolver(uri).then(function (urlOrObj) {
            if (typeof urlOrObj === 'string') {
                return urlOrObj;
            } else {
                return urlOrObj.url;
            }
        });
    } else {
        return Promise.resolve(uri);
    }
}

function DASCoords() {
}

function coordsMatch(c1, c2) {
    return c1.taxon == c2.taxon && c1.auth == c2.auth && c1.version == c2.version;
}

//
// DAS 1.6 entry_points command
//

DASSource.prototype.entryPoints = function(callback) {
    var dasURI = this.dasBaseURI + 'entry_points';
    this.doCrossDomainRequest(dasURI, function(responseXML) {
            if (!responseXML) {
                return callback([]);
            }

                var entryPoints = new Array();
                
                var segs = responseXML.getElementsByTagName('SEGMENT');
                for (var i = 0; i < segs.length; ++i) {
                    var seg = segs[i];
                    var segId = seg.getAttribute('id');
                    
                    var segSize = seg.getAttribute('size');
                    var segMin, segMax;
                    if (segSize) {
                        segMin = 1; segMax = segSize|0;
                    } else {
                        segMin = seg.getAttribute('start');
                        if (segMin) {
                            segMin |= 0;
                        }
                        segMax = seg.getAttribute('stop');
                        if (segMax) {
                            segMax |= 0;
                        }
                    }
                    var segDesc = null;
                    if (seg.firstChild) {
                        segDesc = seg.firstChild.nodeValue;
                    }
                    entryPoints.push(new DASSegment(segId, segMin, segMax, segDesc));
                }          
               callback(entryPoints);
    });         
}

//
// DAS 1.6 sequence command
// Do we need an option to fall back to the dna command?
//

function DASSequence(name, start, end, alpha, seq) {
    this.name = name;
    this.start = start;
    this.end = end;
    this.alphabet = alpha;
    this.seq = seq;
}

DASSource.prototype.sequence = function(segment, callback) {
    var dasURI = this.dasBaseURI + 'sequence?' + segment.toDASQuery();
    this.doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            callback([]);
            return;
        } else {
                var seqs = new Array();
                
                var segs = responseXML.getElementsByTagName('SEQUENCE');
                for (var i = 0; i < segs.length; ++i) {
                    var seg = segs[i];
                    var segId = seg.getAttribute('id');
                    var segMin = seg.getAttribute('start');
                    var segMax = seg.getAttribute('stop');
                    var segAlpha = 'DNA';
                    var segSeq = null;
                    if (seg.firstChild) {
                        var rawSeq = seg.firstChild.nodeValue;
                        segSeq = '';
                        var idx = 0;
                        while (true) {
                            var space = rawSeq.indexOf('\n', idx);
                            if (space >= 0) {
                                segSeq += rawSeq.substring(idx, space).toUpperCase();
                                idx = space + 1;
                            } else {
                                segSeq += rawSeq.substring(idx).toUpperCase();
                                break;
                            }
                        }
                    }
                    seqs.push(new DASSequence(segId, segMin, segMax, segAlpha, segSeq));
                }
                
                callback(seqs);
        }
    });
}

//
// DAS 1.6 features command
//

function DASFeature() {
}

function DASGroup(id) {
    if (id)
        this.id = id;
}

function DASLink(desc, uri) {
    this.desc = desc;
    this.uri = uri;
}

DASSource.prototype.features = function(segment, options, callback) {
    options = options || {};
    var thisB = this;

    var dasURI;
    if (this.features_uri) {
        dasURI = this.features_uri;
    } else {
        var filters = [];

        if (segment) {
            filters.push(segment.toDASQuery());
        } else if (options.group) {
            var g = options.group;
            if (typeof g == 'string') {
                filters.push('group_id=' + g);
            } else {
                for (var gi = 0; gi < g.length; ++gi) {
                    filters.push('group_id=' + g[gi]);
                }
            }
        }

        if (options.adjacent) {
            var adj = options.adjacent;
            if (typeof adj == 'string') {
                adj = [adj];
            }
            for (var ai = 0; ai < adj.length; ++ai) {
                filters.push('adjacent=' + adj[ai]);
            }
        }

        if (options.type) {
            if (typeof options.type == 'string') {
                filters.push('type=' + options.type);
            } else {
                for (var ti = 0; ti < options.type.length; ++ti) {
                    filters.push('type=' + options.type[ti]);
                }
            }
        }
        
        if (options.maxbins) {
            filters.push('maxbins=' + options.maxbins);
        }
        
        if (filters.length > 0) {
            dasURI = this.dasBaseURI + 'features?' + filters.join(';');
        } else {
            callback([], 'No filters specified');
        }
    } 
   

    this.doCrossDomainRequest(dasURI, function(responseXML, req) {
        if (!responseXML) {
            var msg;
            if (req.status == 0) {
                msg = 'server may not support CORS';
            } else {
                msg = 'status=' + req.status;
            }
            callback([], 'Failed request: ' + msg);
            return;
        }
/*      if (req) {
            var caps = req.getResponseHeader('X-DAS-Capabilties');
            if (caps) {
                alert(caps);
            }
        } */

        var features = new Array();
        var segmentMap = {};

        var segs = responseXML.getElementsByTagName('SEGMENT');
        for (var si = 0; si < segs.length; ++si) {
            var segmentXML = segs[si];
            var segmentID = segmentXML.getAttribute('id');
            segmentMap[segmentID] = {
                min: segmentXML.getAttribute('start'),
                max: segmentXML.getAttribute('stop')
            };
            
            var featureXMLs = segmentXML.getElementsByTagName('FEATURE');
            for (var i = 0; i < featureXMLs.length; ++i) {
                var feature = featureXMLs[i];
                var dasFeature = new DASFeature();
                
                dasFeature.segment = segmentID;
                dasFeature.id = feature.getAttribute('id');
                dasFeature.label = feature.getAttribute('label');


/*
                var childNodes = feature.childNodes;
                for (var c = 0; c < childNodes.length; ++c) {
                    var cn = childNodes[c];
                    if (cn.nodeType == Node.ELEMENT_NODE) {
                        var key = cn.tagName;
                        //var val = null;
                        //if (cn.firstChild) {
                        //   val = cn.firstChild.nodeValue;
                        //}
                        dasFeature[key] = 'x';
                    }
                } */


                var spos = elementValue(feature, "START");
                var epos = elementValue(feature, "END");
                if ((spos|0) > (epos|0)) {
                    dasFeature.min = epos|0;
                    dasFeature.max = spos|0;
                } else {
                    dasFeature.min = spos|0;
                    dasFeature.max = epos|0;
                }
                {
                    var tec = feature.getElementsByTagName('TYPE');
                    if (tec.length > 0) {
                        var te = tec[0];
                        if (te.firstChild) {
                            dasFeature.type = te.firstChild.nodeValue;
                        }
                        dasFeature.typeId = te.getAttribute('id');
                        dasFeature.typeCv = te.getAttribute('cvId');
                    }
                }
                dasFeature.type = elementValue(feature, "TYPE");
                if (!dasFeature.type && dasFeature.typeId) {
                    dasFeature.type = dasFeature.typeId; // FIXME?
                }
                
                dasFeature.method = elementValue(feature, "METHOD");
                {
                    var ori = elementValue(feature, "ORIENTATION");
                    if (!ori) {
                        ori = '0';
                    }
                    dasFeature.orientation = ori;
                }
                dasFeature.score = elementValue(feature, "SCORE");
                dasFeature.links = dasLinksOf(feature);
                dasFeature.notes = dasNotesOf(feature);
                
                var groups = feature.getElementsByTagName("GROUP");
                for (var gi  = 0; gi < groups.length; ++gi) {
                    var groupXML = groups[gi];
                    var dasGroup = new DASGroup();
                    dasGroup.type = groupXML.getAttribute('type');
                    dasGroup.id = groupXML.getAttribute('id');
                    dasGroup.links = dasLinksOf(groupXML);
                    dasGroup.notes = dasNotesOf(groupXML);
                    if (!dasFeature.groups) {
                        dasFeature.groups = new Array(dasGroup);
                    } else {
                        dasFeature.groups.push(dasGroup);
                    }
                }

                // Magic notes.  Check with TAD before changing this.
                if (dasFeature.notes) {
                    for (var ni = 0; ni < dasFeature.notes.length; ++ni) {
                        var n = dasFeature.notes[ni];
                        if (n.indexOf('Genename=') == 0) {
                            var gg = new DASGroup();
                            gg.type='gene';
                            gg.id = n.substring(9);
                            if (!dasFeature.groups) {
                                dasFeature.groups = new Array(gg);
                            } else {
                                dasFeature.groups.push(gg);
                            }
                        }
                    }
                }
                
                {
                    var pec = feature.getElementsByTagName('PART');
                    if (pec.length > 0) {
                        var parts = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parts.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parts = parts;
                    }
                }
                {
                    var pec = feature.getElementsByTagName('PARENT');
                    if (pec.length > 0) {
                        var parents = [];
                        for (var pi = 0; pi < pec.length; ++pi) {
                            parents.push(pec[pi].getAttribute('id'));
                        }
                        dasFeature.parents = parents;
                    }
                }
                
                features.push(dasFeature);
            }
        }
                
        callback(features, undefined, segmentMap);
    },
    function (err) {
        callback([], err);
    });
}

function DASAlignment(type) {
    this.type = type;
    this.objects = {};
    this.blocks = [];
}

DASSource.prototype.alignments = function(segment, options, callback) {
    var dasURI = this.dasBaseURI + 'alignment?query=' + segment;
    this.doCrossDomainRequest(dasURI, function(responseXML) {
        if (!responseXML) {
            callback([], 'Failed request ' + dasURI);
            return;
        }

        var alignments = [];
        var aliXMLs = responseXML.getElementsByTagName('alignment');
        for (var ai = 0; ai < aliXMLs.length; ++ai) {
            var aliXML = aliXMLs[ai];
            var ali = new DASAlignment(aliXML.getAttribute('alignType'));
            var objXMLs = aliXML.getElementsByTagName('alignObject');
            for (var oi = 0; oi < objXMLs.length; ++oi) {
                var objXML = objXMLs[oi];
                var obj = {
                    id:          objXML.getAttribute('intObjectId'),
                    accession:   objXML.getAttribute('dbAccessionId'),
                    version:     objXML.getAttribute('objectVersion'),
                    dbSource:    objXML.getAttribute('dbSource'),
                    dbVersion:   objXML.getAttribute('dbVersion')
                };
                ali.objects[obj.id] = obj;
            }
            
            var blockXMLs = aliXML.getElementsByTagName('block');
            for (var bi = 0; bi < blockXMLs.length; ++bi) {
                var blockXML = blockXMLs[bi];
                var block = {
                    order:      blockXML.getAttribute('blockOrder'),
                    segments:   []
                };
                var segXMLs = blockXML.getElementsByTagName('segment');
                for (var si = 0; si < segXMLs.length; ++si) {
                    var segXML = segXMLs[si];
                    var seg = {
                        object:      segXML.getAttribute('intObjectId'),
                        min:         segXML.getAttribute('start'),
                        max:         segXML.getAttribute('end'),
                        strand:      segXML.getAttribute('strand'),
                        cigar:       elementValue(segXML, 'cigar')
                    };
                    block.segments.push(seg);
                }
                ali.blocks.push(block);
            }       
                    
            alignments.push(ali);
        }
        callback(alignments);
    });
}


function DASStylesheet() {
    this.styles = [];
}

DASStylesheet.prototype.pushStyle = function(filters, zoom, style) {
    if (!filters) {
        filters = {type: 'default'};
    }
    var styleHolder = shallowCopy(filters);
    if (zoom) {
        styleHolder.zoom = zoom;
    }
    styleHolder.style = style;
    this.styles.push(styleHolder);
}

function DASStyle() {
}

function parseGradient(grad) {
    var steps = grad.getAttribute('steps');
    if (steps) {
        steps = steps|0;
    } else {
        steps = 50;
    }

    var stops = [];
    var colors = [];
    var se = grad.getElementsByTagName('STOP');
    for (var si = 0; si < se.length; ++si) {
        var stop = se[si];
        stops.push(1.0 * stop.getAttribute('score'));
        colors.push(stop.firstChild.nodeValue);
    }

    return makeColourSteps(steps, stops, colors);
}

DASSource.prototype.stylesheet = function(successCB, failureCB) {
    var dasURI, creds = this.credentials;
    if (this.stylesheet_uri) {
        dasURI = this.stylesheet_uri;
        creds = false;
    } else {
        dasURI = this.dasBaseURI + 'stylesheet';
    }

    this.getURI(dasURI).then(function(dasURI) {
        doCrossDomainRequest(dasURI, function(responseXML) {
            if (!responseXML) {
                if (failureCB) {
                    failureCB();
                } 
                return;
            }
            var stylesheet = new DASStylesheet();
            var typeXMLs = responseXML.getElementsByTagName('TYPE');
            for (var i = 0; i < typeXMLs.length; ++i) {
                var typeStyle = typeXMLs[i];
            
                var filter = {};
                filter.type = typeStyle.getAttribute('id'); // Am I right in thinking that this makes DASSTYLE XML invalid?  Ugh.
                filter.label = typeStyle.getAttribute('label');
                filter.method = typeStyle.getAttribute('method');
                var glyphXMLs = typeStyle.getElementsByTagName('GLYPH');
                for (var gi = 0; gi < glyphXMLs.length; ++gi) {
                    var glyphXML = glyphXMLs[gi];
                    var zoom = glyphXML.getAttribute('zoom');
                    var glyph = childElementOf(glyphXML);
                    var style = new DASStyle();
                    style.glyph = glyph.localName;
                    var child = glyph.firstChild;
                    
                    while (child) {
                        if (child.nodeType == Node.ELEMENT_NODE) {
                            if (child.localName == 'BGGRAD') {
                                style[child.localName] = parseGradient(child);
                            } else {      
                                style[child.localName] = child.firstChild.nodeValue;
                            }
                        }
                        child = child.nextSibling;
                    }
                    stylesheet.pushStyle(filter, zoom, style);
                }
            }
            successCB(stylesheet);
        }, creds);
    }).catch(function(err) {
        console.log(err);
        failureCB();
    });
}

//
// sources command
// 

function DASRegistry(uri, opts)
{
    opts = opts || {};
    this.uri = uri;
    this.opts = opts;   
}

DASRegistry.prototype.sources = function(callback, failure, opts)
{
    if (!opts) {
        opts = {};
    }

    var filters = [];
    if (opts.taxon) {
        filters.push('organism=' + opts.taxon);
    }
    if (opts.auth) {
        filters.push('authority=' + opts.auth);
    }
    if (opts.version) {
        filters.push('version=' + opts.version);
    }
    var quri = this.uri;
    if (filters.length > 0) {
        quri = quri + '?' + filters.join('&');   // '&' as a separator to hack around dasregistry.org bug.
    }

    doCrossDomainRequest(quri, function(responseXML) {
        if (!responseXML && failure) {
            failure();
            return;
        }

        var sources = [];       
        var sourceXMLs = responseXML.getElementsByTagName('SOURCE');
        for (var si = 0; si < sourceXMLs.length; ++si) {
            var sourceXML = sourceXMLs[si];
            var versionXMLs = sourceXML.getElementsByTagName('VERSION');
            if (versionXMLs.length < 1) {
                continue;
            }
            var versionXML = versionXMLs[0];

            var coordXMLs = versionXML.getElementsByTagName('COORDINATES');
            var coords = [];
            for (var ci = 0; ci < coordXMLs.length; ++ci) {
                var coordXML = coordXMLs[ci];
                var coord = new DASCoords();
                coord.auth = coordXML.getAttribute('authority');
                coord.taxon = coordXML.getAttribute('taxid');
                coord.version = coordXML.getAttribute('version');
                coords.push(coord);
            }
            
            var caps = [];
            var capXMLs = versionXML.getElementsByTagName('CAPABILITY');
            var uri;
            for (var ci = 0; ci < capXMLs.length; ++ci) {
                var capXML = capXMLs[ci];
                
                caps.push(capXML.getAttribute('type'));

                if (capXML.getAttribute('type') == 'das1:features') {
                    var fep = capXML.getAttribute('query_uri');
                    uri = fep.substring(0, fep.length - ('features'.length));
                }
            }

            var props = {};
            var propXMLs = versionXML.getElementsByTagName('PROP');
            for (var pi = 0; pi < propXMLs.length; ++pi) {
                pusho(props, propXMLs[pi].getAttribute('name'), propXMLs[pi].getAttribute('value'));
            }
            
            if (uri) {
                var source = new DASSource(uri, {
                    source_uri: sourceXML.getAttribute('uri'),
                    name:  sourceXML.getAttribute('title'),
                    desc:  sourceXML.getAttribute('description'),
                    coords: coords,
                    props: props,
                    capabilities: caps
                });
                sources.push(source);
            }
        }
        
        callback(sources);
    });
}


//
// Utility functions
//

function elementValue(element, tag)
{
    var children = element.getElementsByTagName(tag);
    if (children.length > 0 && children[0].firstChild) {
        var c = children[0];
        if (c.childNodes.length == 1) {
            return c.firstChild.nodeValue;
        } else {
            var s = '';
            for (var ni = 0; ni < c.childNodes.length; ++ni) {
                s += c.childNodes[ni].nodeValue;
            }
            return s;
        }

    } else {
        return null;
    }
}

function childElementOf(element)
{
    if (element.hasChildNodes()) {
        var child = element.firstChild;
        do {
            if (child.nodeType == Node.ELEMENT_NODE) {
                return child;
            } 
            child = child.nextSibling;
        } while (child != null);
    }
    return null;
}


function dasLinksOf(element)
{
    var links = new Array();
    var maybeLinkChilden = element.getElementsByTagName('LINK');
    for (var ci = 0; ci < maybeLinkChilden.length; ++ci) {
        var linkXML = maybeLinkChilden[ci];
        if (linkXML.parentNode == element) {
            links.push(new DASLink(linkXML.firstChild ? linkXML.firstChild.nodeValue : 'Unknown', linkXML.getAttribute('href')));
        }
    }
    
    return links;
}

function dasNotesOf(element)
{
    var notes = [];
    var maybeNotes = element.getElementsByTagName('NOTE');
    for (var ni = 0; ni < maybeNotes.length; ++ni) {
        if (maybeNotes[ni].firstChild) {
            notes.push(maybeNotes[ni].firstChild.nodeValue);
        }
    }
    return notes;
}

function doCrossDomainRequest(url, handler, credentials, custAuth) {
    // TODO: explicit error handlers?

    if (window.XDomainRequest) {
        var req = new XDomainRequest();
        req.onload = function() {
            var dom = new ActiveXObject("Microsoft.XMLDOM");
            dom.async = false;
            dom.loadXML(req.responseText);
            handler(dom);
        }
        req.open("get", url);
        req.send('');
    } else {
        try {
            var req = new XMLHttpRequest();
            var timeout = setTimeout(
                function() {
                    console.log('timing out '  + url);
                    req.abort();
                    handler(null, req);
                },
                5000
            );

            req.timeout = 5000;
            req.ontimeout = function() {
                console.log('timeout on ' + url);
            };

            req.onreadystatechange = function() {
                if (req.readyState == 4) {
                    clearTimeout(timeout);
                    if (req.status >= 200 || req.status == 0) {
                        handler(req.responseXML, req);
                    }
                }
            };
            req.open("get", url, true);
            if (credentials) {
                req.withCredentials = true;
            }
            if (custAuth) {
                req.setRequestHeader('X-DAS-Authorisation', custAuth);
            }
            req.overrideMimeType('text/xml');
            req.setRequestHeader('Accept', 'application/xml,*/*');
            req.send('');
        } catch (e) {
            handler(null, req, e);
        }
    }
}

DASSource.prototype.doCrossDomainRequest = function(url, handler, errHandler) {
    var custAuth;
    if (this.xUser) {
        custAuth = 'Basic ' + btoa(this.xUser + ':' + this.xPass);
    }

    try {
        return doCrossDomainRequest(url, handler, this.credentials, custAuth);
    } catch (err) {
        if (errHandler) {
            errHandler(err);
        } else {
            throw err;
        }
    }
}

function isDasBooleanTrue(s) {
    s = ('' + s).toLowerCase();
    return s==='yes' || s==='true';
}

function isDasBooleanNotFalse(s) {
    if (!s)
        return false;

    s = ('' + s).toLowerCase();
    return s!=='no' || s!=='false';
}

function copyStylesheet(ss) {
    var nss = shallowCopy(ss);
    nss.styles = [];
    for (var si = 0; si < ss.styles.length; ++si) {
        var sh = nss.styles[si] = shallowCopy(ss.styles[si]);
        sh._methodRE = sh._labelRE = sh._typeRE = undefined;
        sh.style = shallowCopy(sh.style);
        sh.style.id = undefined;
        sh.style._gradient = undefined;
    }
    return nss;
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        DASGroup: DASGroup,
        DASFeature: DASFeature,
        DASStylesheet: DASStylesheet,
        DASStyle: DASStyle,
        DASSource: DASSource,
        DASSegment: DASSegment,
        DASRegistry: DASRegistry,
        DASSequence: DASSequence,
        DASLink: DASLink,

        isDasBooleanTrue: isDasBooleanTrue,
        isDasBooleanNotFalse: isDasBooleanNotFalse,
        copyStylesheet: copyStylesheet,
        coordsMatch: coordsMatch
    };
}

},{"./color":4,"./utils":11}],6:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// encode.js: interface for ENCODE DCC services
//

"use strict";

if (typeof(require) !== 'undefined') {
    var Promise = require('es6-promise').Promise;
}

function lookupEncodeURI(uri, json) {
    if (uri.indexOf('?') < 0)
        uri = uri + '?soft=true';

    return new Promise(function(accept, reject) {
        var req = new XMLHttpRequest();
        req.onreadystatechange = function() {
            if (req.readyState == 4) {
                if (req.status >= 300) {
                    reject('Error code ' + req.status);
                } else {
                    var resp = JSON.parse(req.response);
                    accept(json ? resp : resp.location);
                }
            }
        };
    
        req.open('GET', uri, true);
        req.setRequestHeader('Accept', 'application/json');
        req.responseType = 'text';
        req.send('');
    });
}

function EncodeURLHolder(url) {
    this.rawurl = url;
}

EncodeURLHolder.prototype.getURLPromise = function() {
    if (this.urlPromise && this.urlPromiseValidity > Date.now()) {
        return this.urlPromise;
    } else {
        this.urlPromise = lookupEncodeURI(this.rawurl, true).then(function(resp) {
            return resp.location;
        });
        this.urlPromiseValidity = Date.now() + (12 * 3600 * 1000);
        return this.urlPromise;
    }
}

function EncodeFetchable(url, start, end, opts) {
    if (!opts) {
        if (typeof start === 'object') {
            opts = start;
            start = undefined;
        } else {
            opts = {};
        }
    }

    this.url = (typeof url === 'string' ? new EncodeURLHolder(url) : url);
    this.start = start || 0;
    if (end) {
        this.end = end;
    }
    this.opts = opts;
}



EncodeFetchable.prototype.slice = function(s, l) {
    if (s < 0) {
        throw 'Bad slice ' + s;
    }

    var ns = this.start, ne = this.end;
    if (ns && s) {
        ns = ns + s;
    } else {
        ns = s || ns;
    }
    if (l && ns) {
        ne = ns + l - 1;
    } else {
        ne = ne || l - 1;
    }
    return new EncodeFetchable(this.url, ns, ne, this.opts);
}

EncodeFetchable.prototype.fetchAsText = function(callback) {
    var self = this;
    var req = new XMLHttpRequest();
    var length;
    self.url.getURLPromise().then(function(url) {
        req.open('GET', url, true);

        if (self.end) {
            if (self.end - self.start > 100000000) {
                throw 'Monster fetch!';
            }
            req.setRequestHeader('Range', 'bytes=' + self.start + '-' + self.end);
            length = self.end - self.start + 1;
        }

        req.onreadystatechange = function() {
            if (req.readyState == 4) {
                if (req.status == 200 || req.status == 206) {
                    return callback(req.responseText);
                } else {
                    return callback(null);
                }
            }
        };
        if (self.opts.credentials) {
            req.withCredentials = true;
        }
        req.send('');
    }).catch(function(err) {
        console.log(err);
        return callback(null);
    });
}

EncodeFetchable.prototype.salted = function() {
    return this;
}

EncodeFetchable.prototype.fetch = function(callback, attempt, truncatedLength) {
    var self = this;

    attempt = attempt || 1;
    if (attempt > 3) {
        return callback(null);
    }

    self.url.getURLPromise().then(function (url) {
        var req = new XMLHttpRequest();
        var length;
        req.open('GET', url, true);
        req.overrideMimeType('text/plain; charset=x-user-defined');
        if (self.end) {
            if (self.end - self.start > 100000000) {
                throw 'Monster fetch!';
            }
            req.setRequestHeader('Range', 'bytes=' + self.start + '-' + self.end);
            length = self.end - self.start + 1;
        }
        req.responseType = 'arraybuffer';
        req.onreadystatechange = function() {
            if (req.readyState == 4) {
                if (req.status == 200 || req.status == 206) {
                    if (req.response) {
                        var bl = req.response.byteLength;
                        if (length && length != bl && (!truncatedLength || bl != truncatedLength)) {
                            return self.fetch(callback, attempt + 1, bl);
                        } else {
                            return callback(req.response);
                        }
                    } else if (req.mozResponseArrayBuffer) {
                        return callback(req.mozResponseArrayBuffer);
                    } else {
                        var r = req.responseText;
                        if (length && length != r.length && (!truncatedLength || r.length != truncatedLength)) {
                            return self.fetch(callback, attempt + 1, r.length);
                        } else {
                            return callback(bstringToBuffer(req.responseText));
                        }
                    }
                } else {
                    return self.fetch(callback, attempt + 1);
                }
            }
        };
        if (self.opts.credentials) {
            req.withCredentials = true;
        }
        req.send('');
    }).catch(function(err) {
        console.log(err);
    });
}

function bstringToBuffer(result) {
    if (!result) {
        return null;
    }

    var ba = new Uint8Array(result.length);
    for (var i = 0; i < ba.length; ++i) {
        ba[i] = result.charCodeAt(i);
    }
    return ba.buffer;
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        lookupEncodeURI: lookupEncodeURI,
        EncodeFetchable: EncodeFetchable
    };
}

},{"es6-promise":12}],7:[function(require,module,exports){
(function (global){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2014
//
// fetchworker.js
//

"use strict";

var bin = require('./bin');
var bam = require('./bam');
var bigwig = require('./bigwig');
var encode = require('./encode');
var utils = require('./utils');

var Promise = require('es6-promise').Promise;

var connections = {};
var resolveResolvers = {};

var idSeed = 0;

global.newID = function() {
    return 'cn' + (++idSeed);
}

postMessage({tag: 'init'});

self.onmessage = function(event) {
    var d = event.data;
    var command = event.data.command;
    var tag = event.data.tag;

    if (!command) {
        var rr = resolveResolvers[tag];
        if (rr) {
            if (d.err) {
                rr.reject(d.err);
            } else {
                rr.resolve(d.url);
            }
                
            delete resolveResolvers[tag];
        }
    } else if (command === 'connectBAM') {
        var id = newID();
        var resolver;
        if (d.resolver) {
            resolver = proxyResolver(d.resolver);
        }
        var bamF, baiF, indexChunks;
        if (d.blob) {
            bamF = new bin.BlobFetchable(d.blob);
            baiF = new bin.BlobFetchable(d.indexBlob);
        } else {
            bamF = new bin.URLFetchable(d.uri, {credentials: d.credentials, resolver: resolver});
            baiF = new bin.URLFetchable(d.indexUri, {credentials: d.credentials, resolver: resolver});
            indexChunks = d.indexChunks;
        }

        bam.makeBam(bamF, baiF, indexChunks, function(bamObj, err) {
            if (bamObj) {
                connections[id] = new BAMWorkerFetcher(bamObj);
                postMessage({tag: tag, result: id});
            } else {
                postMessage({tag: tag, error: err || "Couldn't fetch BAM"});
            }
        });
    } else if (command === 'connectBBI') {
        var id = newID();
        var resolver;
        if (d.resolver) {
            resolver = proxyResolver(d.resolver);
        }
        var bbi;
        if (d.blob) {
            bbi = new bin.BlobFetchable(d.blob);
        } else if (d.transport == 'encode') {
            bbi = new encode.EncodeFetchable(d.uri, {credentials: d.credentials});
        } else {
            bbi = new bin.URLFetchable(d.uri, {credentials: d.credentials, resolver: resolver});
        }

        bigwig.makeBwg(bbi, function(bwg, err) {
            if (bwg) {
                connections[id] = new BBIWorkerFetcher(bwg);
                postMessage({tag: tag, result: id});
            } else {
                postMessage({tag: tag, error: err || "Couldn't fetch BBI"});
            }
        }, d.uri);
    } else if (command === 'textxhr') {
        utils.textXHR(d.uri, function(resp, err) {
            if (resp) {
                postMessage({tag: tag, result: resp});
            } else {
                postMessage({tag: tag, err: err || "Couldn't fetch resource"});
            }
        });
    } else if (command === 'fetch') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.fetch(d.tag, d.chr, d.min, d.max, d.zoom, d.opts);
    } else if (command === 'leap') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.leap(d.tag, d.chr, d.pos, d.dir);
    } else if (command === 'quantLeap') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.quantLeap(d.tag, d.chr, d.pos, d.dir, d.threshold, d.under);
    } else if (command === 'meta') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.meta(d.tag);
    } else if (command === 'search') {
        var con = connections[event.data.connection];
        if (!con) {
            return postMessage({tag: tag, error: 'No such connection: ' + event.data.connection});
        }

        con.search(d.tag, d.query, d.index);
    } else if (command === 'date') {
        return postMessage({tag: tag, result: Date.now()|0});
    } else {
        postMessage({tag: tag, error: 'Bad command ' + command});
    }
}

function BAMWorkerFetcher(bam) {
    this.bam = bam;
}

BAMWorkerFetcher.prototype.fetch = function(tag, chr, min, max, zoom, opts) {
    opts = opts || {};
    this.bam.fetch(chr, min, max, function(records, err) {
        if (records) {
            postMessage({tag: tag, result: records, time: Date.now()|0});
        } else {
            postMessage({tag: tag, error: err});
        }
    }, opts);
}

function BBIWorkerFetcher(bbi) {
    this.bbi = bbi;
}

BBIWorkerFetcher.prototype.fetch = function(tag, chr, min, max, zoom) {
    if (typeof(zoom) !== 'number')
        zoom = -1;

    var data;
    if (zoom < 0) {
        data = this.bbi.getUnzoomedView();
    } else {
        data = this.bbi.getZoomedView(zoom);
    }

    data.readWigData(chr, min, max, function(features) {
        postMessage({tag: tag, result: features});
    });
}

BBIWorkerFetcher.prototype.meta = function(tag) {
    var scales = [1];
    for (var z = 0; z < this.bbi.zoomLevels.length; ++z) {
        scales.push(this.bbi.zoomLevels[z].reduction);
    }

    var thisB = this;
    var meta = {type: this.bbi.type,
                zoomLevels: scales,
                fieldCount: this.bbi.fieldCount,
                definedFieldCount: this.bbi.definedFieldCount,
                schema: this.bbi.schema};
    if (this.bbi.type === 'bigbed') {
        this.bbi.getExtraIndices(function(ei) {
            if (ei) {
                thisB.extraIndices = ei;
                meta.extraIndices = ei.map(function(i) {return i.field});
            }
            postMessage({tag: tag, result: meta});
        });
    } else {
        postMessage({tag: tag, result: meta});
    }
}

BBIWorkerFetcher.prototype.leap = function(tag, chr, pos, dir) {
    this.bbi.getUnzoomedView().getFirstAdjacent(chr, pos, dir, function(result, err) {
        postMessage({tag: tag, result: result, error: err});
    });
}

BBIWorkerFetcher.prototype.quantLeap = function(tag, chr, pos, dir, threshold, under) {
    this.bbi.thresholdSearch(chr, pos, dir, threshold, function(result, err) {
        postMessage({tag: tag, result: result, error: err});
    });
}

BBIWorkerFetcher.prototype.search = function(tag, query, index) {
    var is = this.extraIndices[0];
    is.lookup(query, function(result, err) {
        postMessage({tag: tag, result: result, error: err});
    });
}

function proxyResolver(tag) {
    return function(url) {
        var lid = newID();
        return new Promise(function (resolve, reject) {
            resolveResolvers[lid] = {resolve: resolve, reject: reject};
            postMessage({tag: lid,
                         cmd: 'resolve',
                         resolver: tag,
                         url: url});
        });
    }
}

}).call(this,typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./bam":1,"./bigwig":2,"./bin":3,"./encode":6,"./utils":11,"es6-promise":12}],8:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// lh3utils.js: common support for lh3's file formats
//

if (typeof(require) !== 'undefined') {
    var jszlib = require('jszlib');
    var jszlib_inflate_buffer = jszlib.inflateBuffer;
    var arrayCopy = jszlib.arrayCopy;
}

function Vob(b, o) {
    this.block = b;
    this.offset = o;
}

Vob.prototype.toString = function() {
    return '' + this.block + ':' + this.offset;
}

function readVob(ba, offset) {
    var block = ((ba[offset+6] & 0xff) * 0x100000000) + ((ba[offset+5] & 0xff) * 0x1000000) + ((ba[offset+4] & 0xff) * 0x10000) + ((ba[offset+3] & 0xff) * 0x100) + ((ba[offset+2] & 0xff));
    var bint = (ba[offset+1] << 8) | (ba[offset]);
    if (block == 0 && bint == 0) {
        return null;  // Should only happen in the linear index?
    } else {
        return new Vob(block, bint);
    }
}

function unbgzf(data, lim) {
    lim = Math.min(lim || 1, data.byteLength - 50);
    var oBlockList = [];
    var ptr = [0];
    var totalSize = 0;

    while (ptr[0] < lim) {
        var ba = new Uint8Array(data, ptr[0], 12); // FIXME is this enough for all credible BGZF block headers?
        var xlen = (ba[11] << 8) | (ba[10]);
        // dlog('xlen[' + (ptr[0]) +']=' + xlen);
        var unc = jszlib_inflate_buffer(data, 12 + xlen + ptr[0], Math.min(65536, data.byteLength - 12 - xlen - ptr[0]), ptr);
        ptr[0] += 8;
        totalSize += unc.byteLength;
        oBlockList.push(unc);
    }

    if (oBlockList.length == 1) {
        return oBlockList[0];
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = new Uint8Array(oBlockList[i]);
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}

function Chunk(minv, maxv) {
    this.minv = minv; this.maxv = maxv;
}


//
// Binning (transliterated from SAM1.3 spec)
//

/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
function reg2bin(beg, end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
var MAX_BIN = (((1<<18)-1)/7);
function reg2bins(beg, end) 
{
    var i = 0, k, list = [];
    --end;
    list.push(0);
    for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) list.push(k);
    for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) list.push(k);
    for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) list.push(k);
    for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) list.push(k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push(k);
    return list;
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        unbgzf: unbgzf,
        readVob: readVob,
        reg2bin: reg2bin,
        reg2bins: reg2bins,
        Chunk: Chunk
    };
}
},{"jszlib":24}],9:[function(require,module,exports){
/*
 * A JavaScript implementation of the Secure Hash Algorithm, SHA-1, as defined
 * in FIPS 180-1
 * Version 2.2 Copyright Paul Johnston 2000 - 2009.
 * Other contributors: Greg Holt, Andrew Kepert, Ydnar, Lostinet
 * Distributed under the BSD License
 * See http://pajhome.org.uk/crypt/md5 for details.
 */

 "use strict";

/*
 * Configurable variables. You may need to tweak these to be compatible with
 * the server-side, but the defaults work in most cases.
 */
var hexcase = 0;  /* hex output format. 0 - lowercase; 1 - uppercase        */
var b64pad  = ""; /* base-64 pad character. "=" for strict RFC compliance   */

/*
 * These are the functions you'll usually want to call
 * They take string arguments and return either hex or base-64 encoded strings
 */
function hex_sha1(s)    { return rstr2hex(rstr_sha1(str2rstr_utf8(s))); }
function b64_sha1(s)    { return rstr2b64(rstr_sha1(str2rstr_utf8(s))); }
function any_sha1(s, e) { return rstr2any(rstr_sha1(str2rstr_utf8(s)), e); }
function hex_hmac_sha1(k, d)
  { return rstr2hex(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d))); }
function b64_hmac_sha1(k, d)
  { return rstr2b64(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d))); }
function any_hmac_sha1(k, d, e)
  { return rstr2any(rstr_hmac_sha1(str2rstr_utf8(k), str2rstr_utf8(d)), e); }

/*
 * Perform a simple self-test to see if the VM is working
 */
function sha1_vm_test()
{
  return hex_sha1("abc").toLowerCase() == "a9993e364706816aba3e25717850c26c9cd0d89d";
}

/*
 * Calculate the SHA1 of a raw string
 */
function rstr_sha1(s)
{
  return binb2rstr(binb_sha1(rstr2binb(s), s.length * 8));
}

/*
 * Calculate the HMAC-SHA1 of a key and some data (raw strings)
 */
function rstr_hmac_sha1(key, data)
{
  var bkey = rstr2binb(key);
  if(bkey.length > 16) bkey = binb_sha1(bkey, key.length * 8);

  var ipad = Array(16), opad = Array(16);
  for(var i = 0; i < 16; i++)
  {
    ipad[i] = bkey[i] ^ 0x36363636;
    opad[i] = bkey[i] ^ 0x5C5C5C5C;
  }

  var hash = binb_sha1(ipad.concat(rstr2binb(data)), 512 + data.length * 8);
  return binb2rstr(binb_sha1(opad.concat(hash), 512 + 160));
}

/*
 * Convert a raw string to a hex string
 */
function rstr2hex(input)
{
  // try { hexcase } catch(e) { hexcase=0; }
  var hex_tab = hexcase ? "0123456789ABCDEF" : "0123456789abcdef";
  var output = "";
  var x;
  for(var i = 0; i < input.length; i++)
  {
    x = input.charCodeAt(i);
    output += hex_tab.charAt((x >>> 4) & 0x0F)
           +  hex_tab.charAt( x        & 0x0F);
  }
  return output;
}

/*
 * Convert a raw string to a base-64 string
 */
function rstr2b64(input)
{
  // try { b64pad } catch(e) { b64pad=''; }
  var tab = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  var output = "";
  var len = input.length;
  for(var i = 0; i < len; i += 3)
  {
    var triplet = (input.charCodeAt(i) << 16)
                | (i + 1 < len ? input.charCodeAt(i+1) << 8 : 0)
                | (i + 2 < len ? input.charCodeAt(i+2)      : 0);
    for(var j = 0; j < 4; j++)
    {
      if(i * 8 + j * 6 > input.length * 8) output += b64pad;
      else output += tab.charAt((triplet >>> 6*(3-j)) & 0x3F);
    }
  }
  return output;
}

/*
 * Convert a raw string to an arbitrary string encoding
 */
function rstr2any(input, encoding)
{
  var divisor = encoding.length;
  var remainders = Array();
  var i, q, x, quotient;

  /* Convert to an array of 16-bit big-endian values, forming the dividend */
  var dividend = Array(Math.ceil(input.length / 2));
  for(i = 0; i < dividend.length; i++)
  {
    dividend[i] = (input.charCodeAt(i * 2) << 8) | input.charCodeAt(i * 2 + 1);
  }

  /*
   * Repeatedly perform a long division. The binary array forms the dividend,
   * the length of the encoding is the divisor. Once computed, the quotient
   * forms the dividend for the next step. We stop when the dividend is zero.
   * All remainders are stored for later use.
   */
  while(dividend.length > 0)
  {
    quotient = Array();
    x = 0;
    for(i = 0; i < dividend.length; i++)
    {
      x = (x << 16) + dividend[i];
      q = Math.floor(x / divisor);
      x -= q * divisor;
      if(quotient.length > 0 || q > 0)
        quotient[quotient.length] = q;
    }
    remainders[remainders.length] = x;
    dividend = quotient;
  }

  /* Convert the remainders to the output string */
  var output = "";
  for(i = remainders.length - 1; i >= 0; i--)
    output += encoding.charAt(remainders[i]);

  /* Append leading zero equivalents */
  var full_length = Math.ceil(input.length * 8 /
                                    (Math.log(encoding.length) / Math.log(2)))
  for(i = output.length; i < full_length; i++)
    output = encoding[0] + output;

  return output;
}

/*
 * Encode a string as utf-8.
 * For efficiency, this assumes the input is valid utf-16.
 */
function str2rstr_utf8(input)
{
  var output = "";
  var i = -1;
  var x, y;

  while(++i < input.length)
  {
    /* Decode utf-16 surrogate pairs */
    x = input.charCodeAt(i);
    y = i + 1 < input.length ? input.charCodeAt(i + 1) : 0;
    if(0xD800 <= x && x <= 0xDBFF && 0xDC00 <= y && y <= 0xDFFF)
    {
      x = 0x10000 + ((x & 0x03FF) << 10) + (y & 0x03FF);
      i++;
    }

    /* Encode output as utf-8 */
    if(x <= 0x7F)
      output += String.fromCharCode(x);
    else if(x <= 0x7FF)
      output += String.fromCharCode(0xC0 | ((x >>> 6 ) & 0x1F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0xFFFF)
      output += String.fromCharCode(0xE0 | ((x >>> 12) & 0x0F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
    else if(x <= 0x1FFFFF)
      output += String.fromCharCode(0xF0 | ((x >>> 18) & 0x07),
                                    0x80 | ((x >>> 12) & 0x3F),
                                    0x80 | ((x >>> 6 ) & 0x3F),
                                    0x80 | ( x         & 0x3F));
  }
  return output;
}

/*
 * Encode a string as utf-16
 */
function str2rstr_utf16le(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode( input.charCodeAt(i)        & 0xFF,
                                  (input.charCodeAt(i) >>> 8) & 0xFF);
  return output;
}

function str2rstr_utf16be(input)
{
  var output = "";
  for(var i = 0; i < input.length; i++)
    output += String.fromCharCode((input.charCodeAt(i) >>> 8) & 0xFF,
                                   input.charCodeAt(i)        & 0xFF);
  return output;
}

/*
 * Convert a raw string to an array of big-endian words
 * Characters >255 have their high-byte silently ignored.
 */
function rstr2binb(input)
{
  var output = Array(input.length >> 2);
  for(var i = 0; i < output.length; i++)
    output[i] = 0;
  for(var i = 0; i < input.length * 8; i += 8)
    output[i>>5] |= (input.charCodeAt(i / 8) & 0xFF) << (24 - i % 32);
  return output;
}

/*
 * Convert an array of big-endian words to a string
 */
function binb2rstr(input)
{
  var output = "";
  for(var i = 0; i < input.length * 32; i += 8)
    output += String.fromCharCode((input[i>>5] >>> (24 - i % 32)) & 0xFF);
  return output;
}

/*
 * Calculate the SHA-1 of an array of big-endian words, and a bit length
 */
function binb_sha1(x, len)
{
  /* append padding */
  x[len >> 5] |= 0x80 << (24 - len % 32);
  x[((len + 64 >> 9) << 4) + 15] = len;

  var w = Array(80);
  var a =  1732584193;
  var b = -271733879;
  var c = -1732584194;
  var d =  271733878;
  var e = -1009589776;

  for(var i = 0; i < x.length; i += 16)
  {
    var olda = a;
    var oldb = b;
    var oldc = c;
    var oldd = d;
    var olde = e;

    for(var j = 0; j < 80; j++)
    {
      if(j < 16) w[j] = x[i + j];
      else w[j] = bit_rol(w[j-3] ^ w[j-8] ^ w[j-14] ^ w[j-16], 1);
      var t = safe_add(safe_add(bit_rol(a, 5), sha1_ft(j, b, c, d)),
                       safe_add(safe_add(e, w[j]), sha1_kt(j)));
      e = d;
      d = c;
      c = bit_rol(b, 30);
      b = a;
      a = t;
    }

    a = safe_add(a, olda);
    b = safe_add(b, oldb);
    c = safe_add(c, oldc);
    d = safe_add(d, oldd);
    e = safe_add(e, olde);
  }
  return Array(a, b, c, d, e);

}

/*
 * Perform the appropriate triplet combination function for the current
 * iteration
 */
function sha1_ft(t, b, c, d)
{
  if(t < 20) return (b & c) | ((~b) & d);
  if(t < 40) return b ^ c ^ d;
  if(t < 60) return (b & c) | (b & d) | (c & d);
  return b ^ c ^ d;
}

/*
 * Determine the appropriate additive constant for the current iteration
 */
function sha1_kt(t)
{
  return (t < 20) ?  1518500249 : (t < 40) ?  1859775393 :
         (t < 60) ? -1894007588 : -899497514;
}

/*
 * Add integers, wrapping at 2^32. This uses 16-bit operations internally
 * to work around bugs in some JS interpreters.
 */
function safe_add(x, y)
{
  var lsw = (x & 0xFFFF) + (y & 0xFFFF);
  var msw = (x >> 16) + (y >> 16) + (lsw >> 16);
  return (msw << 16) | (lsw & 0xFFFF);
}

/*
 * Bitwise rotate a 32-bit number to the left.
 */
function bit_rol(num, cnt)
{
  return (num << cnt) | (num >>> (32 - cnt));
}

if (typeof(module) !== 'undefined') {
  module.exports = {
    b64_sha1: b64_sha1,
    hex_sha1: hex_sha1
  }
}

},{}],10:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// spans.js: JavaScript Intset/Location port.
//

"use strict";


function Range(min, max)
{
    if (typeof(min) != 'number' || typeof(max) != 'number')
        throw 'Bad range ' + min + ',' + max;
    this._min = min;
    this._max = max;
}

Range.prototype.min = function() {
    return this._min;
}

Range.prototype.max = function() {
    return this._max;
}

Range.prototype.contains = function(pos) {
    return pos >= this._min && pos <= this._max;
}

Range.prototype.isContiguous = function() {
    return true;
}

Range.prototype.ranges = function() {
    return [this];
}

Range.prototype._pushRanges = function(ranges) {
    ranges.push(this);
}

Range.prototype.toString = function() {
    return '[' + this._min + '-' + this._max + ']';
}

function _Compound(ranges) {
    // given: a set of unsorted possibly overlapping ranges
    // sort the input ranges
    var sorted = ranges.sort(_rangeOrder);
    // merge overlaps between adjacent ranges
    var merged = [];
    var current = sorted.shift();
    sorted.forEach(function(range) {
        if (range._min <= current._max) {
            if (range._max > current._max) {
                current._max = range._max;
            }
        }
        else {
            merged.push(current);
            current = range;
        }
    });
    merged.push(current);
    this._ranges = merged;
}

_Compound.prototype.min = function() {
    return this._ranges[0].min();
}

_Compound.prototype.max = function() {
    return this._ranges[this._ranges.length - 1].max();
}

// returns the index of the first range that is not less than pos
_Compound.prototype.lower_bound = function(pos) {
    // first check if pos is out of range
    var r = this.ranges();
    if (pos > this.max()) return r.length;
    if (pos < this.min()) return 0;
    // do a binary search
    var a=0, b=r.length - 1;
    while (a <= b) {
        var m = Math.floor((a+b)/2);
        if (pos > r[m]._max) {
            a = m+1;
        }
        else if (pos < r[m]._min) {
            b = m-1;
        }
        else {
            return m;
        }
    }
    return a;
}

_Compound.prototype.contains = function(pos) {
    var lb = this.lower_bound(pos);
    if (lb < this._ranges.length && this._ranges[lb].contains(pos)) {
        return true;
    }
    return false;
}

_Compound.prototype.insertRange = function(range) {
    var lb = this.lower_bound(range._min);
    if (lb === this._ranges.length) { // range follows this
        this._ranges.push(range);
        return;
    }
    
    var r = this.ranges();
    if (range._max < r[lb]._min) { // range preceeds lb
        this._ranges.splice(lb,0,range);
        return;
    }

    // range overlaps lb (at least)
    if (r[lb]._min < range._min) range._min = r[lb]._min;
    var ub = lb+1;
    while (ub < r.length && r[ub]._min <= range._max) {
        ub++;
    }
    ub--;
    // ub is the upper bound of the new range
    if (r[ub]._max > range._max) range._max = r[ub]._max;
    
    // splice range into this._ranges
    this._ranges.splice(lb,ub-lb+1,range);
    return;
}

_Compound.prototype.isContiguous = function() {
    return this._ranges.length > 1;
}

_Compound.prototype.ranges = function() {
    return this._ranges;
}

_Compound.prototype._pushRanges = function(ranges) {
    for (var ri = 0; ri < this._ranges.length; ++ri)
        ranges.push(this._ranges[ri]);
}

_Compound.prototype.toString = function() {
    var s = '';
    for (var r = 0; r < this._ranges.length; ++r) {
        if (r>0) {
            s = s + ',';
        }
        s = s + this._ranges[r].toString();
    }
    return s;
}

function union(s0, s1) {
    if (! (s0 instanceof _Compound)) {
        if (! (s0 instanceof Array))
            s0 = [s0];
        s0 = new _Compound(s0);
    }
    
    if (s1)
        s0.insertRange(s1);

    return s0;
}

function intersection(s0, s1) {
    var r0 = s0.ranges();
    var r1 = s1.ranges();
    var l0 = r0.length, l1 = r1.length;
    var i0 = 0, i1 = 0;
    var or = [];

    while (i0 < l0 && i1 < l1) {
        var s0 = r0[i0], s1 = r1[i1];
        var lapMin = Math.max(s0.min(), s1.min());
        var lapMax = Math.min(s0.max(), s1.max());
        if (lapMax >= lapMin) {
            or.push(new Range(lapMin, lapMax));
        }
        if (s0.max() > s1.max()) {
            ++i1;
        } else {
            ++i0;
        }
    }
    
    if (or.length == 0) {
        return null; // FIXME
    } else if (or.length == 1) {
        return or[0];
    } else {
        return new _Compound(or);
    }
}

function coverage(s) {
    var tot = 0;
    var rl = s.ranges();
    for (var ri = 0; ri < rl.length; ++ri) {
        var r = rl[ri];
        tot += (r.max() - r.min() + 1);
    }
    return tot;
}



function rangeOrder(a, b)
{
    if (a.min() < b.min()) {
        return -1;
    } else if (a.min() > b.min()) {
        return 1;
    } else if (a.max() < b.max()) {
        return -1;
    } else if (b.max() > a.max()) {
        return 1;
    } else {
        return 0;
    }
}

function _rangeOrder(a, b)
{
    if (a._min < b._min) {
        return -1;
    } else if (a._min > b._min) {
        return 1;
    } else if (a._max < b._max) {
        return -1;
    } else if (b._max > a._max) {
        return 1;
    } else {
        return 0;
    }
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        Range: Range,
        union: union,
        intersection: intersection,
        coverage: coverage,
        rangeOver: rangeOrder,
        _rangeOrder: _rangeOrder
    }
}
},{}],11:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2010
//
// utils.js: odds, sods, and ends.
//

"use strict";

if (typeof(require) !== 'undefined') {
    var sha1 = require('./sha1');
    var b64_sha1 = sha1.b64_sha1;
}

var NUM_REGEXP = new RegExp('[0-9]+');

function stringToNumbersArray(str) {
    var nums = new Array();
    var m;
    while (m = NUM_REGEXP.exec(str)) {
        nums.push(m[0]);
        str=str.substring(m.index + (m[0].length));
    }
    return nums;
}

var STRICT_NUM_REGEXP = new RegExp('^[0-9]+$');

function stringToInt(str) {
    str = str.replace(new RegExp(',', 'g'), '');
    if (!STRICT_NUM_REGEXP.test(str)) {
        return null;
    }
    return str|0;
}

function pushnew(a, v) {
    for (var i = 0; i < a.length; ++i) {
        if (a[i] == v) {
            return;
        }
    }
    a.push(v);
}

function pusho(obj, k, v) {
    if (obj[k]) {
        obj[k].push(v);
    } else {
        obj[k] = [v];
    }
}

function pushnewo(obj, k, v) {
    var a = obj[k];
    if (a) {
        for (var i = 0; i < a.length; ++i) {    // indexOf requires JS16 :-(.
            if (a[i] == v) {
                return;
            }
        }
        a.push(v);
    } else {
        obj[k] = [v];
    }
}


function pick(a, b, c, d)
{
    if (a) {
        return a;
    } else if (b) {
        return b;
    } else if (c) {
        return c;
    } else if (d) {
        return d;
    }
}

function pushnew(l, o)
{
    for (var i = 0; i < l.length; ++i) {
        if (l[i] == o) {
            return;
        }
    }
    l.push(o);
}



function arrayIndexOf(a, x) {
    if (!a) {
        return -1;
    }

    for (var i = 0; i < a.length; ++i) {
        if (a[i] === x) {
            return i;
        }
    }
    return -1;
}

function arrayRemove(a, x) {
    var i = arrayIndexOf(a, x);
    if (i >= 0) {
        a.splice(i, 1);
        return true;
    }
    return false;
}

//
// DOM utilities
//


function makeElement(tag, children, attribs, styles)
{
    var ele = document.createElement(tag);
    if (children) {
        if (! (children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (c) {
                if (typeof c == 'string') {
                    c = document.createTextNode(c);
                } else if (typeof c == 'number') {
                    c = document.createTextNode('' + c);
                }
                ele.appendChild(c);
            }
        }
    }
    
    if (attribs) {
        for (var l in attribs) {
            try {
                ele[l] = attribs[l];
            } catch (e) {
                console.log('error setting ' + l);
                throw(e);
            }
        }
    }
    if (styles) {
        for (var l in styles) {
            ele.style[l] = styles[l];
        }
    }
    return ele;
}

function makeElementNS(namespace, tag, children, attribs)
{
    var ele = document.createElementNS(namespace, tag);
    if (children) {
        if (! (children instanceof Array)) {
            children = [children];
        }
        for (var i = 0; i < children.length; ++i) {
            var c = children[i];
            if (typeof c == 'string') {
                c = document.createTextNode(c);
            }
            ele.appendChild(c);
        }
    }
    
    setAttrs(ele, attribs);
    return ele;
}

var attr_name_cache = {};

function setAttr(node, key, value)
{
    var attr = attr_name_cache[key];
    if (!attr) {
        var _attr = '';
        for (var c = 0; c < key.length; ++c) {
            var cc = key.substring(c, c+1);
            var lcc = cc.toLowerCase();
            if (lcc != cc) {
                _attr = _attr + '-' + lcc;
            } else {
                _attr = _attr + cc;
            }
        }
        attr_name_cache[key] = _attr;
        attr = _attr;
    }
    node.setAttribute(attr, value);
}

function setAttrs(node, attribs)
{
    if (attribs) {
        for (var l in attribs) {
            setAttr(node, l, attribs[l]);
        }
    }
}



function removeChildren(node)
{
    if (!node || !node.childNodes) {
        return;
    }

    while (node.childNodes.length > 0) {
        node.removeChild(node.firstChild);
    }
}



//
// WARNING: not for general use!
//

function miniJSONify(o, exc) {
    if (typeof o === 'undefined') {
        return 'undefined';
    } else if (o == null) {
        return 'null';
    } else if (typeof o == 'string') {
        return "'" + o + "'";
    } else if (typeof o == 'number') {
        return "" + o;
    } else if (typeof o == 'boolean') {
        return "" + o;
    } else if (typeof o == 'object') {
        if (o instanceof Array) {
            var s = null;
            for (var i = 0; i < o.length; ++i) {
                s = (s == null ? '' : (s + ', ')) + miniJSONify(o[i], exc);
            }
            return '[' + (s?s:'') + ']';
        } else {
            exc = exc || {};
            var s = null;
            for (var k in o) {
                if (exc[k])
                    continue;
                if (k != undefined && typeof(o[k]) != 'function') {
                    s = (s == null ? '' : (s + ', ')) + k + ': ' + miniJSONify(o[k], exc);
                }
            }
            return '{' + (s?s:'') + '}';
        }
    } else {
        return (typeof o);
    }
}

function shallowCopy(o) {
    var n = {};
    for (var k in o) {
        n[k] = o[k];
    }
    return n;
}

function Observed(x) {
    this.value = x;
    this.listeners = [];
}

Observed.prototype.addListener = function(f) {
    this.listeners.push(f);
}

Observed.prototype.addListenerAndFire = function(f) {
    this.listeners.push(f);
    f(this.value);
}

Observed.prototype.removeListener = function(f) {
    arrayRemove(this.listeners, f);
}

Observed.prototype.get = function() {
    return this.value;
}

Observed.prototype.set = function(x) {
    this.value = x;
    for (var i = 0; i < this.listeners.length; ++i) {
        this.listeners[i](x);
    }
}

function Awaited() {
    this.queue = [];
}

Awaited.prototype.provide = function(x) {
    if (this.res !== undefined) {
        throw "Resource has already been provided.";
    }

    this.res = x;
    for (var i = 0; i < this.queue.length; ++i) {
        this.queue[i](x);
    }
    this.queue = null;   // avoid leaking closures.
}

Awaited.prototype.await = function(f) {
    if (this.res !== undefined) {
        f(this.res);
        return this.res;
    } else {
        this.queue.push(f);
    }
}

var __dalliance_saltSeed = 0;

function saltURL(url) {
    return url + '?salt=' + b64_sha1('' + Date.now() + ',' + (++__dalliance_saltSeed));
}

function textXHR(url, callback, opts) {
    if (opts && opts.salt) 
        url = saltURL(url);

    try {
        var timeout;
        if (opts.timeout) {
            timeout = setTimeout(
                function() {
                    console.log('timing out ' + url);
                    req.abort();
                    return callback(null, 'Timeout');
                },
                opts.timeout
            );
        }

        var req = new XMLHttpRequest();
        req.onreadystatechange = function() {
    	    if (req.readyState == 4) {
                if (timeout)
                    clearTimeout(timeout);
    	        if (req.status < 200 || req.status >= 300) {
    		    callback(null, 'Error code ' + req.status);
    	        } else {
    		    callback(req.responseText);
    	        }
    	    }
        };
        
        req.open('GET', url, true);
        req.responseType = 'text';

        if (opts && opts.credentials) {
            req.withCredentials = true;
        }
        req.send('');
    } catch (e) {
        callback(null, 'Exception ' + e);
    }
}

function relativeURL(base, rel) {
    // FIXME quite naive -- good enough for trackhubs?

    if (rel.indexOf('http:') == 0 || rel.indexOf('https:') == 0) {
        return rel;
    }

    var li = base.lastIndexOf('/');
    if (li >= 0) {
        return base.substr(0, li + 1) + rel;
    } else {
        return rel;
    }
}

var AMINO_ACID_TRANSLATION = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',
    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',
    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',
    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',
    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'TAT': 'Y',
    'TAC': 'Y',
    'TAA': '*',  // stop
    'TAG': '*',  // stop
    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'TGT': 'C',
    'TGC': 'C',
    'TGA': '*',  // stop
    'TGG': 'W',
    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G'
}

function resolveUrlToPage(rel) {
    return makeElement('a', null, {href: rel}).href;
}

//
// Missing APIs
// 

if (!('trim' in String.prototype)) {
    String.prototype.trim = function() {
        return this.replace(/^\s+/, '').replace(/\s+$/, '');
    };
}

if (typeof(module) !== 'undefined') {
    module.exports = {
        textXHR: textXHR,
        relativeURL: relativeURL,
        resolveUrlToPage: resolveUrlToPage,
        shallowCopy: shallowCopy,
        pusho: pusho,
        pushnew: pushnew,
        pushnewo: pushnewo,
        arrayIndexOf: arrayIndexOf,
        pick: pick,

        makeElement: makeElement,
        makeElementNS: makeElementNS,
        removeChildren: removeChildren,

        miniJSONify: miniJSONify,

        Observed: Observed,
        Awaited: Awaited,

        AMINO_ACID_TRANSLATION: AMINO_ACID_TRANSLATION
    }
}

},{"./sha1":9}],12:[function(require,module,exports){
"use strict";
var Promise = require("./promise/promise").Promise;
var polyfill = require("./promise/polyfill").polyfill;
exports.Promise = Promise;
exports.polyfill = polyfill;
},{"./promise/polyfill":17,"./promise/promise":18}],13:[function(require,module,exports){
"use strict";
/* global toString */

var isArray = require("./utils").isArray;
var isFunction = require("./utils").isFunction;

/**
  Returns a promise that is fulfilled when all the given promises have been
  fulfilled, or rejected if any of them become rejected. The return promise
  is fulfilled with an array that gives all the values in the order they were
  passed in the `promises` array argument.

  Example:

  ```javascript
  var promise1 = RSVP.resolve(1);
  var promise2 = RSVP.resolve(2);
  var promise3 = RSVP.resolve(3);
  var promises = [ promise1, promise2, promise3 ];

  RSVP.all(promises).then(function(array){
    // The array here would be [ 1, 2, 3 ];
  });
  ```

  If any of the `promises` given to `RSVP.all` are rejected, the first promise
  that is rejected will be given as an argument to the returned promises's
  rejection handler. For example:

  Example:

  ```javascript
  var promise1 = RSVP.resolve(1);
  var promise2 = RSVP.reject(new Error("2"));
  var promise3 = RSVP.reject(new Error("3"));
  var promises = [ promise1, promise2, promise3 ];

  RSVP.all(promises).then(function(array){
    // Code here never runs because there are rejected promises!
  }, function(error) {
    // error.message === "2"
  });
  ```

  @method all
  @for RSVP
  @param {Array} promises
  @param {String} label
  @return {Promise} promise that is fulfilled when all `promises` have been
  fulfilled, or rejected if any of them become rejected.
*/
function all(promises) {
  /*jshint validthis:true */
  var Promise = this;

  if (!isArray(promises)) {
    throw new TypeError('You must pass an array to all.');
  }

  return new Promise(function(resolve, reject) {
    var results = [], remaining = promises.length,
    promise;

    if (remaining === 0) {
      resolve([]);
    }

    function resolver(index) {
      return function(value) {
        resolveAll(index, value);
      };
    }

    function resolveAll(index, value) {
      results[index] = value;
      if (--remaining === 0) {
        resolve(results);
      }
    }

    for (var i = 0; i < promises.length; i++) {
      promise = promises[i];

      if (promise && isFunction(promise.then)) {
        promise.then(resolver(i), reject);
      } else {
        resolveAll(i, promise);
      }
    }
  });
}

exports.all = all;
},{"./utils":22}],14:[function(require,module,exports){
(function (process,global){
"use strict";
var browserGlobal = (typeof window !== 'undefined') ? window : {};
var BrowserMutationObserver = browserGlobal.MutationObserver || browserGlobal.WebKitMutationObserver;
var local = (typeof global !== 'undefined') ? global : (this === undefined? window:this);

// node
function useNextTick() {
  return function() {
    process.nextTick(flush);
  };
}

function useMutationObserver() {
  var iterations = 0;
  var observer = new BrowserMutationObserver(flush);
  var node = document.createTextNode('');
  observer.observe(node, { characterData: true });

  return function() {
    node.data = (iterations = ++iterations % 2);
  };
}

function useSetTimeout() {
  return function() {
    local.setTimeout(flush, 1);
  };
}

var queue = [];
function flush() {
  for (var i = 0; i < queue.length; i++) {
    var tuple = queue[i];
    var callback = tuple[0], arg = tuple[1];
    callback(arg);
  }
  queue = [];
}

var scheduleFlush;

// Decide what async method to use to triggering processing of queued callbacks:
if (typeof process !== 'undefined' && {}.toString.call(process) === '[object process]') {
  scheduleFlush = useNextTick();
} else if (BrowserMutationObserver) {
  scheduleFlush = useMutationObserver();
} else {
  scheduleFlush = useSetTimeout();
}

function asap(callback, arg) {
  var length = queue.push([callback, arg]);
  if (length === 1) {
    // If length is 1, that means that we need to schedule an async flush.
    // If additional callbacks are queued before the queue is flushed, they
    // will be processed by this flush that we are scheduling.
    scheduleFlush();
  }
}

exports.asap = asap;
}).call(this,require("1YiZ5S"),typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"1YiZ5S":23}],15:[function(require,module,exports){
"use strict";
/**
  `RSVP.Promise.cast` returns the same promise if that promise shares a constructor
  with the promise being casted.

  Example:

  ```javascript
  var promise = RSVP.resolve(1);
  var casted = RSVP.Promise.cast(promise);

  console.log(promise === casted); // true
  ```

  In the case of a promise whose constructor does not match, it is assimilated.
  The resulting promise will fulfill or reject based on the outcome of the
  promise being casted.

  In the case of a non-promise, a promise which will fulfill with that value is
  returned.

  Example:

  ```javascript
  var value = 1; // could be a number, boolean, string, undefined...
  var casted = RSVP.Promise.cast(value);

  console.log(value === casted); // false
  console.log(casted instanceof RSVP.Promise) // true

  casted.then(function(val) {
    val === value // => true
  });
  ```

  `RSVP.Promise.cast` is similar to `RSVP.resolve`, but `RSVP.Promise.cast` differs in the
  following ways:
  * `RSVP.Promise.cast` serves as a memory-efficient way of getting a promise, when you
  have something that could either be a promise or a value. RSVP.resolve
  will have the same effect but will create a new promise wrapper if the
  argument is a promise.
  * `RSVP.Promise.cast` is a way of casting incoming thenables or promise subclasses to
  promises of the exact class specified, so that the resulting object's `then` is
  ensured to have the behavior of the constructor you are calling cast on (i.e., RSVP.Promise).

  @method cast
  @for RSVP
  @param {Object} object to be casted
  @return {Promise} promise that is fulfilled when all properties of `promises`
  have been fulfilled, or rejected if any of them become rejected.
*/


function cast(object) {
  /*jshint validthis:true */
  if (object && typeof object === 'object' && object.constructor === this) {
    return object;
  }

  var Promise = this;

  return new Promise(function(resolve) {
    resolve(object);
  });
}

exports.cast = cast;
},{}],16:[function(require,module,exports){
"use strict";
var config = {
  instrument: false
};

function configure(name, value) {
  if (arguments.length === 2) {
    config[name] = value;
  } else {
    return config[name];
  }
}

exports.config = config;
exports.configure = configure;
},{}],17:[function(require,module,exports){
(function (global){
"use strict";
/*global self*/
var RSVPPromise = require("./promise").Promise;
var isFunction = require("./utils").isFunction;

function polyfill() {
  var local;

  if (typeof global !== 'undefined') {
    local = global;
  } else if (typeof window !== 'undefined' && window.document) {
    local = window;
  } else {
    local = self;
  }

  var es6PromiseSupport = 
    "Promise" in local &&
    // Some of these methods are missing from
    // Firefox/Chrome experimental implementations
    "cast" in local.Promise &&
    "resolve" in local.Promise &&
    "reject" in local.Promise &&
    "all" in local.Promise &&
    "race" in local.Promise &&
    // Older version of the spec had a resolver object
    // as the arg rather than a function
    (function() {
      var resolve;
      new local.Promise(function(r) { resolve = r; });
      return isFunction(resolve);
    }());

  if (!es6PromiseSupport) {
    local.Promise = RSVPPromise;
  }
}

exports.polyfill = polyfill;
}).call(this,typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./promise":18,"./utils":22}],18:[function(require,module,exports){
"use strict";
var config = require("./config").config;
var configure = require("./config").configure;
var objectOrFunction = require("./utils").objectOrFunction;
var isFunction = require("./utils").isFunction;
var now = require("./utils").now;
var cast = require("./cast").cast;
var all = require("./all").all;
var race = require("./race").race;
var staticResolve = require("./resolve").resolve;
var staticReject = require("./reject").reject;
var asap = require("./asap").asap;

var counter = 0;

config.async = asap; // default async is asap;

function Promise(resolver) {
  if (!isFunction(resolver)) {
    throw new TypeError('You must pass a resolver function as the first argument to the promise constructor');
  }

  if (!(this instanceof Promise)) {
    throw new TypeError("Failed to construct 'Promise': Please use the 'new' operator, this object constructor cannot be called as a function.");
  }

  this._subscribers = [];

  invokeResolver(resolver, this);
}

function invokeResolver(resolver, promise) {
  function resolvePromise(value) {
    resolve(promise, value);
  }

  function rejectPromise(reason) {
    reject(promise, reason);
  }

  try {
    resolver(resolvePromise, rejectPromise);
  } catch(e) {
    rejectPromise(e);
  }
}

function invokeCallback(settled, promise, callback, detail) {
  var hasCallback = isFunction(callback),
      value, error, succeeded, failed;

  if (hasCallback) {
    try {
      value = callback(detail);
      succeeded = true;
    } catch(e) {
      failed = true;
      error = e;
    }
  } else {
    value = detail;
    succeeded = true;
  }

  if (handleThenable(promise, value)) {
    return;
  } else if (hasCallback && succeeded) {
    resolve(promise, value);
  } else if (failed) {
    reject(promise, error);
  } else if (settled === FULFILLED) {
    resolve(promise, value);
  } else if (settled === REJECTED) {
    reject(promise, value);
  }
}

var PENDING   = void 0;
var SEALED    = 0;
var FULFILLED = 1;
var REJECTED  = 2;

function subscribe(parent, child, onFulfillment, onRejection) {
  var subscribers = parent._subscribers;
  var length = subscribers.length;

  subscribers[length] = child;
  subscribers[length + FULFILLED] = onFulfillment;
  subscribers[length + REJECTED]  = onRejection;
}

function publish(promise, settled) {
  var child, callback, subscribers = promise._subscribers, detail = promise._detail;

  for (var i = 0; i < subscribers.length; i += 3) {
    child = subscribers[i];
    callback = subscribers[i + settled];

    invokeCallback(settled, child, callback, detail);
  }

  promise._subscribers = null;
}

Promise.prototype = {
  constructor: Promise,

  _state: undefined,
  _detail: undefined,
  _subscribers: undefined,

  then: function(onFulfillment, onRejection) {
    var promise = this;

    var thenPromise = new this.constructor(function() {});

    if (this._state) {
      var callbacks = arguments;
      config.async(function invokePromiseCallback() {
        invokeCallback(promise._state, thenPromise, callbacks[promise._state - 1], promise._detail);
      });
    } else {
      subscribe(this, thenPromise, onFulfillment, onRejection);
    }

    return thenPromise;
  },

  'catch': function(onRejection) {
    return this.then(null, onRejection);
  }
};

Promise.all = all;
Promise.cast = cast;
Promise.race = race;
Promise.resolve = staticResolve;
Promise.reject = staticReject;

function handleThenable(promise, value) {
  var then = null,
  resolved;

  try {
    if (promise === value) {
      throw new TypeError("A promises callback cannot return that same promise.");
    }

    if (objectOrFunction(value)) {
      then = value.then;

      if (isFunction(then)) {
        then.call(value, function(val) {
          if (resolved) { return true; }
          resolved = true;

          if (value !== val) {
            resolve(promise, val);
          } else {
            fulfill(promise, val);
          }
        }, function(val) {
          if (resolved) { return true; }
          resolved = true;

          reject(promise, val);
        });

        return true;
      }
    }
  } catch (error) {
    if (resolved) { return true; }
    reject(promise, error);
    return true;
  }

  return false;
}

function resolve(promise, value) {
  if (promise === value) {
    fulfill(promise, value);
  } else if (!handleThenable(promise, value)) {
    fulfill(promise, value);
  }
}

function fulfill(promise, value) {
  if (promise._state !== PENDING) { return; }
  promise._state = SEALED;
  promise._detail = value;

  config.async(publishFulfillment, promise);
}

function reject(promise, reason) {
  if (promise._state !== PENDING) { return; }
  promise._state = SEALED;
  promise._detail = reason;

  config.async(publishRejection, promise);
}

function publishFulfillment(promise) {
  publish(promise, promise._state = FULFILLED);
}

function publishRejection(promise) {
  publish(promise, promise._state = REJECTED);
}

exports.Promise = Promise;
},{"./all":13,"./asap":14,"./cast":15,"./config":16,"./race":19,"./reject":20,"./resolve":21,"./utils":22}],19:[function(require,module,exports){
"use strict";
/* global toString */
var isArray = require("./utils").isArray;

/**
  `RSVP.race` allows you to watch a series of promises and act as soon as the
  first promise given to the `promises` argument fulfills or rejects.

  Example:

  ```javascript
  var promise1 = new RSVP.Promise(function(resolve, reject){
    setTimeout(function(){
      resolve("promise 1");
    }, 200);
  });

  var promise2 = new RSVP.Promise(function(resolve, reject){
    setTimeout(function(){
      resolve("promise 2");
    }, 100);
  });

  RSVP.race([promise1, promise2]).then(function(result){
    // result === "promise 2" because it was resolved before promise1
    // was resolved.
  });
  ```

  `RSVP.race` is deterministic in that only the state of the first completed
  promise matters. For example, even if other promises given to the `promises`
  array argument are resolved, but the first completed promise has become
  rejected before the other promises became fulfilled, the returned promise
  will become rejected:

  ```javascript
  var promise1 = new RSVP.Promise(function(resolve, reject){
    setTimeout(function(){
      resolve("promise 1");
    }, 200);
  });

  var promise2 = new RSVP.Promise(function(resolve, reject){
    setTimeout(function(){
      reject(new Error("promise 2"));
    }, 100);
  });

  RSVP.race([promise1, promise2]).then(function(result){
    // Code here never runs because there are rejected promises!
  }, function(reason){
    // reason.message === "promise2" because promise 2 became rejected before
    // promise 1 became fulfilled
  });
  ```

  @method race
  @for RSVP
  @param {Array} promises array of promises to observe
  @param {String} label optional string for describing the promise returned.
  Useful for tooling.
  @return {Promise} a promise that becomes fulfilled with the value the first
  completed promises is resolved with if the first completed promise was
  fulfilled, or rejected with the reason that the first completed promise
  was rejected with.
*/
function race(promises) {
  /*jshint validthis:true */
  var Promise = this;

  if (!isArray(promises)) {
    throw new TypeError('You must pass an array to race.');
  }
  return new Promise(function(resolve, reject) {
    var results = [], promise;

    for (var i = 0; i < promises.length; i++) {
      promise = promises[i];

      if (promise && typeof promise.then === 'function') {
        promise.then(resolve, reject);
      } else {
        resolve(promise);
      }
    }
  });
}

exports.race = race;
},{"./utils":22}],20:[function(require,module,exports){
"use strict";
/**
  `RSVP.reject` returns a promise that will become rejected with the passed
  `reason`. `RSVP.reject` is essentially shorthand for the following:

  ```javascript
  var promise = new RSVP.Promise(function(resolve, reject){
    reject(new Error('WHOOPS'));
  });

  promise.then(function(value){
    // Code here doesn't run because the promise is rejected!
  }, function(reason){
    // reason.message === 'WHOOPS'
  });
  ```

  Instead of writing the above, your code now simply becomes the following:

  ```javascript
  var promise = RSVP.reject(new Error('WHOOPS'));

  promise.then(function(value){
    // Code here doesn't run because the promise is rejected!
  }, function(reason){
    // reason.message === 'WHOOPS'
  });
  ```

  @method reject
  @for RSVP
  @param {Any} reason value that the returned promise will be rejected with.
  @param {String} label optional string for identifying the returned promise.
  Useful for tooling.
  @return {Promise} a promise that will become rejected with the given
  `reason`.
*/
function reject(reason) {
  /*jshint validthis:true */
  var Promise = this;

  return new Promise(function (resolve, reject) {
    reject(reason);
  });
}

exports.reject = reject;
},{}],21:[function(require,module,exports){
"use strict";
/**
  `RSVP.resolve` returns a promise that will become fulfilled with the passed
  `value`. `RSVP.resolve` is essentially shorthand for the following:

  ```javascript
  var promise = new RSVP.Promise(function(resolve, reject){
    resolve(1);
  });

  promise.then(function(value){
    // value === 1
  });
  ```

  Instead of writing the above, your code now simply becomes the following:

  ```javascript
  var promise = RSVP.resolve(1);

  promise.then(function(value){
    // value === 1
  });
  ```

  @method resolve
  @for RSVP
  @param {Any} value value that the returned promise will be resolved with
  @param {String} label optional string for identifying the returned promise.
  Useful for tooling.
  @return {Promise} a promise that will become fulfilled with the given
  `value`
*/
function resolve(value) {
  /*jshint validthis:true */
  var Promise = this;
  return new Promise(function(resolve, reject) {
    resolve(value);
  });
}

exports.resolve = resolve;
},{}],22:[function(require,module,exports){
"use strict";
function objectOrFunction(x) {
  return isFunction(x) || (typeof x === "object" && x !== null);
}

function isFunction(x) {
  return typeof x === "function";
}

function isArray(x) {
  return Object.prototype.toString.call(x) === "[object Array]";
}

// Date.now is not available in browsers < IE9
// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Date/now#Compatibility
var now = Date.now || function() { return new Date().getTime(); };


exports.objectOrFunction = objectOrFunction;
exports.isFunction = isFunction;
exports.isArray = isArray;
exports.now = now;
},{}],23:[function(require,module,exports){
// shim for using process in browser

var process = module.exports = {};

process.nextTick = (function () {
    var canSetImmediate = typeof window !== 'undefined'
    && window.setImmediate;
    var canPost = typeof window !== 'undefined'
    && window.postMessage && window.addEventListener
    ;

    if (canSetImmediate) {
        return function (f) { return window.setImmediate(f) };
    }

    if (canPost) {
        var queue = [];
        window.addEventListener('message', function (ev) {
            var source = ev.source;
            if ((source === window || source === null) && ev.data === 'process-tick') {
                ev.stopPropagation();
                if (queue.length > 0) {
                    var fn = queue.shift();
                    fn();
                }
            }
        }, true);

        return function nextTick(fn) {
            queue.push(fn);
            window.postMessage('process-tick', '*');
        };
    }

    return function nextTick(fn) {
        setTimeout(fn, 0);
    };
})();

process.title = 'browser';
process.browser = true;
process.env = {};
process.argv = [];

function noop() {}

process.on = noop;
process.addListener = noop;
process.once = noop;
process.off = noop;
process.removeListener = noop;
process.removeAllListeners = noop;
process.emit = noop;

process.binding = function (name) {
    throw new Error('process.binding is not supported');
}

// TODO(shtylman)
process.cwd = function () { return '/' };
process.chdir = function (dir) {
    throw new Error('process.chdir is not supported');
};

},{}],24:[function(require,module,exports){
/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Javascript ZLib
// By Thomas Down 2010-2011
//
// Based very heavily on portions of jzlib (by ymnk@jcraft.com), who in
// turn credits Jean-loup Gailly and Mark Adler for the original zlib code.
//
// inflate.js: ZLib inflate code
//

//
// Shared constants
//

var MAX_WBITS=15; // 32K LZ77 window
var DEF_WBITS=MAX_WBITS;
var MAX_MEM_LEVEL=9;
var MANY=1440;
var BMAX = 15;

// preset dictionary flag in zlib header
var PRESET_DICT=0x20;

var Z_NO_FLUSH=0;
var Z_PARTIAL_FLUSH=1;
var Z_SYNC_FLUSH=2;
var Z_FULL_FLUSH=3;
var Z_FINISH=4;

var Z_DEFLATED=8;

var Z_OK=0;
var Z_STREAM_END=1;
var Z_NEED_DICT=2;
var Z_ERRNO=-1;
var Z_STREAM_ERROR=-2;
var Z_DATA_ERROR=-3;
var Z_MEM_ERROR=-4;
var Z_BUF_ERROR=-5;
var Z_VERSION_ERROR=-6;

var METHOD=0;   // waiting for method byte
var FLAG=1;     // waiting for flag byte
var DICT4=2;    // four dictionary check bytes to go
var DICT3=3;    // three dictionary check bytes to go
var DICT2=4;    // two dictionary check bytes to go
var DICT1=5;    // one dictionary check byte to go
var DICT0=6;    // waiting for inflateSetDictionary
var BLOCKS=7;   // decompressing blocks
var CHECK4=8;   // four check bytes to go
var CHECK3=9;   // three check bytes to go
var CHECK2=10;  // two check bytes to go
var CHECK1=11;  // one check byte to go
var DONE=12;    // finished check, done
var BAD=13;     // got an error--stay here

var inflate_mask = [0x00000000, 0x00000001, 0x00000003, 0x00000007, 0x0000000f, 0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff, 0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff, 0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff];

var IB_TYPE=0;  // get type bits (3, including end bit)
var IB_LENS=1;  // get lengths for stored
var IB_STORED=2;// processing stored block
var IB_TABLE=3; // get table lengths
var IB_BTREE=4; // get bit lengths tree for a dynamic block
var IB_DTREE=5; // get length, distance trees for a dynamic block
var IB_CODES=6; // processing fixed or dynamic block
var IB_DRY=7;   // output remaining window bytes
var IB_DONE=8;  // finished last block, done
var IB_BAD=9;   // ot a data error--stuck here

var fixed_bl = 9;
var fixed_bd = 5;

var fixed_tl = [
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,192,
    80,7,10, 0,8,96, 0,8,32, 0,9,160,
    0,8,0, 0,8,128, 0,8,64, 0,9,224,
    80,7,6, 0,8,88, 0,8,24, 0,9,144,
    83,7,59, 0,8,120, 0,8,56, 0,9,208,
    81,7,17, 0,8,104, 0,8,40, 0,9,176,
    0,8,8, 0,8,136, 0,8,72, 0,9,240,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,200,
    81,7,13, 0,8,100, 0,8,36, 0,9,168,
    0,8,4, 0,8,132, 0,8,68, 0,9,232,
    80,7,8, 0,8,92, 0,8,28, 0,9,152,
    84,7,83, 0,8,124, 0,8,60, 0,9,216,
    82,7,23, 0,8,108, 0,8,44, 0,9,184,
    0,8,12, 0,8,140, 0,8,76, 0,9,248,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,196,
    81,7,11, 0,8,98, 0,8,34, 0,9,164,
    0,8,2, 0,8,130, 0,8,66, 0,9,228,
    80,7,7, 0,8,90, 0,8,26, 0,9,148,
    84,7,67, 0,8,122, 0,8,58, 0,9,212,
    82,7,19, 0,8,106, 0,8,42, 0,9,180,
    0,8,10, 0,8,138, 0,8,74, 0,9,244,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,204,
    81,7,15, 0,8,102, 0,8,38, 0,9,172,
    0,8,6, 0,8,134, 0,8,70, 0,9,236,
    80,7,9, 0,8,94, 0,8,30, 0,9,156,
    84,7,99, 0,8,126, 0,8,62, 0,9,220,
    82,7,27, 0,8,110, 0,8,46, 0,9,188,
    0,8,14, 0,8,142, 0,8,78, 0,9,252,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,194,
    80,7,10, 0,8,97, 0,8,33, 0,9,162,
    0,8,1, 0,8,129, 0,8,65, 0,9,226,
    80,7,6, 0,8,89, 0,8,25, 0,9,146,
    83,7,59, 0,8,121, 0,8,57, 0,9,210,
    81,7,17, 0,8,105, 0,8,41, 0,9,178,
    0,8,9, 0,8,137, 0,8,73, 0,9,242,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,202,
    81,7,13, 0,8,101, 0,8,37, 0,9,170,
    0,8,5, 0,8,133, 0,8,69, 0,9,234,
    80,7,8, 0,8,93, 0,8,29, 0,9,154,
    84,7,83, 0,8,125, 0,8,61, 0,9,218,
    82,7,23, 0,8,109, 0,8,45, 0,9,186,
    0,8,13, 0,8,141, 0,8,77, 0,9,250,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,198,
    81,7,11, 0,8,99, 0,8,35, 0,9,166,
    0,8,3, 0,8,131, 0,8,67, 0,9,230,
    80,7,7, 0,8,91, 0,8,27, 0,9,150,
    84,7,67, 0,8,123, 0,8,59, 0,9,214,
    82,7,19, 0,8,107, 0,8,43, 0,9,182,
    0,8,11, 0,8,139, 0,8,75, 0,9,246,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,206,
    81,7,15, 0,8,103, 0,8,39, 0,9,174,
    0,8,7, 0,8,135, 0,8,71, 0,9,238,
    80,7,9, 0,8,95, 0,8,31, 0,9,158,
    84,7,99, 0,8,127, 0,8,63, 0,9,222,
    82,7,27, 0,8,111, 0,8,47, 0,9,190,
    0,8,15, 0,8,143, 0,8,79, 0,9,254,
    96,7,256, 0,8,80, 0,8,16, 84,8,115,
    82,7,31, 0,8,112, 0,8,48, 0,9,193,

    80,7,10, 0,8,96, 0,8,32, 0,9,161,
    0,8,0, 0,8,128, 0,8,64, 0,9,225,
    80,7,6, 0,8,88, 0,8,24, 0,9,145,
    83,7,59, 0,8,120, 0,8,56, 0,9,209,
    81,7,17, 0,8,104, 0,8,40, 0,9,177,
    0,8,8, 0,8,136, 0,8,72, 0,9,241,
    80,7,4, 0,8,84, 0,8,20, 85,8,227,
    83,7,43, 0,8,116, 0,8,52, 0,9,201,
    81,7,13, 0,8,100, 0,8,36, 0,9,169,
    0,8,4, 0,8,132, 0,8,68, 0,9,233,
    80,7,8, 0,8,92, 0,8,28, 0,9,153,
    84,7,83, 0,8,124, 0,8,60, 0,9,217,
    82,7,23, 0,8,108, 0,8,44, 0,9,185,
    0,8,12, 0,8,140, 0,8,76, 0,9,249,
    80,7,3, 0,8,82, 0,8,18, 85,8,163,
    83,7,35, 0,8,114, 0,8,50, 0,9,197,
    81,7,11, 0,8,98, 0,8,34, 0,9,165,
    0,8,2, 0,8,130, 0,8,66, 0,9,229,
    80,7,7, 0,8,90, 0,8,26, 0,9,149,
    84,7,67, 0,8,122, 0,8,58, 0,9,213,
    82,7,19, 0,8,106, 0,8,42, 0,9,181,
    0,8,10, 0,8,138, 0,8,74, 0,9,245,
    80,7,5, 0,8,86, 0,8,22, 192,8,0,
    83,7,51, 0,8,118, 0,8,54, 0,9,205,
    81,7,15, 0,8,102, 0,8,38, 0,9,173,
    0,8,6, 0,8,134, 0,8,70, 0,9,237,
    80,7,9, 0,8,94, 0,8,30, 0,9,157,
    84,7,99, 0,8,126, 0,8,62, 0,9,221,
    82,7,27, 0,8,110, 0,8,46, 0,9,189,
    0,8,14, 0,8,142, 0,8,78, 0,9,253,
    96,7,256, 0,8,81, 0,8,17, 85,8,131,
    82,7,31, 0,8,113, 0,8,49, 0,9,195,
    80,7,10, 0,8,97, 0,8,33, 0,9,163,
    0,8,1, 0,8,129, 0,8,65, 0,9,227,
    80,7,6, 0,8,89, 0,8,25, 0,9,147,
    83,7,59, 0,8,121, 0,8,57, 0,9,211,
    81,7,17, 0,8,105, 0,8,41, 0,9,179,
    0,8,9, 0,8,137, 0,8,73, 0,9,243,
    80,7,4, 0,8,85, 0,8,21, 80,8,258,
    83,7,43, 0,8,117, 0,8,53, 0,9,203,
    81,7,13, 0,8,101, 0,8,37, 0,9,171,
    0,8,5, 0,8,133, 0,8,69, 0,9,235,
    80,7,8, 0,8,93, 0,8,29, 0,9,155,
    84,7,83, 0,8,125, 0,8,61, 0,9,219,
    82,7,23, 0,8,109, 0,8,45, 0,9,187,
    0,8,13, 0,8,141, 0,8,77, 0,9,251,
    80,7,3, 0,8,83, 0,8,19, 85,8,195,
    83,7,35, 0,8,115, 0,8,51, 0,9,199,
    81,7,11, 0,8,99, 0,8,35, 0,9,167,
    0,8,3, 0,8,131, 0,8,67, 0,9,231,
    80,7,7, 0,8,91, 0,8,27, 0,9,151,
    84,7,67, 0,8,123, 0,8,59, 0,9,215,
    82,7,19, 0,8,107, 0,8,43, 0,9,183,
    0,8,11, 0,8,139, 0,8,75, 0,9,247,
    80,7,5, 0,8,87, 0,8,23, 192,8,0,
    83,7,51, 0,8,119, 0,8,55, 0,9,207,
    81,7,15, 0,8,103, 0,8,39, 0,9,175,
    0,8,7, 0,8,135, 0,8,71, 0,9,239,
    80,7,9, 0,8,95, 0,8,31, 0,9,159,
    84,7,99, 0,8,127, 0,8,63, 0,9,223,
    82,7,27, 0,8,111, 0,8,47, 0,9,191,
    0,8,15, 0,8,143, 0,8,79, 0,9,255
];
var fixed_td = [
    80,5,1, 87,5,257, 83,5,17, 91,5,4097,
    81,5,5, 89,5,1025, 85,5,65, 93,5,16385,
    80,5,3, 88,5,513, 84,5,33, 92,5,8193,
    82,5,9, 90,5,2049, 86,5,129, 192,5,24577,
    80,5,2, 87,5,385, 83,5,25, 91,5,6145,
    81,5,7, 89,5,1537, 85,5,97, 93,5,24577,
    80,5,4, 88,5,769, 84,5,49, 92,5,12289,
    82,5,13, 90,5,3073, 86,5,193, 192,5,24577
];

  // Tables for deflate from PKZIP's appnote.txt.
  var cplens = [ // Copy lengths for literal codes 257..285
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258, 0, 0
  ];

  // see note #13 above about 258
  var cplext = [ // Extra bits for literal codes 257..285
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0, 112, 112  // 112==invalid
  ];

 var cpdist = [ // Copy offsets for distance codes 0..29
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577
  ];

  var cpdext = [ // Extra bits for distance codes
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11,
        12, 12, 13, 13];

//
// ZStream.java
//

function ZStream() {
}


ZStream.prototype.inflateInit = function(w, nowrap) {
    if (!w) {
	w = DEF_WBITS;
    }
    if (nowrap) {
	nowrap = false;
    }
    this.istate = new Inflate();
    return this.istate.inflateInit(this, nowrap?-w:w);
}

ZStream.prototype.inflate = function(f) {
    if(this.istate==null) return Z_STREAM_ERROR;
    return this.istate.inflate(this, f);
}

ZStream.prototype.inflateEnd = function(){
    if(this.istate==null) return Z_STREAM_ERROR;
    var ret=istate.inflateEnd(this);
    this.istate = null;
    return ret;
}
ZStream.prototype.inflateSync = function(){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSync(this);
}
ZStream.prototype.inflateSetDictionary = function(dictionary, dictLength){
    // if(istate == null) return Z_STREAM_ERROR;
    return istate.inflateSetDictionary(this, dictionary, dictLength);
}

/*

  public int deflateInit(int level){
    return deflateInit(level, MAX_WBITS);
  }
  public int deflateInit(int level, boolean nowrap){
    return deflateInit(level, MAX_WBITS, nowrap);
  }
  public int deflateInit(int level, int bits){
    return deflateInit(level, bits, false);
  }
  public int deflateInit(int level, int bits, boolean nowrap){
    dstate=new Deflate();
    return dstate.deflateInit(this, level, nowrap?-bits:bits);
  }
  public int deflate(int flush){
    if(dstate==null){
      return Z_STREAM_ERROR;
    }
    return dstate.deflate(this, flush);
  }
  public int deflateEnd(){
    if(dstate==null) return Z_STREAM_ERROR;
    int ret=dstate.deflateEnd();
    dstate=null;
    return ret;
  }
  public int deflateParams(int level, int strategy){
    if(dstate==null) return Z_STREAM_ERROR;
    return dstate.deflateParams(this, level, strategy);
  }
  public int deflateSetDictionary (byte[] dictionary, int dictLength){
    if(dstate == null)
      return Z_STREAM_ERROR;
    return dstate.deflateSetDictionary(this, dictionary, dictLength);
  }

*/

/*
  // Flush as much pending output as possible. All deflate() output goes
  // through this function so some applications may wish to modify it
  // to avoid allocating a large strm->next_out buffer and copying into it.
  // (See also read_buf()).
  void flush_pending(){
    int len=dstate.pending;

    if(len>avail_out) len=avail_out;
    if(len==0) return;

    if(dstate.pending_buf.length<=dstate.pending_out ||
       next_out.length<=next_out_index ||
       dstate.pending_buf.length<(dstate.pending_out+len) ||
       next_out.length<(next_out_index+len)){
      System.out.println(dstate.pending_buf.length+", "+dstate.pending_out+
			 ", "+next_out.length+", "+next_out_index+", "+len);
      System.out.println("avail_out="+avail_out);
    }

    System.arraycopy(dstate.pending_buf, dstate.pending_out,
		     next_out, next_out_index, len);

    next_out_index+=len;
    dstate.pending_out+=len;
    total_out+=len;
    avail_out-=len;
    dstate.pending-=len;
    if(dstate.pending==0){
      dstate.pending_out=0;
    }
  }

  // Read a new buffer from the current input stream, update the adler32
  // and total number of bytes read.  All deflate() input goes through
  // this function so some applications may wish to modify it to avoid
  // allocating a large strm->next_in buffer and copying from it.
  // (See also flush_pending()).
  int read_buf(byte[] buf, int start, int size) {
    int len=avail_in;

    if(len>size) len=size;
    if(len==0) return 0;

    avail_in-=len;

    if(dstate.noheader==0) {
      adler=_adler.adler32(adler, next_in, next_in_index, len);
    }
    System.arraycopy(next_in, next_in_index, buf, start, len);
    next_in_index  += len;
    total_in += len;
    return len;
  }

  public void free(){
    next_in=null;
    next_out=null;
    msg=null;
    _adler=null;
  }
}
*/


//
// Inflate.java
//

function Inflate() {
    this.was = [0];
}

Inflate.prototype.inflateReset = function(z) {
    if(z == null || z.istate == null) return Z_STREAM_ERROR;
    
    z.total_in = z.total_out = 0;
    z.msg = null;
    z.istate.mode = z.istate.nowrap!=0 ? BLOCKS : METHOD;
    z.istate.blocks.reset(z, null);
    return Z_OK;
}

Inflate.prototype.inflateEnd = function(z){
    if(this.blocks != null)
      this.blocks.free(z);
    this.blocks=null;
    return Z_OK;
}

Inflate.prototype.inflateInit = function(z, w){
    z.msg = null;
    this.blocks = null;

    // handle undocumented nowrap option (no zlib header or check)
    nowrap = 0;
    if(w < 0){
      w = - w;
      nowrap = 1;
    }

    // set window size
    if(w<8 ||w>15){
      this.inflateEnd(z);
      return Z_STREAM_ERROR;
    }
    this.wbits=w;

    z.istate.blocks=new InfBlocks(z, 
				  z.istate.nowrap!=0 ? null : this,
				  1<<w);

    // reset state
    this.inflateReset(z);
    return Z_OK;
  }

Inflate.prototype.inflate = function(z, f){
    var r, b;

    if(z == null || z.istate == null || z.next_in == null)
      return Z_STREAM_ERROR;
    f = f == Z_FINISH ? Z_BUF_ERROR : Z_OK;
    r = Z_BUF_ERROR;
    while (true){
      switch (z.istate.mode){
      case METHOD:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        if(((z.istate.method = z.next_in[z.next_in_index++])&0xf)!=Z_DEFLATED){
          z.istate.mode = BAD;
          z.msg="unknown compression method";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        if((z.istate.method>>4)+8>z.istate.wbits){
          z.istate.mode = BAD;
          z.msg="invalid window size";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }
        z.istate.mode=FLAG;
      case FLAG:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        b = (z.next_in[z.next_in_index++])&0xff;

        if((((z.istate.method << 8)+b) % 31)!=0){
          z.istate.mode = BAD;
          z.msg = "incorrect header check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        if((b&PRESET_DICT)==0){
          z.istate.mode = BLOCKS;
          break;
        }
        z.istate.mode = DICT4;
      case DICT4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=DICT3;
      case DICT3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode=DICT2;
      case DICT2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode=DICT1;
      case DICT1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need += (z.next_in[z.next_in_index++]&0xff);
        z.adler = z.istate.need;
        z.istate.mode = DICT0;
        return Z_NEED_DICT;
      case DICT0:
        z.istate.mode = BAD;
        z.msg = "need dictionary";
        z.istate.marker = 0;       // can try inflateSync
        return Z_STREAM_ERROR;
      case BLOCKS:

        r = z.istate.blocks.proc(z, r);
        if(r == Z_DATA_ERROR){
          z.istate.mode = BAD;
          z.istate.marker = 0;     // can try inflateSync
          break;
        }
        if(r == Z_OK){
          r = f;
        }
        if(r != Z_STREAM_END){
          return r;
        }
        r = f;
        z.istate.blocks.reset(z, z.istate.was);
        if(z.istate.nowrap!=0){
          z.istate.mode=DONE;
          break;
        }
        z.istate.mode=CHECK4;
      case CHECK4:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need=((z.next_in[z.next_in_index++]&0xff)<<24)&0xff000000;
        z.istate.mode=CHECK3;
      case CHECK3:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<16)&0xff0000;
        z.istate.mode = CHECK2;
      case CHECK2:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=((z.next_in[z.next_in_index++]&0xff)<<8)&0xff00;
        z.istate.mode = CHECK1;
      case CHECK1:

        if(z.avail_in==0)return r;r=f;

        z.avail_in--; z.total_in++;
        z.istate.need+=(z.next_in[z.next_in_index++]&0xff);

        if(((z.istate.was[0])) != ((z.istate.need))){
          z.istate.mode = BAD;
          z.msg = "incorrect data check";
          z.istate.marker = 5;       // can't try inflateSync
          break;
        }

        z.istate.mode = DONE;
      case DONE:
        return Z_STREAM_END;
      case BAD:
        return Z_DATA_ERROR;
      default:
        return Z_STREAM_ERROR;
      }
    }
  }


Inflate.prototype.inflateSetDictionary = function(z,  dictionary, dictLength) {
    var index=0;
    var length = dictLength;
    if(z==null || z.istate == null|| z.istate.mode != DICT0)
      return Z_STREAM_ERROR;

    if(z._adler.adler32(1, dictionary, 0, dictLength)!=z.adler){
      return Z_DATA_ERROR;
    }

    z.adler = z._adler.adler32(0, null, 0, 0);

    if(length >= (1<<z.istate.wbits)){
      length = (1<<z.istate.wbits)-1;
      index=dictLength - length;
    }
    z.istate.blocks.set_dictionary(dictionary, index, length);
    z.istate.mode = BLOCKS;
    return Z_OK;
  }

//  static private byte[] mark = {(byte)0, (byte)0, (byte)0xff, (byte)0xff};
var mark = [0, 0, 255, 255]

Inflate.prototype.inflateSync = function(z){
    var n;       // number of bytes to look at
    var p;       // pointer to bytes
    var m;       // number of marker bytes found in a row
    var r, w;   // temporaries to save total_in and total_out

    // set up
    if(z == null || z.istate == null)
      return Z_STREAM_ERROR;
    if(z.istate.mode != BAD){
      z.istate.mode = BAD;
      z.istate.marker = 0;
    }
    if((n=z.avail_in)==0)
      return Z_BUF_ERROR;
    p=z.next_in_index;
    m=z.istate.marker;

    // search
    while (n!=0 && m < 4){
      if(z.next_in[p] == mark[m]){
        m++;
      }
      else if(z.next_in[p]!=0){
        m = 0;
      }
      else{
        m = 4 - m;
      }
      p++; n--;
    }

    // restore
    z.total_in += p-z.next_in_index;
    z.next_in_index = p;
    z.avail_in = n;
    z.istate.marker = m;

    // return no joy or set up to restart on a new block
    if(m != 4){
      return Z_DATA_ERROR;
    }
    r=z.total_in;  w=z.total_out;
    this.inflateReset(z);
    z.total_in=r;  z.total_out = w;
    z.istate.mode = BLOCKS;
    return Z_OK;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. This function is used by one PPP
  // implementation to provide an additional safety check. PPP uses Z_SYNC_FLUSH
  // but removes the length bytes of the resulting empty stored block. When
  // decompressing, PPP checks that at the end of input packet, inflate is
  // waiting for these length bytes.
Inflate.prototype.inflateSyncPoint = function(z){
    if(z == null || z.istate == null || z.istate.blocks == null)
      return Z_STREAM_ERROR;
    return z.istate.blocks.sync_point();
}


//
// InfBlocks.java
//

var INFBLOCKS_BORDER = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

function InfBlocks(z, checkfn, w) {
    this.hufts=new Int32Array(MANY*3);
    this.window=new Uint8Array(w);
    this.end=w;
    this.checkfn = checkfn;
    this.mode = IB_TYPE;
    this.reset(z, null);

    this.left = 0;            // if STORED, bytes left to copy 

    this.table = 0;           // table lengths (14 bits) 
    this.index = 0;           // index into blens (or border) 
    this.blens = null;         // bit lengths of codes 
    this.bb=new Int32Array(1); // bit length tree depth 
    this.tb=new Int32Array(1); // bit length decoding tree 

    this.codes = new InfCodes();

    this.last = 0;            // true if this block is the last block 

  // mode independent information 
    this.bitk = 0;            // bits in bit buffer 
    this.bitb = 0;            // bit buffer 
    this.read = 0;            // window read pointer 
    this.write = 0;           // window write pointer 
    this.check = 0;          // check on output 

    this.inftree=new InfTree();
}




InfBlocks.prototype.reset = function(z, c){
    if(c) c[0]=this.check;
    if(this.mode==IB_CODES){
      this.codes.free(z);
    }
    this.mode=IB_TYPE;
    this.bitk=0;
    this.bitb=0;
    this.read=this.write=0;

    if(this.checkfn)
      z.adler=this.check=z._adler.adler32(0, null, 0, 0);
  }

 InfBlocks.prototype.proc = function(z, r){
    var t;              // temporary storage
    var b;              // bit buffer
    var k;              // bits in bit buffer
    var p;              // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer

    // copy input/output information to locals (UPDATE macro restores)
    {p=z.next_in_index;n=z.avail_in;b=this.bitb;k=this.bitk;}
    {q=this.write;m=(q<this.read ? this.read-q-1 : this.end-q);}

    // process input based on current state
    while(true){
      switch (this.mode){
      case IB_TYPE:

	while(k<(3)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}
	t = (b & 7);
	this.last = t & 1;

	switch (t >>> 1){
        case 0:                         // stored 
          {b>>>=(3);k-=(3);}
          t = k & 7;                    // go to byte boundary

          {b>>>=(t);k-=(t);}
          this.mode = IB_LENS;                  // get length of stored block
          break;
        case 1:                         // fixed
          {
              var bl=new Int32Array(1);
	      var bd=new Int32Array(1);
              var tl=[];
	      var td=[];

	      inflate_trees_fixed(bl, bd, tl, td, z);
              this.codes.init(bl[0], bd[0], tl[0], 0, td[0], 0, z);
          }

          {b>>>=(3);k-=(3);}

          this.mode = IB_CODES;
          break;
        case 2:                         // dynamic

          {b>>>=(3);k-=(3);}

          this.mode = IB_TABLE;
          break;
        case 3:                         // illegal

          {b>>>=(3);k-=(3);}
          this.mode = BAD;
          z.msg = "invalid block type";
          r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	break;
      case IB_LENS:
	while(k<(32)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	if ((((~b) >>> 16) & 0xffff) != (b & 0xffff)){
	  this.mode = BAD;
	  z.msg = "invalid stored block lengths";
	  r = Z_DATA_ERROR;

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z,r);
	}
	this.left = (b & 0xffff);
	b = k = 0;                       // dump bits
	this.mode = this.left!=0 ? IB_STORED : (this.last!=0 ? IB_DRY : IB_TYPE);
	break;
      case IB_STORED:
	if (n == 0){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	if(m==0){
	  if(q==end&&read!=0){
	    q=0; m=(q<this.read ? this.read-q-1 : this.end-q);
	  }
	  if(m==0){
	    this.write=q; 
	    r=this.inflate_flush(z,r);
	    q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	    if(q==this.end && this.read != 0){
	      q=0; m = (q < this.read ? this.read-q-1 : this.end-q);
	    }
	    if(m==0){
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	t = this.left;
	if(t>n) t = n;
	if(t>m) t = m;
	arrayCopy(z.next_in, p, this.window, q, t);
	p += t;  n -= t;
	q += t;  m -= t;
	if ((this.left -= t) != 0)
	  break;
	this.mode = (this.last != 0 ? IB_DRY : IB_TYPE);
	break;
      case IB_TABLE:

	while(k<(14)){
	  if(n!=0){
	    r=Z_OK;
	  }
	  else{
	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;
	    z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  };
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.table = t = (b & 0x3fff);
	if ((t & 0x1f) > 29 || ((t >> 5) & 0x1f) > 29)
	  {
	    this.mode = IB_BAD;
	    z.msg = "too many length or distance symbols";
	    r = Z_DATA_ERROR;

	    this.bitb=b; this.bitk=k; 
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    this.write=q;
	    return this.inflate_flush(z,r);
	  }
	t = 258 + (t & 0x1f) + ((t >> 5) & 0x1f);
	if(this.blens==null || this.blens.length<t){
	    this.blens=new Int32Array(t);
	}
	else{
	  for(var i=0; i<t; i++){
              this.blens[i]=0;
          }
	}

	{b>>>=(14);k-=(14);}

	this.index = 0;
	mode = IB_BTREE;
      case IB_BTREE:
	while (this.index < 4 + (this.table >>> 10)){
	  while(k<(3)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

	  this.blens[INFBLOCKS_BORDER[this.index++]] = b&7;

	  {b>>>=(3);k-=(3);}
	}

	while(this.index < 19){
	  this.blens[INFBLOCKS_BORDER[this.index++]] = 0;
	}

	this.bb[0] = 7;
	t = this.inftree.inflate_trees_bits(this.blens, this.bb, this.tb, this.hufts, z);
	if (t != Z_OK){
	  r = t;
	  if (r == Z_DATA_ERROR){
	    this.blens=null;
	    this.mode = IB_BAD;
	  }

	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  write=q;
	  return this.inflate_flush(z,r);
	}

	this.index = 0;
	this.mode = IB_DTREE;
      case IB_DTREE:
	while (true){
	  t = this.table;
	  if(!(this.index < 258 + (t & 0x1f) + ((t >> 5) & 0x1f))){
	    break;
	  }

	  var h; //int[]
	  var i, j, c;

	  t = this.bb[0];

	  while(k<(t)){
	    if(n!=0){
	      r=Z_OK;
	    }
	    else{
	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;
	      z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    };
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;
	    k+=8;
	  }

//	  if (this.tb[0]==-1){
//            dlog("null...");
//	  }

	  t=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+1];
	  c=this.hufts[(this.tb[0]+(b & inflate_mask[t]))*3+2];

	  if (c < 16){
	    b>>>=(t);k-=(t);
	    this.blens[this.index++] = c;
	  }
	  else { // c == 16..18
	    i = c == 18 ? 7 : c - 14;
	    j = c == 18 ? 11 : 3;

	    while(k<(t+i)){
	      if(n!=0){
		r=Z_OK;
	      }
	      else{
		this.bitb=b; this.bitk=k; 
		z.avail_in=n;
		z.total_in+=p-z.next_in_index;z.next_in_index=p;
		this.write=q;
		return this.inflate_flush(z,r);
	      };
	      n--;
	      b|=(z.next_in[p++]&0xff)<<k;
	      k+=8;
	    }

	    b>>>=(t);k-=(t);

	    j += (b & inflate_mask[i]);

	    b>>>=(i);k-=(i);

	    i = this.index;
	    t = this.table;
	    if (i + j > 258 + (t & 0x1f) + ((t >> 5) & 0x1f) ||
		(c == 16 && i < 1)){
	      this.blens=null;
	      this.mode = IB_BAD;
	      z.msg = "invalid bit length repeat";
	      r = Z_DATA_ERROR;

	      this.bitb=b; this.bitk=k; 
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      this.write=q;
	      return this.inflate_flush(z,r);
	    }

	    c = c == 16 ? this.blens[i-1] : 0;
	    do{
	      this.blens[i++] = c;
	    }
	    while (--j!=0);
	    this.index = i;
	  }
	}

	this.tb[0]=-1;
	{
	    var bl=new Int32Array(1);
	    var bd=new Int32Array(1);
	    var tl=new Int32Array(1);
	    var td=new Int32Array(1);
	    bl[0] = 9;         // must be <= 9 for lookahead assumptions
	    bd[0] = 6;         // must be <= 9 for lookahead assumptions

	    t = this.table;
	    t = this.inftree.inflate_trees_dynamic(257 + (t & 0x1f), 
					      1 + ((t >> 5) & 0x1f),
					      this.blens, bl, bd, tl, td, this.hufts, z);

	    if (t != Z_OK){
	        if (t == Z_DATA_ERROR){
	            this.blens=null;
	            this.mode = BAD;
	        }
	        r = t;

	        this.bitb=b; this.bitk=k; 
	        z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	        this.write=q;
	        return this.inflate_flush(z,r);
	    }
	    this.codes.init(bl[0], bd[0], this.hufts, tl[0], this.hufts, td[0], z);
	}
	this.mode = IB_CODES;
      case IB_CODES:
	this.bitb=b; this.bitk=k;
	z.avail_in=n; z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;

	if ((r = this.codes.proc(this, z, r)) != Z_STREAM_END){
	  return this.inflate_flush(z, r);
	}
	r = Z_OK;
	this.codes.free(z);

	p=z.next_in_index; n=z.avail_in;b=this.bitb;k=this.bitk;
	q=this.write;m = (q < this.read ? this.read-q-1 : this.end-q);

	if (this.last==0){
	  this.mode = IB_TYPE;
	  break;
	}
	this.mode = IB_DRY;
      case IB_DRY:
	this.write=q; 
	r = this.inflate_flush(z, r); 
	q=this.write; m = (q < this.read ? this.read-q-1 : this.end-q);
	if (this.read != this.write){
	  this.bitb=b; this.bitk=k; 
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  this.write=q;
	  return this.inflate_flush(z, r);
	}
	mode = DONE;
      case IB_DONE:
	r = Z_STREAM_END;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      case IB_BAD:
	r = Z_DATA_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);

      default:
	r = Z_STREAM_ERROR;

	this.bitb=b; this.bitk=k; 
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	this.write=q;
	return this.inflate_flush(z, r);
      }
    }
  }

InfBlocks.prototype.free = function(z){
    this.reset(z, null);
    this.window=null;
    this.hufts=null;
}

InfBlocks.prototype.set_dictionary = function(d, start, n){
    arrayCopy(d, start, window, 0, n);
    this.read = this.write = n;
}

  // Returns true if inflate is currently at the end of a block generated
  // by Z_SYNC_FLUSH or Z_FULL_FLUSH. 
InfBlocks.prototype.sync_point = function(){
    return this.mode == IB_LENS;
}

  // copy as much as possible from the sliding window to the output area
InfBlocks.prototype.inflate_flush = function(z, r){
    var n;
    var p;
    var q;

    // local copies of source and destination pointers
    p = z.next_out_index;
    q = this.read;

    // compute number of bytes to copy as far as end of window
    n = ((q <= this.write ? this.write : this.end) - q);
    if (n > z.avail_out) n = z.avail_out;
    if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

    // update counters
    z.avail_out -= n;
    z.total_out += n;

    // update check information
    if(this.checkfn != null)
      z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

    // copy as far as end of window
    arrayCopy(this.window, q, z.next_out, p, n);
    p += n;
    q += n;

    // see if more to copy at beginning of window
    if (q == this.end){
      // wrap pointers
      q = 0;
      if (this.write == this.end)
        this.write = 0;

      // compute bytes to copy
      n = this.write - q;
      if (n > z.avail_out) n = z.avail_out;
      if (n!=0 && r == Z_BUF_ERROR) r = Z_OK;

      // update counters
      z.avail_out -= n;
      z.total_out += n;

      // update check information
      if(this.checkfn != null)
	z.adler=this.check=z._adler.adler32(this.check, this.window, q, n);

      // copy
      arrayCopy(this.window, q, z.next_out, p, n);
      p += n;
      q += n;
    }

    // update pointers
    z.next_out_index = p;
    this.read = q;

    // done
    return r;
  }

//
// InfCodes.java
//

var IC_START=0;  // x: set up for LEN
var IC_LEN=1;    // i: get length/literal/eob next
var IC_LENEXT=2; // i: getting length extra (have base)
var IC_DIST=3;   // i: get distance next
var IC_DISTEXT=4;// i: getting distance extra
var IC_COPY=5;   // o: copying bytes in window, waiting for space
var IC_LIT=6;    // o: got literal, waiting for output space
var IC_WASH=7;   // o: got eob, possibly still output waiting
var IC_END=8;    // x: got eob and all data flushed
var IC_BADCODE=9;// x: got error

function InfCodes() {
}

InfCodes.prototype.init = function(bl, bd, tl, tl_index, td, td_index, z) {
    this.mode=IC_START;
    this.lbits=bl;
    this.dbits=bd;
    this.ltree=tl;
    this.ltree_index=tl_index;
    this.dtree = td;
    this.dtree_index=td_index;
    this.tree=null;
}

InfCodes.prototype.proc = function(s, z, r){ 
    var j;              // temporary storage
    var t;              // temporary pointer (int[])
    var tindex;         // temporary pointer
    var e;              // extra bits or operation
    var b=0;            // bit buffer
    var k=0;            // bits in bit buffer
    var p=0;            // input data pointer
    var n;              // bytes available there
    var q;              // output window write pointer
    var m;              // bytes to end of window or read pointer
    var f;              // pointer to copy strings from

    // copy input/output information to locals (UPDATE macro restores)
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // process input and output based on current state
    while (true){
      switch (this.mode){
	// waiting for "i:"=input, "o:"=output, "x:"=nothing
      case IC_START:         // x: set up for LEN
	if (m >= 258 && n >= 10){

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  r = this.inflate_fast(this.lbits, this.dbits, 
			   this.ltree, this.ltree_index, 
			   this.dtree, this.dtree_index,
			   s, z);

	  p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
	  q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	  if (r != Z_OK){
	    this.mode = r == Z_STREAM_END ? IC_WASH : IC_BADCODE;
	    break;
	  }
	}
	this.need = this.lbits;
	this.tree = this.ltree;
	this.tree_index=this.ltree_index;

	this.mode = IC_LEN;
      case IC_LEN:           // i: get length/literal/eob next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--;
	  b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b&inflate_mask[j]))*3;

	b>>>=(this.tree[tindex+1]);
	k-=(this.tree[tindex+1]);

	e=this.tree[tindex];

	if(e == 0){               // literal
	  this.lit = this.tree[tindex+2];
	  this.mode = IC_LIT;
	  break;
	}
	if((e & 16)!=0 ){          // length
	  this.get = e & 15;
	  this.len = this.tree[tindex+2];
	  this.mode = IC_LENEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	if ((e & 32)!=0){               // end of block
	  this.mode = IC_WASH;
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid literal/length code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_LENEXT:        // i: getting length extra (have base)
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.len += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.need = this.dbits;
	this.tree = this.dtree;
	this.tree_index = this.dtree_index;
	this.mode = IC_DIST;
      case IC_DIST:          // i: get distance next
	j = this.need;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	tindex=(this.tree_index+(b & inflate_mask[j]))*3;

	b>>=this.tree[tindex+1];
	k-=this.tree[tindex+1];

	e = (this.tree[tindex]);
	if((e & 16)!=0){               // distance
	  this.get = e & 15;
	  this.dist = this.tree[tindex+2];
	  this.mode = IC_DISTEXT;
	  break;
	}
	if ((e & 64) == 0){        // next table
	  this.need = e;
	  this.tree_index = tindex/3 + this.tree[tindex+2];
	  break;
	}
	this.mode = IC_BADCODE;        // invalid code
	z.msg = "invalid distance code";
	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_DISTEXT:       // i: getting distance extra
	j = this.get;

	while(k<(j)){
	  if(n!=0)r=Z_OK;
	  else{

	    s.bitb=b;s.bitk=k;
	    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	    s.write=q;
	    return s.inflate_flush(z,r);
	  }
	  n--; b|=(z.next_in[p++]&0xff)<<k;
	  k+=8;
	}

	this.dist += (b & inflate_mask[j]);

	b>>=j;
	k-=j;

	this.mode = IC_COPY;
      case IC_COPY:          // o: copying bytes in window, waiting for space
        f = q - this.dist;
        while(f < 0){     // modulo window size-"while" instead
          f += s.end;     // of "if" handles invalid distances
	}
	while (this.len!=0){

	  if(m==0){
	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.write=q; r=s.inflate_flush(z,r);
	      q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	      if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}

	      if(m==0){
		s.bitb=b;s.bitk=k;
		z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
		s.write=q;
		return s.inflate_flush(z,r);
	      }  
	    }
	  }

	  s.window[q++]=s.window[f++]; m--;

	  if (f == s.end)
            f = 0;
	  this.len--;
	}
	this.mode = IC_START;
	break;
      case IC_LIT:           // o: got literal, waiting for output space
	if(m==0){
	  if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	  if(m==0){
	    s.write=q; r=s.inflate_flush(z,r);
	    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	    if(q==s.end&&s.read!=0){q=0;m=q<s.read?s.read-q-1:s.end-q;}
	    if(m==0){
	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;
	      return s.inflate_flush(z,r);
	    }
	  }
	}
	r=Z_OK;

	s.window[q++]=this.lit; m--;

	this.mode = IC_START;
	break;
      case IC_WASH:           // o: got eob, possibly more output
	if (k > 7){        // return unused byte, if any
	  k -= 8;
	  n++;
	  p--;             // can always return one
	}

	s.write=q; r=s.inflate_flush(z,r);
	q=s.write;m=q<s.read?s.read-q-1:s.end-q;

	if (s.read != s.write){
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;
	  return s.inflate_flush(z,r);
	}
	this.mode = IC_END;
      case IC_END:
	r = Z_STREAM_END;
	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      case IC_BADCODE:       // x: got error

	r = Z_DATA_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);

      default:
	r = Z_STREAM_ERROR;

	s.bitb=b;s.bitk=k;
	z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	s.write=q;
	return s.inflate_flush(z,r);
      }
    }
  }

InfCodes.prototype.free = function(z){
    //  ZFREE(z, c);
}

  // Called with number of bytes left to write in window at least 258
  // (the maximum string length) and number of input bytes available
  // at least ten.  The ten bytes are six bytes for the longest length/
  // distance pair plus four bytes for overloading the bit buffer.

InfCodes.prototype.inflate_fast = function(bl, bd, tl, tl_index, td, td_index, s, z) {
    var t;                // temporary pointer
    var   tp;             // temporary pointer (int[])
    var tp_index;         // temporary pointer
    var e;                // extra bits or operation
    var b;                // bit buffer
    var k;                // bits in bit buffer
    var p;                // input data pointer
    var n;                // bytes available there
    var q;                // output window write pointer
    var m;                // bytes to end of window or read pointer
    var ml;               // mask for literal/length tree
    var md;               // mask for distance tree
    var c;                // bytes to copy
    var d;                // distance back to copy from
    var r;                // copy source pointer

    var tp_index_t_3;     // (tp_index+t)*3

    // load input, output, bit values
    p=z.next_in_index;n=z.avail_in;b=s.bitb;k=s.bitk;
    q=s.write;m=q<s.read?s.read-q-1:s.end-q;

    // initialize masks
    ml = inflate_mask[bl];
    md = inflate_mask[bd];

    // do until not enough input or output space for fast loop
    do {                          // assume called with m >= 258 && n >= 10
      // get literal/length code
      while(k<(20)){              // max bits for literal/length code
	n--;
	b|=(z.next_in[p++]&0xff)<<k;k+=8;
      }

      t= b&ml;
      tp=tl; 
      tp_index=tl_index;
      tp_index_t_3=(tp_index+t)*3;
      if ((e = tp[tp_index_t_3]) == 0){
	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	s.window[q++] = tp[tp_index_t_3+2];
	m--;
	continue;
      }
      do {

	b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	if((e&16)!=0){
	  e &= 15;
	  c = tp[tp_index_t_3+2] + (b & inflate_mask[e]);

	  b>>=e; k-=e;

	  // decode distance base of block to copy
	  while(k<(15)){           // max bits for distance code
	    n--;
	    b|=(z.next_in[p++]&0xff)<<k;k+=8;
	  }

	  t= b&md;
	  tp=td;
	  tp_index=td_index;
          tp_index_t_3=(tp_index+t)*3;
	  e = tp[tp_index_t_3];

	  do {

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    if((e&16)!=0){
	      // get extra bits to add to distance base
	      e &= 15;
	      while(k<(e)){         // get extra bits (up to 13)
		n--;
		b|=(z.next_in[p++]&0xff)<<k;k+=8;
	      }

	      d = tp[tp_index_t_3+2] + (b&inflate_mask[e]);

	      b>>=(e); k-=(e);

	      // do the copy
	      m -= c;
	      if (q >= d){                // offset before dest
		//  just copy
		r=q-d;
		if(q-r>0 && 2>(q-r)){           
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
		else{
		  s.window[q++]=s.window[r++]; // minimum count is three,
		  s.window[q++]=s.window[r++]; // so unroll loop a little
		  c-=2;
		}
	      }
	      else{                  // else offset after destination
                r=q-d;
                do{
                  r+=s.end;          // force pointer in window
                }while(r<0);         // covers invalid distances
		e=s.end-r;
		if(c>e){             // if source crosses,
		  c-=e;              // wrapped copy
		  if(q-r>0 && e>(q-r)){           
		    do{s.window[q++] = s.window[r++];}
		    while(--e!=0);
		  }
		  else{
		    arrayCopy(s.window, r, s.window, q, e);
		    q+=e; r+=e; e=0;
		  }
		  r = 0;                  // copy rest from start of window
		}

	      }

	      // copy all or what's left
              do{s.window[q++] = s.window[r++];}
		while(--c!=0);
	      break;
	    }
	    else if((e&64)==0){
	      t+=tp[tp_index_t_3+2];
	      t+=(b&inflate_mask[e]);
	      tp_index_t_3=(tp_index+t)*3;
	      e=tp[tp_index_t_3];
	    }
	    else{
	      z.msg = "invalid distance code";

	      c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	      s.bitb=b;s.bitk=k;
	      z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	      s.write=q;

	      return Z_DATA_ERROR;
	    }
	  }
	  while(true);
	  break;
	}

	if((e&64)==0){
	  t+=tp[tp_index_t_3+2];
	  t+=(b&inflate_mask[e]);
	  tp_index_t_3=(tp_index+t)*3;
	  if((e=tp[tp_index_t_3])==0){

	    b>>=(tp[tp_index_t_3+1]); k-=(tp[tp_index_t_3+1]);

	    s.window[q++]=tp[tp_index_t_3+2];
	    m--;
	    break;
	  }
	}
	else if((e&32)!=0){

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;
 
	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_STREAM_END;
	}
	else{
	  z.msg="invalid literal/length code";

	  c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

	  s.bitb=b;s.bitk=k;
	  z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
	  s.write=q;

	  return Z_DATA_ERROR;
	}
      } 
      while(true);
    } 
    while(m>=258 && n>= 10);

    // not enough input or output--restore pointers and return
    c=z.avail_in-n;c=(k>>3)<c?k>>3:c;n+=c;p-=c;k-=c<<3;

    s.bitb=b;s.bitk=k;
    z.avail_in=n;z.total_in+=p-z.next_in_index;z.next_in_index=p;
    s.write=q;

    return Z_OK;
}

//
// InfTree.java
//

function InfTree() {
}

InfTree.prototype.huft_build = function(b, bindex, n, s, d, e, t, m, hp, hn, v) {

    // Given a list of code lengths and a maximum table size, make a set of
    // tables to decode that set of codes.  Return Z_OK on success, Z_BUF_ERROR
    // if the given code set is incomplete (the tables are still built in this
    // case), Z_DATA_ERROR if the input is invalid (an over-subscribed set of
    // lengths), or Z_MEM_ERROR if not enough memory.

    var a;                       // counter for codes of length k
    var f;                       // i repeats in table every f entries
    var g;                       // maximum code length
    var h;                       // table level
    var i;                       // counter, current code
    var j;                       // counter
    var k;                       // number of bits in current code
    var l;                       // bits per table (returned in m)
    var mask;                    // (1 << w) - 1, to avoid cc -O bug on HP
    var p;                       // pointer into c[], b[], or v[]
    var q;                       // points to current table
    var w;                       // bits before this table == (l * h)
    var xp;                      // pointer into x
    var y;                       // number of dummy codes added
    var z;                       // number of entries in current table

    // Generate counts for each bit length

    p = 0; i = n;
    do {
      this.c[b[bindex+p]]++; p++; i--;   // assume all entries <= BMAX
    }while(i!=0);

    if(this.c[0] == n){                // null input--all zero length codes
      t[0] = -1;
      m[0] = 0;
      return Z_OK;
    }

    // Find minimum and maximum length, bound *m by those
    l = m[0];
    for (j = 1; j <= BMAX; j++)
      if(this.c[j]!=0) break;
    k = j;                        // minimum code length
    if(l < j){
      l = j;
    }
    for (i = BMAX; i!=0; i--){
      if(this.c[i]!=0) break;
    }
    g = i;                        // maximum code length
    if(l > i){
      l = i;
    }
    m[0] = l;

    // Adjust last length count to fill out codes, if needed
    for (y = 1 << j; j < i; j++, y <<= 1){
      if ((y -= this.c[j]) < 0){
        return Z_DATA_ERROR;
      }
    }
    if ((y -= this.c[i]) < 0){
      return Z_DATA_ERROR;
    }
    this.c[i] += y;

    // Generate starting offsets into the value table for each length
    this.x[1] = j = 0;
    p = 1;  xp = 2;
    while (--i!=0) {                 // note that i == g from above
      this.x[xp] = (j += this.c[p]);
      xp++;
      p++;
    }

    // Make a table of values in order of bit lengths
    i = 0; p = 0;
    do {
      if ((j = b[bindex+p]) != 0){
        this.v[this.x[j]++] = i;
      }
      p++;
    }
    while (++i < n);
    n = this.x[g];                     // set n to length of v

    // Generate the Huffman codes and for each, make the table entries
    this.x[0] = i = 0;                 // first Huffman code is zero
    p = 0;                        // grab values in bit order
    h = -1;                       // no tables yet--level -1
    w = -l;                       // bits decoded == (l * h)
    this.u[0] = 0;                     // just to keep compilers happy
    q = 0;                        // ditto
    z = 0;                        // ditto

    // go through the bit lengths (k already is bits in shortest code)
    for (; k <= g; k++){
      a = this.c[k];
      while (a--!=0){
	// here i is the Huffman code of length k bits for value *p
	// make tables up to required level
        while (k > w + l){
          h++;
          w += l;                 // previous table always l bits
	  // compute minimum size table less than or equal to l bits
          z = g - w;
          z = (z > l) ? l : z;        // table size upper limit
          if((f=1<<(j=k-w))>a+1){     // try a k-w bit table
                                      // too few codes for k-w bit table
            f -= a + 1;               // deduct codes from patterns left
            xp = k;
            if(j < z){
              while (++j < z){        // try smaller tables up to z bits
                if((f <<= 1) <= this.c[++xp])
                  break;              // enough codes to use up j bits
                f -= this.c[xp];           // else deduct codes from patterns
              }
	    }
          }
          z = 1 << j;                 // table entries for j-bit table

	  // allocate new table
          if (this.hn[0] + z > MANY){       // (note: doesn't matter for fixed)
            return Z_DATA_ERROR;       // overflow of MANY
          }
          this.u[h] = q = /*hp+*/ this.hn[0];   // DEBUG
          this.hn[0] += z;
 
	  // connect to last table, if there is one
	  if(h!=0){
            this.x[h]=i;           // save pattern for backing up
            this.r[0]=j;     // bits in this table
            this.r[1]=l;     // bits to dump before this table
            j=i>>>(w - l);
            this.r[2] = (q - this.u[h-1] - j);               // offset to this table
            arrayCopy(this.r, 0, hp, (this.u[h-1]+j)*3, 3); // connect to last table
          }
          else{
            t[0] = q;               // first table is returned result
	  }
        }

	// set up table entry in r
        this.r[1] = (k - w);
        if (p >= n){
          this.r[0] = 128 + 64;      // out of values--invalid code
	}
        else if (v[p] < s){
          this.r[0] = (this.v[p] < 256 ? 0 : 32 + 64);  // 256 is end-of-block
          this.r[2] = this.v[p++];          // simple code is just the value
        }
        else{
          this.r[0]=(e[this.v[p]-s]+16+64); // non-simple--look up in lists
          this.r[2]=d[this.v[p++] - s];
        }

        // fill code-like entries with r
        f=1<<(k-w);
        for (j=i>>>w;j<z;j+=f){
          arrayCopy(this.r, 0, hp, (q+j)*3, 3);
	}

	// backwards increment the k-bit code i
        for (j = 1 << (k - 1); (i & j)!=0; j >>>= 1){
          i ^= j;
	}
        i ^= j;

	// backup over finished tables
        mask = (1 << w) - 1;      // needed on HP, cc -O bug
        while ((i & mask) != this.x[h]){
          h--;                    // don't need to update q
          w -= l;
          mask = (1 << w) - 1;
        }
      }
    }
    // Return Z_BUF_ERROR if we were given an incomplete table
    return y != 0 && g != 1 ? Z_BUF_ERROR : Z_OK;
}

InfTree.prototype.inflate_trees_bits = function(c, bb, tb, hp, z) {
    var result;
    this.initWorkArea(19);
    this.hn[0]=0;
    result = this.huft_build(c, 0, 19, 19, null, null, tb, bb, hp, this.hn, this.v);

    if(result == Z_DATA_ERROR){
      z.msg = "oversubscribed dynamic bit lengths tree";
    }
    else if(result == Z_BUF_ERROR || bb[0] == 0){
      z.msg = "incomplete dynamic bit lengths tree";
      result = Z_DATA_ERROR;
    }
    return result;
}

InfTree.prototype.inflate_trees_dynamic = function(nl, nd, c, bl, bd, tl, td, hp, z) {
    var result;

    // build literal/length tree
    this.initWorkArea(288);
    this.hn[0]=0;
    result = this.huft_build(c, 0, nl, 257, cplens, cplext, tl, bl, hp, this.hn, this.v);
    if (result != Z_OK || bl[0] == 0){
      if(result == Z_DATA_ERROR){
        z.msg = "oversubscribed literal/length tree";
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "incomplete literal/length tree";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    // build distance tree
    this.initWorkArea(288);
    result = this.huft_build(c, nl, nd, 0, cpdist, cpdext, td, bd, hp, this.hn, this.v);

    if (result != Z_OK || (bd[0] == 0 && nl > 257)){
      if (result == Z_DATA_ERROR){
        z.msg = "oversubscribed distance tree";
      }
      else if (result == Z_BUF_ERROR) {
        z.msg = "incomplete distance tree";
        result = Z_DATA_ERROR;
      }
      else if (result != Z_MEM_ERROR){
        z.msg = "empty distance tree with lengths";
        result = Z_DATA_ERROR;
      }
      return result;
    }

    return Z_OK;
}
/*
  static int inflate_trees_fixed(int[] bl,  //literal desired/actual bit depth
                                 int[] bd,  //distance desired/actual bit depth
                                 int[][] tl,//literal/length tree result
                                 int[][] td,//distance tree result 
                                 ZStream z  //for memory allocation
				 ){

*/

function inflate_trees_fixed(bl, bd, tl, td, z) {
    bl[0]=fixed_bl;
    bd[0]=fixed_bd;
    tl[0]=fixed_tl;
    td[0]=fixed_td;
    return Z_OK;
}

InfTree.prototype.initWorkArea = function(vsize){
    if(this.hn==null){
        this.hn=new Int32Array(1);
        this.v=new Int32Array(vsize);
        this.c=new Int32Array(BMAX+1);
        this.r=new Int32Array(3);
        this.u=new Int32Array(BMAX);
        this.x=new Int32Array(BMAX+1);
    }
    if(this.v.length<vsize){ 
        this.v=new Int32Array(vsize); 
    }
    for(var i=0; i<vsize; i++){this.v[i]=0;}
    for(var i=0; i<BMAX+1; i++){this.c[i]=0;}
    for(var i=0; i<3; i++){this.r[i]=0;}
//  for(int i=0; i<BMAX; i++){u[i]=0;}
    arrayCopy(this.c, 0, this.u, 0, BMAX);
//  for(int i=0; i<BMAX+1; i++){x[i]=0;}
    arrayCopy(this.c, 0, this.x, 0, BMAX+1);
}

var testArray = new Uint8Array(1);
var hasSubarray = (typeof testArray.subarray === 'function');
var hasSlice = false; /* (typeof testArray.slice === 'function'); */ // Chrome slice performance is so dire that we're currently not using it...

function arrayCopy(src, srcOffset, dest, destOffset, count) {
    if (count == 0) {
        return;
    } 
    if (!src) {
        throw "Undef src";
    } else if (!dest) {
        throw "Undef dest";
    }

    if (srcOffset == 0 && count == src.length) {
        arrayCopy_fast(src, dest, destOffset);
    } else if (hasSubarray) {
        arrayCopy_fast(src.subarray(srcOffset, srcOffset + count), dest, destOffset); 
    } else if (src.BYTES_PER_ELEMENT == 1 && count > 100) {
        arrayCopy_fast(new Uint8Array(src.buffer, src.byteOffset + srcOffset, count), dest, destOffset);
    } else { 
        arrayCopy_slow(src, srcOffset, dest, destOffset, count);
    }

}

function arrayCopy_slow(src, srcOffset, dest, destOffset, count) {

    // dlog('_slow call: srcOffset=' + srcOffset + '; destOffset=' + destOffset + '; count=' + count);

     for (var i = 0; i < count; ++i) {
        dest[destOffset + i] = src[srcOffset + i];
    }
}

function arrayCopy_fast(src, dest, destOffset) {
    dest.set(src, destOffset);
}


  // largest prime smaller than 65536
var ADLER_BASE=65521; 
  // NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1
var ADLER_NMAX=5552;

function adler32(adler, /* byte[] */ buf,  index, len){
    if(buf == null){ return 1; }

    var s1=adler&0xffff;
    var s2=(adler>>16)&0xffff;
    var k;

    while(len > 0) {
      k=len<ADLER_NMAX?len:ADLER_NMAX;
      len-=k;
      while(k>=16){
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        s1+=buf[index++]&0xff; s2+=s1;
        k-=16;
      }
      if(k!=0){
        do{
          s1+=buf[index++]&0xff; s2+=s1;
        }
        while(--k!=0);
      }
      s1%=ADLER_BASE;
      s2%=ADLER_BASE;
    }
    return (s2<<16)|s1;
}



function jszlib_inflate_buffer(buffer, start, length, afterUncOffset) {
    if (!start) {
        buffer = new Uint8Array(buffer);
    } else if (!length) {
        buffer = new Uint8Array(buffer, start, buffer.byteLength - start);
    } else {
        buffer = new Uint8Array(buffer, start, length);
    }

    var z = new ZStream();
    z.inflateInit(DEF_WBITS, true);
    z.next_in = buffer;
    z.next_in_index = 0;
    z.avail_in = buffer.length;

    var oBlockList = [];
    var totalSize = 0;
    while (true) {
        var obuf = new Uint8Array(32000);
        z.next_out = obuf;
        z.next_out_index = 0;
        z.avail_out = obuf.length;
        var status = z.inflate(Z_NO_FLUSH);
        if (status != Z_OK && status != Z_STREAM_END && status != Z_BUF_ERROR) {
            throw z.msg;
        }
        if (z.avail_out != 0) {
            var newob = new Uint8Array(obuf.length - z.avail_out);
            arrayCopy(obuf, 0, newob, 0, (obuf.length - z.avail_out));
            obuf = newob;
        }
        oBlockList.push(obuf);
        totalSize += obuf.length;
        if (status == Z_STREAM_END || status == Z_BUF_ERROR) {
            break;
        }
    }

    if (afterUncOffset) {
        afterUncOffset[0] = (start || 0) + z.next_in_index;
    }

    if (oBlockList.length == 1) {
        return oBlockList[0].buffer;
    } else {
        var out = new Uint8Array(totalSize);
        var cursor = 0;
        for (var i = 0; i < oBlockList.length; ++i) {
            var b = oBlockList[i];
            arrayCopy(b, 0, out, cursor, b.length);
            cursor += b.length;
        }
        return out.buffer;
    }
}

if (typeof(module) !== 'undefined') {
  module.exports = {
    inflateBuffer: jszlib_inflate_buffer,
    arrayCopy: arrayCopy
  };
}

},{}]},{},[7])
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9ub2RlX21vZHVsZXMvZ3VscC1icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2UvanMvYmFtLmpzIiwiL1VzZXJzL2pheWhvZGdzb24vZGFsbGlhbmNlL2pzL2JpZ3dpZy5qcyIsIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9qcy9iaW4uanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2UvanMvY29sb3IuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2UvanMvZGFzLmpzIiwiL1VzZXJzL2pheWhvZGdzb24vZGFsbGlhbmNlL2pzL2VuY29kZS5qcyIsIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9qcy9mYWtlX2ViNTk0YWI5LmpzIiwiL1VzZXJzL2pheWhvZGdzb24vZGFsbGlhbmNlL2pzL2xoM3V0aWxzLmpzIiwiL1VzZXJzL2pheWhvZGdzb24vZGFsbGlhbmNlL2pzL3NoYTEuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2UvanMvc3BhbnMuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2UvanMvdXRpbHMuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2Uvbm9kZV9tb2R1bGVzL2VzNi1wcm9taXNlL2Rpc3QvY29tbW9uanMvbWFpbi5qcyIsIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9ub2RlX21vZHVsZXMvZXM2LXByb21pc2UvZGlzdC9jb21tb25qcy9wcm9taXNlL2FsbC5qcyIsIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9ub2RlX21vZHVsZXMvZXM2LXByb21pc2UvZGlzdC9jb21tb25qcy9wcm9taXNlL2FzYXAuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2Uvbm9kZV9tb2R1bGVzL2VzNi1wcm9taXNlL2Rpc3QvY29tbW9uanMvcHJvbWlzZS9jYXN0LmpzIiwiL1VzZXJzL2pheWhvZGdzb24vZGFsbGlhbmNlL25vZGVfbW9kdWxlcy9lczYtcHJvbWlzZS9kaXN0L2NvbW1vbmpzL3Byb21pc2UvY29uZmlnLmpzIiwiL1VzZXJzL2pheWhvZGdzb24vZGFsbGlhbmNlL25vZGVfbW9kdWxlcy9lczYtcHJvbWlzZS9kaXN0L2NvbW1vbmpzL3Byb21pc2UvcG9seWZpbGwuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2Uvbm9kZV9tb2R1bGVzL2VzNi1wcm9taXNlL2Rpc3QvY29tbW9uanMvcHJvbWlzZS9wcm9taXNlLmpzIiwiL1VzZXJzL2pheWhvZGdzb24vZGFsbGlhbmNlL25vZGVfbW9kdWxlcy9lczYtcHJvbWlzZS9kaXN0L2NvbW1vbmpzL3Byb21pc2UvcmFjZS5qcyIsIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9ub2RlX21vZHVsZXMvZXM2LXByb21pc2UvZGlzdC9jb21tb25qcy9wcm9taXNlL3JlamVjdC5qcyIsIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9ub2RlX21vZHVsZXMvZXM2LXByb21pc2UvZGlzdC9jb21tb25qcy9wcm9taXNlL3Jlc29sdmUuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2Uvbm9kZV9tb2R1bGVzL2VzNi1wcm9taXNlL2Rpc3QvY29tbW9uanMvcHJvbWlzZS91dGlscy5qcyIsIi9Vc2Vycy9qYXlob2Rnc29uL2RhbGxpYW5jZS9ub2RlX21vZHVsZXMvZ3VscC1icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9wcm9jZXNzL2Jyb3dzZXIuanMiLCIvVXNlcnMvamF5aG9kZ3Nvbi9kYWxsaWFuY2Uvbm9kZV9tb2R1bGVzL2pzemxpYi9qcy9pbmZsYXRlLmpzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBO0FDQUE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzNoQkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDcGtDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdlRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUhBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUNoMUJBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDN01BO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUM1T0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDNUdBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ25WQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMvUEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDL2VBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDSkE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzVGQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOURBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2xFQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDZEE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN4Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3BOQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3hGQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQzlDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDekNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3JCQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUMvREE7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBIiwiZmlsZSI6ImdlbmVyYXRlZC5qcyIsInNvdXJjZVJvb3QiOiIiLCJzb3VyY2VzQ29udGVudCI6WyIoZnVuY3Rpb24gZSh0LG4scil7ZnVuY3Rpb24gcyhvLHUpe2lmKCFuW29dKXtpZighdFtvXSl7dmFyIGE9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtpZighdSYmYSlyZXR1cm4gYShvLCEwKTtpZihpKXJldHVybiBpKG8sITApO3Rocm93IG5ldyBFcnJvcihcIkNhbm5vdCBmaW5kIG1vZHVsZSAnXCIrbytcIidcIil9dmFyIGY9bltvXT17ZXhwb3J0czp7fX07dFtvXVswXS5jYWxsKGYuZXhwb3J0cyxmdW5jdGlvbihlKXt2YXIgbj10W29dWzFdW2VdO3JldHVybiBzKG4/bjplKX0sZixmLmV4cG9ydHMsZSx0LG4scil9cmV0dXJuIG5bb10uZXhwb3J0c312YXIgaT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2Zvcih2YXIgbz0wO288ci5sZW5ndGg7bysrKXMocltvXSk7cmV0dXJuIHN9KSIsIi8qIC0qLSBtb2RlOiBqYXZhc2NyaXB0OyBjLWJhc2ljLW9mZnNldDogNDsgaW5kZW50LXRhYnMtbW9kZTogbmlsIC0qLSAqL1xuXG4vLyBcbi8vIERhbGxpYW5jZSBHZW5vbWUgRXhwbG9yZXJcbi8vIChjKSBUaG9tYXMgRG93biAyMDA2LTIwMTFcbi8vXG4vLyBiYW0uanM6IGluZGV4ZWQgYmluYXJ5IGFsaWdubWVudHNcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG5pZiAodHlwZW9mKHJlcXVpcmUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIHZhciBzcGFucyA9IHJlcXVpcmUoJy4vc3BhbnMnKTtcbiAgICB2YXIgUmFuZ2UgPSBzcGFucy5SYW5nZTtcbiAgICB2YXIgdW5pb24gPSBzcGFucy51bmlvbjtcbiAgICB2YXIgaW50ZXJzZWN0aW9uID0gc3BhbnMuaW50ZXJzZWN0aW9uO1xuXG4gICAgdmFyIGJpbiA9IHJlcXVpcmUoJy4vYmluJyk7XG4gICAgdmFyIHJlYWRJbnQgPSBiaW4ucmVhZEludDtcbiAgICB2YXIgcmVhZFNob3J0ID0gYmluLnJlYWRTaG9ydDtcbiAgICB2YXIgcmVhZEJ5dGUgPSBiaW4ucmVhZEJ5dGU7XG4gICAgdmFyIHJlYWRJbnQ2NCA9IGJpbi5yZWFkSW50NjQ7XG4gICAgdmFyIHJlYWRGbG9hdCA9IGJpbi5yZWFkRmxvYXQ7XG5cbiAgICB2YXIgbGgzdXRpbHMgPSByZXF1aXJlKCcuL2xoM3V0aWxzJyk7XG4gICAgdmFyIHJlYWRWb2IgPSBsaDN1dGlscy5yZWFkVm9iO1xuICAgIHZhciB1bmJnemYgPSBsaDN1dGlscy51bmJnemY7XG4gICAgdmFyIHJlZzJiaW5zID0gbGgzdXRpbHMucmVnMmJpbnM7XG4gICAgdmFyIENodW5rID0gbGgzdXRpbHMuQ2h1bms7XG59XG5cblxudmFyIEJBTV9NQUdJQyA9IDB4MTRkNDE0MjtcbnZhciBCQUlfTUFHSUMgPSAweDE0OTQxNDI7XG5cbnZhciBCYW1GbGFncyA9IHtcbiAgICBNVUxUSVBMRV9TRUdNRU5UUzogICAgICAgMHgxLFxuICAgIEFMTF9TRUdNRU5UU19BTElHTjogICAgICAweDIsXG4gICAgU0VHTUVOVF9VTk1BUFBFRDogICAgICAgIDB4NCxcbiAgICBORVhUX1NFR01FTlRfVU5NQVBQRUQ6ICAgMHg4LFxuICAgIFJFVkVSU0VfQ09NUExFTUVOVDogICAgICAweDEwLFxuICAgIE5FWFRfUkVWRVJTRV9DT01QTEVNRU5UOiAweDIwLFxuICAgIEZJUlNUX1NFR01FTlQ6ICAgICAgICAgICAweDQwLFxuICAgIExBU1RfU0VHTUVOVDogICAgICAgICAgICAweDgwLFxuICAgIFNFQ09OREFSWV9BTElHTk1FTlQ6ICAgICAweDEwMCxcbiAgICBRQ19GQUlMOiAgICAgICAgICAgICAgICAgMHgyMDAsXG4gICAgRFVQTElDQVRFOiAgICAgICAgICAgICAgIDB4NDAwLFxuICAgIFNVUFBMRU1FTlRBUlk6ICAgICAgICAgICAweDgwMFxufTtcblxuZnVuY3Rpb24gQmFtRmlsZSgpIHtcbn1cblxuXG4vLyBDYWxjdWxhdGUgdGhlIGxlbmd0aCAoaW4gYnl0ZXMpIG9mIHRoZSBCQUkgcmVmIHN0YXJ0aW5nIGF0IG9mZnNldC5cbi8vIFJldHVybnMge25iaW4sIGxlbmd0aCwgbWluQmxvY2tJbmRleH1cbmZ1bmN0aW9uIF9nZXRCYWlSZWZMZW5ndGgodW5jYmEsIG9mZnNldCkge1xuICAgIHZhciBwID0gb2Zmc2V0O1xuICAgIHZhciBuYmluID0gcmVhZEludCh1bmNiYSwgcCk7IHAgKz0gNDtcbiAgICBmb3IgKHZhciBiID0gMDsgYiA8IG5iaW47ICsrYikge1xuICAgICAgICB2YXIgYmluID0gcmVhZEludCh1bmNiYSwgcCk7XG4gICAgICAgIHZhciBuY2huayA9IHJlYWRJbnQodW5jYmEsIHArNCk7XG4gICAgICAgIHAgKz0gOCArIChuY2huayAqIDE2KTtcbiAgICB9XG4gICAgdmFyIG5pbnR2ID0gcmVhZEludCh1bmNiYSwgcCk7IHAgKz0gNDtcblxuICAgIHZhciBtaW5CbG9ja0luZGV4ID0gMTAwMDAwMDAwMDtcbiAgICB2YXIgcSA9IHA7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBuaW50djsgKytpKSB7XG4gICAgICAgIHZhciB2ID0gcmVhZFZvYih1bmNiYSwgcSk7IHEgKz0gODtcbiAgICAgICAgaWYgKHYpIHtcbiAgICAgICAgICAgIHZhciBiaSA9IHYuYmxvY2s7XG4gICAgICAgICAgICBpZiAodi5vZmZzZXQgPiAwKVxuICAgICAgICAgICAgICAgIGJpICs9IDY1NTM2O1xuXG4gICAgICAgICAgICBpZiAoYmkgPCBtaW5CbG9ja0luZGV4KVxuICAgICAgICAgICAgICAgIG1pbkJsb2NrSW5kZXggPSBiaTtcbiAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG4gICAgfVxuICAgIHAgKz0gKG5pbnR2ICogOCk7XG5cbiAgICByZXR1cm4ge1xuICAgICAgICBtaW5CbG9ja0luZGV4OiBtaW5CbG9ja0luZGV4LFxuICAgICAgICBuYmluOiBuYmluLFxuICAgICAgICBsZW5ndGg6IHAgLSBvZmZzZXRcbiAgICB9O1xufVxuXG5cbmZ1bmN0aW9uIG1ha2VCYW0oZGF0YSwgYmFpLCBpbmRleENodW5rcywgY2FsbGJhY2ssIGF0dGVtcHRlZCkge1xuICAgIC8vIERvIGFuIGluaXRpYWwgcHJvYmUgb24gdGhlIEJBTSBmaWxlIHRvIGNhdGNoIGFueSBtaXhlZC1jb250ZW50IGVycm9ycy5cbiAgICBkYXRhLnNsaWNlKDAsIDEwKS5mZXRjaChmdW5jdGlvbihoZWFkZXIpIHtcbiAgICAgICAgaWYgKGhlYWRlcikge1xuICAgICAgICAgICAgcmV0dXJuIG1ha2VCYW0yKGRhdGEsIGJhaSwgaW5kZXhDaHVua3MsIGNhbGxiYWNrLCBhdHRlbXB0ZWQpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsIFwiQ291bGRuJ3QgYWNjZXNzIEJBTS5cIik7XG4gICAgICAgIH1cbiAgICB9LCB7dGltZW91dDogNTAwMH0pO1xufVxuXG5mdW5jdGlvbiBtYWtlQmFtMihkYXRhLCBiYWksIGluZGV4Q2h1bmtzLCBjYWxsYmFjaywgYXR0ZW1wdGVkKSB7XG4gICAgdmFyIGJhbSA9IG5ldyBCYW1GaWxlKCk7XG4gICAgYmFtLmRhdGEgPSBkYXRhO1xuICAgIGJhbS5iYWkgPSBiYWk7XG4gICAgYmFtLmluZGV4Q2h1bmtzID0gaW5kZXhDaHVua3M7XG5cbiAgICB2YXIgbWluQmxvY2tJbmRleCA9IGJhbS5pbmRleENodW5rcyA/IGJhbS5pbmRleENodW5rcy5taW5CbG9ja0luZGV4IDogMTAwMDAwMDAwMDtcblxuICAgIC8vIEZpbGxzIG91dCBiYW0uY2hyVG9JbmRleCBhbmQgYmFtLmluZGV4VG9DaHIgYmFzZWQgb24gdGhlIGZpcnN0IGZldyBieXRlcyBvZiB0aGUgQkFNLlxuICAgIGZ1bmN0aW9uIHBhcnNlQmFtSGVhZGVyKHIpIHtcbiAgICAgICAgaWYgKCFyKSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgXCJDb3VsZG4ndCBhY2Nlc3MgQkFNXCIpO1xuICAgICAgICB9XG5cbiAgICAgICAgdmFyIHVuYyA9IHVuYmd6ZihyLCByLmJ5dGVMZW5ndGgpO1xuICAgICAgICB2YXIgdW5jYmEgPSBuZXcgVWludDhBcnJheSh1bmMpO1xuXG4gICAgICAgIHZhciBtYWdpYyA9IHJlYWRJbnQodW5jYmEsIDApO1xuICAgICAgICBpZiAobWFnaWMgIT0gQkFNX01BR0lDKSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgXCJOb3QgYSBCQU0gZmlsZSwgbWFnaWM9MHhcIiArIG1hZ2ljLnRvU3RyaW5nKDE2KSk7XG4gICAgICAgIH1cbiAgICAgICAgdmFyIGhlYWRMZW4gPSByZWFkSW50KHVuY2JhLCA0KTtcbiAgICAgICAgdmFyIGhlYWRlciA9ICcnO1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGhlYWRMZW47ICsraSkge1xuICAgICAgICAgICAgaGVhZGVyICs9IFN0cmluZy5mcm9tQ2hhckNvZGUodW5jYmFbaSArIDhdKTtcbiAgICAgICAgfVxuXG4gICAgICAgIHZhciBuUmVmID0gcmVhZEludCh1bmNiYSwgaGVhZExlbiArIDgpO1xuICAgICAgICB2YXIgcCA9IGhlYWRMZW4gKyAxMjtcblxuICAgICAgICBiYW0uY2hyVG9JbmRleCA9IHt9O1xuICAgICAgICBiYW0uaW5kZXhUb0NociA9IFtdO1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG5SZWY7ICsraSkge1xuICAgICAgICAgICAgdmFyIGxOYW1lID0gcmVhZEludCh1bmNiYSwgcCk7XG4gICAgICAgICAgICB2YXIgbmFtZSA9ICcnO1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCBsTmFtZS0xOyArK2opIHtcbiAgICAgICAgICAgICAgICBuYW1lICs9IFN0cmluZy5mcm9tQ2hhckNvZGUodW5jYmFbcCArIDQgKyBqXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICB2YXIgbFJlZiA9IHJlYWRJbnQodW5jYmEsIHAgKyBsTmFtZSArIDQpO1xuICAgICAgICAgICAgYmFtLmNoclRvSW5kZXhbbmFtZV0gPSBpO1xuICAgICAgICAgICAgaWYgKG5hbWUuaW5kZXhPZignY2hyJykgPT0gMCkge1xuICAgICAgICAgICAgICAgIGJhbS5jaHJUb0luZGV4W25hbWUuc3Vic3RyaW5nKDMpXSA9IGk7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIGJhbS5jaHJUb0luZGV4WydjaHInICsgbmFtZV0gPSBpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgYmFtLmluZGV4VG9DaHIucHVzaChuYW1lKTtcblxuICAgICAgICAgICAgcCA9IHAgKyA4ICsgbE5hbWU7XG4gICAgICAgIH1cblxuICAgICAgICBpZiAoYmFtLmluZGljZXMpIHtcbiAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhiYW0pO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgZnVuY3Rpb24gcGFyc2VCYWkoaGVhZGVyKSB7XG4gICAgICAgIGlmICghaGVhZGVyKSB7XG4gICAgICAgICAgICByZXR1cm4gXCJDb3VsZG4ndCBhY2Nlc3MgQkFJXCI7XG4gICAgICAgIH1cblxuICAgICAgICB2YXIgdW5jYmEgPSBuZXcgVWludDhBcnJheShoZWFkZXIpO1xuICAgICAgICB2YXIgYmFpTWFnaWMgPSByZWFkSW50KHVuY2JhLCAwKTtcbiAgICAgICAgaWYgKGJhaU1hZ2ljICE9IEJBSV9NQUdJQykge1xuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsICdOb3QgYSBCQUkgZmlsZSwgbWFnaWM9MHgnICsgYmFpTWFnaWMudG9TdHJpbmcoMTYpKTtcbiAgICAgICAgfVxuXG4gICAgICAgIHZhciBucmVmID0gcmVhZEludCh1bmNiYSwgNCk7XG5cbiAgICAgICAgYmFtLmluZGljZXMgPSBbXTtcblxuICAgICAgICB2YXIgcCA9IDg7XG4gICAgICAgIGZvciAodmFyIHJlZiA9IDA7IHJlZiA8IG5yZWY7ICsrcmVmKSB7XG4gICAgICAgICAgICB2YXIgYmxvY2tTdGFydCA9IHA7XG4gICAgICAgICAgICB2YXIgbyA9IF9nZXRCYWlSZWZMZW5ndGgodW5jYmEsIGJsb2NrU3RhcnQpO1xuICAgICAgICAgICAgcCArPSBvLmxlbmd0aDtcblxuICAgICAgICAgICAgbWluQmxvY2tJbmRleCA9IE1hdGgubWluKG8ubWluQmxvY2tJbmRleCwgbWluQmxvY2tJbmRleCk7XG5cbiAgICAgICAgICAgIHZhciBuYmluID0gby5uYmluO1xuXG4gICAgICAgICAgICBpZiAobmJpbiA+IDApIHtcbiAgICAgICAgICAgICAgICBiYW0uaW5kaWNlc1tyZWZdID0gbmV3IFVpbnQ4QXJyYXkoaGVhZGVyLCBibG9ja1N0YXJ0LCBwIC0gYmxvY2tTdGFydCk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICByZXR1cm4gdHJ1ZTtcbiAgICB9XG5cbiAgICBpZiAoIWJhbS5pbmRleENodW5rcykge1xuICAgICAgICBiYW0uYmFpLmZldGNoKGZ1bmN0aW9uKGhlYWRlcikgeyAgIC8vIERvIHdlIHJlYWxseSBuZWVkIHRvIGZldGNoIHRoZSB3aG9sZSB0aGluZz8gOi0oXG4gICAgICAgICAgICB2YXIgcmVzdWx0ID0gcGFyc2VCYWkoaGVhZGVyKTtcbiAgICAgICAgICAgIGlmIChyZXN1bHQgIT09IHRydWUpIHtcbiAgICAgICAgICAgICAgICBpZiAoYmFtLmJhaS51cmwgJiYgdHlwZW9mKGF0dGVtcHRlZCkgPT09IFwidW5kZWZpbmVkXCIpIHtcbiAgICAgICAgICAgICAgICAgICAgLy8gQWxyZWFkeSBhdHRlbXB0ZWQgeC5iYW0uYmFpIG5vdCB0aGVyZSBzbyBub3cgdHJ5aW5nIHguYmFpXG4gICAgICAgICAgICAgICAgICAgIGJhbS5iYWkudXJsID0gYmFtLmRhdGEudXJsLnJlcGxhY2UobmV3IFJlZ0V4cCgnLmJhbSQnKSwgJy5iYWknKTtcbiAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgICAvLyBUcnVlIGxldHMgdXMga25vdyB3ZSBhcmUgbWFraW5nIGEgc2Vjb25kIGF0dGVtcHRcbiAgICAgICAgICAgICAgICAgICAgbWFrZUJhbTIoZGF0YSwgYmFtLmJhaSwgaW5kZXhDaHVua3MsIGNhbGxiYWNrLCB0cnVlKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIC8vIFdlJ3ZlIGF0dGVtcHRlZCB4LmJhbS5iYWkgJiB4LmJhaSBhbmQgbm90aGluZyB3b3JrZWRcbiAgICAgICAgICAgICAgICAgICAgY2FsbGJhY2sobnVsbCwgcmVzdWx0KTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICBiYW0uZGF0YS5zbGljZSgwLCBtaW5CbG9ja0luZGV4KS5mZXRjaChwYXJzZUJhbUhlYWRlcik7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0pOyAgIC8vIFRpbWVvdXQgb24gZmlyc3QgcmVxdWVzdCB0byBjYXRjaCBDaHJvbWUgbWl4ZWQtY29udGVudCBlcnJvci5cbiAgICB9IGVsc2Uge1xuICAgICAgICB2YXIgY2h1bmtzID0gYmFtLmluZGV4Q2h1bmtzLmNodW5rcztcbiAgICAgICAgYmFtLmluZGljZXMgPSBbXVxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGNodW5rcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICBiYW0uaW5kaWNlc1tpXSA9IG51bGw7ICAvLyBUbyBiZSBmaWxsZWQgb3V0IGxhemlseSBhcyBuZWVkZWRcbiAgICAgICAgfVxuICAgICAgICBiYW0uZGF0YS5zbGljZSgwLCBtaW5CbG9ja0luZGV4KS5mZXRjaChwYXJzZUJhbUhlYWRlcik7XG4gICAgfVxufVxuXG5cblxuQmFtRmlsZS5wcm90b3R5cGUuYmxvY2tzRm9yUmFuZ2UgPSBmdW5jdGlvbihyZWZJZCwgbWluLCBtYXgpIHtcbiAgICB2YXIgaW5kZXggPSB0aGlzLmluZGljZXNbcmVmSWRdO1xuICAgIGlmICghaW5kZXgpIHtcbiAgICAgICAgcmV0dXJuIFtdO1xuICAgIH1cblxuICAgIHZhciBpbnRCaW5zTCA9IHJlZzJiaW5zKG1pbiwgbWF4KTtcbiAgICB2YXIgaW50QmlucyA9IFtdO1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgaW50Qmluc0wubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgaW50Qmluc1tpbnRCaW5zTFtpXV0gPSB0cnVlO1xuICAgIH1cbiAgICB2YXIgbGVhZkNodW5rcyA9IFtdLCBvdGhlckNodW5rcyA9IFtdO1xuXG4gICAgdmFyIG5iaW4gPSByZWFkSW50KGluZGV4LCAwKTtcbiAgICB2YXIgcCA9IDQ7XG4gICAgZm9yICh2YXIgYiA9IDA7IGIgPCBuYmluOyArK2IpIHtcbiAgICAgICAgdmFyIGJpbiA9IHJlYWRJbnQoaW5kZXgsIHApO1xuICAgICAgICB2YXIgbmNobmsgPSByZWFkSW50KGluZGV4LCBwKzQpO1xuLy8gICAgICAgIGRsb2coJ2Jpbj0nICsgYmluICsgJzsgbmNobms9JyArIG5jaG5rKTtcbiAgICAgICAgcCArPSA4O1xuICAgICAgICBpZiAoaW50Qmluc1tiaW5dKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBjID0gMDsgYyA8IG5jaG5rOyArK2MpIHtcbiAgICAgICAgICAgICAgICB2YXIgY3MgPSByZWFkVm9iKGluZGV4LCBwKTtcbiAgICAgICAgICAgICAgICB2YXIgY2UgPSByZWFkVm9iKGluZGV4LCBwICsgOCk7XG4gICAgICAgICAgICAgICAgKGJpbiA8IDQ2ODEgPyBvdGhlckNodW5rcyA6IGxlYWZDaHVua3MpLnB1c2gobmV3IENodW5rKGNzLCBjZSkpO1xuICAgICAgICAgICAgICAgIHAgKz0gMTY7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBwICs9ICAobmNobmsgKiAxNik7XG4gICAgICAgIH1cbiAgICB9XG4gICAgLy8gY29uc29sZS5sb2coJ2xlYWZDaHVua3MgPSAnICsgbWluaUpTT05pZnkobGVhZkNodW5rcykpO1xuICAgIC8vIGNvbnNvbGUubG9nKCdvdGhlckNodW5rcyA9ICcgKyBtaW5pSlNPTmlmeShvdGhlckNodW5rcykpO1xuXG4gICAgdmFyIG5pbnR2ID0gcmVhZEludChpbmRleCwgcCk7XG4gICAgdmFyIGxvd2VzdCA9IG51bGw7XG4gICAgdmFyIG1pbkxpbiA9IE1hdGgubWluKG1pbj4+MTQsIG5pbnR2IC0gMSksIG1heExpbiA9IE1hdGgubWluKG1heD4+MTQsIG5pbnR2IC0gMSk7XG4gICAgZm9yICh2YXIgaSA9IG1pbkxpbjsgaSA8PSBtYXhMaW47ICsraSkge1xuICAgICAgICB2YXIgbGIgPSAgcmVhZFZvYihpbmRleCwgcCArIDQgKyAoaSAqIDgpKTtcbiAgICAgICAgaWYgKCFsYikge1xuICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKCFsb3dlc3QgfHwgbGIuYmxvY2sgPCBsb3dlc3QuYmxvY2sgfHwgbGIub2Zmc2V0IDwgbG93ZXN0Lm9mZnNldCkge1xuICAgICAgICAgICAgbG93ZXN0ID0gbGI7XG4gICAgICAgIH1cbiAgICB9XG4gICAgLy8gY29uc29sZS5sb2coJ0xvd2VzdCBMQiA9ICcgKyBsb3dlc3QpO1xuICAgIFxuICAgIHZhciBwcnVuZWRPdGhlckNodW5rcyA9IFtdO1xuICAgIGlmIChsb3dlc3QgIT0gbnVsbCkge1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG90aGVyQ2h1bmtzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICB2YXIgY2huayA9IG90aGVyQ2h1bmtzW2ldO1xuICAgICAgICAgICAgaWYgKGNobmsubWF4di5ibG9jayA+PSBsb3dlc3QuYmxvY2sgJiYgY2huay5tYXh2Lm9mZnNldCA+PSBsb3dlc3Qub2Zmc2V0KSB7XG4gICAgICAgICAgICAgICAgcHJ1bmVkT3RoZXJDaHVua3MucHVzaChjaG5rKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICAvLyBjb25zb2xlLmxvZygncHJ1bmVkT3RoZXJDaHVua3MgPSAnICsgbWluaUpTT05pZnkocHJ1bmVkT3RoZXJDaHVua3MpKTtcbiAgICBvdGhlckNodW5rcyA9IHBydW5lZE90aGVyQ2h1bmtzO1xuXG4gICAgdmFyIGludENodW5rcyA9IFtdO1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgb3RoZXJDaHVua3MubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgaW50Q2h1bmtzLnB1c2gob3RoZXJDaHVua3NbaV0pO1xuICAgIH1cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGxlYWZDaHVua3MubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgaW50Q2h1bmtzLnB1c2gobGVhZkNodW5rc1tpXSk7XG4gICAgfVxuXG4gICAgaW50Q2h1bmtzLnNvcnQoZnVuY3Rpb24oYzAsIGMxKSB7XG4gICAgICAgIHZhciBkaWYgPSBjMC5taW52LmJsb2NrIC0gYzEubWludi5ibG9jaztcbiAgICAgICAgaWYgKGRpZiAhPSAwKSB7XG4gICAgICAgICAgICByZXR1cm4gZGlmO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgcmV0dXJuIGMwLm1pbnYub2Zmc2V0IC0gYzEubWludi5vZmZzZXQ7XG4gICAgICAgIH1cbiAgICB9KTtcbiAgICB2YXIgbWVyZ2VkQ2h1bmtzID0gW107XG4gICAgaWYgKGludENodW5rcy5sZW5ndGggPiAwKSB7XG4gICAgICAgIHZhciBjdXIgPSBpbnRDaHVua3NbMF07XG4gICAgICAgIGZvciAodmFyIGkgPSAxOyBpIDwgaW50Q2h1bmtzLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICB2YXIgbmMgPSBpbnRDaHVua3NbaV07XG4gICAgICAgICAgICBpZiAobmMubWludi5ibG9jayA9PSBjdXIubWF4di5ibG9jayAvKiAmJiBuYy5taW52Lm9mZnNldCA9PSBjdXIubWF4di5vZmZzZXQgKi8pIHsgLy8gbm8gcG9pbnQgc3BsaXR0aW5nIG1pZC1ibG9ja1xuICAgICAgICAgICAgICAgIGN1ciA9IG5ldyBDaHVuayhjdXIubWludiwgbmMubWF4dik7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIG1lcmdlZENodW5rcy5wdXNoKGN1cik7XG4gICAgICAgICAgICAgICAgY3VyID0gbmM7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgbWVyZ2VkQ2h1bmtzLnB1c2goY3VyKTtcbiAgICB9XG4gICAgLy8gY29uc29sZS5sb2coJ21lcmdlZENodW5rcyA9ICcgKyBtaW5pSlNPTmlmeShtZXJnZWRDaHVua3MpKTtcblxuICAgIHJldHVybiBtZXJnZWRDaHVua3M7XG59XG5cbkJhbUZpbGUucHJvdG90eXBlLmZldGNoID0gZnVuY3Rpb24oY2hyLCBtaW4sIG1heCwgY2FsbGJhY2ssIG9wdHMpIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuICAgIG9wdHMgPSBvcHRzIHx8IHt9O1xuXG4gICAgdmFyIGNocklkID0gdGhpcy5jaHJUb0luZGV4W2Nocl07XG4gICAgdmFyIGNodW5rcztcbiAgICBpZiAoY2hySWQgPT09IHVuZGVmaW5lZCkge1xuICAgICAgICBjaHVua3MgPSBbXTtcbiAgICB9IGVsc2Uge1xuICAgICAgICAvLyBGZXRjaCB0aGlzIHBvcnRpb24gb2YgdGhlIEJBSSBpZiBpdCBoYXNuJ3QgYmVlbiBsb2FkZWQgeWV0LlxuICAgICAgICBpZiAodGhpcy5pbmRpY2VzW2NocklkXSA9PT0gbnVsbCAmJiB0aGlzLmluZGV4Q2h1bmtzLmNodW5rc1tjaHJJZF0pIHtcbiAgICAgICAgICAgIHZhciBzdGFydF9zdG9wID0gdGhpcy5pbmRleENodW5rcy5jaHVua3NbY2hySWRdO1xuICAgICAgICAgICAgcmV0dXJuIHRoaXMuYmFpLnNsaWNlKHN0YXJ0X3N0b3BbMF0sIHN0YXJ0X3N0b3BbMV0pLmZldGNoKGZ1bmN0aW9uKGRhdGEpIHtcbiAgICAgICAgICAgICAgICB2YXIgYnVmZmVyID0gbmV3IFVpbnQ4QXJyYXkoZGF0YSk7XG4gICAgICAgICAgICAgICAgdGhpcy5pbmRpY2VzW2NocklkXSA9IGJ1ZmZlcjtcbiAgICAgICAgICAgICAgICByZXR1cm4gdGhpcy5mZXRjaChjaHIsIG1pbiwgbWF4LCBjYWxsYmFjaywgb3B0cyk7XG4gICAgICAgICAgICB9LmJpbmQodGhpcykpO1xuICAgICAgICB9XG5cbiAgICAgICAgY2h1bmtzID0gdGhpcy5ibG9ja3NGb3JSYW5nZShjaHJJZCwgbWluLCBtYXgpO1xuICAgICAgICBpZiAoIWNodW5rcykge1xuICAgICAgICAgICAgY2FsbGJhY2sobnVsbCwgJ0Vycm9yIGluIGluZGV4IGZldGNoJyk7XG4gICAgICAgIH1cbiAgICB9XG4gICAgXG4gICAgdmFyIHJlY29yZHMgPSBbXTtcbiAgICB2YXIgaW5kZXggPSAwO1xuICAgIHZhciBkYXRhO1xuXG4gICAgZnVuY3Rpb24gdHJhbXAoKSB7XG4gICAgICAgIGlmIChpbmRleCA+PSBjaHVua3MubGVuZ3RoKSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2socmVjb3Jkcyk7XG4gICAgICAgIH0gZWxzZSBpZiAoIWRhdGEpIHtcbiAgICAgICAgICAgIHZhciBjID0gY2h1bmtzW2luZGV4XTtcbiAgICAgICAgICAgIHZhciBmZXRjaE1pbiA9IGMubWludi5ibG9jaztcbiAgICAgICAgICAgIHZhciBmZXRjaE1heCA9IGMubWF4di5ibG9jayArICgxPDwxNik7IC8vICpzaWdoKlxuICAgICAgICAgICAgLy8gY29uc29sZS5sb2coJ2ZldGNoaW5nICcgKyBmZXRjaE1pbiArICc6JyArIGZldGNoTWF4KTtcbiAgICAgICAgICAgIHRoaXNCLmRhdGEuc2xpY2UoZmV0Y2hNaW4sIGZldGNoTWF4IC0gZmV0Y2hNaW4pLmZldGNoKGZ1bmN0aW9uKHIpIHtcbiAgICAgICAgICAgICAgICBkYXRhID0gdW5iZ3pmKHIsIGMubWF4di5ibG9jayAtIGMubWludi5ibG9jayArIDEpO1xuICAgICAgICAgICAgICAgIHJldHVybiB0cmFtcCgpO1xuICAgICAgICAgICAgfSk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShkYXRhKTtcbiAgICAgICAgICAgIHZhciBmaW5pc2hlZCA9IHRoaXNCLnJlYWRCYW1SZWNvcmRzKGJhLCBjaHVua3NbaW5kZXhdLm1pbnYub2Zmc2V0LCByZWNvcmRzLCBtaW4sIG1heCwgY2hySWQsIG9wdHMpO1xuICAgICAgICAgICAgZGF0YSA9IG51bGw7XG4gICAgICAgICAgICArK2luZGV4O1xuICAgICAgICAgICAgaWYgKGZpbmlzaGVkKVxuICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhyZWNvcmRzKTtcbiAgICAgICAgICAgIGVsc2VcbiAgICAgICAgICAgICAgICByZXR1cm4gdHJhbXAoKTtcbiAgICAgICAgfVxuICAgIH1cbiAgICB0cmFtcCgpO1xufVxuXG52YXIgU0VRUkVUX0RFQ09ERVIgPSBbJz0nLCAnQScsICdDJywgJ3gnLCAnRycsICd4JywgJ3gnLCAneCcsICdUJywgJ3gnLCAneCcsICd4JywgJ3gnLCAneCcsICd4JywgJ04nXTtcbnZhciBDSUdBUl9ERUNPREVSID0gWydNJywgJ0knLCAnRCcsICdOJywgJ1MnLCAnSCcsICdQJywgJz0nLCAnWCcsICc/JywgJz8nLCAnPycsICc/JywgJz8nLCAnPycsICc/J107XG5cbmZ1bmN0aW9uIEJhbVJlY29yZCgpIHtcbn1cblxuQmFtRmlsZS5wcm90b3R5cGUucmVhZEJhbVJlY29yZHMgPSBmdW5jdGlvbihiYSwgb2Zmc2V0LCBzaW5rLCBtaW4sIG1heCwgY2hySWQsIG9wdHMpIHtcbiAgICB3aGlsZSAodHJ1ZSkge1xuICAgICAgICB2YXIgYmxvY2tTaXplID0gcmVhZEludChiYSwgb2Zmc2V0KTtcbiAgICAgICAgdmFyIGJsb2NrRW5kID0gb2Zmc2V0ICsgYmxvY2tTaXplICsgNDtcbiAgICAgICAgaWYgKGJsb2NrRW5kID49IGJhLmxlbmd0aCkge1xuICAgICAgICAgICAgcmV0dXJuIGZhbHNlO1xuICAgICAgICB9XG5cbiAgICAgICAgdmFyIHJlY29yZCA9IG5ldyBCYW1SZWNvcmQoKTtcblxuICAgICAgICB2YXIgcmVmSUQgPSByZWFkSW50KGJhLCBvZmZzZXQgKyA0KTtcbiAgICAgICAgdmFyIHBvcyA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDgpO1xuICAgICAgICBcbiAgICAgICAgdmFyIGJtbiA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDEyKTtcbiAgICAgICAgdmFyIGJpbiA9IChibW4gJiAweGZmZmYwMDAwKSA+PiAxNjtcbiAgICAgICAgdmFyIG1xID0gKGJtbiAmIDB4ZmYwMCkgPj4gODtcbiAgICAgICAgdmFyIG5sID0gYm1uICYgMHhmZjtcblxuICAgICAgICB2YXIgZmxhZ19uYyA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDE2KTtcbiAgICAgICAgdmFyIGZsYWcgPSAoZmxhZ19uYyAmIDB4ZmZmZjAwMDApID4+IDE2O1xuICAgICAgICB2YXIgbmMgPSBmbGFnX25jICYgMHhmZmZmO1xuICAgIFxuICAgICAgICB2YXIgbHNlcSA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDIwKTtcbiAgICAgICAgXG4gICAgICAgIHZhciBuZXh0UmVmICA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDI0KTtcbiAgICAgICAgdmFyIG5leHRQb3MgPSByZWFkSW50KGJhLCBvZmZzZXQgKyAyOCk7XG4gICAgICAgIFxuICAgICAgICB2YXIgdGxlbiA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDMyKTtcbiAgICBcbiAgICAgICAgcmVjb3JkLnNlZ21lbnQgPSB0aGlzLmluZGV4VG9DaHJbcmVmSURdO1xuICAgICAgICByZWNvcmQuZmxhZyA9IGZsYWc7XG4gICAgICAgIHJlY29yZC5wb3MgPSBwb3M7XG4gICAgICAgIHJlY29yZC5tcSA9IG1xO1xuICAgICAgICBpZiAob3B0cy5saWdodClcbiAgICAgICAgICAgIHJlY29yZC5zZXFMZW5ndGggPSBsc2VxO1xuXG4gICAgICAgIGlmICghb3B0cy5saWdodCkge1xuICAgICAgICAgICAgaWYgKG5leHRSZWYgPj0gMCkge1xuICAgICAgICAgICAgICAgIHJlY29yZC5uZXh0U2VnbWVudCA9IHRoaXMuaW5kZXhUb0NocltuZXh0UmVmXTtcbiAgICAgICAgICAgICAgICByZWNvcmQubmV4dFBvcyA9IG5leHRQb3M7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHZhciByZWFkTmFtZSA9ICcnO1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCBubC0xOyArK2opIHtcbiAgICAgICAgICAgICAgICByZWFkTmFtZSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW29mZnNldCArIDM2ICsgal0pO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmVjb3JkLnJlYWROYW1lID0gcmVhZE5hbWU7XG4gICAgICAgIFxuICAgICAgICAgICAgdmFyIHAgPSBvZmZzZXQgKyAzNiArIG5sO1xuXG4gICAgICAgICAgICB2YXIgY2lnYXIgPSAnJztcbiAgICAgICAgICAgIGZvciAodmFyIGMgPSAwOyBjIDwgbmM7ICsrYykge1xuICAgICAgICAgICAgICAgIHZhciBjaWdvcCA9IHJlYWRJbnQoYmEsIHApO1xuICAgICAgICAgICAgICAgIGNpZ2FyID0gY2lnYXIgKyAoY2lnb3A+PjQpICsgQ0lHQVJfREVDT0RFUltjaWdvcCAmIDB4Zl07XG4gICAgICAgICAgICAgICAgcCArPSA0O1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmVjb3JkLmNpZ2FyID0gY2lnYXI7XG4gICAgICAgIFxuICAgICAgICAgICAgdmFyIHNlcSA9ICcnO1xuICAgICAgICAgICAgdmFyIHNlcUJ5dGVzID0gKGxzZXEgKyAxKSA+PiAxO1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCBzZXFCeXRlczsgKytqKSB7XG4gICAgICAgICAgICAgICAgdmFyIHNiID0gYmFbcCArIGpdO1xuICAgICAgICAgICAgICAgIHNlcSArPSBTRVFSRVRfREVDT0RFUlsoc2IgJiAweGYwKSA+PiA0XTtcbiAgICAgICAgICAgICAgICBpZiAoc2VxLmxlbmd0aCA8IGxzZXEpXG4gICAgICAgICAgICAgICAgICAgIHNlcSArPSBTRVFSRVRfREVDT0RFUlsoc2IgJiAweDBmKV07XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBwICs9IHNlcUJ5dGVzO1xuICAgICAgICAgICAgcmVjb3JkLnNlcSA9IHNlcTtcblxuICAgICAgICAgICAgdmFyIHFzZXEgPSAnJztcbiAgICAgICAgICAgIGZvciAodmFyIGogPSAwOyBqIDwgbHNlcTsgKytqKSB7XG4gICAgICAgICAgICAgICAgcXNlcSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW3AgKyBqXSArIDMzKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHAgKz0gbHNlcTtcbiAgICAgICAgICAgIHJlY29yZC5xdWFscyA9IHFzZXE7XG5cbiAgICAgICAgICAgIHdoaWxlIChwIDwgYmxvY2tFbmQpIHtcbiAgICAgICAgICAgICAgICB2YXIgdGFnID0gU3RyaW5nLmZyb21DaGFyQ29kZShiYVtwXSwgYmFbcCArIDFdKTtcbiAgICAgICAgICAgICAgICB2YXIgdHlwZSA9IFN0cmluZy5mcm9tQ2hhckNvZGUoYmFbcCArIDJdKTtcbiAgICAgICAgICAgICAgICB2YXIgdmFsdWU7XG5cbiAgICAgICAgICAgICAgICBpZiAodHlwZSA9PSAnQScpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFsdWUgPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW3AgKyAzXSk7XG4gICAgICAgICAgICAgICAgICAgIHAgKz0gNDtcbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT0gJ2knIHx8IHR5cGUgPT0gJ0knKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhbHVlID0gcmVhZEludChiYSwgcCArIDMpO1xuICAgICAgICAgICAgICAgICAgICBwICs9IDc7XG4gICAgICAgICAgICAgICAgfSBlbHNlIGlmICh0eXBlID09ICdjJyB8fCB0eXBlID09ICdDJykge1xuICAgICAgICAgICAgICAgICAgICB2YWx1ZSA9IGJhW3AgKyAzXTtcbiAgICAgICAgICAgICAgICAgICAgcCArPSA0O1xuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PSAncycgfHwgdHlwZSA9PSAnUycpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFsdWUgPSByZWFkU2hvcnQoYmEsIHAgKyAzKTtcbiAgICAgICAgICAgICAgICAgICAgcCArPSA1O1xuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PSAnZicpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFsdWUgPSByZWFkRmxvYXQoYmEsIHAgKyAzKTtcbiAgICAgICAgICAgICAgICAgICAgcCArPSA3O1xuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAodHlwZSA9PSAnWicgfHwgdHlwZSA9PSAnSCcpIHtcbiAgICAgICAgICAgICAgICAgICAgcCArPSAzO1xuICAgICAgICAgICAgICAgICAgICB2YWx1ZSA9ICcnO1xuICAgICAgICAgICAgICAgICAgICBmb3IgKDs7KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgY2MgPSBiYVtwKytdO1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGNjID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFsdWUgKz0gU3RyaW5nLmZyb21DaGFyQ29kZShjYyk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKHR5cGUgPT0gJ0InKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBhdHlwZSA9IFN0cmluZy5mcm9tQ2hhckNvZGUoYmFbcCArIDNdKTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGFsZW4gPSByZWFkSW50KGJhLCBwICsgNCk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBlbGVuO1xuICAgICAgICAgICAgICAgICAgICB2YXIgcmVhZGVyO1xuICAgICAgICAgICAgICAgICAgICBpZiAoYXR5cGUgPT0gJ2knIHx8IGF0eXBlID09ICdJJyB8fCBhdHlwZSA9PSAnZicpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGVsZW4gPSA0O1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGF0eXBlID09ICdmJylcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICByZWFkZXIgPSByZWFkRmxvYXQ7XG4gICAgICAgICAgICAgICAgICAgICAgICBlbHNlXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgcmVhZGVyID0gcmVhZEludDtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIGlmIChhdHlwZSA9PSAncycgfHwgYXR5cGUgPT0gJ1MnKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBlbGVuID0gMjtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJlYWRlciA9IHJlYWRTaG9ydDtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIGlmIChhdHlwZSA9PSAnYycgfHwgYXR5cGUgPT0gJ0MnKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBlbGVuID0gMTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJlYWRlciA9IHJlYWRCeXRlO1xuICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgdGhyb3cgJ1Vua25vd24gYXJyYXkgdHlwZSAnICsgYXR5cGU7XG4gICAgICAgICAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgICAgICAgICBwICs9IDg7XG4gICAgICAgICAgICAgICAgICAgIHZhbHVlID0gW107XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgYWxlbjsgKytpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YWx1ZS5wdXNoKHJlYWRlcihiYSwgcCkpO1xuICAgICAgICAgICAgICAgICAgICAgICAgcCArPSBlbGVuO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgdGhyb3cgJ1Vua25vd24gdHlwZSAnKyB0eXBlO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICByZWNvcmRbdGFnXSA9IHZhbHVlO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgaWYgKCFtaW4gfHwgcmVjb3JkLnBvcyA8PSBtYXggJiYgcmVjb3JkLnBvcyArIGxzZXEgPj0gbWluKSB7XG4gICAgICAgICAgICBpZiAoY2hySWQgPT09IHVuZGVmaW5lZCB8fCByZWZJRCA9PSBjaHJJZCkge1xuICAgICAgICAgICAgICAgIHNpbmsucHVzaChyZWNvcmQpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGlmIChyZWNvcmQucG9zID4gbWF4KSB7XG4gICAgICAgICAgICByZXR1cm4gdHJ1ZTtcbiAgICAgICAgfVxuICAgICAgICBvZmZzZXQgPSBibG9ja0VuZDtcbiAgICB9XG5cbiAgICAvLyBFeGl0cyB2aWEgdG9wIG9mIGxvb3AuXG59O1xuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIG1ha2VCYW06IG1ha2VCYW0sXG4gICAgICAgIEJBTV9NQUdJQzogQkFNX01BR0lDLFxuICAgICAgICBCQUlfTUFHSUM6IEJBSV9NQUdJQyxcbiAgICAgICAgQmFtRmxhZ3M6IEJhbUZsYWdzXG4gICAgfTtcbn1cbiIsIi8qIC0qLSBtb2RlOiBqYXZhc2NyaXB0OyBjLWJhc2ljLW9mZnNldDogNDsgaW5kZW50LXRhYnMtbW9kZTogbmlsIC0qLSAqL1xuXG4vLyBcbi8vIERhbGxpYW5jZSBHZW5vbWUgRXhwbG9yZXJcbi8vIChjKSBUaG9tYXMgRG93biAyMDA2LTIwMTBcbi8vXG4vLyBiaWd3aWcuanM6IGluZGV4ZWQgYmluYXJ5IFdJRyAoYW5kIEJFRCkgZmlsZXNcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG5cbmlmICh0eXBlb2YocmVxdWlyZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgdmFyIHNwYW5zID0gcmVxdWlyZSgnLi9zcGFucycpO1xuICAgIHZhciBSYW5nZSA9IHNwYW5zLlJhbmdlO1xuICAgIHZhciB1bmlvbiA9IHNwYW5zLnVuaW9uO1xuICAgIHZhciBpbnRlcnNlY3Rpb24gPSBzcGFucy5pbnRlcnNlY3Rpb247XG5cbiAgICB2YXIgZGFzID0gcmVxdWlyZSgnLi9kYXMnKTtcbiAgICB2YXIgREFTRmVhdHVyZSA9IGRhcy5EQVNGZWF0dXJlO1xuICAgIHZhciBEQVNHcm91cCA9IGRhcy5EQVNHcm91cDtcblxuICAgIHZhciB1dGlscyA9IHJlcXVpcmUoJy4vdXRpbHMnKTtcbiAgICB2YXIgc2hhbGxvd0NvcHkgPSB1dGlscy5zaGFsbG93Q29weTtcblxuICAgIHZhciBiaW4gPSByZXF1aXJlKCcuL2JpbicpO1xuICAgIHZhciByZWFkSW50ID0gYmluLnJlYWRJbnQ7XG5cbiAgICB2YXIganN6bGliID0gcmVxdWlyZSgnanN6bGliJyk7XG4gICAgdmFyIGpzemxpYl9pbmZsYXRlX2J1ZmZlciA9IGpzemxpYi5pbmZsYXRlQnVmZmVyO1xuICAgIHZhciBhcnJheUNvcHkgPSBqc3psaWIuYXJyYXlDb3B5O1xufVxuXG52YXIgQklHX1dJR19NQUdJQyA9IDB4ODg4RkZDMjY7XG52YXIgQklHX1dJR19NQUdJQ19CRSA9IDB4MjZGQzhGODg7XG52YXIgQklHX0JFRF9NQUdJQyA9IDB4ODc4OUYyRUI7XG52YXIgQklHX0JFRF9NQUdJQ19CRSA9IDB4RUJGMjg5ODc7XG5cblxudmFyIEJJR19XSUdfVFlQRV9HUkFQSCA9IDE7XG52YXIgQklHX1dJR19UWVBFX1ZTVEVQID0gMjtcbnZhciBCSUdfV0lHX1RZUEVfRlNURVAgPSAzO1xuICBcbnZhciBNMSA9IDI1NjtcbnZhciBNMiA9IDI1NioyNTY7XG52YXIgTTMgPSAyNTYqMjU2KjI1NjtcbnZhciBNNCA9IDI1NioyNTYqMjU2KjI1NjtcblxudmFyIEJFRF9DT0xPUl9SRUdFWFAgPSBuZXcgUmVnRXhwKFwiXlswLTldKyxbMC05XSssWzAtOV0rXCIpO1xuXG5mdW5jdGlvbiBid2dfcmVhZE9mZnNldChiYSwgbykge1xuICAgIHZhciBvZmZzZXQgPSBiYVtvXSArIGJhW28rMV0qTTEgKyBiYVtvKzJdKk0yICsgYmFbbyszXSpNMyArIGJhW28rNF0qTTQ7XG4gICAgcmV0dXJuIG9mZnNldDtcbn1cblxuZnVuY3Rpb24gQmlnV2lnKCkge1xufVxuXG5CaWdXaWcucHJvdG90eXBlLnJlYWRDaHJvbVRyZWUgPSBmdW5jdGlvbihjYWxsYmFjaykge1xuICAgIHZhciB0aGlzQiA9IHRoaXM7XG4gICAgdGhpcy5jaHJvbXNUb0lEcyA9IHt9O1xuICAgIHRoaXMuaWRzVG9DaHJvbXMgPSB7fTtcbiAgICB0aGlzLm1heElEID0gMDtcblxuICAgIHZhciB1ZG8gPSB0aGlzLnVuem9vbWVkRGF0YU9mZnNldDtcbiAgICB2YXIgZWIgPSAodWRvIC0gdGhpcy5jaHJvbVRyZWVPZmZzZXQpICYgMztcbiAgICB1ZG8gPSB1ZG8gKyA0IC0gZWI7XG5cbiAgICB0aGlzLmRhdGEuc2xpY2UodGhpcy5jaHJvbVRyZWVPZmZzZXQsIHVkbyAtIHRoaXMuY2hyb21UcmVlT2Zmc2V0KS5mZXRjaChmdW5jdGlvbihicHQpIHtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIGJwdE1hZ2ljID0gbGFbMF07XG4gICAgICAgIHZhciBibG9ja1NpemUgPSBsYVsxXTtcbiAgICAgICAgdmFyIGtleVNpemUgPSBsYVsyXTtcbiAgICAgICAgdmFyIHZhbFNpemUgPSBsYVszXTtcbiAgICAgICAgdmFyIGl0ZW1Db3VudCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCAxNik7XG4gICAgICAgIHZhciByb290Tm9kZU9mZnNldCA9IDMyO1xuXG4gICAgICAgIHZhciBicHRSZWFkTm9kZSA9IGZ1bmN0aW9uKG9mZnNldCkge1xuICAgICAgICAgICAgdmFyIG5vZGVUeXBlID0gYmFbb2Zmc2V0XTtcbiAgICAgICAgICAgIHZhciBjbnQgPSBzYVsob2Zmc2V0LzIpICsgMV07XG4gICAgICAgICAgICBvZmZzZXQgKz0gNDtcbiAgICAgICAgICAgIGZvciAodmFyIG4gPSAwOyBuIDwgY250OyArK24pIHtcbiAgICAgICAgICAgICAgICBpZiAobm9kZVR5cGUgPT0gMCkge1xuICAgICAgICAgICAgICAgICAgICBvZmZzZXQgKz0ga2V5U2l6ZTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGNoaWxkT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIG9mZnNldCk7XG4gICAgICAgICAgICAgICAgICAgIG9mZnNldCArPSA4O1xuICAgICAgICAgICAgICAgICAgICBjaGlsZE9mZnNldCAtPSB0aGlzQi5jaHJvbVRyZWVPZmZzZXQ7XG4gICAgICAgICAgICAgICAgICAgIGJwdFJlYWROb2RlKGNoaWxkT2Zmc2V0KTtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICB2YXIga2V5ID0gJyc7XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGtpID0gMDsga2kgPCBrZXlTaXplOyArK2tpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgY2hhckNvZGUgPSBiYVtvZmZzZXQrK107XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoY2hhckNvZGUgIT0gMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGtleSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGNoYXJDb2RlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB2YXIgY2hyb21JZCA9IChiYVtvZmZzZXQrM108PDI0KSB8IChiYVtvZmZzZXQrMl08PDE2KSB8IChiYVtvZmZzZXQrMV08PDgpIHwgKGJhW29mZnNldCswXSk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBjaHJvbVNpemUgPSAoYmFbb2Zmc2V0ICsgN108PDI0KSB8IChiYVtvZmZzZXQrNl08PDE2KSB8IChiYVtvZmZzZXQrNV08PDgpIHwgKGJhW29mZnNldCs0XSk7XG4gICAgICAgICAgICAgICAgICAgIG9mZnNldCArPSA4O1xuXG4gICAgICAgICAgICAgICAgICAgIHRoaXNCLmNocm9tc1RvSURzW2tleV0gPSBjaHJvbUlkO1xuICAgICAgICAgICAgICAgICAgICBpZiAoa2V5LmluZGV4T2YoJ2NocicpID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHRoaXNCLmNocm9tc1RvSURzW2tleS5zdWJzdHIoMyldID0gY2hyb21JZDtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB0aGlzQi5pZHNUb0Nocm9tc1tjaHJvbUlkXSA9IGtleTtcbiAgICAgICAgICAgICAgICAgICAgdGhpc0IubWF4SUQgPSBNYXRoLm1heCh0aGlzQi5tYXhJRCwgY2hyb21JZCk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9O1xuICAgICAgICBicHRSZWFkTm9kZShyb290Tm9kZU9mZnNldCk7XG5cbiAgICAgICAgY2FsbGJhY2sodGhpc0IpO1xuICAgIH0pO1xufVxuXG5mdW5jdGlvbiBCaWdXaWdWaWV3KGJ3ZywgY2lyVHJlZU9mZnNldCwgY2lyVHJlZUxlbmd0aCwgaXNTdW1tYXJ5KSB7XG4gICAgdGhpcy5id2cgPSBid2c7XG4gICAgdGhpcy5jaXJUcmVlT2Zmc2V0ID0gY2lyVHJlZU9mZnNldDtcbiAgICB0aGlzLmNpclRyZWVMZW5ndGggPSBjaXJUcmVlTGVuZ3RoO1xuICAgIHRoaXMuaXNTdW1tYXJ5ID0gaXNTdW1tYXJ5O1xufVxuXG5cblxuQmlnV2lnVmlldy5wcm90b3R5cGUucmVhZFdpZ0RhdGEgPSBmdW5jdGlvbihjaHJOYW1lLCBtaW4sIG1heCwgY2FsbGJhY2spIHtcbiAgICB2YXIgY2hyID0gdGhpcy5id2cuY2hyb21zVG9JRHNbY2hyTmFtZV07XG4gICAgaWYgKGNociA9PT0gdW5kZWZpbmVkKSB7XG4gICAgICAgIC8vIE5vdCBhbiBlcnJvciBiZWNhdXNlIHNvbWUgLmJ3Z3Mgd29uJ3QgaGF2ZSBkYXRhIGZvciBhbGwgY2hyb21vc29tZXMuXG4gICAgICAgIHJldHVybiBjYWxsYmFjayhbXSk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdGhpcy5yZWFkV2lnRGF0YUJ5SWQoY2hyLCBtaW4sIG1heCwgY2FsbGJhY2spO1xuICAgIH1cbn1cblxuQmlnV2lnVmlldy5wcm90b3R5cGUucmVhZFdpZ0RhdGFCeUlkID0gZnVuY3Rpb24oY2hyLCBtaW4sIG1heCwgY2FsbGJhY2spIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuICAgIGlmICghdGhpcy5jaXJIZWFkZXIpIHtcbiAgICAgICAgdGhpcy5id2cuZGF0YS5zbGljZSh0aGlzLmNpclRyZWVPZmZzZXQsIDQ4KS5mZXRjaChmdW5jdGlvbihyZXN1bHQpIHtcbiAgICAgICAgICAgIHRoaXNCLmNpckhlYWRlciA9IHJlc3VsdDtcbiAgICAgICAgICAgIHZhciBsYSA9IG5ldyBJbnQzMkFycmF5KHRoaXNCLmNpckhlYWRlcik7XG4gICAgICAgICAgICB0aGlzQi5jaXJCbG9ja1NpemUgPSBsYVsxXTtcbiAgICAgICAgICAgIHRoaXNCLnJlYWRXaWdEYXRhQnlJZChjaHIsIG1pbiwgbWF4LCBjYWxsYmFjayk7XG4gICAgICAgIH0pO1xuICAgICAgICByZXR1cm47XG4gICAgfVxuXG4gICAgdmFyIGJsb2Nrc1RvRmV0Y2ggPSBbXTtcbiAgICB2YXIgb3V0c3RhbmRpbmcgPSAwO1xuXG4gICAgdmFyIGJlZm9yZUJXRyA9IERhdGUubm93KCk7XG5cbiAgICB2YXIgZmlsdGVyID0gZnVuY3Rpb24oY2hyb21JZCwgZm1pbiwgZm1heCwgdG9rcykge1xuICAgICAgICByZXR1cm4gKChjaHIgPCAwIHx8IGNocm9tSWQgPT0gY2hyKSAmJiBmbWluIDw9IG1heCAmJiBmbWF4ID49IG1pbik7XG4gICAgfVxuXG4gICAgdmFyIGNpckZvYlJlY3VyID0gZnVuY3Rpb24ob2Zmc2V0LCBsZXZlbCkge1xuICAgICAgICBpZiAodGhpc0IuYndnLmluc3RydW1lbnQpXG4gICAgICAgICAgICBjb25zb2xlLmxvZygnbGV2ZWw9JyArIGxldmVsICsgJzsgb2Zmc2V0PScgKyBvZmZzZXQgKyAnOyB0aW1lPScgKyAoRGF0ZS5ub3coKXwwKSk7XG5cbiAgICAgICAgb3V0c3RhbmRpbmcgKz0gb2Zmc2V0Lmxlbmd0aDtcblxuICAgICAgICBpZiAob2Zmc2V0Lmxlbmd0aCA9PSAxICYmIG9mZnNldFswXSAtIHRoaXNCLmNpclRyZWVPZmZzZXQgPT0gNDggJiYgdGhpc0IuY2FjaGVkQ2lyUm9vdCkge1xuICAgICAgICAgICAgY2lyRm9iUmVjdXIyKHRoaXNCLmNhY2hlZENpclJvb3QsIDAsIGxldmVsKTtcbiAgICAgICAgICAgIC0tb3V0c3RhbmRpbmc7XG4gICAgICAgICAgICBpZiAob3V0c3RhbmRpbmcgPT0gMCkge1xuICAgICAgICAgICAgICAgIHRoaXNCLmZldGNoRmVhdHVyZXMoZmlsdGVyLCBibG9ja3NUb0ZldGNoLCBjYWxsYmFjayk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXR1cm47XG4gICAgICAgIH1cblxuICAgICAgICB2YXIgbWF4Q2lyQmxvY2tTcGFuID0gNCArICAodGhpc0IuY2lyQmxvY2tTaXplICogMzIpOyAgIC8vIFVwcGVyIGJvdW5kIG9uIHNpemUsIGJhc2VkIG9uIGEgY29tcGxldGVseSBmdWxsIGxlYWYgbm9kZS5cbiAgICAgICAgdmFyIHNwYW5zO1xuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG9mZnNldC5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgdmFyIGJsb2NrU3BhbiA9IG5ldyBSYW5nZShvZmZzZXRbaV0sIG9mZnNldFtpXSArIG1heENpckJsb2NrU3Bhbik7XG4gICAgICAgICAgICBzcGFucyA9IHNwYW5zID8gdW5pb24oc3BhbnMsIGJsb2NrU3BhbikgOiBibG9ja1NwYW47XG4gICAgICAgIH1cbiAgICAgICAgXG4gICAgICAgIHZhciBmZXRjaFJhbmdlcyA9IHNwYW5zLnJhbmdlcygpO1xuICAgICAgICBmb3IgKHZhciByID0gMDsgciA8IGZldGNoUmFuZ2VzLmxlbmd0aDsgKytyKSB7XG4gICAgICAgICAgICB2YXIgZnIgPSBmZXRjaFJhbmdlc1tyXTtcbiAgICAgICAgICAgIGNpckZvYlN0YXJ0RmV0Y2gob2Zmc2V0LCBmciwgbGV2ZWwpO1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgdmFyIGNpckZvYlN0YXJ0RmV0Y2ggPSBmdW5jdGlvbihvZmZzZXQsIGZyLCBsZXZlbCwgYXR0ZW1wdHMpIHtcbiAgICAgICAgdmFyIGxlbmd0aCA9IGZyLm1heCgpIC0gZnIubWluKCk7XG4gICAgICAgIHRoaXNCLmJ3Zy5kYXRhLnNsaWNlKGZyLm1pbigpLCBmci5tYXgoKSAtIGZyLm1pbigpKS5mZXRjaChmdW5jdGlvbihyZXN1bHRCdWZmZXIpIHtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgb2Zmc2V0Lmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICAgICAgaWYgKGZyLmNvbnRhaW5zKG9mZnNldFtpXSkpIHtcbiAgICAgICAgICAgICAgICAgICAgY2lyRm9iUmVjdXIyKHJlc3VsdEJ1ZmZlciwgb2Zmc2V0W2ldIC0gZnIubWluKCksIGxldmVsKTtcblxuICAgICAgICAgICAgICAgICAgICBpZiAob2Zmc2V0W2ldIC0gdGhpc0IuY2lyVHJlZU9mZnNldCA9PSA0OCAmJiBvZmZzZXRbaV0gLSBmci5taW4oKSA9PSAwKVxuICAgICAgICAgICAgICAgICAgICAgICAgdGhpc0IuY2FjaGVkQ2lyUm9vdCA9IHJlc3VsdEJ1ZmZlcjtcblxuICAgICAgICAgICAgICAgICAgICAtLW91dHN0YW5kaW5nO1xuICAgICAgICAgICAgICAgICAgICBpZiAob3V0c3RhbmRpbmcgPT0gMCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdGhpc0IuZmV0Y2hGZWF0dXJlcyhmaWx0ZXIsIGJsb2Nrc1RvRmV0Y2gsIGNhbGxiYWNrKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgfVxuXG4gICAgdmFyIGNpckZvYlJlY3VyMiA9IGZ1bmN0aW9uKGNpckJsb2NrRGF0YSwgb2Zmc2V0LCBsZXZlbCkge1xuICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShjaXJCbG9ja0RhdGEpO1xuICAgICAgICB2YXIgc2EgPSBuZXcgSW50MTZBcnJheShjaXJCbG9ja0RhdGEpO1xuICAgICAgICB2YXIgbGEgPSBuZXcgSW50MzJBcnJheShjaXJCbG9ja0RhdGEpO1xuXG4gICAgICAgIHZhciBpc0xlYWYgPSBiYVtvZmZzZXRdO1xuICAgICAgICB2YXIgY250ID0gc2Fbb2Zmc2V0LzIgKyAxXTtcbiAgICAgICAgb2Zmc2V0ICs9IDQ7XG5cbiAgICAgICAgaWYgKGlzTGVhZiAhPSAwKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGNudDsgKytpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGxvID0gb2Zmc2V0LzQ7XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0Q2hyb20gPSBsYVtsb107XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0QmFzZSA9IGxhW2xvICsgMV07XG4gICAgICAgICAgICAgICAgdmFyIGVuZENocm9tID0gbGFbbG8gKyAyXTtcbiAgICAgICAgICAgICAgICB2YXIgZW5kQmFzZSA9IGxhW2xvICsgM107XG4gICAgICAgICAgICAgICAgdmFyIGJsb2NrT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIG9mZnNldCsxNik7XG4gICAgICAgICAgICAgICAgdmFyIGJsb2NrU2l6ZSA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQrMjQpO1xuICAgICAgICAgICAgICAgIGlmICgoKGNociA8IDAgfHwgc3RhcnRDaHJvbSA8IGNocikgfHwgKHN0YXJ0Q2hyb20gPT0gY2hyICYmIHN0YXJ0QmFzZSA8PSBtYXgpKSAmJlxuICAgICAgICAgICAgICAgICAgICAoKGNociA8IDAgfHwgZW5kQ2hyb20gICA+IGNocikgfHwgKGVuZENocm9tID09IGNociAmJiBlbmRCYXNlID49IG1pbikpKVxuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgYmxvY2tzVG9GZXRjaC5wdXNoKHtvZmZzZXQ6IGJsb2NrT2Zmc2V0LCBzaXplOiBibG9ja1NpemV9KTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgb2Zmc2V0ICs9IDMyO1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdmFyIHJlY3VyT2Zmc2V0cyA9IFtdO1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBjbnQ7ICsraSkge1xuICAgICAgICAgICAgICAgIHZhciBsbyA9IG9mZnNldC80O1xuICAgICAgICAgICAgICAgIHZhciBzdGFydENocm9tID0gbGFbbG9dO1xuICAgICAgICAgICAgICAgIHZhciBzdGFydEJhc2UgPSBsYVtsbyArIDFdO1xuICAgICAgICAgICAgICAgIHZhciBlbmRDaHJvbSA9IGxhW2xvICsgMl07XG4gICAgICAgICAgICAgICAgdmFyIGVuZEJhc2UgPSBsYVtsbyArIDNdO1xuICAgICAgICAgICAgICAgIHZhciBibG9ja09mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQrMTYpO1xuICAgICAgICAgICAgICAgIGlmICgoY2hyIDwgMCB8fCBzdGFydENocm9tIDwgY2hyIHx8IChzdGFydENocm9tID09IGNociAmJiBzdGFydEJhc2UgPD0gbWF4KSkgJiZcbiAgICAgICAgICAgICAgICAgICAgKGNociA8IDAgfHwgZW5kQ2hyb20gICA+IGNociB8fCAoZW5kQ2hyb20gPT0gY2hyICYmIGVuZEJhc2UgPj0gbWluKSkpXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICByZWN1ck9mZnNldHMucHVzaChibG9ja09mZnNldCk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIG9mZnNldCArPSAyNDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmIChyZWN1ck9mZnNldHMubGVuZ3RoID4gMCkge1xuICAgICAgICAgICAgICAgIGNpckZvYlJlY3VyKHJlY3VyT2Zmc2V0cywgbGV2ZWwgKyAxKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH07XG5cbiAgICBjaXJGb2JSZWN1cihbdGhpc0IuY2lyVHJlZU9mZnNldCArIDQ4XSwgMSk7XG59XG5cblxuQmlnV2lnVmlldy5wcm90b3R5cGUuZmV0Y2hGZWF0dXJlcyA9IGZ1bmN0aW9uKGZpbHRlciwgYmxvY2tzVG9GZXRjaCwgY2FsbGJhY2spIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuXG4gICAgYmxvY2tzVG9GZXRjaC5zb3J0KGZ1bmN0aW9uKGIwLCBiMSkge1xuICAgICAgICByZXR1cm4gKGIwLm9mZnNldHwwKSAtIChiMS5vZmZzZXR8MCk7XG4gICAgfSk7XG5cbiAgICBpZiAoYmxvY2tzVG9GZXRjaC5sZW5ndGggPT0gMCkge1xuICAgICAgICBjYWxsYmFjayhbXSk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIGZlYXR1cmVzID0gW107XG4gICAgICAgIHZhciBjcmVhdGVGZWF0dXJlID0gZnVuY3Rpb24oY2hyLCBmbWluLCBmbWF4LCBvcHRzKSB7XG4gICAgICAgICAgICBpZiAoIW9wdHMpIHtcbiAgICAgICAgICAgICAgICBvcHRzID0ge307XG4gICAgICAgICAgICB9XG4gICAgICAgIFxuICAgICAgICAgICAgdmFyIGYgPSBuZXcgREFTRmVhdHVyZSgpO1xuICAgICAgICAgICAgZi5fY2hyb21JZCA9IGNocjtcbiAgICAgICAgICAgIGYuc2VnbWVudCA9IHRoaXNCLmJ3Zy5pZHNUb0Nocm9tc1tjaHJdO1xuICAgICAgICAgICAgZi5taW4gPSBmbWluO1xuICAgICAgICAgICAgZi5tYXggPSBmbWF4O1xuICAgICAgICAgICAgZi50eXBlID0gdGhpc0IuYndnLnR5cGU7XG4gICAgICAgICAgICBcbiAgICAgICAgICAgIGZvciAodmFyIGsgaW4gb3B0cykge1xuICAgICAgICAgICAgICAgIGZba10gPSBvcHRzW2tdO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgXG4gICAgICAgICAgICBmZWF0dXJlcy5wdXNoKGYpO1xuICAgICAgICB9O1xuXG4gICAgICAgIHZhciB0cmFtcCA9IGZ1bmN0aW9uKCkge1xuICAgICAgICAgICAgaWYgKGJsb2Nrc1RvRmV0Y2gubGVuZ3RoID09IDApIHtcbiAgICAgICAgICAgICAgICB2YXIgYWZ0ZXJCV0cgPSBEYXRlLm5vdygpO1xuICAgICAgICAgICAgICAgIC8vIGRsb2coJ0JXRyBmZXRjaCB0b29rICcgKyAoYWZ0ZXJCV0cgLSBiZWZvcmVCV0cpICsgJ21zJyk7XG4gICAgICAgICAgICAgICAgY2FsbGJhY2soZmVhdHVyZXMpO1xuICAgICAgICAgICAgICAgIHJldHVybjsgIC8vIGp1c3QgaW4gY2FzZS4uLlxuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICB2YXIgYmxvY2sgPSBibG9ja3NUb0ZldGNoWzBdO1xuICAgICAgICAgICAgICAgIGlmIChibG9jay5kYXRhKSB7XG4gICAgICAgICAgICAgICAgICAgIHRoaXNCLnBhcnNlRmVhdHVyZXMoYmxvY2suZGF0YSwgY3JlYXRlRmVhdHVyZSwgZmlsdGVyKTtcbiAgICAgICAgICAgICAgICAgICAgYmxvY2tzVG9GZXRjaC5zcGxpY2UoMCwgMSk7XG4gICAgICAgICAgICAgICAgICAgIHRyYW1wKCk7XG4gICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGZldGNoU3RhcnQgPSBibG9jay5vZmZzZXQ7XG4gICAgICAgICAgICAgICAgICAgIHZhciBmZXRjaFNpemUgPSBibG9jay5zaXplO1xuICAgICAgICAgICAgICAgICAgICB2YXIgYmkgPSAxO1xuICAgICAgICAgICAgICAgICAgICB3aGlsZSAoYmkgPCBibG9ja3NUb0ZldGNoLmxlbmd0aCAmJiBibG9ja3NUb0ZldGNoW2JpXS5vZmZzZXQgPT0gKGZldGNoU3RhcnQgKyBmZXRjaFNpemUpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBmZXRjaFNpemUgKz0gYmxvY2tzVG9GZXRjaFtiaV0uc2l6ZTtcbiAgICAgICAgICAgICAgICAgICAgICAgICsrYmk7XG4gICAgICAgICAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgICAgICAgICB0aGlzQi5id2cuZGF0YS5zbGljZShmZXRjaFN0YXJ0LCBmZXRjaFNpemUpLmZldGNoKGZ1bmN0aW9uKHJlc3VsdCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIG9mZnNldCA9IDA7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgYmkgPSAwO1xuICAgICAgICAgICAgICAgICAgICAgICAgd2hpbGUgKG9mZnNldCA8IGZldGNoU2l6ZSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBmYiA9IGJsb2Nrc1RvRmV0Y2hbYmldO1xuICAgICAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGRhdGE7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYgKHRoaXNCLmJ3Zy51bmNvbXByZXNzQnVmU2l6ZSA+IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgZGF0YSA9IGpzemxpYl9pbmZsYXRlX2J1ZmZlcihyZXN1bHQsIG9mZnNldCArIDIsIGZiLnNpemUgLSAyKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgdG1wID0gbmV3IFVpbnQ4QXJyYXkoZmIuc2l6ZSk7ICAgIC8vIEZJWE1FIGlzIHRoaXMgcmVhbGx5IHRoZSBiZXN0IHdlIGNhbiBkbz9cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgYXJyYXlDb3B5KG5ldyBVaW50OEFycmF5KHJlc3VsdCwgb2Zmc2V0LCBmYi5zaXplKSwgMCwgdG1wLCAwLCBmYi5zaXplKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgZGF0YSA9IHRtcC5idWZmZXI7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGZiLmRhdGEgPSBkYXRhO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIG9mZnNldCArPSBmYi5zaXplO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICsrYmk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB0cmFtcCgpO1xuICAgICAgICAgICAgICAgICAgICB9KTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgdHJhbXAoKTtcbiAgICB9XG59XG5cbkJpZ1dpZ1ZpZXcucHJvdG90eXBlLnBhcnNlRmVhdHVyZXMgPSBmdW5jdGlvbihkYXRhLCBjcmVhdGVGZWF0dXJlLCBmaWx0ZXIpIHtcbiAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShkYXRhKTtcblxuICAgIGlmICh0aGlzLmlzU3VtbWFyeSkge1xuICAgICAgICB2YXIgc2EgPSBuZXcgSW50MTZBcnJheShkYXRhKTtcbiAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoZGF0YSk7XG4gICAgICAgIHZhciBmYSA9IG5ldyBGbG9hdDMyQXJyYXkoZGF0YSk7XG5cbiAgICAgICAgdmFyIGl0ZW1Db3VudCA9IGRhdGEuYnl0ZUxlbmd0aC8zMjtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBpdGVtQ291bnQ7ICsraSkge1xuICAgICAgICAgICAgdmFyIGNocm9tSWQgPSAgIGxhWyhpKjgpXTtcbiAgICAgICAgICAgIHZhciBzdGFydCA9ICAgICBsYVsoaSo4KSsxXTtcbiAgICAgICAgICAgIHZhciBlbmQgPSAgICAgICBsYVsoaSo4KSsyXTtcbiAgICAgICAgICAgIHZhciB2YWxpZENudCA9ICBsYVsoaSo4KSszXTtcbiAgICAgICAgICAgIHZhciBtaW5WYWwgICAgPSBmYVsoaSo4KSs0XTtcbiAgICAgICAgICAgIHZhciBtYXhWYWwgICAgPSBmYVsoaSo4KSs1XTtcbiAgICAgICAgICAgIHZhciBzdW1EYXRhICAgPSBmYVsoaSo4KSs2XTtcbiAgICAgICAgICAgIHZhciBzdW1TcURhdGEgPSBmYVsoaSo4KSs3XTtcbiAgICAgICAgICAgIFxuICAgICAgICAgICAgaWYgKGZpbHRlcihjaHJvbUlkLCBzdGFydCArIDEsIGVuZCkpIHtcbiAgICAgICAgICAgICAgICB2YXIgc3VtbWFyeU9wdHMgPSB7dHlwZTogJ2JpZ3dpZycsIHNjb3JlOiBzdW1EYXRhL3ZhbGlkQ250LCBtYXhTY29yZTogbWF4VmFsfTtcbiAgICAgICAgICAgICAgICBpZiAodGhpcy5id2cudHlwZSA9PSAnYmlnYmVkJykge1xuICAgICAgICAgICAgICAgICAgICBzdW1tYXJ5T3B0cy50eXBlID0gJ2RlbnNpdHknO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBjcmVhdGVGZWF0dXJlKGNocm9tSWQsIHN0YXJ0ICsgMSwgZW5kLCBzdW1tYXJ5T3B0cyk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9IGVsc2UgaWYgKHRoaXMuYndnLnR5cGUgPT0gJ2JpZ3dpZycpIHtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoZGF0YSk7XG4gICAgICAgIHZhciBsYSA9IG5ldyBJbnQzMkFycmF5KGRhdGEpO1xuICAgICAgICB2YXIgZmEgPSBuZXcgRmxvYXQzMkFycmF5KGRhdGEpO1xuXG4gICAgICAgIHZhciBjaHJvbUlkID0gbGFbMF07XG4gICAgICAgIHZhciBibG9ja1N0YXJ0ID0gbGFbMV07XG4gICAgICAgIHZhciBibG9ja0VuZCA9IGxhWzJdO1xuICAgICAgICB2YXIgaXRlbVN0ZXAgPSBsYVszXTtcbiAgICAgICAgdmFyIGl0ZW1TcGFuID0gbGFbNF07XG4gICAgICAgIHZhciBibG9ja1R5cGUgPSBiYVsyMF07XG4gICAgICAgIHZhciBpdGVtQ291bnQgPSBzYVsxMV07XG4gICAgICAgIFxuICAgICAgICBpZiAoYmxvY2tUeXBlID09IEJJR19XSUdfVFlQRV9GU1RFUCkge1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBpdGVtQ291bnQ7ICsraSkge1xuICAgICAgICAgICAgICAgIHZhciBzY29yZSA9IGZhW2kgKyA2XTtcbiAgICAgICAgICAgICAgICB2YXIgZm1pbiA9IGJsb2NrU3RhcnQgKyAoaSppdGVtU3RlcCkgKyAxLCBmbWF4ID0gYmxvY2tTdGFydCArIChpKml0ZW1TdGVwKSArIGl0ZW1TcGFuO1xuICAgICAgICAgICAgICAgIGlmIChmaWx0ZXIoY2hyb21JZCwgZm1pbiwgZm1heCkpXG4gICAgICAgICAgICAgICAgICAgIGNyZWF0ZUZlYXR1cmUoY2hyb21JZCwgZm1pbiwgZm1heCwge3Njb3JlOiBzY29yZX0pO1xuICAgICAgICAgICAgfVxuICAgICAgICB9IGVsc2UgaWYgKGJsb2NrVHlwZSA9PSBCSUdfV0lHX1RZUEVfVlNURVApIHtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgaXRlbUNvdW50OyArK2kpIHtcbiAgICAgICAgICAgICAgICB2YXIgc3RhcnQgPSBsYVsoaSoyKSArIDZdICsgMTtcbiAgICAgICAgICAgICAgICB2YXIgZW5kID0gc3RhcnQgKyBpdGVtU3BhbiAtIDE7XG4gICAgICAgICAgICAgICAgdmFyIHNjb3JlID0gZmFbKGkqMikgKyA3XTtcbiAgICAgICAgICAgICAgICBpZiAoZmlsdGVyKGNocm9tSWQsIHN0YXJ0LCBlbmQpKVxuICAgICAgICAgICAgICAgICAgICBjcmVhdGVGZWF0dXJlKGNocm9tSWQsIHN0YXJ0LCBlbmQsIHtzY29yZTogc2NvcmV9KTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSBlbHNlIGlmIChibG9ja1R5cGUgPT0gQklHX1dJR19UWVBFX0dSQVBIKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGl0ZW1Db3VudDsgKytpKSB7XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0ID0gbGFbKGkqMykgKyA2XSArIDE7XG4gICAgICAgICAgICAgICAgdmFyIGVuZCAgID0gbGFbKGkqMykgKyA3XTtcbiAgICAgICAgICAgICAgICB2YXIgc2NvcmUgPSBmYVsoaSozKSArIDhdO1xuICAgICAgICAgICAgICAgIGlmIChzdGFydCA+IGVuZCkge1xuICAgICAgICAgICAgICAgICAgICBzdGFydCA9IGVuZDtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgaWYgKGZpbHRlcihjaHJvbUlkLCBzdGFydCwgZW5kKSlcbiAgICAgICAgICAgICAgICAgICAgY3JlYXRlRmVhdHVyZShjaHJvbUlkLCBzdGFydCwgZW5kLCB7c2NvcmU6IHNjb3JlfSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBjb25zb2xlLmxvZygnQ3VycmVudGx5IG5vdCBoYW5kbGluZyBid2dUeXBlPScgKyBibG9ja1R5cGUpO1xuICAgICAgICB9XG4gICAgfSBlbHNlIGlmICh0aGlzLmJ3Zy50eXBlID09ICdiaWdiZWQnKSB7XG4gICAgICAgIHZhciBvZmZzZXQgPSAwO1xuICAgICAgICB2YXIgZGZjID0gdGhpcy5id2cuZGVmaW5lZEZpZWxkQ291bnQ7XG4gICAgICAgIHZhciBzY2hlbWEgPSB0aGlzLmJ3Zy5zY2hlbWE7XG5cbiAgICAgICAgd2hpbGUgKG9mZnNldCA8IGJhLmxlbmd0aCkge1xuICAgICAgICAgICAgdmFyIGNocm9tSWQgPSAoYmFbb2Zmc2V0KzNdPDwyNCkgfCAoYmFbb2Zmc2V0KzJdPDwxNikgfCAoYmFbb2Zmc2V0KzFdPDw4KSB8IChiYVtvZmZzZXQrMF0pO1xuICAgICAgICAgICAgdmFyIHN0YXJ0ID0gKGJhW29mZnNldCs3XTw8MjQpIHwgKGJhW29mZnNldCs2XTw8MTYpIHwgKGJhW29mZnNldCs1XTw8OCkgfCAoYmFbb2Zmc2V0KzRdKTtcbiAgICAgICAgICAgIHZhciBlbmQgPSAoYmFbb2Zmc2V0KzExXTw8MjQpIHwgKGJhW29mZnNldCsxMF08PDE2KSB8IChiYVtvZmZzZXQrOV08PDgpIHwgKGJhW29mZnNldCs4XSk7XG4gICAgICAgICAgICBvZmZzZXQgKz0gMTI7XG4gICAgICAgICAgICB2YXIgcmVzdCA9ICcnO1xuICAgICAgICAgICAgd2hpbGUgKHRydWUpIHtcbiAgICAgICAgICAgICAgICB2YXIgY2ggPSBiYVtvZmZzZXQrK107XG4gICAgICAgICAgICAgICAgaWYgKGNoICE9IDApIHtcbiAgICAgICAgICAgICAgICAgICAgcmVzdCArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGNoKTtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICBicmVhaztcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHZhciBmZWF0dXJlT3B0cyA9IHt9O1xuICAgICAgICAgICAgXG4gICAgICAgICAgICB2YXIgYmVkQ29sdW1ucztcbiAgICAgICAgICAgIGlmIChyZXN0Lmxlbmd0aCA+IDApIHtcbiAgICAgICAgICAgICAgICBiZWRDb2x1bW5zID0gcmVzdC5zcGxpdCgnXFx0Jyk7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIGJlZENvbHVtbnMgPSBbXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmIChiZWRDb2x1bW5zLmxlbmd0aCA+IDAgJiYgZGZjID4gMykge1xuICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzLmxhYmVsID0gYmVkQ29sdW1uc1swXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmIChiZWRDb2x1bW5zLmxlbmd0aCA+IDEgJiYgZGZjID4gNCkge1xuICAgICAgICAgICAgICAgIHZhciBzY29yZSA9IHBhcnNlSW50KGJlZENvbHVtbnNbMV0pO1xuICAgICAgICAgICAgICAgIGlmICghaXNOYU4oc2NvcmUpKVxuICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy5zY29yZSA9IHNjb3JlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgaWYgKGJlZENvbHVtbnMubGVuZ3RoID4gMiAmJiBkZmMgPiA1KSB7XG4gICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMub3JpZW50YXRpb24gPSBiZWRDb2x1bW5zWzJdO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgaWYgKGJlZENvbHVtbnMubGVuZ3RoID4gNSAmJiBkZmMgPiA4KSB7XG4gICAgICAgICAgICAgICAgdmFyIGNvbG9yID0gYmVkQ29sdW1uc1s1XTtcbiAgICAgICAgICAgICAgICBpZiAoQkVEX0NPTE9SX1JFR0VYUC50ZXN0KGNvbG9yKSkge1xuICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy5pdGVtUmdiID0gJ3JnYignICsgY29sb3IgKyAnKSc7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBpZiAoYmVkQ29sdW1ucy5sZW5ndGggPiBkZmMtMyAmJiBzY2hlbWEpIHtcbiAgICAgICAgICAgICAgICBmb3IgKHZhciBjb2wgPSBkZmMgLSAzOyBjb2wgPCBiZWRDb2x1bW5zLmxlbmd0aDsgKytjb2wpIHtcbiAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHNbc2NoZW1hLmZpZWxkc1tjb2wrM10ubmFtZV0gPSBiZWRDb2x1bW5zW2NvbF07XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICBpZiAoZmlsdGVyKGNocm9tSWQsIHN0YXJ0ICsgMSwgZW5kLCBiZWRDb2x1bW5zKSkge1xuICAgICAgICAgICAgICAgIGlmIChkZmMgPCAxMikge1xuICAgICAgICAgICAgICAgICAgICBjcmVhdGVGZWF0dXJlKGNocm9tSWQsIHN0YXJ0ICsgMSwgZW5kLCBmZWF0dXJlT3B0cyk7XG4gICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHRoaWNrU3RhcnQgPSBiZWRDb2x1bW5zWzNdfDA7XG4gICAgICAgICAgICAgICAgICAgIHZhciB0aGlja0VuZCAgID0gYmVkQ29sdW1uc1s0XXwwO1xuICAgICAgICAgICAgICAgICAgICB2YXIgYmxvY2tDb3VudCA9IGJlZENvbHVtbnNbNl18MDtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGJsb2NrU2l6ZXMgPSBiZWRDb2x1bW5zWzddLnNwbGl0KCcsJyk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBibG9ja1N0YXJ0cyA9IGJlZENvbHVtbnNbOF0uc3BsaXQoJywnKTtcblxuICAgICAgICAgICAgICAgICAgICBpZiAoZmVhdHVyZU9wdHMuZXhvbkZyYW1lcykge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGV4b25GcmFtZXMgPSBmZWF0dXJlT3B0cy5leG9uRnJhbWVzLnNwbGl0KCcsJyk7XG4gICAgICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy5leG9uRnJhbWVzID0gdW5kZWZpbmVkO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy50eXBlID0gJ3RyYW5zY3JpcHQnXG4gICAgICAgICAgICAgICAgICAgIHZhciBncnAgPSBuZXcgREFTR3JvdXAoKTtcbiAgICAgICAgICAgICAgICAgICAgZm9yICh2YXIgayBpbiBmZWF0dXJlT3B0cykge1xuICAgICAgICAgICAgICAgICAgICAgICAgZ3JwW2tdID0gZmVhdHVyZU9wdHNba107XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgZ3JwLmlkID0gYmVkQ29sdW1uc1swXTtcbiAgICAgICAgICAgICAgICAgICAgZ3JwLnNlZ21lbnQgPSB0aGlzLmJ3Zy5pZHNUb0Nocm9tc1tjaHJvbUlkXTtcbiAgICAgICAgICAgICAgICAgICAgZ3JwLm1pbiA9IHN0YXJ0ICsgMTtcbiAgICAgICAgICAgICAgICAgICAgZ3JwLm1heCA9IGVuZDtcbiAgICAgICAgICAgICAgICAgICAgZ3JwLm5vdGVzID0gW107XG4gICAgICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzLmdyb3VwcyA9IFtncnBdO1xuXG4gICAgICAgICAgICAgICAgICAgIC8vIE1vdmluZyB0b3dhcmRzIHVzaW5nIGJpZ0dlbmVQcmVkIG1vZGVsLCBidXQgd2lsbFxuICAgICAgICAgICAgICAgICAgICAvLyBzdGlsbCBzdXBwb3J0IG9sZCBEYWxsaWFuY2Utc3R5bGUgQkVEMTIrZ2VuZS1uYW1lIGZvciB0aGVcbiAgICAgICAgICAgICAgICAgICAgLy8gZm9yZXNlZWFibGUgZnV0dXJlLlxuICAgICAgICAgICAgICAgICAgICBpZiAoYmVkQ29sdW1ucy5sZW5ndGggPiA5KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgZ2VuZUlkID0gZmVhdHVyZU9wdHMuZ2VuZU5hbWUgfHwgYmVkQ29sdW1uc1s5XTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBnZW5lTmFtZSA9IGdlbmVJZDtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChiZWRDb2x1bW5zLmxlbmd0aCA+IDEwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZ2VuZU5hbWUgPSBiZWRDb2x1bW5zWzEwXTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChmZWF0dXJlT3B0cy5nZW5lTmFtZTIpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZ2VuZU5hbWUgPSBmZWF0dXJlT3B0cy5nZW5lTmFtZTI7XG5cbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBnZyA9IHNoYWxsb3dDb3B5KGdycCk7XG4gICAgICAgICAgICAgICAgICAgICAgICBnZy5pZCA9IGdlbmVJZDtcbiAgICAgICAgICAgICAgICAgICAgICAgIGdnLmxhYmVsID0gZ2VuZU5hbWU7XG4gICAgICAgICAgICAgICAgICAgICAgICBnZy50eXBlID0gJ2dlbmUnO1xuICAgICAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMuZ3JvdXBzLnB1c2goZ2cpO1xuICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgdmFyIHNwYW5MaXN0ID0gW107XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGIgPSAwOyBiIDwgYmxvY2tDb3VudDsgKytiKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgYm1pbiA9IChibG9ja1N0YXJ0c1tiXXwwKSArIHN0YXJ0O1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGJtYXggPSBibWluICsgKGJsb2NrU2l6ZXNbYl18MCk7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgc3BhbiA9IG5ldyBSYW5nZShibWluLCBibWF4KTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHNwYW5MaXN0LnB1c2goc3Bhbik7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgdmFyIHNwYW5zID0gdW5pb24oc3Bhbkxpc3QpO1xuICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgdmFyIHRzTGlzdCA9IHNwYW5zLnJhbmdlcygpO1xuICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBzID0gMDsgcyA8IHRzTGlzdC5sZW5ndGg7ICsrcykge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHRzID0gdHNMaXN0W3NdO1xuICAgICAgICAgICAgICAgICAgICAgICAgY3JlYXRlRmVhdHVyZShjaHJvbUlkLCB0cy5taW4oKSArIDEsIHRzLm1heCgpLCBmZWF0dXJlT3B0cyk7XG4gICAgICAgICAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgICAgICAgICBpZiAodGhpY2tFbmQgPiB0aGlja1N0YXJ0KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgY29kaW5nUmVnaW9uID0gKGZlYXR1cmVPcHRzLm9yaWVudGF0aW9uID09ICcrJykgP1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIG5ldyBSYW5nZSh0aGlja1N0YXJ0LCB0aGlja0VuZCArIDMpIDpcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBuZXcgUmFuZ2UodGhpY2tTdGFydCAtIDMsIHRoaWNrRW5kKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAvLyArLy0gMyB0byBhY2NvdW50IGZvciBzdG9wIGNvZG9uXG5cbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciB0bCA9IGludGVyc2VjdGlvbihzcGFucywgY29kaW5nUmVnaW9uKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmICh0bCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzLnR5cGUgPSAndHJhbnNsYXRpb24nO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciB0bExpc3QgPSB0bC5yYW5nZXMoKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgcmVhZGluZ0ZyYW1lID0gMDtcblxuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciB0bE9mZnNldCA9IDA7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgd2hpbGUgKHRsTGlzdFswXS5taW4oKSA+IHRzTGlzdFt0bE9mZnNldF0ubWF4KCkpXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHRsT2Zmc2V0Kys7XG5cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBzID0gMDsgcyA8IHRsTGlzdC5sZW5ndGg7ICsrcykge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAvLyBSZWNvcmQgcmVhZGluZyBmcmFtZSBmb3IgZXZlcnkgZXhvblxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgaW5kZXggPSBzO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoZmVhdHVyZU9wdHMub3JpZW50YXRpb24gPT0gJy0nKVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaW5kZXggPSB0bExpc3QubGVuZ3RoIC0gcyAtIDE7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciB0cyA9IHRsTGlzdFtpbmRleF07XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGZlYXR1cmVPcHRzLnJlYWRmcmFtZSA9IHJlYWRpbmdGcmFtZTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGV4b25GcmFtZXMpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBicmYgPSBwYXJzZUludChleG9uRnJhbWVzW2luZGV4ICsgdGxPZmZzZXRdKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmICh0eXBlb2YoYnJmKSA9PT0gJ251bWJlcicgJiYgYnJmID49IDAgJiYgYnJmIDw9IDIpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBmZWF0dXJlT3B0cy5yZWFkZnJhbWUgPSBicmY7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgZmVhdHVyZU9wdHMucmVhZGZyYW1lRXhwbGljaXQgPSB0cnVlO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBsZW5ndGggPSB0cy5tYXgoKSAtIHRzLm1pbigpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICByZWFkaW5nRnJhbWUgPSAocmVhZGluZ0ZyYW1lICsgbGVuZ3RoKSAlIDM7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNyZWF0ZUZlYXR1cmUoY2hyb21JZCwgdHMubWluKCkgKyAxLCB0cy5tYXgoKSwgZmVhdHVyZU9wdHMpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH0gZWxzZSB7XG4gICAgICAgIHRocm93IEVycm9yKFwiRG9uJ3Qga25vdyB3aGF0IHRvIGRvIHdpdGggXCIgKyB0aGlzLmJ3Zy50eXBlKTtcbiAgICB9XG59XG5cbi8vXG4vLyBuYXN0eSBjdXQvcGFzdGUsIHNob3VsZCByb2xsIGJhY2sgaW4hXG4vL1xuXG5CaWdXaWdWaWV3LnByb3RvdHlwZS5nZXRGaXJzdEFkamFjZW50ID0gZnVuY3Rpb24oY2hyTmFtZSwgcG9zLCBkaXIsIGNhbGxiYWNrKSB7XG4gICAgdmFyIGNociA9IHRoaXMuYndnLmNocm9tc1RvSURzW2Nock5hbWVdO1xuICAgIGlmIChjaHIgPT09IHVuZGVmaW5lZCkge1xuICAgICAgICAvLyBOb3QgYW4gZXJyb3IgYmVjYXVzZSBzb21lIC5id2dzIHdvbid0IGhhdmUgZGF0YSBmb3IgYWxsIGNocm9tb3NvbWVzLlxuICAgICAgICByZXR1cm4gY2FsbGJhY2soW10pO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHRoaXMuZ2V0Rmlyc3RBZGphY2VudEJ5SWQoY2hyLCBwb3MsIGRpciwgY2FsbGJhY2spO1xuICAgIH1cbn1cblxuQmlnV2lnVmlldy5wcm90b3R5cGUuZ2V0Rmlyc3RBZGphY2VudEJ5SWQgPSBmdW5jdGlvbihjaHIsIHBvcywgZGlyLCBjYWxsYmFjaykge1xuICAgIHZhciB0aGlzQiA9IHRoaXM7XG4gICAgaWYgKCF0aGlzLmNpckhlYWRlcikge1xuICAgICAgICB0aGlzLmJ3Zy5kYXRhLnNsaWNlKHRoaXMuY2lyVHJlZU9mZnNldCwgNDgpLmZldGNoKGZ1bmN0aW9uKHJlc3VsdCkge1xuICAgICAgICAgICAgdGhpc0IuY2lySGVhZGVyID0gcmVzdWx0O1xuICAgICAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkodGhpc0IuY2lySGVhZGVyKTtcbiAgICAgICAgICAgIHRoaXNCLmNpckJsb2NrU2l6ZSA9IGxhWzFdO1xuICAgICAgICAgICAgdGhpc0IuZ2V0Rmlyc3RBZGphY2VudEJ5SWQoY2hyLCBwb3MsIGRpciwgY2FsbGJhY2spO1xuICAgICAgICB9KTtcbiAgICAgICAgcmV0dXJuO1xuICAgIH1cblxuICAgIHZhciBibG9ja1RvRmV0Y2ggPSBudWxsO1xuICAgIHZhciBiZXN0QmxvY2tDaHIgPSAtMTtcbiAgICB2YXIgYmVzdEJsb2NrT2Zmc2V0ID0gLTE7XG5cbiAgICB2YXIgb3V0c3RhbmRpbmcgPSAwO1xuXG4gICAgdmFyIGJlZm9yZUJXRyA9IERhdGUubm93KCk7XG5cbiAgICB2YXIgY2lyRm9iUmVjdXIgPSBmdW5jdGlvbihvZmZzZXQsIGxldmVsKSB7XG4gICAgICAgIG91dHN0YW5kaW5nICs9IG9mZnNldC5sZW5ndGg7XG5cbiAgICAgICAgdmFyIG1heENpckJsb2NrU3BhbiA9IDQgKyAgKHRoaXNCLmNpckJsb2NrU2l6ZSAqIDMyKTsgICAvLyBVcHBlciBib3VuZCBvbiBzaXplLCBiYXNlZCBvbiBhIGNvbXBsZXRlbHkgZnVsbCBsZWFmIG5vZGUuXG4gICAgICAgIHZhciBzcGFucztcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBvZmZzZXQubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgIHZhciBibG9ja1NwYW4gPSBuZXcgUmFuZ2Uob2Zmc2V0W2ldLCBvZmZzZXRbaV0gKyBtYXhDaXJCbG9ja1NwYW4pO1xuICAgICAgICAgICAgc3BhbnMgPSBzcGFucyA/IHVuaW9uKHNwYW5zLCBibG9ja1NwYW4pIDogYmxvY2tTcGFuO1xuICAgICAgICB9XG4gICAgICAgIFxuICAgICAgICB2YXIgZmV0Y2hSYW5nZXMgPSBzcGFucy5yYW5nZXMoKTtcbiAgICAgICAgZm9yICh2YXIgciA9IDA7IHIgPCBmZXRjaFJhbmdlcy5sZW5ndGg7ICsrcikge1xuICAgICAgICAgICAgdmFyIGZyID0gZmV0Y2hSYW5nZXNbcl07XG4gICAgICAgICAgICBjaXJGb2JTdGFydEZldGNoKG9mZnNldCwgZnIsIGxldmVsKTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHZhciBjaXJGb2JTdGFydEZldGNoID0gZnVuY3Rpb24ob2Zmc2V0LCBmciwgbGV2ZWwsIGF0dGVtcHRzKSB7XG4gICAgICAgIHZhciBsZW5ndGggPSBmci5tYXgoKSAtIGZyLm1pbigpO1xuICAgICAgICB0aGlzQi5id2cuZGF0YS5zbGljZShmci5taW4oKSwgZnIubWF4KCkgLSBmci5taW4oKSkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0QnVmZmVyKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IG9mZnNldC5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgICAgIGlmIChmci5jb250YWlucyhvZmZzZXRbaV0pKSB7XG4gICAgICAgICAgICAgICAgICAgIGNpckZvYlJlY3VyMihyZXN1bHRCdWZmZXIsIG9mZnNldFtpXSAtIGZyLm1pbigpLCBsZXZlbCk7XG4gICAgICAgICAgICAgICAgICAgIC0tb3V0c3RhbmRpbmc7XG4gICAgICAgICAgICAgICAgICAgIGlmIChvdXRzdGFuZGluZyA9PSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoIWJsb2NrVG9GZXRjaCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmIChkaXIgPiAwICYmIChjaHIgIT0gMCB8fCBwb3MgPiAwKSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gdGhpc0IuZ2V0Rmlyc3RBZGphY2VudEJ5SWQoMCwgMCwgZGlyLCBjYWxsYmFjayk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIGlmIChkaXIgPCAwICYmIChjaHIgIT0gdGhpc0IuYndnLm1heElEIHx8IHBvcyA8IDEwMDAwMDAwMDApKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiB0aGlzQi5nZXRGaXJzdEFkamFjZW50QnlJZCh0aGlzQi5id2cubWF4SUQsIDEwMDAwMDAwMDAsIGRpciwgY2FsbGJhY2spO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2soW10pO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzQi5mZXRjaEZlYXR1cmVzKGZ1bmN0aW9uKGNocngsIGZtaW4sIGZtYXgsIHRva3MpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gKGRpciA8IDAgJiYgKGNocnggPCBjaHIgfHwgZm1heCA8IHBvcykpIHx8IChkaXIgPiAwICYmIChjaHJ4ID4gY2hyIHx8IGZtaW4gPiBwb3MpKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH0sIFtibG9ja1RvRmV0Y2hdLCBmdW5jdGlvbihmZWF0dXJlcykge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBiZXN0RmVhdHVyZSA9IG51bGw7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGJlc3RDaHIgPSAtMTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgYmVzdFBvcyA9IC0xO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGZpID0gMDsgZmkgPCBmZWF0dXJlcy5sZW5ndGg7ICsrZmkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGYgPSBmZWF0dXJlc1tmaV07XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciBjaHJ4ID0gZi5fY2hyb21JZCwgZm1pbiA9IGYubWluLCBmbWF4ID0gZi5tYXg7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmIChiZXN0RmVhdHVyZSA9PSBudWxsIHx8ICgoZGlyIDwgMCkgJiYgKGNocnggPiBiZXN0Q2hyIHx8IGZtYXggPiBiZXN0UG9zKSkgfHwgKChkaXIgPiAwKSAmJiAoY2hyeCA8IGJlc3RDaHIgfHwgZm1pbiA8IGJlc3RQb3MpKSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgYmVzdEZlYXR1cmUgPSBmO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgYmVzdFBvcyA9IChkaXIgPCAwKSA/IGZtYXggOiBmbWluO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgYmVzdENociA9IGNocng7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoYmVzdEZlYXR1cmUgIT0gbnVsbCkgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhbYmVzdEZlYXR1cmVdKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBlbHNlXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhbXSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9KTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgfVxuXG4gICAgdmFyIGNpckZvYlJlY3VyMiA9IGZ1bmN0aW9uKGNpckJsb2NrRGF0YSwgb2Zmc2V0LCBsZXZlbCkge1xuICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShjaXJCbG9ja0RhdGEpO1xuICAgICAgICB2YXIgc2EgPSBuZXcgSW50MTZBcnJheShjaXJCbG9ja0RhdGEpO1xuICAgICAgICB2YXIgbGEgPSBuZXcgSW50MzJBcnJheShjaXJCbG9ja0RhdGEpO1xuXG4gICAgICAgIHZhciBpc0xlYWYgPSBiYVtvZmZzZXRdO1xuICAgICAgICB2YXIgY250ID0gc2Fbb2Zmc2V0LzIgKyAxXTtcbiAgICAgICAgb2Zmc2V0ICs9IDQ7XG5cbiAgICAgICAgaWYgKGlzTGVhZiAhPSAwKSB7XG4gICAgICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGNudDsgKytpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGxvID0gb2Zmc2V0LzQ7XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0Q2hyb20gPSBsYVtsb107XG4gICAgICAgICAgICAgICAgdmFyIHN0YXJ0QmFzZSA9IGxhW2xvICsgMV07XG4gICAgICAgICAgICAgICAgdmFyIGVuZENocm9tID0gbGFbbG8gKyAyXTtcbiAgICAgICAgICAgICAgICB2YXIgZW5kQmFzZSA9IGxhW2xvICsgM107XG4gICAgICAgICAgICAgICAgdmFyIGJsb2NrT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIG9mZnNldCsxNik7XG4gICAgICAgICAgICAgICAgdmFyIGJsb2NrU2l6ZSA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQrMjQpO1xuICAgICAgICAgICAgICAgIGlmICgoZGlyIDwgMCAmJiAoKHN0YXJ0Q2hyb20gPCBjaHIgfHwgKHN0YXJ0Q2hyb20gPT0gY2hyICYmIHN0YXJ0QmFzZSA8PSBwb3MpKSkpIHx8XG4gICAgICAgICAgICAgICAgICAgIChkaXIgPiAwICYmICgoZW5kQ2hyb20gPiBjaHIgfHwgKGVuZENocm9tID09IGNociAmJiBlbmRCYXNlID49IHBvcykpKSkpXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAvLyBjb25zb2xlLmxvZygnR290IGFuIGludGVyZXN0aW5nIGJsb2NrOiBzdGFydEJhc2U9JyArIHN0YXJ0Q2hyb20gKyAnOicgKyBzdGFydEJhc2UgKyAnOyBlbmRCYXNlPScgKyBlbmRDaHJvbSArICc6JyArIGVuZEJhc2UgKyAnOyBvZmZzZXQ9JyArIGJsb2NrT2Zmc2V0ICsgJzsgc2l6ZT0nICsgYmxvY2tTaXplKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKC9fcmFuZG9tLy5leGVjKHRoaXNCLmJ3Zy5pZHNUb0Nocm9tc1tzdGFydENocm9tXSkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vIGRsb2coJ3NraXBwaW5nIHJhbmRvbTogJyArIHRoaXNCLmJ3Zy5pZHNUb0Nocm9tc1tzdGFydENocm9tXSk7XG4gICAgICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAoYmxvY2tUb0ZldGNoID09IG51bGwgfHwgKChkaXIgPCAwKSAmJiAoZW5kQ2hyb20gPiBiZXN0QmxvY2tDaHIgfHwgKGVuZENocm9tID09IGJlc3RCbG9ja0NociAmJiBlbmRCYXNlID4gYmVzdEJsb2NrT2Zmc2V0KSkgfHxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAoZGlyID4gMCkgJiYgKHN0YXJ0Q2hyb20gPCBiZXN0QmxvY2tDaHIgfHwgKHN0YXJ0Q2hyb20gPT0gYmVzdEJsb2NrQ2hyICYmIHN0YXJ0QmFzZSA8IGJlc3RCbG9ja09mZnNldCkpKSlcbiAgICAgICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICAgICAgLy8gICAgICAgICAgICAgICAgICAgICAgICBkbG9nKCdiZXN0IGlzOiBzdGFydEJhc2U9JyArIHN0YXJ0Q2hyb20gKyAnOicgKyBzdGFydEJhc2UgKyAnOyBlbmRCYXNlPScgKyBlbmRDaHJvbSArICc6JyArIGVuZEJhc2UgKyAnOyBvZmZzZXQ9JyArIGJsb2NrT2Zmc2V0ICsgJzsgc2l6ZT0nICsgYmxvY2tTaXplKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGJsb2NrVG9GZXRjaCA9IHtvZmZzZXQ6IGJsb2NrT2Zmc2V0LCBzaXplOiBibG9ja1NpemV9O1xuICAgICAgICAgICAgICAgICAgICAgICAgYmVzdEJsb2NrT2Zmc2V0ID0gKGRpciA8IDApID8gZW5kQmFzZSA6IHN0YXJ0QmFzZTtcbiAgICAgICAgICAgICAgICAgICAgICAgIGJlc3RCbG9ja0NociA9IChkaXIgPCAwKSA/IGVuZENocm9tIDogc3RhcnRDaHJvbTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBvZmZzZXQgKz0gMzI7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2YXIgYmVzdFJlY3VyID0gLTE7XG4gICAgICAgICAgICB2YXIgYmVzdFBvcyA9IC0xO1xuICAgICAgICAgICAgdmFyIGJlc3RDaHIgPSAtMTtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgY250OyArK2kpIHtcbiAgICAgICAgICAgICAgICB2YXIgbG8gPSBvZmZzZXQvNDtcbiAgICAgICAgICAgICAgICB2YXIgc3RhcnRDaHJvbSA9IGxhW2xvXTtcbiAgICAgICAgICAgICAgICB2YXIgc3RhcnRCYXNlID0gbGFbbG8gKyAxXTtcbiAgICAgICAgICAgICAgICB2YXIgZW5kQ2hyb20gPSBsYVtsbyArIDJdO1xuICAgICAgICAgICAgICAgIHZhciBlbmRCYXNlID0gbGFbbG8gKyAzXTtcbiAgICAgICAgICAgICAgICB2YXIgYmxvY2tPZmZzZXQgPSAobGFbbG8gKyA0XTw8MzIpIHwgKGxhW2xvICsgNV0pO1xuICAgICAgICAgICAgICAgIGlmICgoZGlyIDwgMCAmJiAoKHN0YXJ0Q2hyb20gPCBjaHIgfHwgKHN0YXJ0Q2hyb20gPT0gY2hyICYmIHN0YXJ0QmFzZSA8PSBwb3MpKSAmJlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgKGVuZENocm9tICAgPj0gY2hyKSkpIHx8XG4gICAgICAgICAgICAgICAgICAgICAoZGlyID4gMCAmJiAoKGVuZENocm9tID4gY2hyIHx8IChlbmRDaHJvbSA9PSBjaHIgJiYgZW5kQmFzZSA+PSBwb3MpKSAmJlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIChzdGFydENocm9tIDw9IGNocikpKSlcbiAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgIGlmIChiZXN0UmVjdXIgPCAwIHx8IGVuZEJhc2UgPiBiZXN0UG9zKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBiZXN0UmVjdXIgPSBibG9ja09mZnNldDtcbiAgICAgICAgICAgICAgICAgICAgICAgIGJlc3RQb3MgPSAoZGlyIDwgMCkgPyBlbmRCYXNlIDogc3RhcnRCYXNlO1xuICAgICAgICAgICAgICAgICAgICAgICAgYmVzdENociA9IChkaXIgPCAwKSA/IGVuZENocm9tIDogc3RhcnRDaHJvbTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBvZmZzZXQgKz0gMjQ7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBpZiAoYmVzdFJlY3VyID49IDApIHtcbiAgICAgICAgICAgICAgICBjaXJGb2JSZWN1cihbYmVzdFJlY3VyXSwgbGV2ZWwgKyAxKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH07XG4gICAgXG5cbiAgICBjaXJGb2JSZWN1cihbdGhpc0IuY2lyVHJlZU9mZnNldCArIDQ4XSwgMSk7XG59XG5cbkJpZ1dpZy5wcm90b3R5cGUucmVhZFdpZ0RhdGEgPSBmdW5jdGlvbihjaHJOYW1lLCBtaW4sIG1heCwgY2FsbGJhY2spIHtcbiAgICB0aGlzLmdldFVuem9vbWVkVmlldygpLnJlYWRXaWdEYXRhKGNock5hbWUsIG1pbiwgbWF4LCBjYWxsYmFjayk7XG59XG5cbkJpZ1dpZy5wcm90b3R5cGUuZ2V0VW56b29tZWRWaWV3ID0gZnVuY3Rpb24oKSB7XG4gICAgaWYgKCF0aGlzLnVuem9vbWVkVmlldykge1xuICAgICAgICB2YXIgY2lyTGVuID0gNDAwMDtcbiAgICAgICAgdmFyIG56bCA9IHRoaXMuem9vbUxldmVsc1swXTtcbiAgICAgICAgaWYgKG56bCkge1xuICAgICAgICAgICAgY2lyTGVuID0gdGhpcy56b29tTGV2ZWxzWzBdLmRhdGFPZmZzZXQgLSB0aGlzLnVuem9vbWVkSW5kZXhPZmZzZXQ7XG4gICAgICAgIH1cbiAgICAgICAgdGhpcy51bnpvb21lZFZpZXcgPSBuZXcgQmlnV2lnVmlldyh0aGlzLCB0aGlzLnVuem9vbWVkSW5kZXhPZmZzZXQsIGNpckxlbiwgZmFsc2UpO1xuICAgIH1cbiAgICByZXR1cm4gdGhpcy51bnpvb21lZFZpZXc7XG59XG5cbkJpZ1dpZy5wcm90b3R5cGUuZ2V0Wm9vbWVkVmlldyA9IGZ1bmN0aW9uKHopIHtcbiAgICB2YXIgemggPSB0aGlzLnpvb21MZXZlbHNbel07XG4gICAgaWYgKCF6aC52aWV3KSB7XG4gICAgICAgIHpoLnZpZXcgPSBuZXcgQmlnV2lnVmlldyh0aGlzLCB6aC5pbmRleE9mZnNldCwgLyogdGhpcy56b29tTGV2ZWxzW3ogKyAxXS5kYXRhT2Zmc2V0IC0gemguaW5kZXhPZmZzZXQgKi8gNDAwMCwgdHJ1ZSk7XG4gICAgfVxuICAgIHJldHVybiB6aC52aWV3O1xufVxuXG5mdW5jdGlvbiBtYWtlQndnKGRhdGEsIGNhbGxiYWNrLCBuYW1lKSB7XG4gICAgdmFyIGJ3ZyA9IG5ldyBCaWdXaWcoKTtcbiAgICBid2cuZGF0YSA9IGRhdGE7XG4gICAgYndnLm5hbWUgPSBuYW1lO1xuICAgIGJ3Zy5kYXRhLnNsaWNlKDAsIDUxMikuc2FsdGVkKCkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0KSB7XG4gICAgICAgIGlmICghcmVzdWx0KSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgXCJDb3VsZG4ndCBmZXRjaCBmaWxlXCIpO1xuICAgICAgICB9XG5cbiAgICAgICAgdmFyIGhlYWRlciA9IHJlc3VsdDtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoaGVhZGVyKTtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoaGVhZGVyKTtcbiAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoaGVhZGVyKTtcbiAgICAgICAgdmFyIG1hZ2ljID0gYmFbMF0gKyAoTTEgKiBiYVsxXSkgKyAoTTIgKiBiYVsyXSkgKyAoTTMgKiBiYVszXSk7XG4gICAgICAgIGlmIChtYWdpYyA9PSBCSUdfV0lHX01BR0lDKSB7XG4gICAgICAgICAgICBid2cudHlwZSA9ICdiaWd3aWcnO1xuICAgICAgICB9IGVsc2UgaWYgKG1hZ2ljID09IEJJR19CRURfTUFHSUMpIHtcbiAgICAgICAgICAgIGJ3Zy50eXBlID0gJ2JpZ2JlZCc7XG4gICAgICAgIH0gZWxzZSBpZiAobWFnaWMgPT0gQklHX1dJR19NQUdJQ19CRSB8fCBtYWdpYyA9PSBCSUdfQkVEX01BR0lDX0JFKSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCwgXCJDdXJyZW50bHkgZG9uJ3Qgc3VwcG9ydCBiaWctZW5kaWFuIEJCSSBmaWxlc1wiKTtcbiAgICAgICAgICAgIFxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsIFwiTm90IGEgc3VwcG9ydGVkIGZvcm1hdCwgbWFnaWM9MHhcIiArIG1hZ2ljLnRvU3RyaW5nKDE2KSk7XG4gICAgICAgICAgICBcbiAgICAgICAgfVxuXG4gICAgICAgIGJ3Zy52ZXJzaW9uID0gc2FbMl07ICAgICAgICAgICAgIC8vIDRcbiAgICAgICAgYndnLm51bVpvb21MZXZlbHMgPSBzYVszXTsgICAgICAgLy8gNlxuICAgICAgICBid2cuY2hyb21UcmVlT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIDgpO1xuICAgICAgICBid2cudW56b29tZWREYXRhT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIDE2KTtcbiAgICAgICAgYndnLnVuem9vbWVkSW5kZXhPZmZzZXQgPSBid2dfcmVhZE9mZnNldChiYSwgMjQpO1xuICAgICAgICBid2cuZmllbGRDb3VudCA9IHNhWzE2XTsgICAgICAgICAvLyAzMlxuICAgICAgICBid2cuZGVmaW5lZEZpZWxkQ291bnQgPSBzYVsxN107ICAvLyAzNFxuICAgICAgICBid2cuYXNPZmZzZXQgPSBid2dfcmVhZE9mZnNldChiYSwgMzYpO1xuICAgICAgICBid2cudG90YWxTdW1tYXJ5T2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIDQ0KTtcbiAgICAgICAgYndnLnVuY29tcHJlc3NCdWZTaXplID0gbGFbMTNdOyAgLy8gNTJcbiAgICAgICAgYndnLmV4dEhlYWRlck9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCA1Nik7XG5cbiAgICAgICAgYndnLnpvb21MZXZlbHMgPSBbXTtcbiAgICAgICAgZm9yICh2YXIgemwgPSAwOyB6bCA8IGJ3Zy5udW1ab29tTGV2ZWxzOyArK3psKSB7XG4gICAgICAgICAgICB2YXIgemxSZWR1Y3Rpb24gPSBsYVt6bCo2ICsgMTZdXG4gICAgICAgICAgICB2YXIgemxEYXRhID0gYndnX3JlYWRPZmZzZXQoYmEsIHpsKjI0ICsgNzIpO1xuICAgICAgICAgICAgdmFyIHpsSW5kZXggPSBid2dfcmVhZE9mZnNldChiYSwgemwqMjQgKyA4MCk7XG4gICAgICAgICAgICBid2cuem9vbUxldmVscy5wdXNoKHtyZWR1Y3Rpb246IHpsUmVkdWN0aW9uLCBkYXRhT2Zmc2V0OiB6bERhdGEsIGluZGV4T2Zmc2V0OiB6bEluZGV4fSk7XG4gICAgICAgIH1cblxuICAgICAgICBid2cucmVhZENocm9tVHJlZShmdW5jdGlvbigpIHtcbiAgICAgICAgICAgIGJ3Zy5nZXRBdXRvU1FMKGZ1bmN0aW9uKGFzKSB7XG4gICAgICAgICAgICAgICAgYndnLnNjaGVtYSA9IGFzO1xuICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhid2cpO1xuICAgICAgICAgICAgfSk7XG4gICAgICAgIH0pO1xuICAgIH0sIHt0aW1lb3V0OiA1MDAwfSk7ICAgIC8vIFBvdGVudGlhbCB0aW1lb3V0IG9uIGZpcnN0IHJlcXVlc3QgdG8gY2F0Y2ggbWl4ZWQtY29udGVudCBlcnJvcnMgb25cbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAvLyBDaHJvbWl1bS5cbn1cblxuXG5CaWdXaWcucHJvdG90eXBlLl90c0ZldGNoID0gZnVuY3Rpb24oem9vbSwgY2hyLCBtaW4sIG1heCwgY2FsbGJhY2spIHtcbiAgICB2YXIgYndnID0gdGhpcztcbiAgICBpZiAoem9vbSA+PSB0aGlzLnpvb21MZXZlbHMubGVuZ3RoIC0gMSkge1xuICAgICAgICBpZiAoIXRoaXMudG9wTGV2ZWxSZWR1Y3Rpb25DYWNoZSkge1xuICAgICAgICAgICAgdGhpcy5nZXRab29tZWRWaWV3KHRoaXMuem9vbUxldmVscy5sZW5ndGggLSAxKS5yZWFkV2lnRGF0YUJ5SWQoLTEsIDAsIDMwMDAwMDAwMCwgZnVuY3Rpb24oZmVhdHMpIHtcbiAgICAgICAgICAgICAgICBid2cudG9wTGV2ZWxSZWR1Y3Rpb25DYWNoZSA9IGZlYXRzO1xuICAgICAgICAgICAgICAgIHJldHVybiBid2cuX3RzRmV0Y2goem9vbSwgY2hyLCBtaW4sIG1heCwgY2FsbGJhY2spO1xuICAgICAgICAgICAgfSk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2YXIgZiA9IFtdO1xuICAgICAgICAgICAgdmFyIGMgPSB0aGlzLnRvcExldmVsUmVkdWN0aW9uQ2FjaGU7XG4gICAgICAgICAgICBmb3IgKHZhciBmaSA9IDA7IGZpIDwgYy5sZW5ndGg7ICsrZmkpIHtcbiAgICAgICAgICAgICAgICBpZiAoY1tmaV0uX2Nocm9tSWQgPT0gY2hyKSB7XG4gICAgICAgICAgICAgICAgICAgIGYucHVzaChjW2ZpXSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGYpO1xuICAgICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIHZpZXc7XG4gICAgICAgIGlmICh6b29tIDwgMCkge1xuICAgICAgICAgICAgdmlldyA9IHRoaXMuZ2V0VW56b29tZWRWaWV3KCk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB2aWV3ID0gdGhpcy5nZXRab29tZWRWaWV3KHpvb20pO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiB2aWV3LnJlYWRXaWdEYXRhQnlJZChjaHIsIG1pbiwgbWF4LCBjYWxsYmFjayk7XG4gICAgfVxufVxuXG5CaWdXaWcucHJvdG90eXBlLnRocmVzaG9sZFNlYXJjaCA9IGZ1bmN0aW9uKGNock5hbWUsIHJlZmVyZW5jZVBvaW50LCBkaXIsIHRocmVzaG9sZCwgY2FsbGJhY2spIHtcbiAgICBkaXIgPSAoZGlyPDApID8gLTEgOiAxO1xuICAgIHZhciBid2cgPSB0aGlzO1xuICAgIHZhciBpbml0aWFsQ2hyID0gdGhpcy5jaHJvbXNUb0lEc1tjaHJOYW1lXTtcbiAgICB2YXIgY2FuZGlkYXRlcyA9IFt7Y2hyT3JkOiAwLCBjaHI6IGluaXRpYWxDaHIsIHpvb206IGJ3Zy56b29tTGV2ZWxzLmxlbmd0aCAtIDQsIG1pbjogMCwgbWF4OiAzMDAwMDAwMDAsIGZyb21SZWY6IHRydWV9XVxuICAgIGZvciAodmFyIGkgPSAxOyBpIDw9IHRoaXMubWF4SUQgKyAxOyArK2kpIHtcbiAgICAgICAgdmFyIGNocklkID0gKGluaXRpYWxDaHIgKyAoZGlyKmkpKSAlICh0aGlzLm1heElEICsgMSk7XG4gICAgICAgIGlmIChjaHJJZCA8IDApIFxuICAgICAgICAgICAgY2hySWQgKz0gKHRoaXMubWF4SUQgKyAxKTtcbiAgICAgICAgY2FuZGlkYXRlcy5wdXNoKHtjaHJPcmQ6IGksIGNocjogY2hySWQsIHpvb206IGJ3Zy56b29tTGV2ZWxzLmxlbmd0aCAtIDEsIG1pbjogMCwgbWF4OiAzMDAwMDAwMDB9KVxuICAgIH1cbiAgICAgICBcbiAgICBmdW5jdGlvbiBmYlRocmVzaG9sZFNlYXJjaFJlY3VyKCkge1xuICAgIFx0aWYgKGNhbmRpZGF0ZXMubGVuZ3RoID09IDApIHtcbiAgICBcdCAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgXHR9XG4gICAgXHRjYW5kaWRhdGVzLnNvcnQoZnVuY3Rpb24oYzEsIGMyKSB7XG4gICAgXHQgICAgdmFyIGQgPSBjMS56b29tIC0gYzIuem9vbTtcbiAgICBcdCAgICBpZiAoZCAhPSAwKVxuICAgIFx0XHQgICAgcmV0dXJuIGQ7XG5cbiAgICAgICAgICAgIGQgPSBjMS5jaHJPcmQgLSBjMi5jaHJPcmQ7XG4gICAgICAgICAgICBpZiAoZCAhPSAwKVxuICAgICAgICAgICAgICAgIHJldHVybiBkO1xuICAgIFx0ICAgIGVsc2VcbiAgICBcdFx0ICAgIHJldHVybiBjMS5taW4gLSBjMi5taW4gKiBkaXI7XG4gICAgXHR9KTtcblxuXHQgICAgdmFyIGNhbmRpZGF0ZSA9IGNhbmRpZGF0ZXMuc3BsaWNlKDAsIDEpWzBdO1xuICAgICAgICBid2cuX3RzRmV0Y2goY2FuZGlkYXRlLnpvb20sIGNhbmRpZGF0ZS5jaHIsIGNhbmRpZGF0ZS5taW4sIGNhbmRpZGF0ZS5tYXgsIGZ1bmN0aW9uKGZlYXRzKSB7XG4gICAgICAgICAgICB2YXIgcnAgPSBkaXIgPiAwID8gMCA6IDMwMDAwMDAwMDtcbiAgICAgICAgICAgIGlmIChjYW5kaWRhdGUuZnJvbVJlZilcbiAgICAgICAgICAgICAgICBycCA9IHJlZmVyZW5jZVBvaW50O1xuICAgICAgICAgICAgXG4gICAgICAgICAgICBmb3IgKHZhciBmaSA9IDA7IGZpIDwgZmVhdHMubGVuZ3RoOyArK2ZpKSB7XG4gICAgXHQgICAgICAgIHZhciBmID0gZmVhdHNbZmldO1xuICAgICAgICAgICAgICAgIHZhciBzY29yZTtcbiAgICAgICAgICAgICAgICBpZiAoZi5tYXhTY29yZSAhPSB1bmRlZmluZWQpXG4gICAgICAgICAgICAgICAgICAgIHNjb3JlID0gZi5tYXhTY29yZTtcbiAgICAgICAgICAgICAgICBlbHNlXG4gICAgICAgICAgICAgICAgICAgIHNjb3JlID0gZi5zY29yZTtcblxuICAgICAgICAgICAgICAgIGlmIChkaXIgPiAwKSB7XG4gICAgXHQgICAgICAgICAgICBpZiAoc2NvcmUgPiB0aHJlc2hvbGQpIHtcbiAgICAgICAgXHRcdCAgICAgICAgaWYgKGNhbmRpZGF0ZS56b29tIDwgMCkge1xuICAgICAgICBcdFx0ICAgICAgICAgICAgaWYgKGYubWluID4gcnApXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhmKTtcbiAgICAgICAgXHRcdCAgICAgICAgfSBlbHNlIGlmIChmLm1heCA+IHJwKSB7XG4gICAgICAgIFx0XHQgICAgICAgICAgICBjYW5kaWRhdGVzLnB1c2goe2NocjogY2FuZGlkYXRlLmNociwgY2hyT3JkOiBjYW5kaWRhdGUuY2hyT3JkLCB6b29tOiBjYW5kaWRhdGUuem9vbSAtIDIsIG1pbjogZi5taW4sIG1heDogZi5tYXgsIGZyb21SZWY6IGNhbmRpZGF0ZS5mcm9tUmVmfSk7XG4gICAgICAgIFx0XHQgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIGlmIChzY29yZSA+IHRocmVzaG9sZCkge1xuICAgICAgICAgICAgXHRcdCAgICBpZiAoY2FuZGlkYXRlLnpvb20gPCAwKSB7XG4gICAgICAgICAgICAgICAgXHQgICAgICAgIGlmIChmLm1heCA8IHJwKVxuICAgICAgICAgICAgICAgIFx0XHRcdCAgICByZXR1cm4gY2FsbGJhY2soZik7XG4gICAgICAgICAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKGYubWluIDwgcnApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBjYW5kaWRhdGVzLnB1c2goe2NocjogY2FuZGlkYXRlLmNociwgY2hyT3JkOiBjYW5kaWRhdGUuY2hyT3JkLCB6b29tOiBjYW5kaWRhdGUuem9vbSAtIDIsIG1pbjogZi5taW4sIG1heDogZi5tYXgsIGZyb21SZWY6IGNhbmRpZGF0ZS5mcm9tUmVmfSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgXHQgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgIFx0ICAgIH1cbiAgICAgICAgICAgIGZiVGhyZXNob2xkU2VhcmNoUmVjdXIoKTtcbiAgICAgICAgfSk7XG4gICAgfVxuICAgIFxuICAgIGZiVGhyZXNob2xkU2VhcmNoUmVjdXIoKTtcbn1cblxuQmlnV2lnLnByb3RvdHlwZS5nZXRBdXRvU1FMID0gZnVuY3Rpb24oY2FsbGJhY2spIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuICAgIGlmICghdGhpcy5hc09mZnNldClcbiAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwpO1xuXG5cbiAgICB0aGlzLmRhdGEuc2xpY2UodGhpcy5hc09mZnNldCwgMjA0OCkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0KSB7XG4gICAgICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KHJlc3VsdCk7XG4gICAgICAgIHZhciBzID0gJyc7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgYmEubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgIGlmIChiYVtpXSA9PSAwKVxuICAgICAgICAgICAgICAgIGJyZWFrO1xuICAgICAgICAgICAgcyArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGJhW2ldKTtcbiAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgLyogXG4gICAgICAgICAqIFF1aWNrJ24nZGlydHkgYXR0ZW1wdCB0byBwYXJzZSBhdXRvU3FsIGZvcm1hdC5cbiAgICAgICAgICogU2VlOiBodHRwOi8vd3d3LmxpbnV4am91cm5hbC5jb20vZmlsZXMvbGludXhqb3VybmFsLmNvbS9saW51eGpvdXJuYWwvYXJ0aWNsZXMvMDU5LzU5NDkvNTk0OWwyLmh0bWxcbiAgICAgICAgICovXG5cbiAgICAgICAgdmFyIGhlYWRlcl9yZSA9IC8oXFx3KylcXHMrKFxcdyspXFxzKyhcIihbXlwiXSspXCIpP1xccytcXChcXHMqLztcbiAgICAgICAgdmFyIGZpZWxkX3JlID0gLyhbXFx3XFxbXFxdXSspXFxzKyhcXHcrKVxccyo7XFxzKihcIihbXlwiXSspXCIpP1xccyovZztcblxuICAgICAgICB2YXIgaGVhZGVyTWF0Y2ggPSBoZWFkZXJfcmUuZXhlYyhzKTtcbiAgICAgICAgaWYgKGhlYWRlck1hdGNoKSB7XG4gICAgICAgICAgICB2YXIgYXMgPSB7XG4gICAgICAgICAgICAgICAgZGVjbFR5cGU6IGhlYWRlck1hdGNoWzFdLFxuICAgICAgICAgICAgICAgIG5hbWU6IGhlYWRlck1hdGNoWzJdLFxuICAgICAgICAgICAgICAgIGNvbW1lbnQ6IGhlYWRlck1hdGNoWzRdLFxuXG4gICAgICAgICAgICAgICAgZmllbGRzOiBbXVxuICAgICAgICAgICAgfTtcblxuICAgICAgICAgICAgcyA9IHMuc3Vic3RyaW5nKGhlYWRlck1hdGNoWzBdKTtcbiAgICAgICAgICAgIGZvciAodmFyIG0gPSBmaWVsZF9yZS5leGVjKHMpOyBtICE9IG51bGw7IG0gPSBmaWVsZF9yZS5leGVjKHMpKSB7XG4gICAgICAgICAgICAgICAgYXMuZmllbGRzLnB1c2goe3R5cGU6IG1bMV0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgIG5hbWU6IG1bMl0sXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgIGNvbW1lbnQ6IG1bNF19KTtcbiAgICAgICAgICAgIH1cblxuICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGFzKTtcbiAgICAgICAgfVxuICAgIH0pO1xufVxuXG5CaWdXaWcucHJvdG90eXBlLmdldEV4dHJhSW5kaWNlcyA9IGZ1bmN0aW9uKGNhbGxiYWNrKSB7XG4gICAgdmFyIHRoaXNCID0gdGhpcztcbiAgICBpZiAodGhpcy52ZXJzaW9uIDwgNCB8fCB0aGlzLmV4dEhlYWRlck9mZnNldCA9PSAwIHx8IHRoaXMudHlwZSAhPSAnYmlnYmVkJykge1xuICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdGhpcy5kYXRhLnNsaWNlKHRoaXMuZXh0SGVhZGVyT2Zmc2V0LCA2NCkuZmV0Y2goZnVuY3Rpb24ocmVzdWx0KSB7XG4gICAgICAgICAgICBpZiAoIXJlc3VsdCkge1xuICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsLCBcIkNvdWxkbid0IGZldGNoIGV4dGVuc2lvbiBoZWFkZXJcIik7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KHJlc3VsdCk7XG4gICAgICAgICAgICB2YXIgc2EgPSBuZXcgSW50MTZBcnJheShyZXN1bHQpO1xuICAgICAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkocmVzdWx0KTtcbiAgICAgICAgICAgIFxuICAgICAgICAgICAgdmFyIGV4dEhlYWRlclNpemUgPSBzYVswXTtcbiAgICAgICAgICAgIHZhciBleHRyYUluZGV4Q291bnQgPSBzYVsxXTtcbiAgICAgICAgICAgIHZhciBleHRyYUluZGV4TGlzdE9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCA0KTtcblxuICAgICAgICAgICAgaWYgKGV4dHJhSW5kZXhDb3VudCA9PSAwKSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwpO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAvLyBGSVhNRSAyMGJ5dGUgcmVjb3JkcyBvbmx5IG1ha2Ugc2Vuc2UgZm9yIHNpbmdsZS1maWVsZCBpbmRpY2VzLlxuICAgICAgICAgICAgLy8gUmlnaHQgbm93LCB0aGVzZSBzZWVtIHRvIGJlIHRoZSBvbmx5IHRoaW5ncyBhcm91bmQsIGJ1dCB0aGUgZm9ybWF0XG4gICAgICAgICAgICAvLyBpcyBhY3R1YWxseSBtb3JlIGdlbmVyYWwuXG4gICAgICAgICAgICB0aGlzQi5kYXRhLnNsaWNlKGV4dHJhSW5kZXhMaXN0T2Zmc2V0LCBleHRyYUluZGV4Q291bnQgKiAyMCkuZmV0Y2goZnVuY3Rpb24oZWlsKSB7XG4gICAgICAgICAgICAgICAgaWYgKCFlaWwpIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsIFwiQ291bGRuJ3QgZmV0Y2ggaW5kZXggaW5mb1wiKTtcbiAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShlaWwpO1xuICAgICAgICAgICAgICAgIHZhciBzYSA9IG5ldyBJbnQxNkFycmF5KGVpbCk7XG4gICAgICAgICAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoZWlsKTtcblxuICAgICAgICAgICAgICAgIHZhciBpbmRpY2VzID0gW107XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgaWkgPSAwOyBpaSA8IGV4dHJhSW5kZXhDb3VudDsgKytpaSkge1xuICAgICAgICAgICAgICAgICAgICB2YXIgZWlUeXBlID0gc2FbaWkqMTBdO1xuICAgICAgICAgICAgICAgICAgICB2YXIgZWlGaWVsZENvdW50ID0gc2FbaWkqMTAgKyAxXTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGVpT2Zmc2V0ID0gYndnX3JlYWRPZmZzZXQoYmEsIGlpKjIwICsgNCk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBlaUZpZWxkID0gc2FbaWkqMTAgKyA4XVxuICAgICAgICAgICAgICAgICAgICB2YXIgaW5kZXggPSBuZXcgQkJJRXh0cmFJbmRleCh0aGlzQiwgZWlUeXBlLCBlaUZpZWxkQ291bnQsIGVpT2Zmc2V0LCBlaUZpZWxkKTtcbiAgICAgICAgICAgICAgICAgICAgaW5kaWNlcy5wdXNoKGluZGV4KTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgY2FsbGJhY2soaW5kaWNlcyk7XG4gICAgICAgICAgICB9KTtcbiAgICAgICAgfSk7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBCQklFeHRyYUluZGV4KGJiaSwgdHlwZSwgZmllbGRDb3VudCwgb2Zmc2V0LCBmaWVsZCkge1xuICAgIHRoaXMuYmJpID0gYmJpO1xuICAgIHRoaXMudHlwZSA9IHR5cGU7XG4gICAgdGhpcy5maWVsZENvdW50ID0gZmllbGRDb3VudDtcbiAgICB0aGlzLm9mZnNldCA9IG9mZnNldDtcbiAgICB0aGlzLmZpZWxkID0gZmllbGQ7XG59XG5cbkJCSUV4dHJhSW5kZXgucHJvdG90eXBlLmxvb2t1cCA9IGZ1bmN0aW9uKG5hbWUsIGNhbGxiYWNrKSB7XG4gICAgdmFyIHRoaXNCID0gdGhpcztcblxuICAgIHRoaXMuYmJpLmRhdGEuc2xpY2UodGhpcy5vZmZzZXQsIDMyKS5mZXRjaChmdW5jdGlvbihicHQpIHtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIHNhID0gbmV3IEludDE2QXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIGxhID0gbmV3IEludDMyQXJyYXkoYnB0KTtcbiAgICAgICAgdmFyIGJwdE1hZ2ljID0gbGFbMF07XG4gICAgICAgIHZhciBibG9ja1NpemUgPSBsYVsxXTtcbiAgICAgICAgdmFyIGtleVNpemUgPSBsYVsyXTtcbiAgICAgICAgdmFyIHZhbFNpemUgPSBsYVszXTtcbiAgICAgICAgdmFyIGl0ZW1Db3VudCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCAxNik7XG4gICAgICAgIHZhciByb290Tm9kZU9mZnNldCA9IDMyO1xuXG4gICAgICAgIGZ1bmN0aW9uIGJwdFJlYWROb2RlKG5vZGVPZmZzZXQpIHtcbiAgICAgICAgICAgIHRoaXNCLmJiaS5kYXRhLnNsaWNlKG5vZGVPZmZzZXQsIDQgKyAoYmxvY2tTaXplICogKGtleVNpemUgKyB2YWxTaXplKSkpLmZldGNoKGZ1bmN0aW9uKG5vZGUpIHtcbiAgICAgICAgICAgICAgICB2YXIgYmEgPSBuZXcgVWludDhBcnJheShub2RlKTtcbiAgICAgICAgICAgICAgICB2YXIgc2EgPSBuZXcgVWludDE2QXJyYXkobm9kZSk7XG4gICAgICAgICAgICAgICAgdmFyIGxhID0gbmV3IFVpbnQzMkFycmF5KG5vZGUpO1xuXG4gICAgICAgICAgICAgICAgdmFyIG5vZGVUeXBlID0gYmFbMF07XG4gICAgICAgICAgICAgICAgdmFyIGNudCA9IHNhWzFdO1xuXG4gICAgICAgICAgICAgICAgdmFyIG9mZnNldCA9IDQ7XG4gICAgICAgICAgICAgICAgaWYgKG5vZGVUeXBlID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGxhc3RDaGlsZE9mZnNldCA9IG51bGw7XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIG4gPSAwOyBuIDwgY250OyArK24pIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBrZXkgPSAnJztcbiAgICAgICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGtpID0gMDsga2kgPCBrZXlTaXplOyArK2tpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGNoYXJDb2RlID0gYmFbb2Zmc2V0KytdO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmIChjaGFyQ29kZSAhPSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGtleSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGNoYXJDb2RlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG5cbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBjaGlsZE9mZnNldCA9IGJ3Z19yZWFkT2Zmc2V0KGJhLCBvZmZzZXQpO1xuICAgICAgICAgICAgICAgICAgICAgICAgb2Zmc2V0ICs9IDg7XG4gICAgICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChuYW1lLmxvY2FsZUNvbXBhcmUoa2V5KSA8IDAgJiYgbGFzdENoaWxkT2Zmc2V0KSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgYnB0UmVhZE5vZGUobGFzdENoaWxkT2Zmc2V0KTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm47XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBsYXN0Q2hpbGRPZmZzZXQgPSBjaGlsZE9mZnNldDtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBicHRSZWFkTm9kZShsYXN0Q2hpbGRPZmZzZXQpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIGZvciAodmFyIG4gPSAwOyBuIDwgY250OyArK24pIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBrZXkgPSAnJztcbiAgICAgICAgICAgICAgICAgICAgICAgIGZvciAodmFyIGtpID0gMDsga2kgPCBrZXlTaXplOyArK2tpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGNoYXJDb2RlID0gYmFbb2Zmc2V0KytdO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGlmIChjaGFyQ29kZSAhPSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGtleSArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKGNoYXJDb2RlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vIFNwZWNpZmljIGZvciBFSSBjYXNlLlxuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGtleSA9PSBuYW1lKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHN0YXJ0ID0gYndnX3JlYWRPZmZzZXQoYmEsIG9mZnNldCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGxlbmd0aCA9IHJlYWRJbnQoYmEsIG9mZnNldCArIDgpO1xuXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIHRoaXNCLmJiaS5nZXRVbnpvb21lZFZpZXcoKS5mZXRjaEZlYXR1cmVzKFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBmdW5jdGlvbihjaHIsIG1pbiwgbWF4LCB0b2tzKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAodG9rcyAmJiB0b2tzLmxlbmd0aCA+IHRoaXNCLmZpZWxkIC0gMylcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gdG9rc1t0aGlzQi5maWVsZCAtIDNdID09IG5hbWU7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0sIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBbe29mZnNldDogc3RhcnQsIHNpemU6IGxlbmd0aH1dLCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY2FsbGJhY2spO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgb2Zmc2V0ICs9IHZhbFNpemU7XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKFtdKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9KTtcbiAgICAgICAgfVxuXG4gICAgICAgIGJwdFJlYWROb2RlKHRoaXNCLm9mZnNldCArIHJvb3ROb2RlT2Zmc2V0KTtcbiAgICB9KTtcbn1cblxuaWYgKHR5cGVvZihtb2R1bGUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIG1vZHVsZS5leHBvcnRzID0ge1xuICAgICAgICBtYWtlQndnOiBtYWtlQndnLFxuICAgICAgICBCSUdfQkVEX01BR0lDOiBCSUdfQkVEX01BR0lDLFxuICAgICAgICBCSUdfV0lHX01BR0lDOiBCSUdfV0lHX01BR0lDXG4gICAgfVxufVxuIiwiLyogLSotIG1vZGU6IGphdmFzY3JpcHQ7IGMtYmFzaWMtb2Zmc2V0OiA0OyBpbmRlbnQtdGFicy1tb2RlOiBuaWwgLSotICovXG5cbi8vIFxuLy8gRGFsbGlhbmNlIEdlbm9tZSBFeHBsb3JlclxuLy8gKGMpIFRob21hcyBEb3duIDIwMDYtMjAxMVxuLy9cbi8vIGJpbi5qcyBnZW5lcmFsIGJpbmFyeSBkYXRhIHN1cHBvcnRcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG5pZiAodHlwZW9mKHJlcXVpcmUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIHZhciB1dGlscyA9IHJlcXVpcmUoJy4vdXRpbHMnKTtcbiAgICB2YXIgc2hhbGxvd0NvcHkgPSB1dGlscy5zaGFsbG93Q29weTtcblxuICAgIHZhciBzaGExID0gcmVxdWlyZSgnLi9zaGExJyk7XG4gICAgdmFyIGI2NF9zaGExID0gc2hhMS5iNjRfc2hhMTtcblxuICAgIHZhciBQcm9taXNlID0gcmVxdWlyZSgnZXM2LXByb21pc2UnKS5Qcm9taXNlO1xufVxuXG5mdW5jdGlvbiBCbG9iRmV0Y2hhYmxlKGIpIHtcbiAgICB0aGlzLmJsb2IgPSBiO1xufVxuXG5CbG9iRmV0Y2hhYmxlLnByb3RvdHlwZS5zbGljZSA9IGZ1bmN0aW9uKHN0YXJ0LCBsZW5ndGgpIHtcbiAgICB2YXIgYjtcblxuICAgIGlmICh0aGlzLmJsb2Iuc2xpY2UpIHtcbiAgICAgICAgaWYgKGxlbmd0aCkge1xuICAgICAgICAgICAgYiA9IHRoaXMuYmxvYi5zbGljZShzdGFydCwgc3RhcnQgKyBsZW5ndGgpO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgYiA9IHRoaXMuYmxvYi5zbGljZShzdGFydCk7XG4gICAgICAgIH1cbiAgICB9IGVsc2Uge1xuICAgICAgICBpZiAobGVuZ3RoKSB7XG4gICAgICAgICAgICBiID0gdGhpcy5ibG9iLndlYmtpdFNsaWNlKHN0YXJ0LCBzdGFydCArIGxlbmd0aCk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBiID0gdGhpcy5ibG9iLndlYmtpdFNsaWNlKHN0YXJ0KTtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gbmV3IEJsb2JGZXRjaGFibGUoYik7XG59XG5cbkJsb2JGZXRjaGFibGUucHJvdG90eXBlLnNhbHRlZCA9IGZ1bmN0aW9uKCkge3JldHVybiB0aGlzO31cblxuaWYgKHR5cGVvZihGaWxlUmVhZGVyKSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgICAvLyBjb25zb2xlLmxvZygnZGVmaW5pbmcgYXN5bmMgQmxvYkZldGNoYWJsZS5mZXRjaCcpO1xuXG4gICAgQmxvYkZldGNoYWJsZS5wcm90b3R5cGUuZmV0Y2ggPSBmdW5jdGlvbihjYWxsYmFjaykge1xuICAgICAgICB2YXIgcmVhZGVyID0gbmV3IEZpbGVSZWFkZXIoKTtcbiAgICAgICAgcmVhZGVyLm9ubG9hZGVuZCA9IGZ1bmN0aW9uKGV2KSB7XG4gICAgICAgICAgICBjYWxsYmFjayhic3RyaW5nVG9CdWZmZXIocmVhZGVyLnJlc3VsdCkpO1xuICAgICAgICB9O1xuICAgICAgICByZWFkZXIucmVhZEFzQmluYXJ5U3RyaW5nKHRoaXMuYmxvYik7XG4gICAgfVxuXG59IGVsc2Uge1xuICAgIC8vIGlmIChjb25zb2xlICYmIGNvbnNvbGUubG9nKVxuICAgIC8vICAgIGNvbnNvbGUubG9nKCdkZWZpbmluZyBzeW5jIEJsb2JGZXRjaGFibGUuZmV0Y2gnKTtcblxuICAgIEJsb2JGZXRjaGFibGUucHJvdG90eXBlLmZldGNoID0gZnVuY3Rpb24oY2FsbGJhY2spIHtcbiAgICAgICAgdmFyIHJlYWRlciA9IG5ldyBGaWxlUmVhZGVyU3luYygpO1xuICAgICAgICB0cnkge1xuICAgICAgICAgICAgdmFyIHJlcyA9IHJlYWRlci5yZWFkQXNBcnJheUJ1ZmZlcih0aGlzLmJsb2IpO1xuICAgICAgICAgICAgY2FsbGJhY2socmVzKTtcbiAgICAgICAgfSBjYXRjaCAoZSkge1xuICAgICAgICAgICAgY2FsbGJhY2sobnVsbCwgZSk7XG4gICAgICAgIH1cbiAgICB9XG59XG5cbmZ1bmN0aW9uIFVSTEZldGNoYWJsZSh1cmwsIHN0YXJ0LCBlbmQsIG9wdHMpIHtcbiAgICBpZiAoIW9wdHMpIHtcbiAgICAgICAgaWYgKHR5cGVvZiBzdGFydCA9PT0gJ29iamVjdCcpIHtcbiAgICAgICAgICAgIG9wdHMgPSBzdGFydDtcbiAgICAgICAgICAgIHN0YXJ0ID0gdW5kZWZpbmVkO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgb3B0cyA9IHt9O1xuICAgICAgICB9XG4gICAgfVxuXG4gICAgdGhpcy51cmwgPSB1cmw7XG4gICAgdGhpcy5zdGFydCA9IHN0YXJ0IHx8IDA7XG4gICAgaWYgKGVuZCkge1xuICAgICAgICB0aGlzLmVuZCA9IGVuZDtcbiAgICB9XG4gICAgdGhpcy5vcHRzID0gb3B0cztcbn1cblxuVVJMRmV0Y2hhYmxlLnByb3RvdHlwZS5zbGljZSA9IGZ1bmN0aW9uKHMsIGwpIHtcbiAgICBpZiAocyA8IDApIHtcbiAgICAgICAgdGhyb3cgJ0JhZCBzbGljZSAnICsgcztcbiAgICB9XG5cbiAgICB2YXIgbnMgPSB0aGlzLnN0YXJ0LCBuZSA9IHRoaXMuZW5kO1xuICAgIGlmIChucyAmJiBzKSB7XG4gICAgICAgIG5zID0gbnMgKyBzO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIG5zID0gcyB8fCBucztcbiAgICB9XG4gICAgaWYgKGwgJiYgbnMpIHtcbiAgICAgICAgbmUgPSBucyArIGwgLSAxO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIG5lID0gbmUgfHwgbCAtIDE7XG4gICAgfVxuICAgIHJldHVybiBuZXcgVVJMRmV0Y2hhYmxlKHRoaXMudXJsLCBucywgbmUsIHRoaXMub3B0cyk7XG59XG5cbnZhciBzZWVkPTA7XG52YXIgaXNTYWZhcmkgPSBuYXZpZ2F0b3IudXNlckFnZW50LmluZGV4T2YoJ1NhZmFyaScpID49IDAgJiYgbmF2aWdhdG9yLnVzZXJBZ2VudC5pbmRleE9mKCdDaHJvbWUnKSA8IDAgO1xuXG5VUkxGZXRjaGFibGUucHJvdG90eXBlLmZldGNoQXNUZXh0ID0gZnVuY3Rpb24oY2FsbGJhY2spIHtcbiAgICB0cnkge1xuICAgICAgICB2YXIgcmVxID0gbmV3IFhNTEh0dHBSZXF1ZXN0KCk7XG4gICAgICAgIHZhciBsZW5ndGg7XG4gICAgICAgIHZhciB1cmwgPSB0aGlzLnVybDtcbiAgICAgICAgaWYgKChpc1NhZmFyaSB8fCB0aGlzLm9wdHMuc2FsdCkgJiYgdXJsLmluZGV4T2YoJz8nKSA8IDApIHtcbiAgICAgICAgICAgIHVybCA9IHVybCArICc/c2FsdD0nICsgYjY0X3NoYTEoJycgKyBEYXRlLm5vdygpICsgJywnICsgKCsrc2VlZCkpO1xuICAgICAgICB9XG4gICAgICAgIHJlcS5vcGVuKCdHRVQnLCB1cmwsIHRydWUpO1xuXG4gICAgICAgIGlmICh0aGlzLmVuZCkge1xuICAgICAgICAgICAgaWYgKHRoaXMuZW5kIC0gdGhpcy5zdGFydCA+IDEwMDAwMDAwMCkge1xuICAgICAgICAgICAgICAgIHRocm93ICdNb25zdGVyIGZldGNoISc7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXEuc2V0UmVxdWVzdEhlYWRlcignUmFuZ2UnLCAnYnl0ZXM9JyArIHRoaXMuc3RhcnQgKyAnLScgKyB0aGlzLmVuZCk7XG4gICAgICAgICAgICBsZW5ndGggPSB0aGlzLmVuZCAtIHRoaXMuc3RhcnQgKyAxO1xuICAgICAgICB9XG5cbiAgICAgICAgcmVxLm9ucmVhZHlzdGF0ZWNoYW5nZSA9IGZ1bmN0aW9uKCkge1xuICAgICAgICAgICAgaWYgKHJlcS5yZWFkeVN0YXRlID09IDQpIHtcbiAgICAgICAgICAgICAgICBpZiAocmVxLnN0YXR1cyA9PSAyMDAgfHwgcmVxLnN0YXR1cyA9PSAyMDYpIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKHJlcS5yZXNwb25zZVRleHQpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH07XG4gICAgICAgIGlmICh0aGlzLm9wdHMuY3JlZGVudGlhbHMpIHtcbiAgICAgICAgICAgIHJlcS53aXRoQ3JlZGVudGlhbHMgPSB0cnVlO1xuICAgICAgICB9XG4gICAgICAgIHJlcS5zZW5kKCcnKTtcbiAgICB9IGNhdGNoIChlKSB7XG4gICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsKTtcbiAgICB9XG59XG5cblVSTEZldGNoYWJsZS5wcm90b3R5cGUuc2FsdGVkID0gZnVuY3Rpb24oKSB7XG4gICAgdmFyIG8gPSBzaGFsbG93Q29weSh0aGlzLm9wdHMpO1xuICAgIG8uc2FsdCA9IHRydWU7XG4gICAgcmV0dXJuIG5ldyBVUkxGZXRjaGFibGUodGhpcy51cmwsIHRoaXMuc3RhcnQsIHRoaXMuZW5kLCBvKTtcbn1cblxuVVJMRmV0Y2hhYmxlLnByb3RvdHlwZS5nZXRVUkwgPSBmdW5jdGlvbigpIHtcbiAgICBpZiAodGhpcy5vcHRzLnJlc29sdmVyKSB7XG4gICAgICAgIHJldHVybiB0aGlzLm9wdHMucmVzb2x2ZXIodGhpcy51cmwpLnRoZW4oZnVuY3Rpb24gKHVybE9yT2JqKSB7XG4gICAgICAgICAgICBpZiAodHlwZW9mIHVybE9yT2JqID09PSAnc3RyaW5nJykge1xuICAgICAgICAgICAgICAgIHJldHVybiB1cmxPck9iajtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIHVybE9yT2JqLnVybDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuIFByb21pc2UucmVzb2x2ZSh0aGlzLnVybCk7XG4gICAgfVxufVxuXG5VUkxGZXRjaGFibGUucHJvdG90eXBlLmZldGNoID0gZnVuY3Rpb24oY2FsbGJhY2ssIG9wdHMpIHtcbiAgICB2YXIgdGhpc0IgPSB0aGlzO1xuIFxuICAgIG9wdHMgPSBvcHRzIHx8IHt9O1xuICAgIHZhciBhdHRlbXB0ID0gb3B0cy5hdHRlbXB0IHx8IDE7XG4gICAgdmFyIHRydW5jYXRlZExlbmd0aCA9IG9wdHMudHJ1bmNhdGVkTGVuZ3RoO1xuICAgIGlmIChhdHRlbXB0ID4gMykge1xuICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgfVxuXG4gICAgdGhpcy5nZXRVUkwoKS50aGVuKGZ1bmN0aW9uKHVybCkge1xuICAgICAgICB0cnkge1xuICAgICAgICAgICAgdmFyIHRpbWVvdXQ7XG4gICAgICAgICAgICBpZiAob3B0cy50aW1lb3V0ICYmICF0aGlzQi5vcHRzLmNyZWRlbnRpYWxzKSB7XG4gICAgICAgICAgICAgICAgdGltZW91dCA9IHNldFRpbWVvdXQoXG4gICAgICAgICAgICAgICAgICAgIGZ1bmN0aW9uKCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgY29uc29sZS5sb2coJ3RpbWluZyBvdXQgJyArIHVybCk7XG4gICAgICAgICAgICAgICAgICAgICAgICByZXEuYWJvcnQoKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhudWxsLCAnVGltZW91dCcpO1xuICAgICAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgICAgICBvcHRzLnRpbWVvdXRcbiAgICAgICAgICAgICAgICApO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgXG4gICAgICAgICAgICB2YXIgcmVxID0gbmV3IFhNTEh0dHBSZXF1ZXN0KCk7XG4gICAgICAgICAgICB2YXIgbGVuZ3RoO1xuICAgICAgICAgICAgaWYgKChpc1NhZmFyaSB8fCB0aGlzQi5vcHRzLnNhbHQpICYmIHVybC5pbmRleE9mKCc/JykgPCAwKSB7XG4gICAgICAgICAgICAgICAgdXJsID0gdXJsICsgJz9zYWx0PScgKyBiNjRfc2hhMSgnJyArIERhdGUubm93KCkgKyAnLCcgKyAoKytzZWVkKSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICByZXEub3BlbignR0VUJywgdXJsLCB0cnVlKTtcbiAgICAgICAgICAgIHJlcS5vdmVycmlkZU1pbWVUeXBlKCd0ZXh0L3BsYWluOyBjaGFyc2V0PXgtdXNlci1kZWZpbmVkJyk7XG4gICAgICAgICAgICBpZiAodGhpc0IuZW5kKSB7XG4gICAgICAgICAgICAgICAgaWYgKHRoaXNCLmVuZCAtIHRoaXNCLnN0YXJ0ID4gMTAwMDAwMDAwKSB7XG4gICAgICAgICAgICAgICAgICAgIHRocm93ICdNb25zdGVyIGZldGNoISc7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHJlcS5zZXRSZXF1ZXN0SGVhZGVyKCdSYW5nZScsICdieXRlcz0nICsgdGhpc0Iuc3RhcnQgKyAnLScgKyB0aGlzQi5lbmQpO1xuICAgICAgICAgICAgICAgIGxlbmd0aCA9IHRoaXNCLmVuZCAtIHRoaXNCLnN0YXJ0ICsgMTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJlcS5yZXNwb25zZVR5cGUgPSAnYXJyYXlidWZmZXInO1xuICAgICAgICAgICAgcmVxLm9ucmVhZHlzdGF0ZWNoYW5nZSA9IGZ1bmN0aW9uKCkge1xuICAgICAgICAgICAgICAgIGlmIChyZXEucmVhZHlTdGF0ZSA9PSA0KSB7XG4gICAgICAgICAgICAgICAgICAgIGlmICh0aW1lb3V0KVxuICAgICAgICAgICAgICAgICAgICAgICAgY2xlYXJUaW1lb3V0KHRpbWVvdXQpO1xuICAgICAgICAgICAgICAgICAgICBpZiAocmVxLnN0YXR1cyA9PSAyMDAgfHwgcmVxLnN0YXR1cyA9PSAyMDYpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChyZXEucmVzcG9uc2UpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB2YXIgYmwgPSByZXEucmVzcG9uc2UuYnl0ZUxlbmd0aDtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAobGVuZ3RoICYmIGxlbmd0aCAhPSBibCAmJiAoIXRydW5jYXRlZExlbmd0aCB8fCBibCAhPSB0cnVuY2F0ZWRMZW5ndGgpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiB0aGlzQi5mZXRjaChjYWxsYmFjaywge2F0dGVtcHQ6IGF0dGVtcHQgKyAxLCB0cnVuY2F0ZWRMZW5ndGg6IGJsfSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKHJlcS5yZXNwb25zZSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIGlmIChyZXEubW96UmVzcG9uc2VBcnJheUJ1ZmZlcikge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBjYWxsYmFjayhyZXEubW96UmVzcG9uc2VBcnJheUJ1ZmZlcik7XG4gICAgICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHZhciByID0gcmVxLnJlc3BvbnNlVGV4dDtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAobGVuZ3RoICYmIGxlbmd0aCAhPSByLmxlbmd0aCAmJiAoIXRydW5jYXRlZExlbmd0aCB8fCByLmxlbmd0aCAhPSB0cnVuY2F0ZWRMZW5ndGgpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiB0aGlzQi5mZXRjaChjYWxsYmFjaywge2F0dGVtcHQ6IGF0dGVtcHQgKyAxLCB0cnVuY2F0ZWRMZW5ndGg6IHIubGVuZ3RofSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGJzdHJpbmdUb0J1ZmZlcihyZXEucmVzcG9uc2VUZXh0KSk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIHRoaXNCLmZldGNoKGNhbGxiYWNrLCB7YXR0ZW1wdDogYXR0ZW1wdCArIDF9KTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH07XG4gICAgICAgICAgICBpZiAodGhpc0Iub3B0cy5jcmVkZW50aWFscykge1xuICAgICAgICAgICAgICAgIHJlcS53aXRoQ3JlZGVudGlhbHMgPSB0cnVlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmVxLnNlbmQoJycpO1xuICAgICAgICB9IGNhdGNoIChlKSB7XG4gICAgICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgICAgIH1cbiAgICB9KS5jYXRjaChmdW5jdGlvbihlcnIpIHtcbiAgICAgICAgY29uc29sZS5sb2coZXJyKTtcbiAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsIGVycik7XG4gICAgfSk7XG59XG4gICAgICAgICAgICAgICAgICAgICAgIFxuZnVuY3Rpb24gYnN0cmluZ1RvQnVmZmVyKHJlc3VsdCkge1xuICAgIGlmICghcmVzdWx0KSB7XG4gICAgICAgIHJldHVybiBudWxsO1xuICAgIH1cblxuICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KHJlc3VsdC5sZW5ndGgpO1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgYmEubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgYmFbaV0gPSByZXN1bHQuY2hhckNvZGVBdChpKTtcbiAgICB9XG4gICAgcmV0dXJuIGJhLmJ1ZmZlcjtcbn1cblxuLy8gUmVhZCBmcm9tIFVpbnQ4QXJyYXlcblxuKGZ1bmN0aW9uKGdsb2JhbCkge1xuICAgIHZhciBjb252ZXJ0QnVmZmVyID0gbmV3IEFycmF5QnVmZmVyKDgpO1xuICAgIHZhciBiYSA9IG5ldyBVaW50OEFycmF5KGNvbnZlcnRCdWZmZXIpO1xuICAgIHZhciBmYSA9IG5ldyBGbG9hdDMyQXJyYXkoY29udmVydEJ1ZmZlcik7XG5cblxuICAgIGdsb2JhbC5yZWFkRmxvYXQgPSBmdW5jdGlvbihidWYsIG9mZnNldCkge1xuICAgICAgICBiYVswXSA9IGJ1ZltvZmZzZXRdO1xuICAgICAgICBiYVsxXSA9IGJ1ZltvZmZzZXQrMV07XG4gICAgICAgIGJhWzJdID0gYnVmW29mZnNldCsyXTtcbiAgICAgICAgYmFbM10gPSBidWZbb2Zmc2V0KzNdO1xuICAgICAgICByZXR1cm4gZmFbMF07XG4gICAgfTtcbiB9KHRoaXMpKTtcblxuZnVuY3Rpb24gcmVhZEludDY0KGJhLCBvZmZzZXQpIHtcbiAgICByZXR1cm4gKGJhW29mZnNldCArIDddIDw8IDI0KSB8IChiYVtvZmZzZXQgKyA2XSA8PCAxNikgfCAoYmFbb2Zmc2V0ICsgNV0gPDwgOCkgfCAoYmFbb2Zmc2V0ICsgNF0pO1xufVxuXG5mdW5jdGlvbiByZWFkSW50KGJhLCBvZmZzZXQpIHtcbiAgICByZXR1cm4gKGJhW29mZnNldCArIDNdIDw8IDI0KSB8IChiYVtvZmZzZXQgKyAyXSA8PCAxNikgfCAoYmFbb2Zmc2V0ICsgMV0gPDwgOCkgfCAoYmFbb2Zmc2V0XSk7XG59XG5cbmZ1bmN0aW9uIHJlYWRTaG9ydChiYSwgb2Zmc2V0KSB7XG4gICAgcmV0dXJuIChiYVtvZmZzZXQgKyAxXSA8PCA4KSB8IChiYVtvZmZzZXRdKTtcbn1cblxuZnVuY3Rpb24gcmVhZEJ5dGUoYmEsIG9mZnNldCkge1xuICAgIHJldHVybiBiYVtvZmZzZXRdO1xufVxuXG5mdW5jdGlvbiByZWFkSW50QkUoYmEsIG9mZnNldCkge1xuICAgIHJldHVybiAoYmFbb2Zmc2V0XSA8PCAyNCkgfCAoYmFbb2Zmc2V0ICsgMV0gPDwgMTYpIHwgKGJhW29mZnNldCArIDJdIDw8IDgpIHwgKGJhW29mZnNldCArIDNdKTtcbn1cblxuLy8gRXhwb3J0cyBpZiB3ZSBhcmUgYmVpbmcgdXNlZCBhcyBhIG1vZHVsZVxuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIEJsb2JGZXRjaGFibGU6IEJsb2JGZXRjaGFibGUsXG4gICAgICAgIFVSTEZldGNoYWJsZTogVVJMRmV0Y2hhYmxlLFxuXG4gICAgICAgIHJlYWRJbnQ6IHJlYWRJbnQsXG4gICAgICAgIHJlYWRJbnRCRTogcmVhZEludEJFLFxuICAgICAgICByZWFkSW50NjQ6IHJlYWRJbnQ2NCxcbiAgICAgICAgcmVhZFNob3J0OiByZWFkU2hvcnQsXG4gICAgICAgIHJlYWRCeXRlOiByZWFkQnl0ZSxcbiAgICAgICAgcmVhZEZsb2F0OiB0aGlzLnJlYWRGbG9hdFxuICAgIH1cbn1cbiIsIi8qIC0qLSBtb2RlOiBqYXZhc2NyaXB0OyBjLWJhc2ljLW9mZnNldDogNDsgaW5kZW50LXRhYnMtbW9kZTogbmlsIC0qLSAqL1xuXG4vLyBcbi8vIERhbGxpYW5jZSBHZW5vbWUgRXhwbG9yZXJcbi8vIChjKSBUaG9tYXMgRG93biAyMDA2LTIwMTBcbi8vXG4vLyBjb2xvci5qc1xuLy9cblxuXCJ1c2Ugc3RyaWN0XCI7XG5cbmZ1bmN0aW9uIERDb2xvdXIocmVkLCBncmVlbiwgYmx1ZSwgbmFtZSkge1xuICAgIHRoaXMucmVkID0gcmVkfDA7XG4gICAgdGhpcy5ncmVlbiA9IGdyZWVufDA7XG4gICAgdGhpcy5ibHVlID0gYmx1ZXwwO1xuICAgIGlmIChuYW1lKSB7XG4gICAgICAgIHRoaXMubmFtZSA9IG5hbWU7XG4gICAgfVxufVxuXG5EQ29sb3VyLnByb3RvdHlwZS50b1N2Z1N0cmluZyA9IGZ1bmN0aW9uKCkge1xuICAgIGlmICghdGhpcy5uYW1lKSB7XG4gICAgICAgIHRoaXMubmFtZSA9IFwicmdiKFwiICsgdGhpcy5yZWQgKyBcIixcIiArIHRoaXMuZ3JlZW4gKyBcIixcIiArIHRoaXMuYmx1ZSArIFwiKVwiO1xuICAgIH1cblxuICAgIHJldHVybiB0aGlzLm5hbWU7XG59XG5cbmZ1bmN0aW9uIGhleDIoeCkge1xuICAgIHZhciB5ID0gJzAwJyArIHgudG9TdHJpbmcoMTYpO1xuICAgIHJldHVybiB5LnN1YnN0cmluZyh5Lmxlbmd0aCAtIDIpO1xufVxuXG5EQ29sb3VyLnByb3RvdHlwZS50b0hleFN0cmluZyA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiAnIycgKyBoZXgyKHRoaXMucmVkKSArIGhleDIodGhpcy5ncmVlbikgKyBoZXgyKHRoaXMuYmx1ZSk7XG59XG5cbnZhciBwYWxldHRlID0ge1xuICAgIHJlZDogbmV3IERDb2xvdXIoMjU1LCAwLCAwLCAncmVkJyksXG4gICAgZ3JlZW46IG5ldyBEQ29sb3VyKDAsIDI1NSwgMCwgJ2dyZWVuJyksXG4gICAgYmx1ZTogbmV3IERDb2xvdXIoMCwgMCwgMjU1LCAnYmx1ZScpLFxuICAgIHllbGxvdzogbmV3IERDb2xvdXIoMjU1LCAyNTUsIDAsICd5ZWxsb3cnKSxcbiAgICB3aGl0ZTogbmV3IERDb2xvdXIoMjU1LCAyNTUsIDI1NSwgJ3doaXRlJyksXG4gICAgYmxhY2s6IG5ldyBEQ29sb3VyKDAsIDAsIDAsICdibGFjaycpLFxuICAgIGdyYXk6IG5ldyBEQ29sb3VyKDE4MCwgMTgwLCAxODAsICdncmF5JyksXG4gICAgZ3JleTogbmV3IERDb2xvdXIoMTgwLCAxODAsIDE4MCwgJ2dyZXknKSxcbiAgICBsaWdodHNreWJsdWU6IG5ldyBEQ29sb3VyKDEzNSwgMjA2LCAyNTAsICdsaWdodHNreWJsdWUnKSxcbiAgICBsaWdodHNhbG1vbjogbmV3IERDb2xvdXIoMjU1LCAxNjAsIDEyMiwgJ2xpZ2h0c2FsbW9uJyksXG4gICAgaG90cGluazogbmV3IERDb2xvdXIoMjU1LCAxMDUsIDE4MCwgJ2hvdHBpbmsnKVxufTtcblxudmFyIENPTE9SX1JFID0gbmV3IFJlZ0V4cCgnXiMoWzAtOUEtRmEtZl17Mn0pKFswLTlBLUZhLWZdezJ9KShbMC05QS1GYS1mXXsyfSkkJyk7XG52YXIgQ1NTX0NPTE9SX1JFID0gL3JnYlxcKChbMC05XSspLChbMC05XSspLChbMC05XSspXFwpL1xuXG5mdW5jdGlvbiBkYXNDb2xvdXJGb3JOYW1lKG5hbWUpIHtcbiAgICB2YXIgYyA9IHBhbGV0dGVbbmFtZV07XG4gICAgaWYgKCFjKSB7XG4gICAgICAgIHZhciBtYXRjaCA9IENPTE9SX1JFLmV4ZWMobmFtZSk7XG4gICAgICAgIGlmIChtYXRjaCkge1xuICAgICAgICAgICAgYyA9IG5ldyBEQ29sb3VyKCgnMHgnICsgbWF0Y2hbMV0pfDAsICgnMHgnICsgbWF0Y2hbMl0pfDAsICgnMHgnICsgbWF0Y2hbM10pfDAsIG5hbWUpO1xuICAgICAgICAgICAgcGFsZXR0ZVtuYW1lXSA9IGM7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgXHQgICAgbWF0Y2ggPSBDU1NfQ09MT1JfUkUuZXhlYyhuYW1lKTtcbiAgICBcdCAgICBpZiAobWF0Y2gpIHtcbiAgICAgICAgXHRcdGMgPSBuZXcgRENvbG91cihtYXRjaFsxXXwwLCBtYXRjaFsyXXwwLCBtYXRjaFszXXwwLCBuYW1lKTtcbiAgICAgICAgXHRcdHBhbGV0dGVbbmFtZV0gPSBjO1xuXHQgICAgICAgfSBlbHNlIHtcblx0XHQgICAgICBjb25zb2xlLmxvZyhcImNvdWxkbid0IGhhbmRsZSBjb2xvcjogXCIgKyBuYW1lKTtcblx0XHQgICAgICBjID0gcGFsZXR0ZS5ibGFjaztcblx0XHQgICAgICBwYWxldHRlW25hbWVdID0gYztcblx0ICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gYztcbn1cblxuZnVuY3Rpb24gbWFrZUNvbG91clN0ZXBzKHN0ZXBzLCBzdG9wcywgY29sb3Vycykge1xuICAgIHZhciBkY29sb3VycyA9IFtdO1xuICAgIGZvciAodmFyIGNpID0gMDsgY2kgPCBjb2xvdXJzLmxlbmd0aDsgKytjaSkge1xuICAgICAgICBkY29sb3Vycy5wdXNoKGRhc0NvbG91ckZvck5hbWUoY29sb3Vyc1tjaV0pKTtcbiAgICB9XG5cbiAgICB2YXIgZ3JhZCA9IFtdO1xuICBTVEVQX0xPT1A6XG4gICAgZm9yICh2YXIgc2kgPSAwOyBzaSA8IHN0ZXBzOyArK3NpKSB7XG4gICAgICAgIHZhciBycyA9ICgxLjAgKiBzaSkgLyAoc3RlcHMtMSk7XG4gICAgICAgIHZhciBzY29yZSA9IHN0b3BzWzBdICsgKHN0b3BzW3N0b3BzLmxlbmd0aCAtMV0gLSBzdG9wc1swXSkgKiBycztcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBzdG9wcy5sZW5ndGggLSAxOyArK2kpIHtcbiAgICAgICAgICAgIGlmIChzY29yZSA+PSBzdG9wc1tpXSAmJiBzY29yZSA8PSBzdG9wc1tpKzFdKSB7XG4gICAgICAgICAgICAgICAgdmFyIGZyYWMgPSAoc2NvcmUgLSBzdG9wc1tpXSkgLyAoc3RvcHNbaSsxXSAtIHN0b3BzW2ldKTtcbiAgICAgICAgICAgICAgICB2YXIgY2EgPSBkY29sb3Vyc1tpXTtcbiAgICAgICAgICAgICAgICB2YXIgY2IgPSBkY29sb3Vyc1tpKzFdO1xuXG4gICAgICAgICAgICAgICAgdmFyIGZpbGwgPSBuZXcgRENvbG91cihcbiAgICAgICAgICAgICAgICAgICAgKChjYS5yZWQgKiAoMS4wIC0gZnJhYykpICsgKGNiLnJlZCAqIGZyYWMpKXwwLFxuICAgICAgICAgICAgICAgICAgICAoKGNhLmdyZWVuICogKDEuMCAtIGZyYWMpKSArIChjYi5ncmVlbiAqIGZyYWMpKXwwLFxuICAgICAgICAgICAgICAgICAgICAoKGNhLmJsdWUgKiAoMS4wIC0gZnJhYykpICsgKGNiLmJsdWUgKiBmcmFjKSl8MFxuICAgICAgICAgICAgICAgICkudG9TdmdTdHJpbmcoKTtcbiAgICAgICAgICAgICAgICBncmFkLnB1c2goZmlsbCk7XG5cbiAgICAgICAgICAgICAgICBjb250aW51ZSBTVEVQX0xPT1A7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgdGhyb3cgJ0JhZCBzdGVwJztcbiAgICB9XG5cbiAgICByZXR1cm4gZ3JhZDtcbn1cblxuZnVuY3Rpb24gbWFrZUdyYWRpZW50KHN0ZXBzLCBjb2xvcjEsIGNvbG9yMiwgY29sb3IzKSB7XG4gICAgaWYgKGNvbG9yMykge1xuICAgICAgICByZXR1cm4gbWFrZUNvbG91clN0ZXBzKHN0ZXBzLCBbMCwgMC41LCAxXSwgW2NvbG9yMSwgY29sb3IyLCBjb2xvcjNdKTtcbiAgICB9IGVsc2Uge1xuICAgICAgICByZXR1cm4gbWFrZUNvbG91clN0ZXBzKHN0ZXBzLCBbMCwgMV0sIFtjb2xvcjEsIGNvbG9yMl0pO1xuICAgIH1cbn1cblxuaWYgKHR5cGVvZihtb2R1bGUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIG1vZHVsZS5leHBvcnRzID0ge1xuICAgICAgICBtYWtlQ29sb3VyU3RlcHM6IG1ha2VDb2xvdXJTdGVwcyxcbiAgICAgICAgbWFrZUdyYWRpZW50OiBtYWtlR3JhZGllbnQsXG4gICAgICAgIGRhc0NvbG91ckZvck5hbWU6IGRhc0NvbG91ckZvck5hbWVcbiAgICB9O1xufVxuIiwiLyogLSotIG1vZGU6IGphdmFzY3JpcHQ7IGMtYmFzaWMtb2Zmc2V0OiA0OyBpbmRlbnQtdGFicy1tb2RlOiBuaWwgLSotICovXG5cbi8vIFxuLy8gRGFsbGlhbmNlIEdlbm9tZSBFeHBsb3JlclxuLy8gKGMpIFRob21hcyBEb3duIDIwMDYtMjAxMFxuLy9cbi8vIGRhcy5qczogcXVlcmllcyBhbmQgbG93LWxldmVsIGRhdGEgbW9kZWwuXG4vL1xuXG5cInVzZSBzdHJpY3RcIjtcblxuaWYgKHR5cGVvZihyZXF1aXJlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgICB2YXIgdXRpbHMgPSByZXF1aXJlKCcuL3V0aWxzJyk7XG4gICAgdmFyIHNoYWxsb3dDb3B5ID0gdXRpbHMuc2hhbGxvd0NvcHk7XG4gICAgdmFyIHB1c2hvID0gdXRpbHMucHVzaG87XG5cbiAgICB2YXIgY29sb3IgPSByZXF1aXJlKCcuL2NvbG9yJyk7XG4gICAgdmFyIG1ha2VDb2xvdXJTdGVwcyA9IGNvbG9yLm1ha2VDb2xvdXJTdGVwcztcbn1cblxudmFyIGRhc0xpYkVycm9ySGFuZGxlciA9IGZ1bmN0aW9uKGVyck1zZykge1xuICAgIGFsZXJ0KGVyck1zZyk7XG59XG52YXIgZGFzTGliUmVxdWVzdFF1ZXVlID0gbmV3IEFycmF5KCk7XG5cbmZ1bmN0aW9uIERBU1NlZ21lbnQobmFtZSwgc3RhcnQsIGVuZCwgZGVzY3JpcHRpb24pIHtcbiAgICB0aGlzLm5hbWUgPSBuYW1lO1xuICAgIHRoaXMuc3RhcnQgPSBzdGFydDtcbiAgICB0aGlzLmVuZCA9IGVuZDtcbiAgICB0aGlzLmRlc2NyaXB0aW9uID0gZGVzY3JpcHRpb247XG59XG5EQVNTZWdtZW50LnByb3RvdHlwZS50b1N0cmluZyA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLm5hbWUgKyAnOicgKyB0aGlzLnN0YXJ0ICsgJy4uJyArIHRoaXMuZW5kO1xufTtcbkRBU1NlZ21lbnQucHJvdG90eXBlLmlzQm91bmRlZCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLnN0YXJ0ICYmIHRoaXMuZW5kO1xufVxuREFTU2VnbWVudC5wcm90b3R5cGUudG9EQVNRdWVyeSA9IGZ1bmN0aW9uKCkge1xuICAgIHZhciBxID0gJ3NlZ21lbnQ9JyArIHRoaXMubmFtZTtcbiAgICBpZiAodGhpcy5zdGFydCAmJiB0aGlzLmVuZCkge1xuICAgICAgICBxICs9ICgnOicgKyB0aGlzLnN0YXJ0ICsgJywnICsgdGhpcy5lbmQpO1xuICAgIH1cbiAgICByZXR1cm4gcTtcbn1cblxuXG5mdW5jdGlvbiBEQVNTb3VyY2UoYTEsIGEyKSB7XG4gICAgdmFyIG9wdGlvbnM7XG4gICAgaWYgKHR5cGVvZiBhMSA9PSAnc3RyaW5nJykge1xuICAgICAgICB0aGlzLnVyaSA9IGExO1xuICAgICAgICBvcHRpb25zID0gYTIgfHwge307XG4gICAgfSBlbHNlIHtcbiAgICAgICAgb3B0aW9ucyA9IGExIHx8IHt9O1xuICAgIH1cbiAgICBmb3IgKHZhciBrIGluIG9wdGlvbnMpIHtcbiAgICAgICAgdGhpc1trXSA9IG9wdGlvbnNba107XG4gICAgfVxuXG4gICAgaWYgKCF0aGlzLmNvb3Jkcykge1xuICAgICAgICB0aGlzLmNvb3JkcyA9IFtdO1xuICAgIH1cbiAgICBpZiAoIXRoaXMucHJvcHMpIHtcbiAgICAgICAgdGhpcy5wcm9wcyA9IHt9O1xuICAgIH1cblxuICAgIHRoaXMuZGFzQmFzZVVSSSA9IHRoaXMudXJpO1xuICAgIGlmICh0aGlzLmRhc0Jhc2VVUkkgJiYgdGhpcy5kYXNCYXNlVVJJLnN1YnN0cih0aGlzLnVyaS5sZW5ndGggLSAxKSAhPSAnLycpIHtcbiAgICAgICAgdGhpcy5kYXNCYXNlVVJJID0gdGhpcy5kYXNCYXNlVVJJICsgJy8nO1xuICAgIH1cbn1cblxuREFTU291cmNlLnByb3RvdHlwZS5nZXRVUkkgPSBmdW5jdGlvbih1cmkpIHtcbiAgICBpZiAodGhpcy5yZXNvbHZlcikge1xuICAgICAgICByZXR1cm4gdGhpcy5yZXNvbHZlcih1cmkpLnRoZW4oZnVuY3Rpb24gKHVybE9yT2JqKSB7XG4gICAgICAgICAgICBpZiAodHlwZW9mIHVybE9yT2JqID09PSAnc3RyaW5nJykge1xuICAgICAgICAgICAgICAgIHJldHVybiB1cmxPck9iajtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIHVybE9yT2JqLnVybDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuIFByb21pc2UucmVzb2x2ZSh1cmkpO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gREFTQ29vcmRzKCkge1xufVxuXG5mdW5jdGlvbiBjb29yZHNNYXRjaChjMSwgYzIpIHtcbiAgICByZXR1cm4gYzEudGF4b24gPT0gYzIudGF4b24gJiYgYzEuYXV0aCA9PSBjMi5hdXRoICYmIGMxLnZlcnNpb24gPT0gYzIudmVyc2lvbjtcbn1cblxuLy9cbi8vIERBUyAxLjYgZW50cnlfcG9pbnRzIGNvbW1hbmRcbi8vXG5cbkRBU1NvdXJjZS5wcm90b3R5cGUuZW50cnlQb2ludHMgPSBmdW5jdGlvbihjYWxsYmFjaykge1xuICAgIHZhciBkYXNVUkkgPSB0aGlzLmRhc0Jhc2VVUkkgKyAnZW50cnlfcG9pbnRzJztcbiAgICB0aGlzLmRvQ3Jvc3NEb21haW5SZXF1ZXN0KGRhc1VSSSwgZnVuY3Rpb24ocmVzcG9uc2VYTUwpIHtcbiAgICAgICAgICAgIGlmICghcmVzcG9uc2VYTUwpIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2soW10pO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgdmFyIGVudHJ5UG9pbnRzID0gbmV3IEFycmF5KCk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgdmFyIHNlZ3MgPSByZXNwb25zZVhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnU0VHTUVOVCcpO1xuICAgICAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgc2Vncy5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnID0gc2Vnc1tpXTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZ0lkID0gc2VnLmdldEF0dHJpYnV0ZSgnaWQnKTtcbiAgICAgICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWdTaXplID0gc2VnLmdldEF0dHJpYnV0ZSgnc2l6ZScpO1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnTWluLCBzZWdNYXg7XG4gICAgICAgICAgICAgICAgICAgIGlmIChzZWdTaXplKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBzZWdNaW4gPSAxOyBzZWdNYXggPSBzZWdTaXplfDA7XG4gICAgICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBzZWdNaW4gPSBzZWcuZ2V0QXR0cmlidXRlKCdzdGFydCcpO1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKHNlZ01pbikge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHNlZ01pbiB8PSAwO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgc2VnTWF4ID0gc2VnLmdldEF0dHJpYnV0ZSgnc3RvcCcpO1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKHNlZ01heCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHNlZ01heCB8PSAwO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWdEZXNjID0gbnVsbDtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHNlZy5maXJzdENoaWxkKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBzZWdEZXNjID0gc2VnLmZpcnN0Q2hpbGQubm9kZVZhbHVlO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIGVudHJ5UG9pbnRzLnB1c2gobmV3IERBU1NlZ21lbnQoc2VnSWQsIHNlZ01pbiwgc2VnTWF4LCBzZWdEZXNjKSk7XG4gICAgICAgICAgICAgICAgfSAgICAgICAgICBcbiAgICAgICAgICAgICAgIGNhbGxiYWNrKGVudHJ5UG9pbnRzKTtcbiAgICB9KTsgICAgICAgICBcbn1cblxuLy9cbi8vIERBUyAxLjYgc2VxdWVuY2UgY29tbWFuZFxuLy8gRG8gd2UgbmVlZCBhbiBvcHRpb24gdG8gZmFsbCBiYWNrIHRvIHRoZSBkbmEgY29tbWFuZD9cbi8vXG5cbmZ1bmN0aW9uIERBU1NlcXVlbmNlKG5hbWUsIHN0YXJ0LCBlbmQsIGFscGhhLCBzZXEpIHtcbiAgICB0aGlzLm5hbWUgPSBuYW1lO1xuICAgIHRoaXMuc3RhcnQgPSBzdGFydDtcbiAgICB0aGlzLmVuZCA9IGVuZDtcbiAgICB0aGlzLmFscGhhYmV0ID0gYWxwaGE7XG4gICAgdGhpcy5zZXEgPSBzZXE7XG59XG5cbkRBU1NvdXJjZS5wcm90b3R5cGUuc2VxdWVuY2UgPSBmdW5jdGlvbihzZWdtZW50LCBjYWxsYmFjaykge1xuICAgIHZhciBkYXNVUkkgPSB0aGlzLmRhc0Jhc2VVUkkgKyAnc2VxdWVuY2U/JyArIHNlZ21lbnQudG9EQVNRdWVyeSgpO1xuICAgIHRoaXMuZG9Dcm9zc0RvbWFpblJlcXVlc3QoZGFzVVJJLCBmdW5jdGlvbihyZXNwb25zZVhNTCkge1xuICAgICAgICBpZiAoIXJlc3BvbnNlWE1MKSB7XG4gICAgICAgICAgICBjYWxsYmFjayhbXSk7XG4gICAgICAgICAgICByZXR1cm47XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgdmFyIHNlcXMgPSBuZXcgQXJyYXkoKTtcbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICB2YXIgc2VncyA9IHJlc3BvbnNlWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdTRVFVRU5DRScpO1xuICAgICAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgc2Vncy5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnID0gc2Vnc1tpXTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZ0lkID0gc2VnLmdldEF0dHJpYnV0ZSgnaWQnKTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHNlZ01pbiA9IHNlZy5nZXRBdHRyaWJ1dGUoJ3N0YXJ0Jyk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWdNYXggPSBzZWcuZ2V0QXR0cmlidXRlKCdzdG9wJyk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWdBbHBoYSA9ICdETkEnO1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnU2VxID0gbnVsbDtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHNlZy5maXJzdENoaWxkKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgcmF3U2VxID0gc2VnLmZpcnN0Q2hpbGQubm9kZVZhbHVlO1xuICAgICAgICAgICAgICAgICAgICAgICAgc2VnU2VxID0gJyc7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgaWR4ID0gMDtcbiAgICAgICAgICAgICAgICAgICAgICAgIHdoaWxlICh0cnVlKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIHNwYWNlID0gcmF3U2VxLmluZGV4T2YoJ1xcbicsIGlkeCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgaWYgKHNwYWNlID49IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgc2VnU2VxICs9IHJhd1NlcS5zdWJzdHJpbmcoaWR4LCBzcGFjZSkudG9VcHBlckNhc2UoKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaWR4ID0gc3BhY2UgKyAxO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHNlZ1NlcSArPSByYXdTZXEuc3Vic3RyaW5nKGlkeCkudG9VcHBlckNhc2UoKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIHNlcXMucHVzaChuZXcgREFTU2VxdWVuY2Uoc2VnSWQsIHNlZ01pbiwgc2VnTWF4LCBzZWdBbHBoYSwgc2VnU2VxKSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIGNhbGxiYWNrKHNlcXMpO1xuICAgICAgICB9XG4gICAgfSk7XG59XG5cbi8vXG4vLyBEQVMgMS42IGZlYXR1cmVzIGNvbW1hbmRcbi8vXG5cbmZ1bmN0aW9uIERBU0ZlYXR1cmUoKSB7XG59XG5cbmZ1bmN0aW9uIERBU0dyb3VwKGlkKSB7XG4gICAgaWYgKGlkKVxuICAgICAgICB0aGlzLmlkID0gaWQ7XG59XG5cbmZ1bmN0aW9uIERBU0xpbmsoZGVzYywgdXJpKSB7XG4gICAgdGhpcy5kZXNjID0gZGVzYztcbiAgICB0aGlzLnVyaSA9IHVyaTtcbn1cblxuREFTU291cmNlLnByb3RvdHlwZS5mZWF0dXJlcyA9IGZ1bmN0aW9uKHNlZ21lbnQsIG9wdGlvbnMsIGNhbGxiYWNrKSB7XG4gICAgb3B0aW9ucyA9IG9wdGlvbnMgfHwge307XG4gICAgdmFyIHRoaXNCID0gdGhpcztcblxuICAgIHZhciBkYXNVUkk7XG4gICAgaWYgKHRoaXMuZmVhdHVyZXNfdXJpKSB7XG4gICAgICAgIGRhc1VSSSA9IHRoaXMuZmVhdHVyZXNfdXJpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHZhciBmaWx0ZXJzID0gW107XG5cbiAgICAgICAgaWYgKHNlZ21lbnQpIHtcbiAgICAgICAgICAgIGZpbHRlcnMucHVzaChzZWdtZW50LnRvREFTUXVlcnkoKSk7XG4gICAgICAgIH0gZWxzZSBpZiAob3B0aW9ucy5ncm91cCkge1xuICAgICAgICAgICAgdmFyIGcgPSBvcHRpb25zLmdyb3VwO1xuICAgICAgICAgICAgaWYgKHR5cGVvZiBnID09ICdzdHJpbmcnKSB7XG4gICAgICAgICAgICAgICAgZmlsdGVycy5wdXNoKCdncm91cF9pZD0nICsgZyk7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIGZvciAodmFyIGdpID0gMDsgZ2kgPCBnLmxlbmd0aDsgKytnaSkge1xuICAgICAgICAgICAgICAgICAgICBmaWx0ZXJzLnB1c2goJ2dyb3VwX2lkPScgKyBnW2dpXSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgaWYgKG9wdGlvbnMuYWRqYWNlbnQpIHtcbiAgICAgICAgICAgIHZhciBhZGogPSBvcHRpb25zLmFkamFjZW50O1xuICAgICAgICAgICAgaWYgKHR5cGVvZiBhZGogPT0gJ3N0cmluZycpIHtcbiAgICAgICAgICAgICAgICBhZGogPSBbYWRqXTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGZvciAodmFyIGFpID0gMDsgYWkgPCBhZGoubGVuZ3RoOyArK2FpKSB7XG4gICAgICAgICAgICAgICAgZmlsdGVycy5wdXNoKCdhZGphY2VudD0nICsgYWRqW2FpXSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICBpZiAob3B0aW9ucy50eXBlKSB7XG4gICAgICAgICAgICBpZiAodHlwZW9mIG9wdGlvbnMudHlwZSA9PSAnc3RyaW5nJykge1xuICAgICAgICAgICAgICAgIGZpbHRlcnMucHVzaCgndHlwZT0nICsgb3B0aW9ucy50eXBlKTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgdGkgPSAwOyB0aSA8IG9wdGlvbnMudHlwZS5sZW5ndGg7ICsrdGkpIHtcbiAgICAgICAgICAgICAgICAgICAgZmlsdGVycy5wdXNoKCd0eXBlPScgKyBvcHRpb25zLnR5cGVbdGldKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgXG4gICAgICAgIGlmIChvcHRpb25zLm1heGJpbnMpIHtcbiAgICAgICAgICAgIGZpbHRlcnMucHVzaCgnbWF4Ymlucz0nICsgb3B0aW9ucy5tYXhiaW5zKTtcbiAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgaWYgKGZpbHRlcnMubGVuZ3RoID4gMCkge1xuICAgICAgICAgICAgZGFzVVJJID0gdGhpcy5kYXNCYXNlVVJJICsgJ2ZlYXR1cmVzPycgKyBmaWx0ZXJzLmpvaW4oJzsnKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGNhbGxiYWNrKFtdLCAnTm8gZmlsdGVycyBzcGVjaWZpZWQnKTtcbiAgICAgICAgfVxuICAgIH0gXG4gICBcblxuICAgIHRoaXMuZG9Dcm9zc0RvbWFpblJlcXVlc3QoZGFzVVJJLCBmdW5jdGlvbihyZXNwb25zZVhNTCwgcmVxKSB7XG4gICAgICAgIGlmICghcmVzcG9uc2VYTUwpIHtcbiAgICAgICAgICAgIHZhciBtc2c7XG4gICAgICAgICAgICBpZiAocmVxLnN0YXR1cyA9PSAwKSB7XG4gICAgICAgICAgICAgICAgbXNnID0gJ3NlcnZlciBtYXkgbm90IHN1cHBvcnQgQ09SUyc7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIG1zZyA9ICdzdGF0dXM9JyArIHJlcS5zdGF0dXM7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBjYWxsYmFjayhbXSwgJ0ZhaWxlZCByZXF1ZXN0OiAnICsgbXNnKTtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuLyogICAgICBpZiAocmVxKSB7XG4gICAgICAgICAgICB2YXIgY2FwcyA9IHJlcS5nZXRSZXNwb25zZUhlYWRlcignWC1EQVMtQ2FwYWJpbHRpZXMnKTtcbiAgICAgICAgICAgIGlmIChjYXBzKSB7XG4gICAgICAgICAgICAgICAgYWxlcnQoY2Fwcyk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0gKi9cblxuICAgICAgICB2YXIgZmVhdHVyZXMgPSBuZXcgQXJyYXkoKTtcbiAgICAgICAgdmFyIHNlZ21lbnRNYXAgPSB7fTtcblxuICAgICAgICB2YXIgc2VncyA9IHJlc3BvbnNlWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdTRUdNRU5UJyk7XG4gICAgICAgIGZvciAodmFyIHNpID0gMDsgc2kgPCBzZWdzLmxlbmd0aDsgKytzaSkge1xuICAgICAgICAgICAgdmFyIHNlZ21lbnRYTUwgPSBzZWdzW3NpXTtcbiAgICAgICAgICAgIHZhciBzZWdtZW50SUQgPSBzZWdtZW50WE1MLmdldEF0dHJpYnV0ZSgnaWQnKTtcbiAgICAgICAgICAgIHNlZ21lbnRNYXBbc2VnbWVudElEXSA9IHtcbiAgICAgICAgICAgICAgICBtaW46IHNlZ21lbnRYTUwuZ2V0QXR0cmlidXRlKCdzdGFydCcpLFxuICAgICAgICAgICAgICAgIG1heDogc2VnbWVudFhNTC5nZXRBdHRyaWJ1dGUoJ3N0b3AnKVxuICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIFxuICAgICAgICAgICAgdmFyIGZlYXR1cmVYTUxzID0gc2VnbWVudFhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnRkVBVFVSRScpO1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBmZWF0dXJlWE1Mcy5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgICAgIHZhciBmZWF0dXJlID0gZmVhdHVyZVhNTHNbaV07XG4gICAgICAgICAgICAgICAgdmFyIGRhc0ZlYXR1cmUgPSBuZXcgREFTRmVhdHVyZSgpO1xuICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUuc2VnbWVudCA9IHNlZ21lbnRJRDtcbiAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLmlkID0gZmVhdHVyZS5nZXRBdHRyaWJ1dGUoJ2lkJyk7XG4gICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5sYWJlbCA9IGZlYXR1cmUuZ2V0QXR0cmlidXRlKCdsYWJlbCcpO1xuXG5cbi8qXG4gICAgICAgICAgICAgICAgdmFyIGNoaWxkTm9kZXMgPSBmZWF0dXJlLmNoaWxkTm9kZXM7XG4gICAgICAgICAgICAgICAgZm9yICh2YXIgYyA9IDA7IGMgPCBjaGlsZE5vZGVzLmxlbmd0aDsgKytjKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBjbiA9IGNoaWxkTm9kZXNbY107XG4gICAgICAgICAgICAgICAgICAgIGlmIChjbi5ub2RlVHlwZSA9PSBOb2RlLkVMRU1FTlRfTk9ERSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGtleSA9IGNuLnRhZ05hbWU7XG4gICAgICAgICAgICAgICAgICAgICAgICAvL3ZhciB2YWwgPSBudWxsO1xuICAgICAgICAgICAgICAgICAgICAgICAgLy9pZiAoY24uZmlyc3RDaGlsZCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgLy8gICB2YWwgPSBjbi5maXJzdENoaWxkLm5vZGVWYWx1ZTtcbiAgICAgICAgICAgICAgICAgICAgICAgIC8vfVxuICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZVtrZXldID0gJ3gnO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfSAqL1xuXG5cbiAgICAgICAgICAgICAgICB2YXIgc3BvcyA9IGVsZW1lbnRWYWx1ZShmZWF0dXJlLCBcIlNUQVJUXCIpO1xuICAgICAgICAgICAgICAgIHZhciBlcG9zID0gZWxlbWVudFZhbHVlKGZlYXR1cmUsIFwiRU5EXCIpO1xuICAgICAgICAgICAgICAgIGlmICgoc3Bvc3wwKSA+IChlcG9zfDApKSB7XG4gICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUubWluID0gZXBvc3wwO1xuICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLm1heCA9IHNwb3N8MDtcbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLm1pbiA9IHNwb3N8MDtcbiAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5tYXggPSBlcG9zfDA7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHRlYyA9IGZlYXR1cmUuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1RZUEUnKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHRlYy5sZW5ndGggPiAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgdGUgPSB0ZWNbMF07XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAodGUuZmlyc3RDaGlsZCkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUudHlwZSA9IHRlLmZpcnN0Q2hpbGQubm9kZVZhbHVlO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS50eXBlSWQgPSB0ZS5nZXRBdHRyaWJ1dGUoJ2lkJyk7XG4gICAgICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLnR5cGVDdiA9IHRlLmdldEF0dHJpYnV0ZSgnY3ZJZCcpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUudHlwZSA9IGVsZW1lbnRWYWx1ZShmZWF0dXJlLCBcIlRZUEVcIik7XG4gICAgICAgICAgICAgICAgaWYgKCFkYXNGZWF0dXJlLnR5cGUgJiYgZGFzRmVhdHVyZS50eXBlSWQpIHtcbiAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS50eXBlID0gZGFzRmVhdHVyZS50eXBlSWQ7IC8vIEZJWE1FP1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLm1ldGhvZCA9IGVsZW1lbnRWYWx1ZShmZWF0dXJlLCBcIk1FVEhPRFwiKTtcbiAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBvcmkgPSBlbGVtZW50VmFsdWUoZmVhdHVyZSwgXCJPUklFTlRBVElPTlwiKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCFvcmkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIG9yaSA9ICcwJztcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLm9yaWVudGF0aW9uID0gb3JpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLnNjb3JlID0gZWxlbWVudFZhbHVlKGZlYXR1cmUsIFwiU0NPUkVcIik7XG4gICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5saW5rcyA9IGRhc0xpbmtzT2YoZmVhdHVyZSk7XG4gICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5ub3RlcyA9IGRhc05vdGVzT2YoZmVhdHVyZSk7XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAgdmFyIGdyb3VwcyA9IGZlYXR1cmUuZ2V0RWxlbWVudHNCeVRhZ05hbWUoXCJHUk9VUFwiKTtcbiAgICAgICAgICAgICAgICBmb3IgKHZhciBnaSAgPSAwOyBnaSA8IGdyb3Vwcy5sZW5ndGg7ICsrZ2kpIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGdyb3VwWE1MID0gZ3JvdXBzW2dpXTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIGRhc0dyb3VwID0gbmV3IERBU0dyb3VwKCk7XG4gICAgICAgICAgICAgICAgICAgIGRhc0dyb3VwLnR5cGUgPSBncm91cFhNTC5nZXRBdHRyaWJ1dGUoJ3R5cGUnKTtcbiAgICAgICAgICAgICAgICAgICAgZGFzR3JvdXAuaWQgPSBncm91cFhNTC5nZXRBdHRyaWJ1dGUoJ2lkJyk7XG4gICAgICAgICAgICAgICAgICAgIGRhc0dyb3VwLmxpbmtzID0gZGFzTGlua3NPZihncm91cFhNTCk7XG4gICAgICAgICAgICAgICAgICAgIGRhc0dyb3VwLm5vdGVzID0gZGFzTm90ZXNPZihncm91cFhNTCk7XG4gICAgICAgICAgICAgICAgICAgIGlmICghZGFzRmVhdHVyZS5ncm91cHMpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUuZ3JvdXBzID0gbmV3IEFycmF5KGRhc0dyb3VwKTtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUuZ3JvdXBzLnB1c2goZGFzR3JvdXApO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICAgICAgLy8gTWFnaWMgbm90ZXMuICBDaGVjayB3aXRoIFRBRCBiZWZvcmUgY2hhbmdpbmcgdGhpcy5cbiAgICAgICAgICAgICAgICBpZiAoZGFzRmVhdHVyZS5ub3Rlcykge1xuICAgICAgICAgICAgICAgICAgICBmb3IgKHZhciBuaSA9IDA7IG5pIDwgZGFzRmVhdHVyZS5ub3Rlcy5sZW5ndGg7ICsrbmkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBuID0gZGFzRmVhdHVyZS5ub3Rlc1tuaV07XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAobi5pbmRleE9mKCdHZW5lbmFtZT0nKSA9PSAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgdmFyIGdnID0gbmV3IERBU0dyb3VwKCk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZ2cudHlwZT0nZ2VuZSc7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgZ2cuaWQgPSBuLnN1YnN0cmluZyg5KTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoIWRhc0ZlYXR1cmUuZ3JvdXBzKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUuZ3JvdXBzID0gbmV3IEFycmF5KGdnKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBkYXNGZWF0dXJlLmdyb3Vwcy5wdXNoKGdnKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgXG4gICAgICAgICAgICAgICAge1xuICAgICAgICAgICAgICAgICAgICB2YXIgcGVjID0gZmVhdHVyZS5nZXRFbGVtZW50c0J5VGFnTmFtZSgnUEFSVCcpO1xuICAgICAgICAgICAgICAgICAgICBpZiAocGVjLmxlbmd0aCA+IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBwYXJ0cyA9IFtdO1xuICAgICAgICAgICAgICAgICAgICAgICAgZm9yICh2YXIgcGkgPSAwOyBwaSA8IHBlYy5sZW5ndGg7ICsrcGkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBwYXJ0cy5wdXNoKHBlY1twaV0uZ2V0QXR0cmlidXRlKCdpZCcpKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgICAgIGRhc0ZlYXR1cmUucGFydHMgPSBwYXJ0cztcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBwZWMgPSBmZWF0dXJlLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdQQVJFTlQnKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKHBlYy5sZW5ndGggPiAwKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICB2YXIgcGFyZW50cyA9IFtdO1xuICAgICAgICAgICAgICAgICAgICAgICAgZm9yICh2YXIgcGkgPSAwOyBwaSA8IHBlYy5sZW5ndGg7ICsrcGkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBwYXJlbnRzLnB1c2gocGVjW3BpXS5nZXRBdHRyaWJ1dGUoJ2lkJykpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgZGFzRmVhdHVyZS5wYXJlbnRzID0gcGFyZW50cztcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICBmZWF0dXJlcy5wdXNoKGRhc0ZlYXR1cmUpO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgICAgICAgICAgXG4gICAgICAgIGNhbGxiYWNrKGZlYXR1cmVzLCB1bmRlZmluZWQsIHNlZ21lbnRNYXApO1xuICAgIH0sXG4gICAgZnVuY3Rpb24gKGVycikge1xuICAgICAgICBjYWxsYmFjayhbXSwgZXJyKTtcbiAgICB9KTtcbn1cblxuZnVuY3Rpb24gREFTQWxpZ25tZW50KHR5cGUpIHtcbiAgICB0aGlzLnR5cGUgPSB0eXBlO1xuICAgIHRoaXMub2JqZWN0cyA9IHt9O1xuICAgIHRoaXMuYmxvY2tzID0gW107XG59XG5cbkRBU1NvdXJjZS5wcm90b3R5cGUuYWxpZ25tZW50cyA9IGZ1bmN0aW9uKHNlZ21lbnQsIG9wdGlvbnMsIGNhbGxiYWNrKSB7XG4gICAgdmFyIGRhc1VSSSA9IHRoaXMuZGFzQmFzZVVSSSArICdhbGlnbm1lbnQ/cXVlcnk9JyArIHNlZ21lbnQ7XG4gICAgdGhpcy5kb0Nyb3NzRG9tYWluUmVxdWVzdChkYXNVUkksIGZ1bmN0aW9uKHJlc3BvbnNlWE1MKSB7XG4gICAgICAgIGlmICghcmVzcG9uc2VYTUwpIHtcbiAgICAgICAgICAgIGNhbGxiYWNrKFtdLCAnRmFpbGVkIHJlcXVlc3QgJyArIGRhc1VSSSk7XG4gICAgICAgICAgICByZXR1cm47XG4gICAgICAgIH1cblxuICAgICAgICB2YXIgYWxpZ25tZW50cyA9IFtdO1xuICAgICAgICB2YXIgYWxpWE1McyA9IHJlc3BvbnNlWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdhbGlnbm1lbnQnKTtcbiAgICAgICAgZm9yICh2YXIgYWkgPSAwOyBhaSA8IGFsaVhNTHMubGVuZ3RoOyArK2FpKSB7XG4gICAgICAgICAgICB2YXIgYWxpWE1MID0gYWxpWE1Mc1thaV07XG4gICAgICAgICAgICB2YXIgYWxpID0gbmV3IERBU0FsaWdubWVudChhbGlYTUwuZ2V0QXR0cmlidXRlKCdhbGlnblR5cGUnKSk7XG4gICAgICAgICAgICB2YXIgb2JqWE1McyA9IGFsaVhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnYWxpZ25PYmplY3QnKTtcbiAgICAgICAgICAgIGZvciAodmFyIG9pID0gMDsgb2kgPCBvYmpYTUxzLmxlbmd0aDsgKytvaSkge1xuICAgICAgICAgICAgICAgIHZhciBvYmpYTUwgPSBvYmpYTUxzW29pXTtcbiAgICAgICAgICAgICAgICB2YXIgb2JqID0ge1xuICAgICAgICAgICAgICAgICAgICBpZDogICAgICAgICAgb2JqWE1MLmdldEF0dHJpYnV0ZSgnaW50T2JqZWN0SWQnKSxcbiAgICAgICAgICAgICAgICAgICAgYWNjZXNzaW9uOiAgIG9ialhNTC5nZXRBdHRyaWJ1dGUoJ2RiQWNjZXNzaW9uSWQnKSxcbiAgICAgICAgICAgICAgICAgICAgdmVyc2lvbjogICAgIG9ialhNTC5nZXRBdHRyaWJ1dGUoJ29iamVjdFZlcnNpb24nKSxcbiAgICAgICAgICAgICAgICAgICAgZGJTb3VyY2U6ICAgIG9ialhNTC5nZXRBdHRyaWJ1dGUoJ2RiU291cmNlJyksXG4gICAgICAgICAgICAgICAgICAgIGRiVmVyc2lvbjogICBvYmpYTUwuZ2V0QXR0cmlidXRlKCdkYlZlcnNpb24nKVxuICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICAgICAgYWxpLm9iamVjdHNbb2JqLmlkXSA9IG9iajtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIFxuICAgICAgICAgICAgdmFyIGJsb2NrWE1McyA9IGFsaVhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnYmxvY2snKTtcbiAgICAgICAgICAgIGZvciAodmFyIGJpID0gMDsgYmkgPCBibG9ja1hNTHMubGVuZ3RoOyArK2JpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGJsb2NrWE1MID0gYmxvY2tYTUxzW2JpXTtcbiAgICAgICAgICAgICAgICB2YXIgYmxvY2sgPSB7XG4gICAgICAgICAgICAgICAgICAgIG9yZGVyOiAgICAgIGJsb2NrWE1MLmdldEF0dHJpYnV0ZSgnYmxvY2tPcmRlcicpLFxuICAgICAgICAgICAgICAgICAgICBzZWdtZW50czogICBbXVxuICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICAgICAgdmFyIHNlZ1hNTHMgPSBibG9ja1hNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnc2VnbWVudCcpO1xuICAgICAgICAgICAgICAgIGZvciAodmFyIHNpID0gMDsgc2kgPCBzZWdYTUxzLmxlbmd0aDsgKytzaSkge1xuICAgICAgICAgICAgICAgICAgICB2YXIgc2VnWE1MID0gc2VnWE1Mc1tzaV07XG4gICAgICAgICAgICAgICAgICAgIHZhciBzZWcgPSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBvYmplY3Q6ICAgICAgc2VnWE1MLmdldEF0dHJpYnV0ZSgnaW50T2JqZWN0SWQnKSxcbiAgICAgICAgICAgICAgICAgICAgICAgIG1pbjogICAgICAgICBzZWdYTUwuZ2V0QXR0cmlidXRlKCdzdGFydCcpLFxuICAgICAgICAgICAgICAgICAgICAgICAgbWF4OiAgICAgICAgIHNlZ1hNTC5nZXRBdHRyaWJ1dGUoJ2VuZCcpLFxuICAgICAgICAgICAgICAgICAgICAgICAgc3RyYW5kOiAgICAgIHNlZ1hNTC5nZXRBdHRyaWJ1dGUoJ3N0cmFuZCcpLFxuICAgICAgICAgICAgICAgICAgICAgICAgY2lnYXI6ICAgICAgIGVsZW1lbnRWYWx1ZShzZWdYTUwsICdjaWdhcicpXG4gICAgICAgICAgICAgICAgICAgIH07XG4gICAgICAgICAgICAgICAgICAgIGJsb2NrLnNlZ21lbnRzLnB1c2goc2VnKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgYWxpLmJsb2Nrcy5wdXNoKGJsb2NrKTtcbiAgICAgICAgICAgIH0gICAgICAgXG4gICAgICAgICAgICAgICAgICAgIFxuICAgICAgICAgICAgYWxpZ25tZW50cy5wdXNoKGFsaSk7XG4gICAgICAgIH1cbiAgICAgICAgY2FsbGJhY2soYWxpZ25tZW50cyk7XG4gICAgfSk7XG59XG5cblxuZnVuY3Rpb24gREFTU3R5bGVzaGVldCgpIHtcbiAgICB0aGlzLnN0eWxlcyA9IFtdO1xufVxuXG5EQVNTdHlsZXNoZWV0LnByb3RvdHlwZS5wdXNoU3R5bGUgPSBmdW5jdGlvbihmaWx0ZXJzLCB6b29tLCBzdHlsZSkge1xuICAgIGlmICghZmlsdGVycykge1xuICAgICAgICBmaWx0ZXJzID0ge3R5cGU6ICdkZWZhdWx0J307XG4gICAgfVxuICAgIHZhciBzdHlsZUhvbGRlciA9IHNoYWxsb3dDb3B5KGZpbHRlcnMpO1xuICAgIGlmICh6b29tKSB7XG4gICAgICAgIHN0eWxlSG9sZGVyLnpvb20gPSB6b29tO1xuICAgIH1cbiAgICBzdHlsZUhvbGRlci5zdHlsZSA9IHN0eWxlO1xuICAgIHRoaXMuc3R5bGVzLnB1c2goc3R5bGVIb2xkZXIpO1xufVxuXG5mdW5jdGlvbiBEQVNTdHlsZSgpIHtcbn1cblxuZnVuY3Rpb24gcGFyc2VHcmFkaWVudChncmFkKSB7XG4gICAgdmFyIHN0ZXBzID0gZ3JhZC5nZXRBdHRyaWJ1dGUoJ3N0ZXBzJyk7XG4gICAgaWYgKHN0ZXBzKSB7XG4gICAgICAgIHN0ZXBzID0gc3RlcHN8MDtcbiAgICB9IGVsc2Uge1xuICAgICAgICBzdGVwcyA9IDUwO1xuICAgIH1cblxuICAgIHZhciBzdG9wcyA9IFtdO1xuICAgIHZhciBjb2xvcnMgPSBbXTtcbiAgICB2YXIgc2UgPSBncmFkLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdTVE9QJyk7XG4gICAgZm9yICh2YXIgc2kgPSAwOyBzaSA8IHNlLmxlbmd0aDsgKytzaSkge1xuICAgICAgICB2YXIgc3RvcCA9IHNlW3NpXTtcbiAgICAgICAgc3RvcHMucHVzaCgxLjAgKiBzdG9wLmdldEF0dHJpYnV0ZSgnc2NvcmUnKSk7XG4gICAgICAgIGNvbG9ycy5wdXNoKHN0b3AuZmlyc3RDaGlsZC5ub2RlVmFsdWUpO1xuICAgIH1cblxuICAgIHJldHVybiBtYWtlQ29sb3VyU3RlcHMoc3RlcHMsIHN0b3BzLCBjb2xvcnMpO1xufVxuXG5EQVNTb3VyY2UucHJvdG90eXBlLnN0eWxlc2hlZXQgPSBmdW5jdGlvbihzdWNjZXNzQ0IsIGZhaWx1cmVDQikge1xuICAgIHZhciBkYXNVUkksIGNyZWRzID0gdGhpcy5jcmVkZW50aWFscztcbiAgICBpZiAodGhpcy5zdHlsZXNoZWV0X3VyaSkge1xuICAgICAgICBkYXNVUkkgPSB0aGlzLnN0eWxlc2hlZXRfdXJpO1xuICAgICAgICBjcmVkcyA9IGZhbHNlO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIGRhc1VSSSA9IHRoaXMuZGFzQmFzZVVSSSArICdzdHlsZXNoZWV0JztcbiAgICB9XG5cbiAgICB0aGlzLmdldFVSSShkYXNVUkkpLnRoZW4oZnVuY3Rpb24oZGFzVVJJKSB7XG4gICAgICAgIGRvQ3Jvc3NEb21haW5SZXF1ZXN0KGRhc1VSSSwgZnVuY3Rpb24ocmVzcG9uc2VYTUwpIHtcbiAgICAgICAgICAgIGlmICghcmVzcG9uc2VYTUwpIHtcbiAgICAgICAgICAgICAgICBpZiAoZmFpbHVyZUNCKSB7XG4gICAgICAgICAgICAgICAgICAgIGZhaWx1cmVDQigpO1xuICAgICAgICAgICAgICAgIH0gXG4gICAgICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgdmFyIHN0eWxlc2hlZXQgPSBuZXcgREFTU3R5bGVzaGVldCgpO1xuICAgICAgICAgICAgdmFyIHR5cGVYTUxzID0gcmVzcG9uc2VYTUwuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1RZUEUnKTtcbiAgICAgICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgdHlwZVhNTHMubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgICAgICAgICB2YXIgdHlwZVN0eWxlID0gdHlwZVhNTHNbaV07XG4gICAgICAgICAgICBcbiAgICAgICAgICAgICAgICB2YXIgZmlsdGVyID0ge307XG4gICAgICAgICAgICAgICAgZmlsdGVyLnR5cGUgPSB0eXBlU3R5bGUuZ2V0QXR0cmlidXRlKCdpZCcpOyAvLyBBbSBJIHJpZ2h0IGluIHRoaW5raW5nIHRoYXQgdGhpcyBtYWtlcyBEQVNTVFlMRSBYTUwgaW52YWxpZD8gIFVnaC5cbiAgICAgICAgICAgICAgICBmaWx0ZXIubGFiZWwgPSB0eXBlU3R5bGUuZ2V0QXR0cmlidXRlKCdsYWJlbCcpO1xuICAgICAgICAgICAgICAgIGZpbHRlci5tZXRob2QgPSB0eXBlU3R5bGUuZ2V0QXR0cmlidXRlKCdtZXRob2QnKTtcbiAgICAgICAgICAgICAgICB2YXIgZ2x5cGhYTUxzID0gdHlwZVN0eWxlLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdHTFlQSCcpO1xuICAgICAgICAgICAgICAgIGZvciAodmFyIGdpID0gMDsgZ2kgPCBnbHlwaFhNTHMubGVuZ3RoOyArK2dpKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBnbHlwaFhNTCA9IGdseXBoWE1Mc1tnaV07XG4gICAgICAgICAgICAgICAgICAgIHZhciB6b29tID0gZ2x5cGhYTUwuZ2V0QXR0cmlidXRlKCd6b29tJyk7XG4gICAgICAgICAgICAgICAgICAgIHZhciBnbHlwaCA9IGNoaWxkRWxlbWVudE9mKGdseXBoWE1MKTtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHN0eWxlID0gbmV3IERBU1N0eWxlKCk7XG4gICAgICAgICAgICAgICAgICAgIHN0eWxlLmdseXBoID0gZ2x5cGgubG9jYWxOYW1lO1xuICAgICAgICAgICAgICAgICAgICB2YXIgY2hpbGQgPSBnbHlwaC5maXJzdENoaWxkO1xuICAgICAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICAgICAgd2hpbGUgKGNoaWxkKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICBpZiAoY2hpbGQubm9kZVR5cGUgPT0gTm9kZS5FTEVNRU5UX05PREUpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBpZiAoY2hpbGQubG9jYWxOYW1lID09ICdCR0dSQUQnKSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIHN0eWxlW2NoaWxkLmxvY2FsTmFtZV0gPSBwYXJzZUdyYWRpZW50KGNoaWxkKTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICB9IGVsc2UgeyAgICAgIFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBzdHlsZVtjaGlsZC5sb2NhbE5hbWVdID0gY2hpbGQuZmlyc3RDaGlsZC5ub2RlVmFsdWU7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICAgICAgY2hpbGQgPSBjaGlsZC5uZXh0U2libGluZztcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgICBzdHlsZXNoZWV0LnB1c2hTdHlsZShmaWx0ZXIsIHpvb20sIHN0eWxlKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBzdWNjZXNzQ0Ioc3R5bGVzaGVldCk7XG4gICAgICAgIH0sIGNyZWRzKTtcbiAgICB9KS5jYXRjaChmdW5jdGlvbihlcnIpIHtcbiAgICAgICAgY29uc29sZS5sb2coZXJyKTtcbiAgICAgICAgZmFpbHVyZUNCKCk7XG4gICAgfSk7XG59XG5cbi8vXG4vLyBzb3VyY2VzIGNvbW1hbmRcbi8vIFxuXG5mdW5jdGlvbiBEQVNSZWdpc3RyeSh1cmksIG9wdHMpXG57XG4gICAgb3B0cyA9IG9wdHMgfHwge307XG4gICAgdGhpcy51cmkgPSB1cmk7XG4gICAgdGhpcy5vcHRzID0gb3B0czsgICBcbn1cblxuREFTUmVnaXN0cnkucHJvdG90eXBlLnNvdXJjZXMgPSBmdW5jdGlvbihjYWxsYmFjaywgZmFpbHVyZSwgb3B0cylcbntcbiAgICBpZiAoIW9wdHMpIHtcbiAgICAgICAgb3B0cyA9IHt9O1xuICAgIH1cblxuICAgIHZhciBmaWx0ZXJzID0gW107XG4gICAgaWYgKG9wdHMudGF4b24pIHtcbiAgICAgICAgZmlsdGVycy5wdXNoKCdvcmdhbmlzbT0nICsgb3B0cy50YXhvbik7XG4gICAgfVxuICAgIGlmIChvcHRzLmF1dGgpIHtcbiAgICAgICAgZmlsdGVycy5wdXNoKCdhdXRob3JpdHk9JyArIG9wdHMuYXV0aCk7XG4gICAgfVxuICAgIGlmIChvcHRzLnZlcnNpb24pIHtcbiAgICAgICAgZmlsdGVycy5wdXNoKCd2ZXJzaW9uPScgKyBvcHRzLnZlcnNpb24pO1xuICAgIH1cbiAgICB2YXIgcXVyaSA9IHRoaXMudXJpO1xuICAgIGlmIChmaWx0ZXJzLmxlbmd0aCA+IDApIHtcbiAgICAgICAgcXVyaSA9IHF1cmkgKyAnPycgKyBmaWx0ZXJzLmpvaW4oJyYnKTsgICAvLyAnJicgYXMgYSBzZXBhcmF0b3IgdG8gaGFjayBhcm91bmQgZGFzcmVnaXN0cnkub3JnIGJ1Zy5cbiAgICB9XG5cbiAgICBkb0Nyb3NzRG9tYWluUmVxdWVzdChxdXJpLCBmdW5jdGlvbihyZXNwb25zZVhNTCkge1xuICAgICAgICBpZiAoIXJlc3BvbnNlWE1MICYmIGZhaWx1cmUpIHtcbiAgICAgICAgICAgIGZhaWx1cmUoKTtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuXG4gICAgICAgIHZhciBzb3VyY2VzID0gW107ICAgICAgIFxuICAgICAgICB2YXIgc291cmNlWE1McyA9IHJlc3BvbnNlWE1MLmdldEVsZW1lbnRzQnlUYWdOYW1lKCdTT1VSQ0UnKTtcbiAgICAgICAgZm9yICh2YXIgc2kgPSAwOyBzaSA8IHNvdXJjZVhNTHMubGVuZ3RoOyArK3NpKSB7XG4gICAgICAgICAgICB2YXIgc291cmNlWE1MID0gc291cmNlWE1Mc1tzaV07XG4gICAgICAgICAgICB2YXIgdmVyc2lvblhNTHMgPSBzb3VyY2VYTUwuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1ZFUlNJT04nKTtcbiAgICAgICAgICAgIGlmICh2ZXJzaW9uWE1Mcy5sZW5ndGggPCAxKSB7XG4gICAgICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICB2YXIgdmVyc2lvblhNTCA9IHZlcnNpb25YTUxzWzBdO1xuXG4gICAgICAgICAgICB2YXIgY29vcmRYTUxzID0gdmVyc2lvblhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnQ09PUkRJTkFURVMnKTtcbiAgICAgICAgICAgIHZhciBjb29yZHMgPSBbXTtcbiAgICAgICAgICAgIGZvciAodmFyIGNpID0gMDsgY2kgPCBjb29yZFhNTHMubGVuZ3RoOyArK2NpKSB7XG4gICAgICAgICAgICAgICAgdmFyIGNvb3JkWE1MID0gY29vcmRYTUxzW2NpXTtcbiAgICAgICAgICAgICAgICB2YXIgY29vcmQgPSBuZXcgREFTQ29vcmRzKCk7XG4gICAgICAgICAgICAgICAgY29vcmQuYXV0aCA9IGNvb3JkWE1MLmdldEF0dHJpYnV0ZSgnYXV0aG9yaXR5Jyk7XG4gICAgICAgICAgICAgICAgY29vcmQudGF4b24gPSBjb29yZFhNTC5nZXRBdHRyaWJ1dGUoJ3RheGlkJyk7XG4gICAgICAgICAgICAgICAgY29vcmQudmVyc2lvbiA9IGNvb3JkWE1MLmdldEF0dHJpYnV0ZSgndmVyc2lvbicpO1xuICAgICAgICAgICAgICAgIGNvb3Jkcy5wdXNoKGNvb3JkKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIFxuICAgICAgICAgICAgdmFyIGNhcHMgPSBbXTtcbiAgICAgICAgICAgIHZhciBjYXBYTUxzID0gdmVyc2lvblhNTC5nZXRFbGVtZW50c0J5VGFnTmFtZSgnQ0FQQUJJTElUWScpO1xuICAgICAgICAgICAgdmFyIHVyaTtcbiAgICAgICAgICAgIGZvciAodmFyIGNpID0gMDsgY2kgPCBjYXBYTUxzLmxlbmd0aDsgKytjaSkge1xuICAgICAgICAgICAgICAgIHZhciBjYXBYTUwgPSBjYXBYTUxzW2NpXTtcbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgICAgICBjYXBzLnB1c2goY2FwWE1MLmdldEF0dHJpYnV0ZSgndHlwZScpKTtcblxuICAgICAgICAgICAgICAgIGlmIChjYXBYTUwuZ2V0QXR0cmlidXRlKCd0eXBlJykgPT0gJ2RhczE6ZmVhdHVyZXMnKSB7XG4gICAgICAgICAgICAgICAgICAgIHZhciBmZXAgPSBjYXBYTUwuZ2V0QXR0cmlidXRlKCdxdWVyeV91cmknKTtcbiAgICAgICAgICAgICAgICAgICAgdXJpID0gZmVwLnN1YnN0cmluZygwLCBmZXAubGVuZ3RoIC0gKCdmZWF0dXJlcycubGVuZ3RoKSk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICB2YXIgcHJvcHMgPSB7fTtcbiAgICAgICAgICAgIHZhciBwcm9wWE1McyA9IHZlcnNpb25YTUwuZ2V0RWxlbWVudHNCeVRhZ05hbWUoJ1BST1AnKTtcbiAgICAgICAgICAgIGZvciAodmFyIHBpID0gMDsgcGkgPCBwcm9wWE1Mcy5sZW5ndGg7ICsrcGkpIHtcbiAgICAgICAgICAgICAgICBwdXNobyhwcm9wcywgcHJvcFhNTHNbcGldLmdldEF0dHJpYnV0ZSgnbmFtZScpLCBwcm9wWE1Mc1twaV0uZ2V0QXR0cmlidXRlKCd2YWx1ZScpKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIFxuICAgICAgICAgICAgaWYgKHVyaSkge1xuICAgICAgICAgICAgICAgIHZhciBzb3VyY2UgPSBuZXcgREFTU291cmNlKHVyaSwge1xuICAgICAgICAgICAgICAgICAgICBzb3VyY2VfdXJpOiBzb3VyY2VYTUwuZ2V0QXR0cmlidXRlKCd1cmknKSxcbiAgICAgICAgICAgICAgICAgICAgbmFtZTogIHNvdXJjZVhNTC5nZXRBdHRyaWJ1dGUoJ3RpdGxlJyksXG4gICAgICAgICAgICAgICAgICAgIGRlc2M6ICBzb3VyY2VYTUwuZ2V0QXR0cmlidXRlKCdkZXNjcmlwdGlvbicpLFxuICAgICAgICAgICAgICAgICAgICBjb29yZHM6IGNvb3JkcyxcbiAgICAgICAgICAgICAgICAgICAgcHJvcHM6IHByb3BzLFxuICAgICAgICAgICAgICAgICAgICBjYXBhYmlsaXRpZXM6IGNhcHNcbiAgICAgICAgICAgICAgICB9KTtcbiAgICAgICAgICAgICAgICBzb3VyY2VzLnB1c2goc291cmNlKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBcbiAgICAgICAgY2FsbGJhY2soc291cmNlcyk7XG4gICAgfSk7XG59XG5cblxuLy9cbi8vIFV0aWxpdHkgZnVuY3Rpb25zXG4vL1xuXG5mdW5jdGlvbiBlbGVtZW50VmFsdWUoZWxlbWVudCwgdGFnKVxue1xuICAgIHZhciBjaGlsZHJlbiA9IGVsZW1lbnQuZ2V0RWxlbWVudHNCeVRhZ05hbWUodGFnKTtcbiAgICBpZiAoY2hpbGRyZW4ubGVuZ3RoID4gMCAmJiBjaGlsZHJlblswXS5maXJzdENoaWxkKSB7XG4gICAgICAgIHZhciBjID0gY2hpbGRyZW5bMF07XG4gICAgICAgIGlmIChjLmNoaWxkTm9kZXMubGVuZ3RoID09IDEpIHtcbiAgICAgICAgICAgIHJldHVybiBjLmZpcnN0Q2hpbGQubm9kZVZhbHVlO1xuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgdmFyIHMgPSAnJztcbiAgICAgICAgICAgIGZvciAodmFyIG5pID0gMDsgbmkgPCBjLmNoaWxkTm9kZXMubGVuZ3RoOyArK25pKSB7XG4gICAgICAgICAgICAgICAgcyArPSBjLmNoaWxkTm9kZXNbbmldLm5vZGVWYWx1ZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybiBzO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2Uge1xuICAgICAgICByZXR1cm4gbnVsbDtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIGNoaWxkRWxlbWVudE9mKGVsZW1lbnQpXG57XG4gICAgaWYgKGVsZW1lbnQuaGFzQ2hpbGROb2RlcygpKSB7XG4gICAgICAgIHZhciBjaGlsZCA9IGVsZW1lbnQuZmlyc3RDaGlsZDtcbiAgICAgICAgZG8ge1xuICAgICAgICAgICAgaWYgKGNoaWxkLm5vZGVUeXBlID09IE5vZGUuRUxFTUVOVF9OT0RFKSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIGNoaWxkO1xuICAgICAgICAgICAgfSBcbiAgICAgICAgICAgIGNoaWxkID0gY2hpbGQubmV4dFNpYmxpbmc7XG4gICAgICAgIH0gd2hpbGUgKGNoaWxkICE9IG51bGwpO1xuICAgIH1cbiAgICByZXR1cm4gbnVsbDtcbn1cblxuXG5mdW5jdGlvbiBkYXNMaW5rc09mKGVsZW1lbnQpXG57XG4gICAgdmFyIGxpbmtzID0gbmV3IEFycmF5KCk7XG4gICAgdmFyIG1heWJlTGlua0NoaWxkZW4gPSBlbGVtZW50LmdldEVsZW1lbnRzQnlUYWdOYW1lKCdMSU5LJyk7XG4gICAgZm9yICh2YXIgY2kgPSAwOyBjaSA8IG1heWJlTGlua0NoaWxkZW4ubGVuZ3RoOyArK2NpKSB7XG4gICAgICAgIHZhciBsaW5rWE1MID0gbWF5YmVMaW5rQ2hpbGRlbltjaV07XG4gICAgICAgIGlmIChsaW5rWE1MLnBhcmVudE5vZGUgPT0gZWxlbWVudCkge1xuICAgICAgICAgICAgbGlua3MucHVzaChuZXcgREFTTGluayhsaW5rWE1MLmZpcnN0Q2hpbGQgPyBsaW5rWE1MLmZpcnN0Q2hpbGQubm9kZVZhbHVlIDogJ1Vua25vd24nLCBsaW5rWE1MLmdldEF0dHJpYnV0ZSgnaHJlZicpKSk7XG4gICAgICAgIH1cbiAgICB9XG4gICAgXG4gICAgcmV0dXJuIGxpbmtzO1xufVxuXG5mdW5jdGlvbiBkYXNOb3Rlc09mKGVsZW1lbnQpXG57XG4gICAgdmFyIG5vdGVzID0gW107XG4gICAgdmFyIG1heWJlTm90ZXMgPSBlbGVtZW50LmdldEVsZW1lbnRzQnlUYWdOYW1lKCdOT1RFJyk7XG4gICAgZm9yICh2YXIgbmkgPSAwOyBuaSA8IG1heWJlTm90ZXMubGVuZ3RoOyArK25pKSB7XG4gICAgICAgIGlmIChtYXliZU5vdGVzW25pXS5maXJzdENoaWxkKSB7XG4gICAgICAgICAgICBub3Rlcy5wdXNoKG1heWJlTm90ZXNbbmldLmZpcnN0Q2hpbGQubm9kZVZhbHVlKTtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gbm90ZXM7XG59XG5cbmZ1bmN0aW9uIGRvQ3Jvc3NEb21haW5SZXF1ZXN0KHVybCwgaGFuZGxlciwgY3JlZGVudGlhbHMsIGN1c3RBdXRoKSB7XG4gICAgLy8gVE9ETzogZXhwbGljaXQgZXJyb3IgaGFuZGxlcnM/XG5cbiAgICBpZiAod2luZG93LlhEb21haW5SZXF1ZXN0KSB7XG4gICAgICAgIHZhciByZXEgPSBuZXcgWERvbWFpblJlcXVlc3QoKTtcbiAgICAgICAgcmVxLm9ubG9hZCA9IGZ1bmN0aW9uKCkge1xuICAgICAgICAgICAgdmFyIGRvbSA9IG5ldyBBY3RpdmVYT2JqZWN0KFwiTWljcm9zb2Z0LlhNTERPTVwiKTtcbiAgICAgICAgICAgIGRvbS5hc3luYyA9IGZhbHNlO1xuICAgICAgICAgICAgZG9tLmxvYWRYTUwocmVxLnJlc3BvbnNlVGV4dCk7XG4gICAgICAgICAgICBoYW5kbGVyKGRvbSk7XG4gICAgICAgIH1cbiAgICAgICAgcmVxLm9wZW4oXCJnZXRcIiwgdXJsKTtcbiAgICAgICAgcmVxLnNlbmQoJycpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHRyeSB7XG4gICAgICAgICAgICB2YXIgcmVxID0gbmV3IFhNTEh0dHBSZXF1ZXN0KCk7XG4gICAgICAgICAgICB2YXIgdGltZW91dCA9IHNldFRpbWVvdXQoXG4gICAgICAgICAgICAgICAgZnVuY3Rpb24oKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUubG9nKCd0aW1pbmcgb3V0ICcgICsgdXJsKTtcbiAgICAgICAgICAgICAgICAgICAgcmVxLmFib3J0KCk7XG4gICAgICAgICAgICAgICAgICAgIGhhbmRsZXIobnVsbCwgcmVxKTtcbiAgICAgICAgICAgICAgICB9LFxuICAgICAgICAgICAgICAgIDUwMDBcbiAgICAgICAgICAgICk7XG5cbiAgICAgICAgICAgIHJlcS50aW1lb3V0ID0gNTAwMDtcbiAgICAgICAgICAgIHJlcS5vbnRpbWVvdXQgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgICAgICAgICBjb25zb2xlLmxvZygndGltZW91dCBvbiAnICsgdXJsKTtcbiAgICAgICAgICAgIH07XG5cbiAgICAgICAgICAgIHJlcS5vbnJlYWR5c3RhdGVjaGFuZ2UgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgICAgICAgICBpZiAocmVxLnJlYWR5U3RhdGUgPT0gNCkge1xuICAgICAgICAgICAgICAgICAgICBjbGVhclRpbWVvdXQodGltZW91dCk7XG4gICAgICAgICAgICAgICAgICAgIGlmIChyZXEuc3RhdHVzID49IDIwMCB8fCByZXEuc3RhdHVzID09IDApIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIGhhbmRsZXIocmVxLnJlc3BvbnNlWE1MLCByZXEpO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfTtcbiAgICAgICAgICAgIHJlcS5vcGVuKFwiZ2V0XCIsIHVybCwgdHJ1ZSk7XG4gICAgICAgICAgICBpZiAoY3JlZGVudGlhbHMpIHtcbiAgICAgICAgICAgICAgICByZXEud2l0aENyZWRlbnRpYWxzID0gdHJ1ZTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIGlmIChjdXN0QXV0aCkge1xuICAgICAgICAgICAgICAgIHJlcS5zZXRSZXF1ZXN0SGVhZGVyKCdYLURBUy1BdXRob3Jpc2F0aW9uJywgY3VzdEF1dGgpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmVxLm92ZXJyaWRlTWltZVR5cGUoJ3RleHQveG1sJyk7XG4gICAgICAgICAgICByZXEuc2V0UmVxdWVzdEhlYWRlcignQWNjZXB0JywgJ2FwcGxpY2F0aW9uL3htbCwqLyonKTtcbiAgICAgICAgICAgIHJlcS5zZW5kKCcnKTtcbiAgICAgICAgfSBjYXRjaCAoZSkge1xuICAgICAgICAgICAgaGFuZGxlcihudWxsLCByZXEsIGUpO1xuICAgICAgICB9XG4gICAgfVxufVxuXG5EQVNTb3VyY2UucHJvdG90eXBlLmRvQ3Jvc3NEb21haW5SZXF1ZXN0ID0gZnVuY3Rpb24odXJsLCBoYW5kbGVyLCBlcnJIYW5kbGVyKSB7XG4gICAgdmFyIGN1c3RBdXRoO1xuICAgIGlmICh0aGlzLnhVc2VyKSB7XG4gICAgICAgIGN1c3RBdXRoID0gJ0Jhc2ljICcgKyBidG9hKHRoaXMueFVzZXIgKyAnOicgKyB0aGlzLnhQYXNzKTtcbiAgICB9XG5cbiAgICB0cnkge1xuICAgICAgICByZXR1cm4gZG9Dcm9zc0RvbWFpblJlcXVlc3QodXJsLCBoYW5kbGVyLCB0aGlzLmNyZWRlbnRpYWxzLCBjdXN0QXV0aCk7XG4gICAgfSBjYXRjaCAoZXJyKSB7XG4gICAgICAgIGlmIChlcnJIYW5kbGVyKSB7XG4gICAgICAgICAgICBlcnJIYW5kbGVyKGVycik7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICB0aHJvdyBlcnI7XG4gICAgICAgIH1cbiAgICB9XG59XG5cbmZ1bmN0aW9uIGlzRGFzQm9vbGVhblRydWUocykge1xuICAgIHMgPSAoJycgKyBzKS50b0xvd2VyQ2FzZSgpO1xuICAgIHJldHVybiBzPT09J3llcycgfHwgcz09PSd0cnVlJztcbn1cblxuZnVuY3Rpb24gaXNEYXNCb29sZWFuTm90RmFsc2Uocykge1xuICAgIGlmICghcylcbiAgICAgICAgcmV0dXJuIGZhbHNlO1xuXG4gICAgcyA9ICgnJyArIHMpLnRvTG93ZXJDYXNlKCk7XG4gICAgcmV0dXJuIHMhPT0nbm8nIHx8IHMhPT0nZmFsc2UnO1xufVxuXG5mdW5jdGlvbiBjb3B5U3R5bGVzaGVldChzcykge1xuICAgIHZhciBuc3MgPSBzaGFsbG93Q29weShzcyk7XG4gICAgbnNzLnN0eWxlcyA9IFtdO1xuICAgIGZvciAodmFyIHNpID0gMDsgc2kgPCBzcy5zdHlsZXMubGVuZ3RoOyArK3NpKSB7XG4gICAgICAgIHZhciBzaCA9IG5zcy5zdHlsZXNbc2ldID0gc2hhbGxvd0NvcHkoc3Muc3R5bGVzW3NpXSk7XG4gICAgICAgIHNoLl9tZXRob2RSRSA9IHNoLl9sYWJlbFJFID0gc2guX3R5cGVSRSA9IHVuZGVmaW5lZDtcbiAgICAgICAgc2guc3R5bGUgPSBzaGFsbG93Q29weShzaC5zdHlsZSk7XG4gICAgICAgIHNoLnN0eWxlLmlkID0gdW5kZWZpbmVkO1xuICAgICAgICBzaC5zdHlsZS5fZ3JhZGllbnQgPSB1bmRlZmluZWQ7XG4gICAgfVxuICAgIHJldHVybiBuc3M7XG59XG5cbmlmICh0eXBlb2YobW9kdWxlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgICBtb2R1bGUuZXhwb3J0cyA9IHtcbiAgICAgICAgREFTR3JvdXA6IERBU0dyb3VwLFxuICAgICAgICBEQVNGZWF0dXJlOiBEQVNGZWF0dXJlLFxuICAgICAgICBEQVNTdHlsZXNoZWV0OiBEQVNTdHlsZXNoZWV0LFxuICAgICAgICBEQVNTdHlsZTogREFTU3R5bGUsXG4gICAgICAgIERBU1NvdXJjZTogREFTU291cmNlLFxuICAgICAgICBEQVNTZWdtZW50OiBEQVNTZWdtZW50LFxuICAgICAgICBEQVNSZWdpc3RyeTogREFTUmVnaXN0cnksXG4gICAgICAgIERBU1NlcXVlbmNlOiBEQVNTZXF1ZW5jZSxcbiAgICAgICAgREFTTGluazogREFTTGluayxcblxuICAgICAgICBpc0Rhc0Jvb2xlYW5UcnVlOiBpc0Rhc0Jvb2xlYW5UcnVlLFxuICAgICAgICBpc0Rhc0Jvb2xlYW5Ob3RGYWxzZTogaXNEYXNCb29sZWFuTm90RmFsc2UsXG4gICAgICAgIGNvcHlTdHlsZXNoZWV0OiBjb3B5U3R5bGVzaGVldCxcbiAgICAgICAgY29vcmRzTWF0Y2g6IGNvb3Jkc01hdGNoXG4gICAgfTtcbn1cbiIsIi8qIC0qLSBtb2RlOiBqYXZhc2NyaXB0OyBjLWJhc2ljLW9mZnNldDogNDsgaW5kZW50LXRhYnMtbW9kZTogbmlsIC0qLSAqL1xuXG4vLyBcbi8vIERhbGxpYW5jZSBHZW5vbWUgRXhwbG9yZXJcbi8vIChjKSBUaG9tYXMgRG93biAyMDA2LTIwMTRcbi8vXG4vLyBlbmNvZGUuanM6IGludGVyZmFjZSBmb3IgRU5DT0RFIERDQyBzZXJ2aWNlc1xuLy9cblxuXCJ1c2Ugc3RyaWN0XCI7XG5cbmlmICh0eXBlb2YocmVxdWlyZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgdmFyIFByb21pc2UgPSByZXF1aXJlKCdlczYtcHJvbWlzZScpLlByb21pc2U7XG59XG5cbmZ1bmN0aW9uIGxvb2t1cEVuY29kZVVSSSh1cmksIGpzb24pIHtcbiAgICBpZiAodXJpLmluZGV4T2YoJz8nKSA8IDApXG4gICAgICAgIHVyaSA9IHVyaSArICc/c29mdD10cnVlJztcblxuICAgIHJldHVybiBuZXcgUHJvbWlzZShmdW5jdGlvbihhY2NlcHQsIHJlamVjdCkge1xuICAgICAgICB2YXIgcmVxID0gbmV3IFhNTEh0dHBSZXF1ZXN0KCk7XG4gICAgICAgIHJlcS5vbnJlYWR5c3RhdGVjaGFuZ2UgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgICAgIGlmIChyZXEucmVhZHlTdGF0ZSA9PSA0KSB7XG4gICAgICAgICAgICAgICAgaWYgKHJlcS5zdGF0dXMgPj0gMzAwKSB7XG4gICAgICAgICAgICAgICAgICAgIHJlamVjdCgnRXJyb3IgY29kZSAnICsgcmVxLnN0YXR1cyk7XG4gICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgdmFyIHJlc3AgPSBKU09OLnBhcnNlKHJlcS5yZXNwb25zZSk7XG4gICAgICAgICAgICAgICAgICAgIGFjY2VwdChqc29uID8gcmVzcCA6IHJlc3AubG9jYXRpb24pO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfTtcbiAgICBcbiAgICAgICAgcmVxLm9wZW4oJ0dFVCcsIHVyaSwgdHJ1ZSk7XG4gICAgICAgIHJlcS5zZXRSZXF1ZXN0SGVhZGVyKCdBY2NlcHQnLCAnYXBwbGljYXRpb24vanNvbicpO1xuICAgICAgICByZXEucmVzcG9uc2VUeXBlID0gJ3RleHQnO1xuICAgICAgICByZXEuc2VuZCgnJyk7XG4gICAgfSk7XG59XG5cbmZ1bmN0aW9uIEVuY29kZVVSTEhvbGRlcih1cmwpIHtcbiAgICB0aGlzLnJhd3VybCA9IHVybDtcbn1cblxuRW5jb2RlVVJMSG9sZGVyLnByb3RvdHlwZS5nZXRVUkxQcm9taXNlID0gZnVuY3Rpb24oKSB7XG4gICAgaWYgKHRoaXMudXJsUHJvbWlzZSAmJiB0aGlzLnVybFByb21pc2VWYWxpZGl0eSA+IERhdGUubm93KCkpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMudXJsUHJvbWlzZTtcbiAgICB9IGVsc2Uge1xuICAgICAgICB0aGlzLnVybFByb21pc2UgPSBsb29rdXBFbmNvZGVVUkkodGhpcy5yYXd1cmwsIHRydWUpLnRoZW4oZnVuY3Rpb24ocmVzcCkge1xuICAgICAgICAgICAgcmV0dXJuIHJlc3AubG9jYXRpb247XG4gICAgICAgIH0pO1xuICAgICAgICB0aGlzLnVybFByb21pc2VWYWxpZGl0eSA9IERhdGUubm93KCkgKyAoMTIgKiAzNjAwICogMTAwMCk7XG4gICAgICAgIHJldHVybiB0aGlzLnVybFByb21pc2U7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBFbmNvZGVGZXRjaGFibGUodXJsLCBzdGFydCwgZW5kLCBvcHRzKSB7XG4gICAgaWYgKCFvcHRzKSB7XG4gICAgICAgIGlmICh0eXBlb2Ygc3RhcnQgPT09ICdvYmplY3QnKSB7XG4gICAgICAgICAgICBvcHRzID0gc3RhcnQ7XG4gICAgICAgICAgICBzdGFydCA9IHVuZGVmaW5lZDtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIG9wdHMgPSB7fTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHRoaXMudXJsID0gKHR5cGVvZiB1cmwgPT09ICdzdHJpbmcnID8gbmV3IEVuY29kZVVSTEhvbGRlcih1cmwpIDogdXJsKTtcbiAgICB0aGlzLnN0YXJ0ID0gc3RhcnQgfHwgMDtcbiAgICBpZiAoZW5kKSB7XG4gICAgICAgIHRoaXMuZW5kID0gZW5kO1xuICAgIH1cbiAgICB0aGlzLm9wdHMgPSBvcHRzO1xufVxuXG5cblxuRW5jb2RlRmV0Y2hhYmxlLnByb3RvdHlwZS5zbGljZSA9IGZ1bmN0aW9uKHMsIGwpIHtcbiAgICBpZiAocyA8IDApIHtcbiAgICAgICAgdGhyb3cgJ0JhZCBzbGljZSAnICsgcztcbiAgICB9XG5cbiAgICB2YXIgbnMgPSB0aGlzLnN0YXJ0LCBuZSA9IHRoaXMuZW5kO1xuICAgIGlmIChucyAmJiBzKSB7XG4gICAgICAgIG5zID0gbnMgKyBzO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIG5zID0gcyB8fCBucztcbiAgICB9XG4gICAgaWYgKGwgJiYgbnMpIHtcbiAgICAgICAgbmUgPSBucyArIGwgLSAxO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIG5lID0gbmUgfHwgbCAtIDE7XG4gICAgfVxuICAgIHJldHVybiBuZXcgRW5jb2RlRmV0Y2hhYmxlKHRoaXMudXJsLCBucywgbmUsIHRoaXMub3B0cyk7XG59XG5cbkVuY29kZUZldGNoYWJsZS5wcm90b3R5cGUuZmV0Y2hBc1RleHQgPSBmdW5jdGlvbihjYWxsYmFjaykge1xuICAgIHZhciBzZWxmID0gdGhpcztcbiAgICB2YXIgcmVxID0gbmV3IFhNTEh0dHBSZXF1ZXN0KCk7XG4gICAgdmFyIGxlbmd0aDtcbiAgICBzZWxmLnVybC5nZXRVUkxQcm9taXNlKCkudGhlbihmdW5jdGlvbih1cmwpIHtcbiAgICAgICAgcmVxLm9wZW4oJ0dFVCcsIHVybCwgdHJ1ZSk7XG5cbiAgICAgICAgaWYgKHNlbGYuZW5kKSB7XG4gICAgICAgICAgICBpZiAoc2VsZi5lbmQgLSBzZWxmLnN0YXJ0ID4gMTAwMDAwMDAwKSB7XG4gICAgICAgICAgICAgICAgdGhyb3cgJ01vbnN0ZXIgZmV0Y2ghJztcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJlcS5zZXRSZXF1ZXN0SGVhZGVyKCdSYW5nZScsICdieXRlcz0nICsgc2VsZi5zdGFydCArICctJyArIHNlbGYuZW5kKTtcbiAgICAgICAgICAgIGxlbmd0aCA9IHNlbGYuZW5kIC0gc2VsZi5zdGFydCArIDE7XG4gICAgICAgIH1cblxuICAgICAgICByZXEub25yZWFkeXN0YXRlY2hhbmdlID0gZnVuY3Rpb24oKSB7XG4gICAgICAgICAgICBpZiAocmVxLnJlYWR5U3RhdGUgPT0gNCkge1xuICAgICAgICAgICAgICAgIGlmIChyZXEuc3RhdHVzID09IDIwMCB8fCByZXEuc3RhdHVzID09IDIwNikge1xuICAgICAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2socmVxLnJlc3BvbnNlVGV4dCk7XG4gICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfTtcbiAgICAgICAgaWYgKHNlbGYub3B0cy5jcmVkZW50aWFscykge1xuICAgICAgICAgICAgcmVxLndpdGhDcmVkZW50aWFscyA9IHRydWU7XG4gICAgICAgIH1cbiAgICAgICAgcmVxLnNlbmQoJycpO1xuICAgIH0pLmNhdGNoKGZ1bmN0aW9uKGVycikge1xuICAgICAgICBjb25zb2xlLmxvZyhlcnIpO1xuICAgICAgICByZXR1cm4gY2FsbGJhY2sobnVsbCk7XG4gICAgfSk7XG59XG5cbkVuY29kZUZldGNoYWJsZS5wcm90b3R5cGUuc2FsdGVkID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXM7XG59XG5cbkVuY29kZUZldGNoYWJsZS5wcm90b3R5cGUuZmV0Y2ggPSBmdW5jdGlvbihjYWxsYmFjaywgYXR0ZW1wdCwgdHJ1bmNhdGVkTGVuZ3RoKSB7XG4gICAgdmFyIHNlbGYgPSB0aGlzO1xuXG4gICAgYXR0ZW1wdCA9IGF0dGVtcHQgfHwgMTtcbiAgICBpZiAoYXR0ZW1wdCA+IDMpIHtcbiAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwpO1xuICAgIH1cblxuICAgIHNlbGYudXJsLmdldFVSTFByb21pc2UoKS50aGVuKGZ1bmN0aW9uICh1cmwpIHtcbiAgICAgICAgdmFyIHJlcSA9IG5ldyBYTUxIdHRwUmVxdWVzdCgpO1xuICAgICAgICB2YXIgbGVuZ3RoO1xuICAgICAgICByZXEub3BlbignR0VUJywgdXJsLCB0cnVlKTtcbiAgICAgICAgcmVxLm92ZXJyaWRlTWltZVR5cGUoJ3RleHQvcGxhaW47IGNoYXJzZXQ9eC11c2VyLWRlZmluZWQnKTtcbiAgICAgICAgaWYgKHNlbGYuZW5kKSB7XG4gICAgICAgICAgICBpZiAoc2VsZi5lbmQgLSBzZWxmLnN0YXJ0ID4gMTAwMDAwMDAwKSB7XG4gICAgICAgICAgICAgICAgdGhyb3cgJ01vbnN0ZXIgZmV0Y2ghJztcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJlcS5zZXRSZXF1ZXN0SGVhZGVyKCdSYW5nZScsICdieXRlcz0nICsgc2VsZi5zdGFydCArICctJyArIHNlbGYuZW5kKTtcbiAgICAgICAgICAgIGxlbmd0aCA9IHNlbGYuZW5kIC0gc2VsZi5zdGFydCArIDE7XG4gICAgICAgIH1cbiAgICAgICAgcmVxLnJlc3BvbnNlVHlwZSA9ICdhcnJheWJ1ZmZlcic7XG4gICAgICAgIHJlcS5vbnJlYWR5c3RhdGVjaGFuZ2UgPSBmdW5jdGlvbigpIHtcbiAgICAgICAgICAgIGlmIChyZXEucmVhZHlTdGF0ZSA9PSA0KSB7XG4gICAgICAgICAgICAgICAgaWYgKHJlcS5zdGF0dXMgPT0gMjAwIHx8IHJlcS5zdGF0dXMgPT0gMjA2KSB7XG4gICAgICAgICAgICAgICAgICAgIGlmIChyZXEucmVzcG9uc2UpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciBibCA9IHJlcS5yZXNwb25zZS5ieXRlTGVuZ3RoO1xuICAgICAgICAgICAgICAgICAgICAgICAgaWYgKGxlbmd0aCAmJiBsZW5ndGggIT0gYmwgJiYgKCF0cnVuY2F0ZWRMZW5ndGggfHwgYmwgIT0gdHJ1bmNhdGVkTGVuZ3RoKSkge1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgIHJldHVybiBzZWxmLmZldGNoKGNhbGxiYWNrLCBhdHRlbXB0ICsgMSwgYmwpO1xuICAgICAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gY2FsbGJhY2socmVxLnJlc3BvbnNlKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIGlmIChyZXEubW96UmVzcG9uc2VBcnJheUJ1ZmZlcikge1xuICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKHJlcS5tb3pSZXNwb25zZUFycmF5QnVmZmVyKTtcbiAgICAgICAgICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICAgICAgICAgIHZhciByID0gcmVxLnJlc3BvbnNlVGV4dDtcbiAgICAgICAgICAgICAgICAgICAgICAgIGlmIChsZW5ndGggJiYgbGVuZ3RoICE9IHIubGVuZ3RoICYmICghdHJ1bmNhdGVkTGVuZ3RoIHx8IHIubGVuZ3RoICE9IHRydW5jYXRlZExlbmd0aCkpIHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gc2VsZi5mZXRjaChjYWxsYmFjaywgYXR0ZW1wdCArIDEsIHIubGVuZ3RoKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKGJzdHJpbmdUb0J1ZmZlcihyZXEucmVzcG9uc2VUZXh0KSk7XG4gICAgICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgICAgICByZXR1cm4gc2VsZi5mZXRjaChjYWxsYmFjaywgYXR0ZW1wdCArIDEpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfTtcbiAgICAgICAgaWYgKHNlbGYub3B0cy5jcmVkZW50aWFscykge1xuICAgICAgICAgICAgcmVxLndpdGhDcmVkZW50aWFscyA9IHRydWU7XG4gICAgICAgIH1cbiAgICAgICAgcmVxLnNlbmQoJycpO1xuICAgIH0pLmNhdGNoKGZ1bmN0aW9uKGVycikge1xuICAgICAgICBjb25zb2xlLmxvZyhlcnIpO1xuICAgIH0pO1xufVxuXG5mdW5jdGlvbiBic3RyaW5nVG9CdWZmZXIocmVzdWx0KSB7XG4gICAgaWYgKCFyZXN1bHQpIHtcbiAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgfVxuXG4gICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkocmVzdWx0Lmxlbmd0aCk7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBiYS5sZW5ndGg7ICsraSkge1xuICAgICAgICBiYVtpXSA9IHJlc3VsdC5jaGFyQ29kZUF0KGkpO1xuICAgIH1cbiAgICByZXR1cm4gYmEuYnVmZmVyO1xufVxuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIGxvb2t1cEVuY29kZVVSSTogbG9va3VwRW5jb2RlVVJJLFxuICAgICAgICBFbmNvZGVGZXRjaGFibGU6IEVuY29kZUZldGNoYWJsZVxuICAgIH07XG59XG4iLCIoZnVuY3Rpb24gKGdsb2JhbCl7XG4vKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBEYWxsaWFuY2UgR2Vub21lIEV4cGxvcmVyXG4vLyAoYykgVGhvbWFzIERvd24gMjAwNi0yMDE0XG4vL1xuLy8gZmV0Y2h3b3JrZXIuanNcbi8vXG5cblwidXNlIHN0cmljdFwiO1xuXG52YXIgYmluID0gcmVxdWlyZSgnLi9iaW4nKTtcbnZhciBiYW0gPSByZXF1aXJlKCcuL2JhbScpO1xudmFyIGJpZ3dpZyA9IHJlcXVpcmUoJy4vYmlnd2lnJyk7XG52YXIgZW5jb2RlID0gcmVxdWlyZSgnLi9lbmNvZGUnKTtcbnZhciB1dGlscyA9IHJlcXVpcmUoJy4vdXRpbHMnKTtcblxudmFyIFByb21pc2UgPSByZXF1aXJlKCdlczYtcHJvbWlzZScpLlByb21pc2U7XG5cbnZhciBjb25uZWN0aW9ucyA9IHt9O1xudmFyIHJlc29sdmVSZXNvbHZlcnMgPSB7fTtcblxudmFyIGlkU2VlZCA9IDA7XG5cbmdsb2JhbC5uZXdJRCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiAnY24nICsgKCsraWRTZWVkKTtcbn1cblxucG9zdE1lc3NhZ2Uoe3RhZzogJ2luaXQnfSk7XG5cbnNlbGYub25tZXNzYWdlID0gZnVuY3Rpb24oZXZlbnQpIHtcbiAgICB2YXIgZCA9IGV2ZW50LmRhdGE7XG4gICAgdmFyIGNvbW1hbmQgPSBldmVudC5kYXRhLmNvbW1hbmQ7XG4gICAgdmFyIHRhZyA9IGV2ZW50LmRhdGEudGFnO1xuXG4gICAgaWYgKCFjb21tYW5kKSB7XG4gICAgICAgIHZhciByciA9IHJlc29sdmVSZXNvbHZlcnNbdGFnXTtcbiAgICAgICAgaWYgKHJyKSB7XG4gICAgICAgICAgICBpZiAoZC5lcnIpIHtcbiAgICAgICAgICAgICAgICByci5yZWplY3QoZC5lcnIpO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICByci5yZXNvbHZlKGQudXJsKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICBcbiAgICAgICAgICAgIGRlbGV0ZSByZXNvbHZlUmVzb2x2ZXJzW3RhZ107XG4gICAgICAgIH1cbiAgICB9IGVsc2UgaWYgKGNvbW1hbmQgPT09ICdjb25uZWN0QkFNJykge1xuICAgICAgICB2YXIgaWQgPSBuZXdJRCgpO1xuICAgICAgICB2YXIgcmVzb2x2ZXI7XG4gICAgICAgIGlmIChkLnJlc29sdmVyKSB7XG4gICAgICAgICAgICByZXNvbHZlciA9IHByb3h5UmVzb2x2ZXIoZC5yZXNvbHZlcik7XG4gICAgICAgIH1cbiAgICAgICAgdmFyIGJhbUYsIGJhaUYsIGluZGV4Q2h1bmtzO1xuICAgICAgICBpZiAoZC5ibG9iKSB7XG4gICAgICAgICAgICBiYW1GID0gbmV3IGJpbi5CbG9iRmV0Y2hhYmxlKGQuYmxvYik7XG4gICAgICAgICAgICBiYWlGID0gbmV3IGJpbi5CbG9iRmV0Y2hhYmxlKGQuaW5kZXhCbG9iKTtcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGJhbUYgPSBuZXcgYmluLlVSTEZldGNoYWJsZShkLnVyaSwge2NyZWRlbnRpYWxzOiBkLmNyZWRlbnRpYWxzLCByZXNvbHZlcjogcmVzb2x2ZXJ9KTtcbiAgICAgICAgICAgIGJhaUYgPSBuZXcgYmluLlVSTEZldGNoYWJsZShkLmluZGV4VXJpLCB7Y3JlZGVudGlhbHM6IGQuY3JlZGVudGlhbHMsIHJlc29sdmVyOiByZXNvbHZlcn0pO1xuICAgICAgICAgICAgaW5kZXhDaHVua3MgPSBkLmluZGV4Q2h1bmtzO1xuICAgICAgICB9XG5cbiAgICAgICAgYmFtLm1ha2VCYW0oYmFtRiwgYmFpRiwgaW5kZXhDaHVua3MsIGZ1bmN0aW9uKGJhbU9iaiwgZXJyKSB7XG4gICAgICAgICAgICBpZiAoYmFtT2JqKSB7XG4gICAgICAgICAgICAgICAgY29ubmVjdGlvbnNbaWRdID0gbmV3IEJBTVdvcmtlckZldGNoZXIoYmFtT2JqKTtcbiAgICAgICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogaWR9KTtcbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCBlcnJvcjogZXJyIHx8IFwiQ291bGRuJ3QgZmV0Y2ggQkFNXCJ9KTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfSk7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAnY29ubmVjdEJCSScpIHtcbiAgICAgICAgdmFyIGlkID0gbmV3SUQoKTtcbiAgICAgICAgdmFyIHJlc29sdmVyO1xuICAgICAgICBpZiAoZC5yZXNvbHZlcikge1xuICAgICAgICAgICAgcmVzb2x2ZXIgPSBwcm94eVJlc29sdmVyKGQucmVzb2x2ZXIpO1xuICAgICAgICB9XG4gICAgICAgIHZhciBiYmk7XG4gICAgICAgIGlmIChkLmJsb2IpIHtcbiAgICAgICAgICAgIGJiaSA9IG5ldyBiaW4uQmxvYkZldGNoYWJsZShkLmJsb2IpO1xuICAgICAgICB9IGVsc2UgaWYgKGQudHJhbnNwb3J0ID09ICdlbmNvZGUnKSB7XG4gICAgICAgICAgICBiYmkgPSBuZXcgZW5jb2RlLkVuY29kZUZldGNoYWJsZShkLnVyaSwge2NyZWRlbnRpYWxzOiBkLmNyZWRlbnRpYWxzfSk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBiYmkgPSBuZXcgYmluLlVSTEZldGNoYWJsZShkLnVyaSwge2NyZWRlbnRpYWxzOiBkLmNyZWRlbnRpYWxzLCByZXNvbHZlcjogcmVzb2x2ZXJ9KTtcbiAgICAgICAgfVxuXG4gICAgICAgIGJpZ3dpZy5tYWtlQndnKGJiaSwgZnVuY3Rpb24oYndnLCBlcnIpIHtcbiAgICAgICAgICAgIGlmIChid2cpIHtcbiAgICAgICAgICAgICAgICBjb25uZWN0aW9uc1tpZF0gPSBuZXcgQkJJV29ya2VyRmV0Y2hlcihid2cpO1xuICAgICAgICAgICAgICAgIHBvc3RNZXNzYWdlKHt0YWc6IHRhZywgcmVzdWx0OiBpZH0pO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiBlcnIgfHwgXCJDb3VsZG4ndCBmZXRjaCBCQklcIn0pO1xuICAgICAgICAgICAgfVxuICAgICAgICB9LCBkLnVyaSk7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAndGV4dHhocicpIHtcbiAgICAgICAgdXRpbHMudGV4dFhIUihkLnVyaSwgZnVuY3Rpb24ocmVzcCwgZXJyKSB7XG4gICAgICAgICAgICBpZiAocmVzcCkge1xuICAgICAgICAgICAgICAgIHBvc3RNZXNzYWdlKHt0YWc6IHRhZywgcmVzdWx0OiByZXNwfSk7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIHBvc3RNZXNzYWdlKHt0YWc6IHRhZywgZXJyOiBlcnIgfHwgXCJDb3VsZG4ndCBmZXRjaCByZXNvdXJjZVwifSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH0pO1xuICAgIH0gZWxzZSBpZiAoY29tbWFuZCA9PT0gJ2ZldGNoJykge1xuICAgICAgICB2YXIgY29uID0gY29ubmVjdGlvbnNbZXZlbnQuZGF0YS5jb25uZWN0aW9uXTtcbiAgICAgICAgaWYgKCFjb24pIHtcbiAgICAgICAgICAgIHJldHVybiBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnTm8gc3VjaCBjb25uZWN0aW9uOiAnICsgZXZlbnQuZGF0YS5jb25uZWN0aW9ufSk7XG4gICAgICAgIH1cblxuICAgICAgICBjb24uZmV0Y2goZC50YWcsIGQuY2hyLCBkLm1pbiwgZC5tYXgsIGQuem9vbSwgZC5vcHRzKTtcbiAgICB9IGVsc2UgaWYgKGNvbW1hbmQgPT09ICdsZWFwJykge1xuICAgICAgICB2YXIgY29uID0gY29ubmVjdGlvbnNbZXZlbnQuZGF0YS5jb25uZWN0aW9uXTtcbiAgICAgICAgaWYgKCFjb24pIHtcbiAgICAgICAgICAgIHJldHVybiBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnTm8gc3VjaCBjb25uZWN0aW9uOiAnICsgZXZlbnQuZGF0YS5jb25uZWN0aW9ufSk7XG4gICAgICAgIH1cblxuICAgICAgICBjb24ubGVhcChkLnRhZywgZC5jaHIsIGQucG9zLCBkLmRpcik7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAncXVhbnRMZWFwJykge1xuICAgICAgICB2YXIgY29uID0gY29ubmVjdGlvbnNbZXZlbnQuZGF0YS5jb25uZWN0aW9uXTtcbiAgICAgICAgaWYgKCFjb24pIHtcbiAgICAgICAgICAgIHJldHVybiBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnTm8gc3VjaCBjb25uZWN0aW9uOiAnICsgZXZlbnQuZGF0YS5jb25uZWN0aW9ufSk7XG4gICAgICAgIH1cblxuICAgICAgICBjb24ucXVhbnRMZWFwKGQudGFnLCBkLmNociwgZC5wb3MsIGQuZGlyLCBkLnRocmVzaG9sZCwgZC51bmRlcik7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAnbWV0YScpIHtcbiAgICAgICAgdmFyIGNvbiA9IGNvbm5lY3Rpb25zW2V2ZW50LmRhdGEuY29ubmVjdGlvbl07XG4gICAgICAgIGlmICghY29uKSB7XG4gICAgICAgICAgICByZXR1cm4gcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCBlcnJvcjogJ05vIHN1Y2ggY29ubmVjdGlvbjogJyArIGV2ZW50LmRhdGEuY29ubmVjdGlvbn0pO1xuICAgICAgICB9XG5cbiAgICAgICAgY29uLm1ldGEoZC50YWcpO1xuICAgIH0gZWxzZSBpZiAoY29tbWFuZCA9PT0gJ3NlYXJjaCcpIHtcbiAgICAgICAgdmFyIGNvbiA9IGNvbm5lY3Rpb25zW2V2ZW50LmRhdGEuY29ubmVjdGlvbl07XG4gICAgICAgIGlmICghY29uKSB7XG4gICAgICAgICAgICByZXR1cm4gcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCBlcnJvcjogJ05vIHN1Y2ggY29ubmVjdGlvbjogJyArIGV2ZW50LmRhdGEuY29ubmVjdGlvbn0pO1xuICAgICAgICB9XG5cbiAgICAgICAgY29uLnNlYXJjaChkLnRhZywgZC5xdWVyeSwgZC5pbmRleCk7XG4gICAgfSBlbHNlIGlmIChjb21tYW5kID09PSAnZGF0ZScpIHtcbiAgICAgICAgcmV0dXJuIHBvc3RNZXNzYWdlKHt0YWc6IHRhZywgcmVzdWx0OiBEYXRlLm5vdygpfDB9KTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiAnQmFkIGNvbW1hbmQgJyArIGNvbW1hbmR9KTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIEJBTVdvcmtlckZldGNoZXIoYmFtKSB7XG4gICAgdGhpcy5iYW0gPSBiYW07XG59XG5cbkJBTVdvcmtlckZldGNoZXIucHJvdG90eXBlLmZldGNoID0gZnVuY3Rpb24odGFnLCBjaHIsIG1pbiwgbWF4LCB6b29tLCBvcHRzKSB7XG4gICAgb3B0cyA9IG9wdHMgfHwge307XG4gICAgdGhpcy5iYW0uZmV0Y2goY2hyLCBtaW4sIG1heCwgZnVuY3Rpb24ocmVjb3JkcywgZXJyKSB7XG4gICAgICAgIGlmIChyZWNvcmRzKSB7XG4gICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogcmVjb3JkcywgdGltZTogRGF0ZS5ub3coKXwwfSk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIGVycm9yOiBlcnJ9KTtcbiAgICAgICAgfVxuICAgIH0sIG9wdHMpO1xufVxuXG5mdW5jdGlvbiBCQklXb3JrZXJGZXRjaGVyKGJiaSkge1xuICAgIHRoaXMuYmJpID0gYmJpO1xufVxuXG5CQklXb3JrZXJGZXRjaGVyLnByb3RvdHlwZS5mZXRjaCA9IGZ1bmN0aW9uKHRhZywgY2hyLCBtaW4sIG1heCwgem9vbSkge1xuICAgIGlmICh0eXBlb2Yoem9vbSkgIT09ICdudW1iZXInKVxuICAgICAgICB6b29tID0gLTE7XG5cbiAgICB2YXIgZGF0YTtcbiAgICBpZiAoem9vbSA8IDApIHtcbiAgICAgICAgZGF0YSA9IHRoaXMuYmJpLmdldFVuem9vbWVkVmlldygpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIGRhdGEgPSB0aGlzLmJiaS5nZXRab29tZWRWaWV3KHpvb20pO1xuICAgIH1cblxuICAgIGRhdGEucmVhZFdpZ0RhdGEoY2hyLCBtaW4sIG1heCwgZnVuY3Rpb24oZmVhdHVyZXMpIHtcbiAgICAgICAgcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCByZXN1bHQ6IGZlYXR1cmVzfSk7XG4gICAgfSk7XG59XG5cbkJCSVdvcmtlckZldGNoZXIucHJvdG90eXBlLm1ldGEgPSBmdW5jdGlvbih0YWcpIHtcbiAgICB2YXIgc2NhbGVzID0gWzFdO1xuICAgIGZvciAodmFyIHogPSAwOyB6IDwgdGhpcy5iYmkuem9vbUxldmVscy5sZW5ndGg7ICsreikge1xuICAgICAgICBzY2FsZXMucHVzaCh0aGlzLmJiaS56b29tTGV2ZWxzW3pdLnJlZHVjdGlvbik7XG4gICAgfVxuXG4gICAgdmFyIHRoaXNCID0gdGhpcztcbiAgICB2YXIgbWV0YSA9IHt0eXBlOiB0aGlzLmJiaS50eXBlLFxuICAgICAgICAgICAgICAgIHpvb21MZXZlbHM6IHNjYWxlcyxcbiAgICAgICAgICAgICAgICBmaWVsZENvdW50OiB0aGlzLmJiaS5maWVsZENvdW50LFxuICAgICAgICAgICAgICAgIGRlZmluZWRGaWVsZENvdW50OiB0aGlzLmJiaS5kZWZpbmVkRmllbGRDb3VudCxcbiAgICAgICAgICAgICAgICBzY2hlbWE6IHRoaXMuYmJpLnNjaGVtYX07XG4gICAgaWYgKHRoaXMuYmJpLnR5cGUgPT09ICdiaWdiZWQnKSB7XG4gICAgICAgIHRoaXMuYmJpLmdldEV4dHJhSW5kaWNlcyhmdW5jdGlvbihlaSkge1xuICAgICAgICAgICAgaWYgKGVpKSB7XG4gICAgICAgICAgICAgICAgdGhpc0IuZXh0cmFJbmRpY2VzID0gZWk7XG4gICAgICAgICAgICAgICAgbWV0YS5leHRyYUluZGljZXMgPSBlaS5tYXAoZnVuY3Rpb24oaSkge3JldHVybiBpLmZpZWxkfSk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogbWV0YX0pO1xuICAgICAgICB9KTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogbWV0YX0pO1xuICAgIH1cbn1cblxuQkJJV29ya2VyRmV0Y2hlci5wcm90b3R5cGUubGVhcCA9IGZ1bmN0aW9uKHRhZywgY2hyLCBwb3MsIGRpcikge1xuICAgIHRoaXMuYmJpLmdldFVuem9vbWVkVmlldygpLmdldEZpcnN0QWRqYWNlbnQoY2hyLCBwb3MsIGRpciwgZnVuY3Rpb24ocmVzdWx0LCBlcnIpIHtcbiAgICAgICAgcG9zdE1lc3NhZ2Uoe3RhZzogdGFnLCByZXN1bHQ6IHJlc3VsdCwgZXJyb3I6IGVycn0pO1xuICAgIH0pO1xufVxuXG5CQklXb3JrZXJGZXRjaGVyLnByb3RvdHlwZS5xdWFudExlYXAgPSBmdW5jdGlvbih0YWcsIGNociwgcG9zLCBkaXIsIHRocmVzaG9sZCwgdW5kZXIpIHtcbiAgICB0aGlzLmJiaS50aHJlc2hvbGRTZWFyY2goY2hyLCBwb3MsIGRpciwgdGhyZXNob2xkLCBmdW5jdGlvbihyZXN1bHQsIGVycikge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogcmVzdWx0LCBlcnJvcjogZXJyfSk7XG4gICAgfSk7XG59XG5cbkJCSVdvcmtlckZldGNoZXIucHJvdG90eXBlLnNlYXJjaCA9IGZ1bmN0aW9uKHRhZywgcXVlcnksIGluZGV4KSB7XG4gICAgdmFyIGlzID0gdGhpcy5leHRyYUluZGljZXNbMF07XG4gICAgaXMubG9va3VwKHF1ZXJ5LCBmdW5jdGlvbihyZXN1bHQsIGVycikge1xuICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiB0YWcsIHJlc3VsdDogcmVzdWx0LCBlcnJvcjogZXJyfSk7XG4gICAgfSk7XG59XG5cbmZ1bmN0aW9uIHByb3h5UmVzb2x2ZXIodGFnKSB7XG4gICAgcmV0dXJuIGZ1bmN0aW9uKHVybCkge1xuICAgICAgICB2YXIgbGlkID0gbmV3SUQoKTtcbiAgICAgICAgcmV0dXJuIG5ldyBQcm9taXNlKGZ1bmN0aW9uIChyZXNvbHZlLCByZWplY3QpIHtcbiAgICAgICAgICAgIHJlc29sdmVSZXNvbHZlcnNbbGlkXSA9IHtyZXNvbHZlOiByZXNvbHZlLCByZWplY3Q6IHJlamVjdH07XG4gICAgICAgICAgICBwb3N0TWVzc2FnZSh7dGFnOiBsaWQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgY21kOiAncmVzb2x2ZScsXG4gICAgICAgICAgICAgICAgICAgICAgICAgcmVzb2x2ZXI6IHRhZyxcbiAgICAgICAgICAgICAgICAgICAgICAgICB1cmw6IHVybH0pO1xuICAgICAgICB9KTtcbiAgICB9XG59XG5cbn0pLmNhbGwodGhpcyx0eXBlb2Ygc2VsZiAhPT0gXCJ1bmRlZmluZWRcIiA/IHNlbGYgOiB0eXBlb2Ygd2luZG93ICE9PSBcInVuZGVmaW5lZFwiID8gd2luZG93IDoge30pIiwiLyogLSotIG1vZGU6IGphdmFzY3JpcHQ7IGMtYmFzaWMtb2Zmc2V0OiA0OyBpbmRlbnQtdGFicy1tb2RlOiBuaWwgLSotICovXG5cbi8vIFxuLy8gRGFsbGlhbmNlIEdlbm9tZSBFeHBsb3JlclxuLy8gKGMpIFRob21hcyBEb3duIDIwMDYtMjAxMVxuLy9cbi8vIGxoM3V0aWxzLmpzOiBjb21tb24gc3VwcG9ydCBmb3IgbGgzJ3MgZmlsZSBmb3JtYXRzXG4vL1xuXG5pZiAodHlwZW9mKHJlcXVpcmUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIHZhciBqc3psaWIgPSByZXF1aXJlKCdqc3psaWInKTtcbiAgICB2YXIganN6bGliX2luZmxhdGVfYnVmZmVyID0ganN6bGliLmluZmxhdGVCdWZmZXI7XG4gICAgdmFyIGFycmF5Q29weSA9IGpzemxpYi5hcnJheUNvcHk7XG59XG5cbmZ1bmN0aW9uIFZvYihiLCBvKSB7XG4gICAgdGhpcy5ibG9jayA9IGI7XG4gICAgdGhpcy5vZmZzZXQgPSBvO1xufVxuXG5Wb2IucHJvdG90eXBlLnRvU3RyaW5nID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuICcnICsgdGhpcy5ibG9jayArICc6JyArIHRoaXMub2Zmc2V0O1xufVxuXG5mdW5jdGlvbiByZWFkVm9iKGJhLCBvZmZzZXQpIHtcbiAgICB2YXIgYmxvY2sgPSAoKGJhW29mZnNldCs2XSAmIDB4ZmYpICogMHgxMDAwMDAwMDApICsgKChiYVtvZmZzZXQrNV0gJiAweGZmKSAqIDB4MTAwMDAwMCkgKyAoKGJhW29mZnNldCs0XSAmIDB4ZmYpICogMHgxMDAwMCkgKyAoKGJhW29mZnNldCszXSAmIDB4ZmYpICogMHgxMDApICsgKChiYVtvZmZzZXQrMl0gJiAweGZmKSk7XG4gICAgdmFyIGJpbnQgPSAoYmFbb2Zmc2V0KzFdIDw8IDgpIHwgKGJhW29mZnNldF0pO1xuICAgIGlmIChibG9jayA9PSAwICYmIGJpbnQgPT0gMCkge1xuICAgICAgICByZXR1cm4gbnVsbDsgIC8vIFNob3VsZCBvbmx5IGhhcHBlbiBpbiB0aGUgbGluZWFyIGluZGV4P1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHJldHVybiBuZXcgVm9iKGJsb2NrLCBiaW50KTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIHVuYmd6ZihkYXRhLCBsaW0pIHtcbiAgICBsaW0gPSBNYXRoLm1pbihsaW0gfHwgMSwgZGF0YS5ieXRlTGVuZ3RoIC0gNTApO1xuICAgIHZhciBvQmxvY2tMaXN0ID0gW107XG4gICAgdmFyIHB0ciA9IFswXTtcbiAgICB2YXIgdG90YWxTaXplID0gMDtcblxuICAgIHdoaWxlIChwdHJbMF0gPCBsaW0pIHtcbiAgICAgICAgdmFyIGJhID0gbmV3IFVpbnQ4QXJyYXkoZGF0YSwgcHRyWzBdLCAxMik7IC8vIEZJWE1FIGlzIHRoaXMgZW5vdWdoIGZvciBhbGwgY3JlZGlibGUgQkdaRiBibG9jayBoZWFkZXJzP1xuICAgICAgICB2YXIgeGxlbiA9IChiYVsxMV0gPDwgOCkgfCAoYmFbMTBdKTtcbiAgICAgICAgLy8gZGxvZygneGxlblsnICsgKHB0clswXSkgKyddPScgKyB4bGVuKTtcbiAgICAgICAgdmFyIHVuYyA9IGpzemxpYl9pbmZsYXRlX2J1ZmZlcihkYXRhLCAxMiArIHhsZW4gKyBwdHJbMF0sIE1hdGgubWluKDY1NTM2LCBkYXRhLmJ5dGVMZW5ndGggLSAxMiAtIHhsZW4gLSBwdHJbMF0pLCBwdHIpO1xuICAgICAgICBwdHJbMF0gKz0gODtcbiAgICAgICAgdG90YWxTaXplICs9IHVuYy5ieXRlTGVuZ3RoO1xuICAgICAgICBvQmxvY2tMaXN0LnB1c2godW5jKTtcbiAgICB9XG5cbiAgICBpZiAob0Jsb2NrTGlzdC5sZW5ndGggPT0gMSkge1xuICAgICAgICByZXR1cm4gb0Jsb2NrTGlzdFswXTtcbiAgICB9IGVsc2Uge1xuICAgICAgICB2YXIgb3V0ID0gbmV3IFVpbnQ4QXJyYXkodG90YWxTaXplKTtcbiAgICAgICAgdmFyIGN1cnNvciA9IDA7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgb0Jsb2NrTGlzdC5sZW5ndGg7ICsraSkge1xuICAgICAgICAgICAgdmFyIGIgPSBuZXcgVWludDhBcnJheShvQmxvY2tMaXN0W2ldKTtcbiAgICAgICAgICAgIGFycmF5Q29weShiLCAwLCBvdXQsIGN1cnNvciwgYi5sZW5ndGgpO1xuICAgICAgICAgICAgY3Vyc29yICs9IGIubGVuZ3RoO1xuICAgICAgICB9XG4gICAgICAgIHJldHVybiBvdXQuYnVmZmVyO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gQ2h1bmsobWludiwgbWF4dikge1xuICAgIHRoaXMubWludiA9IG1pbnY7IHRoaXMubWF4diA9IG1heHY7XG59XG5cblxuLy9cbi8vIEJpbm5pbmcgKHRyYW5zbGl0ZXJhdGVkIGZyb20gU0FNMS4zIHNwZWMpXG4vL1xuXG4vKiBjYWxjdWxhdGUgYmluIGdpdmVuIGFuIGFsaWdubWVudCBjb3ZlcmluZyBbYmVnLGVuZCkgKHplcm8tYmFzZWQsIGhhbGYtY2xvc2UtaGFsZi1vcGVuKSAqL1xuZnVuY3Rpb24gcmVnMmJpbihiZWcsIGVuZClcbntcbiAgICAtLWVuZDtcbiAgICBpZiAoYmVnPj4xNCA9PSBlbmQ+PjE0KSByZXR1cm4gKCgxPDwxNSktMSkvNyArIChiZWc+PjE0KTtcbiAgICBpZiAoYmVnPj4xNyA9PSBlbmQ+PjE3KSByZXR1cm4gKCgxPDwxMiktMSkvNyArIChiZWc+PjE3KTtcbiAgICBpZiAoYmVnPj4yMCA9PSBlbmQ+PjIwKSByZXR1cm4gKCgxPDw5KS0xKS83ICsgKGJlZz4+MjApO1xuICAgIGlmIChiZWc+PjIzID09IGVuZD4+MjMpIHJldHVybiAoKDE8PDYpLTEpLzcgKyAoYmVnPj4yMyk7XG4gICAgaWYgKGJlZz4+MjYgPT0gZW5kPj4yNikgcmV0dXJuICgoMTw8MyktMSkvNyArIChiZWc+PjI2KTtcbiAgICByZXR1cm4gMDtcbn1cblxuLyogY2FsY3VsYXRlIHRoZSBsaXN0IG9mIGJpbnMgdGhhdCBtYXkgb3ZlcmxhcCB3aXRoIHJlZ2lvbiBbYmVnLGVuZCkgKHplcm8tYmFzZWQpICovXG52YXIgTUFYX0JJTiA9ICgoKDE8PDE4KS0xKS83KTtcbmZ1bmN0aW9uIHJlZzJiaW5zKGJlZywgZW5kKSBcbntcbiAgICB2YXIgaSA9IDAsIGssIGxpc3QgPSBbXTtcbiAgICAtLWVuZDtcbiAgICBsaXN0LnB1c2goMCk7XG4gICAgZm9yIChrID0gMSArIChiZWc+PjI2KTsgayA8PSAxICsgKGVuZD4+MjYpOyArK2spIGxpc3QucHVzaChrKTtcbiAgICBmb3IgKGsgPSA5ICsgKGJlZz4+MjMpOyBrIDw9IDkgKyAoZW5kPj4yMyk7ICsraykgbGlzdC5wdXNoKGspO1xuICAgIGZvciAoayA9IDczICsgKGJlZz4+MjApOyBrIDw9IDczICsgKGVuZD4+MjApOyArK2spIGxpc3QucHVzaChrKTtcbiAgICBmb3IgKGsgPSA1ODUgKyAoYmVnPj4xNyk7IGsgPD0gNTg1ICsgKGVuZD4+MTcpOyArK2spIGxpc3QucHVzaChrKTtcbiAgICBmb3IgKGsgPSA0NjgxICsgKGJlZz4+MTQpOyBrIDw9IDQ2ODEgKyAoZW5kPj4xNCk7ICsraykgbGlzdC5wdXNoKGspO1xuICAgIHJldHVybiBsaXN0O1xufVxuXG5pZiAodHlwZW9mKG1vZHVsZSkgIT09ICd1bmRlZmluZWQnKSB7XG4gICAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgICAgIHVuYmd6ZjogdW5iZ3pmLFxuICAgICAgICByZWFkVm9iOiByZWFkVm9iLFxuICAgICAgICByZWcyYmluOiByZWcyYmluLFxuICAgICAgICByZWcyYmluczogcmVnMmJpbnMsXG4gICAgICAgIENodW5rOiBDaHVua1xuICAgIH07XG59IiwiLypcclxuICogQSBKYXZhU2NyaXB0IGltcGxlbWVudGF0aW9uIG9mIHRoZSBTZWN1cmUgSGFzaCBBbGdvcml0aG0sIFNIQS0xLCBhcyBkZWZpbmVkXHJcbiAqIGluIEZJUFMgMTgwLTFcclxuICogVmVyc2lvbiAyLjIgQ29weXJpZ2h0IFBhdWwgSm9obnN0b24gMjAwMCAtIDIwMDkuXHJcbiAqIE90aGVyIGNvbnRyaWJ1dG9yczogR3JlZyBIb2x0LCBBbmRyZXcgS2VwZXJ0LCBZZG5hciwgTG9zdGluZXRcclxuICogRGlzdHJpYnV0ZWQgdW5kZXIgdGhlIEJTRCBMaWNlbnNlXHJcbiAqIFNlZSBodHRwOi8vcGFqaG9tZS5vcmcudWsvY3J5cHQvbWQ1IGZvciBkZXRhaWxzLlxyXG4gKi9cclxuXHJcbiBcInVzZSBzdHJpY3RcIjtcclxuXHJcbi8qXHJcbiAqIENvbmZpZ3VyYWJsZSB2YXJpYWJsZXMuIFlvdSBtYXkgbmVlZCB0byB0d2VhayB0aGVzZSB0byBiZSBjb21wYXRpYmxlIHdpdGhcclxuICogdGhlIHNlcnZlci1zaWRlLCBidXQgdGhlIGRlZmF1bHRzIHdvcmsgaW4gbW9zdCBjYXNlcy5cclxuICovXHJcbnZhciBoZXhjYXNlID0gMDsgIC8qIGhleCBvdXRwdXQgZm9ybWF0LiAwIC0gbG93ZXJjYXNlOyAxIC0gdXBwZXJjYXNlICAgICAgICAqL1xyXG52YXIgYjY0cGFkICA9IFwiXCI7IC8qIGJhc2UtNjQgcGFkIGNoYXJhY3Rlci4gXCI9XCIgZm9yIHN0cmljdCBSRkMgY29tcGxpYW5jZSAgICovXHJcblxyXG4vKlxyXG4gKiBUaGVzZSBhcmUgdGhlIGZ1bmN0aW9ucyB5b3UnbGwgdXN1YWxseSB3YW50IHRvIGNhbGxcclxuICogVGhleSB0YWtlIHN0cmluZyBhcmd1bWVudHMgYW5kIHJldHVybiBlaXRoZXIgaGV4IG9yIGJhc2UtNjQgZW5jb2RlZCBzdHJpbmdzXHJcbiAqL1xyXG5mdW5jdGlvbiBoZXhfc2hhMShzKSAgICB7IHJldHVybiByc3RyMmhleChyc3RyX3NoYTEoc3RyMnJzdHJfdXRmOChzKSkpOyB9XHJcbmZ1bmN0aW9uIGI2NF9zaGExKHMpICAgIHsgcmV0dXJuIHJzdHIyYjY0KHJzdHJfc2hhMShzdHIycnN0cl91dGY4KHMpKSk7IH1cclxuZnVuY3Rpb24gYW55X3NoYTEocywgZSkgeyByZXR1cm4gcnN0cjJhbnkocnN0cl9zaGExKHN0cjJyc3RyX3V0ZjgocykpLCBlKTsgfVxyXG5mdW5jdGlvbiBoZXhfaG1hY19zaGExKGssIGQpXHJcbiAgeyByZXR1cm4gcnN0cjJoZXgocnN0cl9obWFjX3NoYTEoc3RyMnJzdHJfdXRmOChrKSwgc3RyMnJzdHJfdXRmOChkKSkpOyB9XHJcbmZ1bmN0aW9uIGI2NF9obWFjX3NoYTEoaywgZClcclxuICB7IHJldHVybiByc3RyMmI2NChyc3RyX2htYWNfc2hhMShzdHIycnN0cl91dGY4KGspLCBzdHIycnN0cl91dGY4KGQpKSk7IH1cclxuZnVuY3Rpb24gYW55X2htYWNfc2hhMShrLCBkLCBlKVxyXG4gIHsgcmV0dXJuIHJzdHIyYW55KHJzdHJfaG1hY19zaGExKHN0cjJyc3RyX3V0ZjgoayksIHN0cjJyc3RyX3V0ZjgoZCkpLCBlKTsgfVxyXG5cclxuLypcclxuICogUGVyZm9ybSBhIHNpbXBsZSBzZWxmLXRlc3QgdG8gc2VlIGlmIHRoZSBWTSBpcyB3b3JraW5nXHJcbiAqL1xyXG5mdW5jdGlvbiBzaGExX3ZtX3Rlc3QoKVxyXG57XHJcbiAgcmV0dXJuIGhleF9zaGExKFwiYWJjXCIpLnRvTG93ZXJDYXNlKCkgPT0gXCJhOTk5M2UzNjQ3MDY4MTZhYmEzZTI1NzE3ODUwYzI2YzljZDBkODlkXCI7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENhbGN1bGF0ZSB0aGUgU0hBMSBvZiBhIHJhdyBzdHJpbmdcclxuICovXHJcbmZ1bmN0aW9uIHJzdHJfc2hhMShzKVxyXG57XHJcbiAgcmV0dXJuIGJpbmIycnN0cihiaW5iX3NoYTEocnN0cjJiaW5iKHMpLCBzLmxlbmd0aCAqIDgpKTtcclxufVxyXG5cclxuLypcclxuICogQ2FsY3VsYXRlIHRoZSBITUFDLVNIQTEgb2YgYSBrZXkgYW5kIHNvbWUgZGF0YSAocmF3IHN0cmluZ3MpXHJcbiAqL1xyXG5mdW5jdGlvbiByc3RyX2htYWNfc2hhMShrZXksIGRhdGEpXHJcbntcclxuICB2YXIgYmtleSA9IHJzdHIyYmluYihrZXkpO1xyXG4gIGlmKGJrZXkubGVuZ3RoID4gMTYpIGJrZXkgPSBiaW5iX3NoYTEoYmtleSwga2V5Lmxlbmd0aCAqIDgpO1xyXG5cclxuICB2YXIgaXBhZCA9IEFycmF5KDE2KSwgb3BhZCA9IEFycmF5KDE2KTtcclxuICBmb3IodmFyIGkgPSAwOyBpIDwgMTY7IGkrKylcclxuICB7XHJcbiAgICBpcGFkW2ldID0gYmtleVtpXSBeIDB4MzYzNjM2MzY7XHJcbiAgICBvcGFkW2ldID0gYmtleVtpXSBeIDB4NUM1QzVDNUM7XHJcbiAgfVxyXG5cclxuICB2YXIgaGFzaCA9IGJpbmJfc2hhMShpcGFkLmNvbmNhdChyc3RyMmJpbmIoZGF0YSkpLCA1MTIgKyBkYXRhLmxlbmd0aCAqIDgpO1xyXG4gIHJldHVybiBiaW5iMnJzdHIoYmluYl9zaGExKG9wYWQuY29uY2F0KGhhc2gpLCA1MTIgKyAxNjApKTtcclxufVxyXG5cclxuLypcclxuICogQ29udmVydCBhIHJhdyBzdHJpbmcgdG8gYSBoZXggc3RyaW5nXHJcbiAqL1xyXG5mdW5jdGlvbiByc3RyMmhleChpbnB1dClcclxue1xyXG4gIC8vIHRyeSB7IGhleGNhc2UgfSBjYXRjaChlKSB7IGhleGNhc2U9MDsgfVxyXG4gIHZhciBoZXhfdGFiID0gaGV4Y2FzZSA/IFwiMDEyMzQ1Njc4OUFCQ0RFRlwiIDogXCIwMTIzNDU2Nzg5YWJjZGVmXCI7XHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgdmFyIHg7XHJcbiAgZm9yKHZhciBpID0gMDsgaSA8IGlucHV0Lmxlbmd0aDsgaSsrKVxyXG4gIHtcclxuICAgIHggPSBpbnB1dC5jaGFyQ29kZUF0KGkpO1xyXG4gICAgb3V0cHV0ICs9IGhleF90YWIuY2hhckF0KCh4ID4+PiA0KSAmIDB4MEYpXHJcbiAgICAgICAgICAgKyAgaGV4X3RhYi5jaGFyQXQoIHggICAgICAgICYgMHgwRik7XHJcbiAgfVxyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENvbnZlcnQgYSByYXcgc3RyaW5nIHRvIGEgYmFzZS02NCBzdHJpbmdcclxuICovXHJcbmZ1bmN0aW9uIHJzdHIyYjY0KGlucHV0KVxyXG57XHJcbiAgLy8gdHJ5IHsgYjY0cGFkIH0gY2F0Y2goZSkgeyBiNjRwYWQ9Jyc7IH1cclxuICB2YXIgdGFiID0gXCJBQkNERUZHSElKS0xNTk9QUVJTVFVWV1hZWmFiY2RlZmdoaWprbG1ub3BxcnN0dXZ3eHl6MDEyMzQ1Njc4OSsvXCI7XHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgdmFyIGxlbiA9IGlucHV0Lmxlbmd0aDtcclxuICBmb3IodmFyIGkgPSAwOyBpIDwgbGVuOyBpICs9IDMpXHJcbiAge1xyXG4gICAgdmFyIHRyaXBsZXQgPSAoaW5wdXQuY2hhckNvZGVBdChpKSA8PCAxNilcclxuICAgICAgICAgICAgICAgIHwgKGkgKyAxIDwgbGVuID8gaW5wdXQuY2hhckNvZGVBdChpKzEpIDw8IDggOiAwKVxyXG4gICAgICAgICAgICAgICAgfCAoaSArIDIgPCBsZW4gPyBpbnB1dC5jaGFyQ29kZUF0KGkrMikgICAgICA6IDApO1xyXG4gICAgZm9yKHZhciBqID0gMDsgaiA8IDQ7IGorKylcclxuICAgIHtcclxuICAgICAgaWYoaSAqIDggKyBqICogNiA+IGlucHV0Lmxlbmd0aCAqIDgpIG91dHB1dCArPSBiNjRwYWQ7XHJcbiAgICAgIGVsc2Ugb3V0cHV0ICs9IHRhYi5jaGFyQXQoKHRyaXBsZXQgPj4+IDYqKDMtaikpICYgMHgzRik7XHJcbiAgICB9XHJcbiAgfVxyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENvbnZlcnQgYSByYXcgc3RyaW5nIHRvIGFuIGFyYml0cmFyeSBzdHJpbmcgZW5jb2RpbmdcclxuICovXHJcbmZ1bmN0aW9uIHJzdHIyYW55KGlucHV0LCBlbmNvZGluZylcclxue1xyXG4gIHZhciBkaXZpc29yID0gZW5jb2RpbmcubGVuZ3RoO1xyXG4gIHZhciByZW1haW5kZXJzID0gQXJyYXkoKTtcclxuICB2YXIgaSwgcSwgeCwgcXVvdGllbnQ7XHJcblxyXG4gIC8qIENvbnZlcnQgdG8gYW4gYXJyYXkgb2YgMTYtYml0IGJpZy1lbmRpYW4gdmFsdWVzLCBmb3JtaW5nIHRoZSBkaXZpZGVuZCAqL1xyXG4gIHZhciBkaXZpZGVuZCA9IEFycmF5KE1hdGguY2VpbChpbnB1dC5sZW5ndGggLyAyKSk7XHJcbiAgZm9yKGkgPSAwOyBpIDwgZGl2aWRlbmQubGVuZ3RoOyBpKyspXHJcbiAge1xyXG4gICAgZGl2aWRlbmRbaV0gPSAoaW5wdXQuY2hhckNvZGVBdChpICogMikgPDwgOCkgfCBpbnB1dC5jaGFyQ29kZUF0KGkgKiAyICsgMSk7XHJcbiAgfVxyXG5cclxuICAvKlxyXG4gICAqIFJlcGVhdGVkbHkgcGVyZm9ybSBhIGxvbmcgZGl2aXNpb24uIFRoZSBiaW5hcnkgYXJyYXkgZm9ybXMgdGhlIGRpdmlkZW5kLFxyXG4gICAqIHRoZSBsZW5ndGggb2YgdGhlIGVuY29kaW5nIGlzIHRoZSBkaXZpc29yLiBPbmNlIGNvbXB1dGVkLCB0aGUgcXVvdGllbnRcclxuICAgKiBmb3JtcyB0aGUgZGl2aWRlbmQgZm9yIHRoZSBuZXh0IHN0ZXAuIFdlIHN0b3Agd2hlbiB0aGUgZGl2aWRlbmQgaXMgemVyby5cclxuICAgKiBBbGwgcmVtYWluZGVycyBhcmUgc3RvcmVkIGZvciBsYXRlciB1c2UuXHJcbiAgICovXHJcbiAgd2hpbGUoZGl2aWRlbmQubGVuZ3RoID4gMClcclxuICB7XHJcbiAgICBxdW90aWVudCA9IEFycmF5KCk7XHJcbiAgICB4ID0gMDtcclxuICAgIGZvcihpID0gMDsgaSA8IGRpdmlkZW5kLmxlbmd0aDsgaSsrKVxyXG4gICAge1xyXG4gICAgICB4ID0gKHggPDwgMTYpICsgZGl2aWRlbmRbaV07XHJcbiAgICAgIHEgPSBNYXRoLmZsb29yKHggLyBkaXZpc29yKTtcclxuICAgICAgeCAtPSBxICogZGl2aXNvcjtcclxuICAgICAgaWYocXVvdGllbnQubGVuZ3RoID4gMCB8fCBxID4gMClcclxuICAgICAgICBxdW90aWVudFtxdW90aWVudC5sZW5ndGhdID0gcTtcclxuICAgIH1cclxuICAgIHJlbWFpbmRlcnNbcmVtYWluZGVycy5sZW5ndGhdID0geDtcclxuICAgIGRpdmlkZW5kID0gcXVvdGllbnQ7XHJcbiAgfVxyXG5cclxuICAvKiBDb252ZXJ0IHRoZSByZW1haW5kZXJzIHRvIHRoZSBvdXRwdXQgc3RyaW5nICovXHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgZm9yKGkgPSByZW1haW5kZXJzLmxlbmd0aCAtIDE7IGkgPj0gMDsgaS0tKVxyXG4gICAgb3V0cHV0ICs9IGVuY29kaW5nLmNoYXJBdChyZW1haW5kZXJzW2ldKTtcclxuXHJcbiAgLyogQXBwZW5kIGxlYWRpbmcgemVybyBlcXVpdmFsZW50cyAqL1xyXG4gIHZhciBmdWxsX2xlbmd0aCA9IE1hdGguY2VpbChpbnB1dC5sZW5ndGggKiA4IC9cclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgKE1hdGgubG9nKGVuY29kaW5nLmxlbmd0aCkgLyBNYXRoLmxvZygyKSkpXHJcbiAgZm9yKGkgPSBvdXRwdXQubGVuZ3RoOyBpIDwgZnVsbF9sZW5ndGg7IGkrKylcclxuICAgIG91dHB1dCA9IGVuY29kaW5nWzBdICsgb3V0cHV0O1xyXG5cclxuICByZXR1cm4gb3V0cHV0O1xyXG59XHJcblxyXG4vKlxyXG4gKiBFbmNvZGUgYSBzdHJpbmcgYXMgdXRmLTguXHJcbiAqIEZvciBlZmZpY2llbmN5LCB0aGlzIGFzc3VtZXMgdGhlIGlucHV0IGlzIHZhbGlkIHV0Zi0xNi5cclxuICovXHJcbmZ1bmN0aW9uIHN0cjJyc3RyX3V0ZjgoaW5wdXQpXHJcbntcclxuICB2YXIgb3V0cHV0ID0gXCJcIjtcclxuICB2YXIgaSA9IC0xO1xyXG4gIHZhciB4LCB5O1xyXG5cclxuICB3aGlsZSgrK2kgPCBpbnB1dC5sZW5ndGgpXHJcbiAge1xyXG4gICAgLyogRGVjb2RlIHV0Zi0xNiBzdXJyb2dhdGUgcGFpcnMgKi9cclxuICAgIHggPSBpbnB1dC5jaGFyQ29kZUF0KGkpO1xyXG4gICAgeSA9IGkgKyAxIDwgaW5wdXQubGVuZ3RoID8gaW5wdXQuY2hhckNvZGVBdChpICsgMSkgOiAwO1xyXG4gICAgaWYoMHhEODAwIDw9IHggJiYgeCA8PSAweERCRkYgJiYgMHhEQzAwIDw9IHkgJiYgeSA8PSAweERGRkYpXHJcbiAgICB7XHJcbiAgICAgIHggPSAweDEwMDAwICsgKCh4ICYgMHgwM0ZGKSA8PCAxMCkgKyAoeSAmIDB4MDNGRik7XHJcbiAgICAgIGkrKztcclxuICAgIH1cclxuXHJcbiAgICAvKiBFbmNvZGUgb3V0cHV0IGFzIHV0Zi04ICovXHJcbiAgICBpZih4IDw9IDB4N0YpXHJcbiAgICAgIG91dHB1dCArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKHgpO1xyXG4gICAgZWxzZSBpZih4IDw9IDB4N0ZGKVxyXG4gICAgICBvdXRwdXQgKz0gU3RyaW5nLmZyb21DaGFyQ29kZSgweEMwIHwgKCh4ID4+PiA2ICkgJiAweDFGKSxcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgMHg4MCB8ICggeCAgICAgICAgICYgMHgzRikpO1xyXG4gICAgZWxzZSBpZih4IDw9IDB4RkZGRilcclxuICAgICAgb3V0cHV0ICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoMHhFMCB8ICgoeCA+Pj4gMTIpICYgMHgwRiksXHJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIDB4ODAgfCAoKHggPj4+IDYgKSAmIDB4M0YpLFxyXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAweDgwIHwgKCB4ICAgICAgICAgJiAweDNGKSk7XHJcbiAgICBlbHNlIGlmKHggPD0gMHgxRkZGRkYpXHJcbiAgICAgIG91dHB1dCArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKDB4RjAgfCAoKHggPj4+IDE4KSAmIDB4MDcpLFxyXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAweDgwIHwgKCh4ID4+PiAxMikgJiAweDNGKSxcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgMHg4MCB8ICgoeCA+Pj4gNiApICYgMHgzRiksXHJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIDB4ODAgfCAoIHggICAgICAgICAmIDB4M0YpKTtcclxuICB9XHJcbiAgcmV0dXJuIG91dHB1dDtcclxufVxyXG5cclxuLypcclxuICogRW5jb2RlIGEgc3RyaW5nIGFzIHV0Zi0xNlxyXG4gKi9cclxuZnVuY3Rpb24gc3RyMnJzdHJfdXRmMTZsZShpbnB1dClcclxue1xyXG4gIHZhciBvdXRwdXQgPSBcIlwiO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBpbnB1dC5sZW5ndGg7IGkrKylcclxuICAgIG91dHB1dCArPSBTdHJpbmcuZnJvbUNoYXJDb2RlKCBpbnB1dC5jaGFyQ29kZUF0KGkpICAgICAgICAmIDB4RkYsXHJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAoaW5wdXQuY2hhckNvZGVBdChpKSA+Pj4gOCkgJiAweEZGKTtcclxuICByZXR1cm4gb3V0cHV0O1xyXG59XHJcblxyXG5mdW5jdGlvbiBzdHIycnN0cl91dGYxNmJlKGlucHV0KVxyXG57XHJcbiAgdmFyIG91dHB1dCA9IFwiXCI7XHJcbiAgZm9yKHZhciBpID0gMDsgaSA8IGlucHV0Lmxlbmd0aDsgaSsrKVxyXG4gICAgb3V0cHV0ICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoKGlucHV0LmNoYXJDb2RlQXQoaSkgPj4+IDgpICYgMHhGRixcclxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbnB1dC5jaGFyQ29kZUF0KGkpICAgICAgICAmIDB4RkYpO1xyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENvbnZlcnQgYSByYXcgc3RyaW5nIHRvIGFuIGFycmF5IG9mIGJpZy1lbmRpYW4gd29yZHNcclxuICogQ2hhcmFjdGVycyA+MjU1IGhhdmUgdGhlaXIgaGlnaC1ieXRlIHNpbGVudGx5IGlnbm9yZWQuXHJcbiAqL1xyXG5mdW5jdGlvbiByc3RyMmJpbmIoaW5wdXQpXHJcbntcclxuICB2YXIgb3V0cHV0ID0gQXJyYXkoaW5wdXQubGVuZ3RoID4+IDIpO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBvdXRwdXQubGVuZ3RoOyBpKyspXHJcbiAgICBvdXRwdXRbaV0gPSAwO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBpbnB1dC5sZW5ndGggKiA4OyBpICs9IDgpXHJcbiAgICBvdXRwdXRbaT4+NV0gfD0gKGlucHV0LmNoYXJDb2RlQXQoaSAvIDgpICYgMHhGRikgPDwgKDI0IC0gaSAlIDMyKTtcclxuICByZXR1cm4gb3V0cHV0O1xyXG59XHJcblxyXG4vKlxyXG4gKiBDb252ZXJ0IGFuIGFycmF5IG9mIGJpZy1lbmRpYW4gd29yZHMgdG8gYSBzdHJpbmdcclxuICovXHJcbmZ1bmN0aW9uIGJpbmIycnN0cihpbnB1dClcclxue1xyXG4gIHZhciBvdXRwdXQgPSBcIlwiO1xyXG4gIGZvcih2YXIgaSA9IDA7IGkgPCBpbnB1dC5sZW5ndGggKiAzMjsgaSArPSA4KVxyXG4gICAgb3V0cHV0ICs9IFN0cmluZy5mcm9tQ2hhckNvZGUoKGlucHV0W2k+PjVdID4+PiAoMjQgLSBpICUgMzIpKSAmIDB4RkYpO1xyXG4gIHJldHVybiBvdXRwdXQ7XHJcbn1cclxuXHJcbi8qXHJcbiAqIENhbGN1bGF0ZSB0aGUgU0hBLTEgb2YgYW4gYXJyYXkgb2YgYmlnLWVuZGlhbiB3b3JkcywgYW5kIGEgYml0IGxlbmd0aFxyXG4gKi9cclxuZnVuY3Rpb24gYmluYl9zaGExKHgsIGxlbilcclxue1xyXG4gIC8qIGFwcGVuZCBwYWRkaW5nICovXHJcbiAgeFtsZW4gPj4gNV0gfD0gMHg4MCA8PCAoMjQgLSBsZW4gJSAzMik7XHJcbiAgeFsoKGxlbiArIDY0ID4+IDkpIDw8IDQpICsgMTVdID0gbGVuO1xyXG5cclxuICB2YXIgdyA9IEFycmF5KDgwKTtcclxuICB2YXIgYSA9ICAxNzMyNTg0MTkzO1xyXG4gIHZhciBiID0gLTI3MTczMzg3OTtcclxuICB2YXIgYyA9IC0xNzMyNTg0MTk0O1xyXG4gIHZhciBkID0gIDI3MTczMzg3ODtcclxuICB2YXIgZSA9IC0xMDA5NTg5Nzc2O1xyXG5cclxuICBmb3IodmFyIGkgPSAwOyBpIDwgeC5sZW5ndGg7IGkgKz0gMTYpXHJcbiAge1xyXG4gICAgdmFyIG9sZGEgPSBhO1xyXG4gICAgdmFyIG9sZGIgPSBiO1xyXG4gICAgdmFyIG9sZGMgPSBjO1xyXG4gICAgdmFyIG9sZGQgPSBkO1xyXG4gICAgdmFyIG9sZGUgPSBlO1xyXG5cclxuICAgIGZvcih2YXIgaiA9IDA7IGogPCA4MDsgaisrKVxyXG4gICAge1xyXG4gICAgICBpZihqIDwgMTYpIHdbal0gPSB4W2kgKyBqXTtcclxuICAgICAgZWxzZSB3W2pdID0gYml0X3JvbCh3W2otM10gXiB3W2otOF0gXiB3W2otMTRdIF4gd1tqLTE2XSwgMSk7XHJcbiAgICAgIHZhciB0ID0gc2FmZV9hZGQoc2FmZV9hZGQoYml0X3JvbChhLCA1KSwgc2hhMV9mdChqLCBiLCBjLCBkKSksXHJcbiAgICAgICAgICAgICAgICAgICAgICAgc2FmZV9hZGQoc2FmZV9hZGQoZSwgd1tqXSksIHNoYTFfa3QoaikpKTtcclxuICAgICAgZSA9IGQ7XHJcbiAgICAgIGQgPSBjO1xyXG4gICAgICBjID0gYml0X3JvbChiLCAzMCk7XHJcbiAgICAgIGIgPSBhO1xyXG4gICAgICBhID0gdDtcclxuICAgIH1cclxuXHJcbiAgICBhID0gc2FmZV9hZGQoYSwgb2xkYSk7XHJcbiAgICBiID0gc2FmZV9hZGQoYiwgb2xkYik7XHJcbiAgICBjID0gc2FmZV9hZGQoYywgb2xkYyk7XHJcbiAgICBkID0gc2FmZV9hZGQoZCwgb2xkZCk7XHJcbiAgICBlID0gc2FmZV9hZGQoZSwgb2xkZSk7XHJcbiAgfVxyXG4gIHJldHVybiBBcnJheShhLCBiLCBjLCBkLCBlKTtcclxuXHJcbn1cclxuXHJcbi8qXHJcbiAqIFBlcmZvcm0gdGhlIGFwcHJvcHJpYXRlIHRyaXBsZXQgY29tYmluYXRpb24gZnVuY3Rpb24gZm9yIHRoZSBjdXJyZW50XHJcbiAqIGl0ZXJhdGlvblxyXG4gKi9cclxuZnVuY3Rpb24gc2hhMV9mdCh0LCBiLCBjLCBkKVxyXG57XHJcbiAgaWYodCA8IDIwKSByZXR1cm4gKGIgJiBjKSB8ICgofmIpICYgZCk7XHJcbiAgaWYodCA8IDQwKSByZXR1cm4gYiBeIGMgXiBkO1xyXG4gIGlmKHQgPCA2MCkgcmV0dXJuIChiICYgYykgfCAoYiAmIGQpIHwgKGMgJiBkKTtcclxuICByZXR1cm4gYiBeIGMgXiBkO1xyXG59XHJcblxyXG4vKlxyXG4gKiBEZXRlcm1pbmUgdGhlIGFwcHJvcHJpYXRlIGFkZGl0aXZlIGNvbnN0YW50IGZvciB0aGUgY3VycmVudCBpdGVyYXRpb25cclxuICovXHJcbmZ1bmN0aW9uIHNoYTFfa3QodClcclxue1xyXG4gIHJldHVybiAodCA8IDIwKSA/ICAxNTE4NTAwMjQ5IDogKHQgPCA0MCkgPyAgMTg1OTc3NTM5MyA6XHJcbiAgICAgICAgICh0IDwgNjApID8gLTE4OTQwMDc1ODggOiAtODk5NDk3NTE0O1xyXG59XHJcblxyXG4vKlxyXG4gKiBBZGQgaW50ZWdlcnMsIHdyYXBwaW5nIGF0IDJeMzIuIFRoaXMgdXNlcyAxNi1iaXQgb3BlcmF0aW9ucyBpbnRlcm5hbGx5XHJcbiAqIHRvIHdvcmsgYXJvdW5kIGJ1Z3MgaW4gc29tZSBKUyBpbnRlcnByZXRlcnMuXHJcbiAqL1xyXG5mdW5jdGlvbiBzYWZlX2FkZCh4LCB5KVxyXG57XHJcbiAgdmFyIGxzdyA9ICh4ICYgMHhGRkZGKSArICh5ICYgMHhGRkZGKTtcclxuICB2YXIgbXN3ID0gKHggPj4gMTYpICsgKHkgPj4gMTYpICsgKGxzdyA+PiAxNik7XHJcbiAgcmV0dXJuIChtc3cgPDwgMTYpIHwgKGxzdyAmIDB4RkZGRik7XHJcbn1cclxuXHJcbi8qXHJcbiAqIEJpdHdpc2Ugcm90YXRlIGEgMzItYml0IG51bWJlciB0byB0aGUgbGVmdC5cclxuICovXHJcbmZ1bmN0aW9uIGJpdF9yb2wobnVtLCBjbnQpXHJcbntcclxuICByZXR1cm4gKG51bSA8PCBjbnQpIHwgKG51bSA+Pj4gKDMyIC0gY250KSk7XHJcbn1cclxuXHJcbmlmICh0eXBlb2YobW9kdWxlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcclxuICBtb2R1bGUuZXhwb3J0cyA9IHtcclxuICAgIGI2NF9zaGExOiBiNjRfc2hhMSxcclxuICAgIGhleF9zaGExOiBoZXhfc2hhMVxyXG4gIH1cclxufVxyXG4iLCIvKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBEYWxsaWFuY2UgR2Vub21lIEV4cGxvcmVyXG4vLyAoYykgVGhvbWFzIERvd24gMjAwNi0yMDEwXG4vL1xuLy8gc3BhbnMuanM6IEphdmFTY3JpcHQgSW50c2V0L0xvY2F0aW9uIHBvcnQuXG4vL1xuXG5cInVzZSBzdHJpY3RcIjtcblxuXG5mdW5jdGlvbiBSYW5nZShtaW4sIG1heClcbntcbiAgICBpZiAodHlwZW9mKG1pbikgIT0gJ251bWJlcicgfHwgdHlwZW9mKG1heCkgIT0gJ251bWJlcicpXG4gICAgICAgIHRocm93ICdCYWQgcmFuZ2UgJyArIG1pbiArICcsJyArIG1heDtcbiAgICB0aGlzLl9taW4gPSBtaW47XG4gICAgdGhpcy5fbWF4ID0gbWF4O1xufVxuXG5SYW5nZS5wcm90b3R5cGUubWluID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX21pbjtcbn1cblxuUmFuZ2UucHJvdG90eXBlLm1heCA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9tYXg7XG59XG5cblJhbmdlLnByb3RvdHlwZS5jb250YWlucyA9IGZ1bmN0aW9uKHBvcykge1xuICAgIHJldHVybiBwb3MgPj0gdGhpcy5fbWluICYmIHBvcyA8PSB0aGlzLl9tYXg7XG59XG5cblJhbmdlLnByb3RvdHlwZS5pc0NvbnRpZ3VvdXMgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdHJ1ZTtcbn1cblxuUmFuZ2UucHJvdG90eXBlLnJhbmdlcyA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiBbdGhpc107XG59XG5cblJhbmdlLnByb3RvdHlwZS5fcHVzaFJhbmdlcyA9IGZ1bmN0aW9uKHJhbmdlcykge1xuICAgIHJhbmdlcy5wdXNoKHRoaXMpO1xufVxuXG5SYW5nZS5wcm90b3R5cGUudG9TdHJpbmcgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gJ1snICsgdGhpcy5fbWluICsgJy0nICsgdGhpcy5fbWF4ICsgJ10nO1xufVxuXG5mdW5jdGlvbiBfQ29tcG91bmQocmFuZ2VzKSB7XG4gICAgLy8gZ2l2ZW46IGEgc2V0IG9mIHVuc29ydGVkIHBvc3NpYmx5IG92ZXJsYXBwaW5nIHJhbmdlc1xuICAgIC8vIHNvcnQgdGhlIGlucHV0IHJhbmdlc1xuICAgIHZhciBzb3J0ZWQgPSByYW5nZXMuc29ydChfcmFuZ2VPcmRlcik7XG4gICAgLy8gbWVyZ2Ugb3ZlcmxhcHMgYmV0d2VlbiBhZGphY2VudCByYW5nZXNcbiAgICB2YXIgbWVyZ2VkID0gW107XG4gICAgdmFyIGN1cnJlbnQgPSBzb3J0ZWQuc2hpZnQoKTtcbiAgICBzb3J0ZWQuZm9yRWFjaChmdW5jdGlvbihyYW5nZSkge1xuICAgICAgICBpZiAocmFuZ2UuX21pbiA8PSBjdXJyZW50Ll9tYXgpIHtcbiAgICAgICAgICAgIGlmIChyYW5nZS5fbWF4ID4gY3VycmVudC5fbWF4KSB7XG4gICAgICAgICAgICAgICAgY3VycmVudC5fbWF4ID0gcmFuZ2UuX21heDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBlbHNlIHtcbiAgICAgICAgICAgIG1lcmdlZC5wdXNoKGN1cnJlbnQpO1xuICAgICAgICAgICAgY3VycmVudCA9IHJhbmdlO1xuICAgICAgICB9XG4gICAgfSk7XG4gICAgbWVyZ2VkLnB1c2goY3VycmVudCk7XG4gICAgdGhpcy5fcmFuZ2VzID0gbWVyZ2VkO1xufVxuXG5fQ29tcG91bmQucHJvdG90eXBlLm1pbiA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9yYW5nZXNbMF0ubWluKCk7XG59XG5cbl9Db21wb3VuZC5wcm90b3R5cGUubWF4ID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMuX3Jhbmdlc1t0aGlzLl9yYW5nZXMubGVuZ3RoIC0gMV0ubWF4KCk7XG59XG5cbi8vIHJldHVybnMgdGhlIGluZGV4IG9mIHRoZSBmaXJzdCByYW5nZSB0aGF0IGlzIG5vdCBsZXNzIHRoYW4gcG9zXG5fQ29tcG91bmQucHJvdG90eXBlLmxvd2VyX2JvdW5kID0gZnVuY3Rpb24ocG9zKSB7XG4gICAgLy8gZmlyc3QgY2hlY2sgaWYgcG9zIGlzIG91dCBvZiByYW5nZVxuICAgIHZhciByID0gdGhpcy5yYW5nZXMoKTtcbiAgICBpZiAocG9zID4gdGhpcy5tYXgoKSkgcmV0dXJuIHIubGVuZ3RoO1xuICAgIGlmIChwb3MgPCB0aGlzLm1pbigpKSByZXR1cm4gMDtcbiAgICAvLyBkbyBhIGJpbmFyeSBzZWFyY2hcbiAgICB2YXIgYT0wLCBiPXIubGVuZ3RoIC0gMTtcbiAgICB3aGlsZSAoYSA8PSBiKSB7XG4gICAgICAgIHZhciBtID0gTWF0aC5mbG9vcigoYStiKS8yKTtcbiAgICAgICAgaWYgKHBvcyA+IHJbbV0uX21heCkge1xuICAgICAgICAgICAgYSA9IG0rMTtcbiAgICAgICAgfVxuICAgICAgICBlbHNlIGlmIChwb3MgPCByW21dLl9taW4pIHtcbiAgICAgICAgICAgIGIgPSBtLTE7XG4gICAgICAgIH1cbiAgICAgICAgZWxzZSB7XG4gICAgICAgICAgICByZXR1cm4gbTtcbiAgICAgICAgfVxuICAgIH1cbiAgICByZXR1cm4gYTtcbn1cblxuX0NvbXBvdW5kLnByb3RvdHlwZS5jb250YWlucyA9IGZ1bmN0aW9uKHBvcykge1xuICAgIHZhciBsYiA9IHRoaXMubG93ZXJfYm91bmQocG9zKTtcbiAgICBpZiAobGIgPCB0aGlzLl9yYW5nZXMubGVuZ3RoICYmIHRoaXMuX3Jhbmdlc1tsYl0uY29udGFpbnMocG9zKSkge1xuICAgICAgICByZXR1cm4gdHJ1ZTtcbiAgICB9XG4gICAgcmV0dXJuIGZhbHNlO1xufVxuXG5fQ29tcG91bmQucHJvdG90eXBlLmluc2VydFJhbmdlID0gZnVuY3Rpb24ocmFuZ2UpIHtcbiAgICB2YXIgbGIgPSB0aGlzLmxvd2VyX2JvdW5kKHJhbmdlLl9taW4pO1xuICAgIGlmIChsYiA9PT0gdGhpcy5fcmFuZ2VzLmxlbmd0aCkgeyAvLyByYW5nZSBmb2xsb3dzIHRoaXNcbiAgICAgICAgdGhpcy5fcmFuZ2VzLnB1c2gocmFuZ2UpO1xuICAgICAgICByZXR1cm47XG4gICAgfVxuICAgIFxuICAgIHZhciByID0gdGhpcy5yYW5nZXMoKTtcbiAgICBpZiAocmFuZ2UuX21heCA8IHJbbGJdLl9taW4pIHsgLy8gcmFuZ2UgcHJlY2VlZHMgbGJcbiAgICAgICAgdGhpcy5fcmFuZ2VzLnNwbGljZShsYiwwLHJhbmdlKTtcbiAgICAgICAgcmV0dXJuO1xuICAgIH1cblxuICAgIC8vIHJhbmdlIG92ZXJsYXBzIGxiIChhdCBsZWFzdClcbiAgICBpZiAocltsYl0uX21pbiA8IHJhbmdlLl9taW4pIHJhbmdlLl9taW4gPSByW2xiXS5fbWluO1xuICAgIHZhciB1YiA9IGxiKzE7XG4gICAgd2hpbGUgKHViIDwgci5sZW5ndGggJiYgclt1Yl0uX21pbiA8PSByYW5nZS5fbWF4KSB7XG4gICAgICAgIHViKys7XG4gICAgfVxuICAgIHViLS07XG4gICAgLy8gdWIgaXMgdGhlIHVwcGVyIGJvdW5kIG9mIHRoZSBuZXcgcmFuZ2VcbiAgICBpZiAoclt1Yl0uX21heCA+IHJhbmdlLl9tYXgpIHJhbmdlLl9tYXggPSByW3ViXS5fbWF4O1xuICAgIFxuICAgIC8vIHNwbGljZSByYW5nZSBpbnRvIHRoaXMuX3Jhbmdlc1xuICAgIHRoaXMuX3Jhbmdlcy5zcGxpY2UobGIsdWItbGIrMSxyYW5nZSk7XG4gICAgcmV0dXJuO1xufVxuXG5fQ29tcG91bmQucHJvdG90eXBlLmlzQ29udGlndW91cyA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiB0aGlzLl9yYW5nZXMubGVuZ3RoID4gMTtcbn1cblxuX0NvbXBvdW5kLnByb3RvdHlwZS5yYW5nZXMgPSBmdW5jdGlvbigpIHtcbiAgICByZXR1cm4gdGhpcy5fcmFuZ2VzO1xufVxuXG5fQ29tcG91bmQucHJvdG90eXBlLl9wdXNoUmFuZ2VzID0gZnVuY3Rpb24ocmFuZ2VzKSB7XG4gICAgZm9yICh2YXIgcmkgPSAwOyByaSA8IHRoaXMuX3Jhbmdlcy5sZW5ndGg7ICsrcmkpXG4gICAgICAgIHJhbmdlcy5wdXNoKHRoaXMuX3Jhbmdlc1tyaV0pO1xufVxuXG5fQ29tcG91bmQucHJvdG90eXBlLnRvU3RyaW5nID0gZnVuY3Rpb24oKSB7XG4gICAgdmFyIHMgPSAnJztcbiAgICBmb3IgKHZhciByID0gMDsgciA8IHRoaXMuX3Jhbmdlcy5sZW5ndGg7ICsrcikge1xuICAgICAgICBpZiAocj4wKSB7XG4gICAgICAgICAgICBzID0gcyArICcsJztcbiAgICAgICAgfVxuICAgICAgICBzID0gcyArIHRoaXMuX3Jhbmdlc1tyXS50b1N0cmluZygpO1xuICAgIH1cbiAgICByZXR1cm4gcztcbn1cblxuZnVuY3Rpb24gdW5pb24oczAsIHMxKSB7XG4gICAgaWYgKCEgKHMwIGluc3RhbmNlb2YgX0NvbXBvdW5kKSkge1xuICAgICAgICBpZiAoISAoczAgaW5zdGFuY2VvZiBBcnJheSkpXG4gICAgICAgICAgICBzMCA9IFtzMF07XG4gICAgICAgIHMwID0gbmV3IF9Db21wb3VuZChzMCk7XG4gICAgfVxuICAgIFxuICAgIGlmIChzMSlcbiAgICAgICAgczAuaW5zZXJ0UmFuZ2UoczEpO1xuXG4gICAgcmV0dXJuIHMwO1xufVxuXG5mdW5jdGlvbiBpbnRlcnNlY3Rpb24oczAsIHMxKSB7XG4gICAgdmFyIHIwID0gczAucmFuZ2VzKCk7XG4gICAgdmFyIHIxID0gczEucmFuZ2VzKCk7XG4gICAgdmFyIGwwID0gcjAubGVuZ3RoLCBsMSA9IHIxLmxlbmd0aDtcbiAgICB2YXIgaTAgPSAwLCBpMSA9IDA7XG4gICAgdmFyIG9yID0gW107XG5cbiAgICB3aGlsZSAoaTAgPCBsMCAmJiBpMSA8IGwxKSB7XG4gICAgICAgIHZhciBzMCA9IHIwW2kwXSwgczEgPSByMVtpMV07XG4gICAgICAgIHZhciBsYXBNaW4gPSBNYXRoLm1heChzMC5taW4oKSwgczEubWluKCkpO1xuICAgICAgICB2YXIgbGFwTWF4ID0gTWF0aC5taW4oczAubWF4KCksIHMxLm1heCgpKTtcbiAgICAgICAgaWYgKGxhcE1heCA+PSBsYXBNaW4pIHtcbiAgICAgICAgICAgIG9yLnB1c2gobmV3IFJhbmdlKGxhcE1pbiwgbGFwTWF4KSk7XG4gICAgICAgIH1cbiAgICAgICAgaWYgKHMwLm1heCgpID4gczEubWF4KCkpIHtcbiAgICAgICAgICAgICsraTE7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICArK2kwO1xuICAgICAgICB9XG4gICAgfVxuICAgIFxuICAgIGlmIChvci5sZW5ndGggPT0gMCkge1xuICAgICAgICByZXR1cm4gbnVsbDsgLy8gRklYTUVcbiAgICB9IGVsc2UgaWYgKG9yLmxlbmd0aCA9PSAxKSB7XG4gICAgICAgIHJldHVybiBvclswXTtcbiAgICB9IGVsc2Uge1xuICAgICAgICByZXR1cm4gbmV3IF9Db21wb3VuZChvcik7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBjb3ZlcmFnZShzKSB7XG4gICAgdmFyIHRvdCA9IDA7XG4gICAgdmFyIHJsID0gcy5yYW5nZXMoKTtcbiAgICBmb3IgKHZhciByaSA9IDA7IHJpIDwgcmwubGVuZ3RoOyArK3JpKSB7XG4gICAgICAgIHZhciByID0gcmxbcmldO1xuICAgICAgICB0b3QgKz0gKHIubWF4KCkgLSByLm1pbigpICsgMSk7XG4gICAgfVxuICAgIHJldHVybiB0b3Q7XG59XG5cblxuXG5mdW5jdGlvbiByYW5nZU9yZGVyKGEsIGIpXG57XG4gICAgaWYgKGEubWluKCkgPCBiLm1pbigpKSB7XG4gICAgICAgIHJldHVybiAtMTtcbiAgICB9IGVsc2UgaWYgKGEubWluKCkgPiBiLm1pbigpKSB7XG4gICAgICAgIHJldHVybiAxO1xuICAgIH0gZWxzZSBpZiAoYS5tYXgoKSA8IGIubWF4KCkpIHtcbiAgICAgICAgcmV0dXJuIC0xO1xuICAgIH0gZWxzZSBpZiAoYi5tYXgoKSA+IGEubWF4KCkpIHtcbiAgICAgICAgcmV0dXJuIDE7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuIDA7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBfcmFuZ2VPcmRlcihhLCBiKVxue1xuICAgIGlmIChhLl9taW4gPCBiLl9taW4pIHtcbiAgICAgICAgcmV0dXJuIC0xO1xuICAgIH0gZWxzZSBpZiAoYS5fbWluID4gYi5fbWluKSB7XG4gICAgICAgIHJldHVybiAxO1xuICAgIH0gZWxzZSBpZiAoYS5fbWF4IDwgYi5fbWF4KSB7XG4gICAgICAgIHJldHVybiAtMTtcbiAgICB9IGVsc2UgaWYgKGIuX21heCA+IGEuX21heCkge1xuICAgICAgICByZXR1cm4gMTtcbiAgICB9IGVsc2Uge1xuICAgICAgICByZXR1cm4gMDtcbiAgICB9XG59XG5cbmlmICh0eXBlb2YobW9kdWxlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgICBtb2R1bGUuZXhwb3J0cyA9IHtcbiAgICAgICAgUmFuZ2U6IFJhbmdlLFxuICAgICAgICB1bmlvbjogdW5pb24sXG4gICAgICAgIGludGVyc2VjdGlvbjogaW50ZXJzZWN0aW9uLFxuICAgICAgICBjb3ZlcmFnZTogY292ZXJhZ2UsXG4gICAgICAgIHJhbmdlT3ZlcjogcmFuZ2VPcmRlcixcbiAgICAgICAgX3JhbmdlT3JkZXI6IF9yYW5nZU9yZGVyXG4gICAgfVxufSIsIi8qIC0qLSBtb2RlOiBqYXZhc2NyaXB0OyBjLWJhc2ljLW9mZnNldDogNDsgaW5kZW50LXRhYnMtbW9kZTogbmlsIC0qLSAqL1xuXG4vLyBcbi8vIERhbGxpYW5jZSBHZW5vbWUgRXhwbG9yZXJcbi8vIChjKSBUaG9tYXMgRG93biAyMDA2LTIwMTBcbi8vXG4vLyB1dGlscy5qczogb2Rkcywgc29kcywgYW5kIGVuZHMuXG4vL1xuXG5cInVzZSBzdHJpY3RcIjtcblxuaWYgKHR5cGVvZihyZXF1aXJlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgICB2YXIgc2hhMSA9IHJlcXVpcmUoJy4vc2hhMScpO1xuICAgIHZhciBiNjRfc2hhMSA9IHNoYTEuYjY0X3NoYTE7XG59XG5cbnZhciBOVU1fUkVHRVhQID0gbmV3IFJlZ0V4cCgnWzAtOV0rJyk7XG5cbmZ1bmN0aW9uIHN0cmluZ1RvTnVtYmVyc0FycmF5KHN0cikge1xuICAgIHZhciBudW1zID0gbmV3IEFycmF5KCk7XG4gICAgdmFyIG07XG4gICAgd2hpbGUgKG0gPSBOVU1fUkVHRVhQLmV4ZWMoc3RyKSkge1xuICAgICAgICBudW1zLnB1c2gobVswXSk7XG4gICAgICAgIHN0cj1zdHIuc3Vic3RyaW5nKG0uaW5kZXggKyAobVswXS5sZW5ndGgpKTtcbiAgICB9XG4gICAgcmV0dXJuIG51bXM7XG59XG5cbnZhciBTVFJJQ1RfTlVNX1JFR0VYUCA9IG5ldyBSZWdFeHAoJ15bMC05XSskJyk7XG5cbmZ1bmN0aW9uIHN0cmluZ1RvSW50KHN0cikge1xuICAgIHN0ciA9IHN0ci5yZXBsYWNlKG5ldyBSZWdFeHAoJywnLCAnZycpLCAnJyk7XG4gICAgaWYgKCFTVFJJQ1RfTlVNX1JFR0VYUC50ZXN0KHN0cikpIHtcbiAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgfVxuICAgIHJldHVybiBzdHJ8MDtcbn1cblxuZnVuY3Rpb24gcHVzaG5ldyhhLCB2KSB7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBhLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIGlmIChhW2ldID09IHYpIHtcbiAgICAgICAgICAgIHJldHVybjtcbiAgICAgICAgfVxuICAgIH1cbiAgICBhLnB1c2godik7XG59XG5cbmZ1bmN0aW9uIHB1c2hvKG9iaiwgaywgdikge1xuICAgIGlmIChvYmpba10pIHtcbiAgICAgICAgb2JqW2tdLnB1c2godik7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgb2JqW2tdID0gW3ZdO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gcHVzaG5ld28ob2JqLCBrLCB2KSB7XG4gICAgdmFyIGEgPSBvYmpba107XG4gICAgaWYgKGEpIHtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBhLmxlbmd0aDsgKytpKSB7ICAgIC8vIGluZGV4T2YgcmVxdWlyZXMgSlMxNiA6LSguXG4gICAgICAgICAgICBpZiAoYVtpXSA9PSB2KSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG4gICAgICAgIGEucHVzaCh2KTtcbiAgICB9IGVsc2Uge1xuICAgICAgICBvYmpba10gPSBbdl07XG4gICAgfVxufVxuXG5cbmZ1bmN0aW9uIHBpY2soYSwgYiwgYywgZClcbntcbiAgICBpZiAoYSkge1xuICAgICAgICByZXR1cm4gYTtcbiAgICB9IGVsc2UgaWYgKGIpIHtcbiAgICAgICAgcmV0dXJuIGI7XG4gICAgfSBlbHNlIGlmIChjKSB7XG4gICAgICAgIHJldHVybiBjO1xuICAgIH0gZWxzZSBpZiAoZCkge1xuICAgICAgICByZXR1cm4gZDtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIHB1c2huZXcobCwgbylcbntcbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGwubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgaWYgKGxbaV0gPT0gbykge1xuICAgICAgICAgICAgcmV0dXJuO1xuICAgICAgICB9XG4gICAgfVxuICAgIGwucHVzaChvKTtcbn1cblxuXG5cbmZ1bmN0aW9uIGFycmF5SW5kZXhPZihhLCB4KSB7XG4gICAgaWYgKCFhKSB7XG4gICAgICAgIHJldHVybiAtMTtcbiAgICB9XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGEubGVuZ3RoOyArK2kpIHtcbiAgICAgICAgaWYgKGFbaV0gPT09IHgpIHtcbiAgICAgICAgICAgIHJldHVybiBpO1xuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiAtMTtcbn1cblxuZnVuY3Rpb24gYXJyYXlSZW1vdmUoYSwgeCkge1xuICAgIHZhciBpID0gYXJyYXlJbmRleE9mKGEsIHgpO1xuICAgIGlmIChpID49IDApIHtcbiAgICAgICAgYS5zcGxpY2UoaSwgMSk7XG4gICAgICAgIHJldHVybiB0cnVlO1xuICAgIH1cbiAgICByZXR1cm4gZmFsc2U7XG59XG5cbi8vXG4vLyBET00gdXRpbGl0aWVzXG4vL1xuXG5cbmZ1bmN0aW9uIG1ha2VFbGVtZW50KHRhZywgY2hpbGRyZW4sIGF0dHJpYnMsIHN0eWxlcylcbntcbiAgICB2YXIgZWxlID0gZG9jdW1lbnQuY3JlYXRlRWxlbWVudCh0YWcpO1xuICAgIGlmIChjaGlsZHJlbikge1xuICAgICAgICBpZiAoISAoY2hpbGRyZW4gaW5zdGFuY2VvZiBBcnJheSkpIHtcbiAgICAgICAgICAgIGNoaWxkcmVuID0gW2NoaWxkcmVuXTtcbiAgICAgICAgfVxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGNoaWxkcmVuLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICB2YXIgYyA9IGNoaWxkcmVuW2ldO1xuICAgICAgICAgICAgaWYgKGMpIHtcbiAgICAgICAgICAgICAgICBpZiAodHlwZW9mIGMgPT0gJ3N0cmluZycpIHtcbiAgICAgICAgICAgICAgICAgICAgYyA9IGRvY3VtZW50LmNyZWF0ZVRleHROb2RlKGMpO1xuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAodHlwZW9mIGMgPT0gJ251bWJlcicpIHtcbiAgICAgICAgICAgICAgICAgICAgYyA9IGRvY3VtZW50LmNyZWF0ZVRleHROb2RlKCcnICsgYyk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIGVsZS5hcHBlbmRDaGlsZChjKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgIH1cbiAgICBcbiAgICBpZiAoYXR0cmlicykge1xuICAgICAgICBmb3IgKHZhciBsIGluIGF0dHJpYnMpIHtcbiAgICAgICAgICAgIHRyeSB7XG4gICAgICAgICAgICAgICAgZWxlW2xdID0gYXR0cmlic1tsXTtcbiAgICAgICAgICAgIH0gY2F0Y2ggKGUpIHtcbiAgICAgICAgICAgICAgICBjb25zb2xlLmxvZygnZXJyb3Igc2V0dGluZyAnICsgbCk7XG4gICAgICAgICAgICAgICAgdGhyb3coZSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG4gICAgaWYgKHN0eWxlcykge1xuICAgICAgICBmb3IgKHZhciBsIGluIHN0eWxlcykge1xuICAgICAgICAgICAgZWxlLnN0eWxlW2xdID0gc3R5bGVzW2xdO1xuICAgICAgICB9XG4gICAgfVxuICAgIHJldHVybiBlbGU7XG59XG5cbmZ1bmN0aW9uIG1ha2VFbGVtZW50TlMobmFtZXNwYWNlLCB0YWcsIGNoaWxkcmVuLCBhdHRyaWJzKVxue1xuICAgIHZhciBlbGUgPSBkb2N1bWVudC5jcmVhdGVFbGVtZW50TlMobmFtZXNwYWNlLCB0YWcpO1xuICAgIGlmIChjaGlsZHJlbikge1xuICAgICAgICBpZiAoISAoY2hpbGRyZW4gaW5zdGFuY2VvZiBBcnJheSkpIHtcbiAgICAgICAgICAgIGNoaWxkcmVuID0gW2NoaWxkcmVuXTtcbiAgICAgICAgfVxuICAgICAgICBmb3IgKHZhciBpID0gMDsgaSA8IGNoaWxkcmVuLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICB2YXIgYyA9IGNoaWxkcmVuW2ldO1xuICAgICAgICAgICAgaWYgKHR5cGVvZiBjID09ICdzdHJpbmcnKSB7XG4gICAgICAgICAgICAgICAgYyA9IGRvY3VtZW50LmNyZWF0ZVRleHROb2RlKGMpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgZWxlLmFwcGVuZENoaWxkKGMpO1xuICAgICAgICB9XG4gICAgfVxuICAgIFxuICAgIHNldEF0dHJzKGVsZSwgYXR0cmlicyk7XG4gICAgcmV0dXJuIGVsZTtcbn1cblxudmFyIGF0dHJfbmFtZV9jYWNoZSA9IHt9O1xuXG5mdW5jdGlvbiBzZXRBdHRyKG5vZGUsIGtleSwgdmFsdWUpXG57XG4gICAgdmFyIGF0dHIgPSBhdHRyX25hbWVfY2FjaGVba2V5XTtcbiAgICBpZiAoIWF0dHIpIHtcbiAgICAgICAgdmFyIF9hdHRyID0gJyc7XG4gICAgICAgIGZvciAodmFyIGMgPSAwOyBjIDwga2V5Lmxlbmd0aDsgKytjKSB7XG4gICAgICAgICAgICB2YXIgY2MgPSBrZXkuc3Vic3RyaW5nKGMsIGMrMSk7XG4gICAgICAgICAgICB2YXIgbGNjID0gY2MudG9Mb3dlckNhc2UoKTtcbiAgICAgICAgICAgIGlmIChsY2MgIT0gY2MpIHtcbiAgICAgICAgICAgICAgICBfYXR0ciA9IF9hdHRyICsgJy0nICsgbGNjO1xuICAgICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgICAgICBfYXR0ciA9IF9hdHRyICsgY2M7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICAgICAgYXR0cl9uYW1lX2NhY2hlW2tleV0gPSBfYXR0cjtcbiAgICAgICAgYXR0ciA9IF9hdHRyO1xuICAgIH1cbiAgICBub2RlLnNldEF0dHJpYnV0ZShhdHRyLCB2YWx1ZSk7XG59XG5cbmZ1bmN0aW9uIHNldEF0dHJzKG5vZGUsIGF0dHJpYnMpXG57XG4gICAgaWYgKGF0dHJpYnMpIHtcbiAgICAgICAgZm9yICh2YXIgbCBpbiBhdHRyaWJzKSB7XG4gICAgICAgICAgICBzZXRBdHRyKG5vZGUsIGwsIGF0dHJpYnNbbF0pO1xuICAgICAgICB9XG4gICAgfVxufVxuXG5cblxuZnVuY3Rpb24gcmVtb3ZlQ2hpbGRyZW4obm9kZSlcbntcbiAgICBpZiAoIW5vZGUgfHwgIW5vZGUuY2hpbGROb2Rlcykge1xuICAgICAgICByZXR1cm47XG4gICAgfVxuXG4gICAgd2hpbGUgKG5vZGUuY2hpbGROb2Rlcy5sZW5ndGggPiAwKSB7XG4gICAgICAgIG5vZGUucmVtb3ZlQ2hpbGQobm9kZS5maXJzdENoaWxkKTtcbiAgICB9XG59XG5cblxuXG4vL1xuLy8gV0FSTklORzogbm90IGZvciBnZW5lcmFsIHVzZSFcbi8vXG5cbmZ1bmN0aW9uIG1pbmlKU09OaWZ5KG8sIGV4Yykge1xuICAgIGlmICh0eXBlb2YgbyA9PT0gJ3VuZGVmaW5lZCcpIHtcbiAgICAgICAgcmV0dXJuICd1bmRlZmluZWQnO1xuICAgIH0gZWxzZSBpZiAobyA9PSBudWxsKSB7XG4gICAgICAgIHJldHVybiAnbnVsbCc7XG4gICAgfSBlbHNlIGlmICh0eXBlb2YgbyA9PSAnc3RyaW5nJykge1xuICAgICAgICByZXR1cm4gXCInXCIgKyBvICsgXCInXCI7XG4gICAgfSBlbHNlIGlmICh0eXBlb2YgbyA9PSAnbnVtYmVyJykge1xuICAgICAgICByZXR1cm4gXCJcIiArIG87XG4gICAgfSBlbHNlIGlmICh0eXBlb2YgbyA9PSAnYm9vbGVhbicpIHtcbiAgICAgICAgcmV0dXJuIFwiXCIgKyBvO1xuICAgIH0gZWxzZSBpZiAodHlwZW9mIG8gPT0gJ29iamVjdCcpIHtcbiAgICAgICAgaWYgKG8gaW5zdGFuY2VvZiBBcnJheSkge1xuICAgICAgICAgICAgdmFyIHMgPSBudWxsO1xuICAgICAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBvLmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICAgICAgcyA9IChzID09IG51bGwgPyAnJyA6IChzICsgJywgJykpICsgbWluaUpTT05pZnkob1tpXSwgZXhjKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgICAgIHJldHVybiAnWycgKyAocz9zOicnKSArICddJztcbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGV4YyA9IGV4YyB8fCB7fTtcbiAgICAgICAgICAgIHZhciBzID0gbnVsbDtcbiAgICAgICAgICAgIGZvciAodmFyIGsgaW4gbykge1xuICAgICAgICAgICAgICAgIGlmIChleGNba10pXG4gICAgICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgICAgIGlmIChrICE9IHVuZGVmaW5lZCAmJiB0eXBlb2Yob1trXSkgIT0gJ2Z1bmN0aW9uJykge1xuICAgICAgICAgICAgICAgICAgICBzID0gKHMgPT0gbnVsbCA/ICcnIDogKHMgKyAnLCAnKSkgKyBrICsgJzogJyArIG1pbmlKU09OaWZ5KG9ba10sIGV4Yyk7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgcmV0dXJuICd7JyArIChzP3M6JycpICsgJ30nO1xuICAgICAgICB9XG4gICAgfSBlbHNlIHtcbiAgICAgICAgcmV0dXJuICh0eXBlb2Ygbyk7XG4gICAgfVxufVxuXG5mdW5jdGlvbiBzaGFsbG93Q29weShvKSB7XG4gICAgdmFyIG4gPSB7fTtcbiAgICBmb3IgKHZhciBrIGluIG8pIHtcbiAgICAgICAgbltrXSA9IG9ba107XG4gICAgfVxuICAgIHJldHVybiBuO1xufVxuXG5mdW5jdGlvbiBPYnNlcnZlZCh4KSB7XG4gICAgdGhpcy52YWx1ZSA9IHg7XG4gICAgdGhpcy5saXN0ZW5lcnMgPSBbXTtcbn1cblxuT2JzZXJ2ZWQucHJvdG90eXBlLmFkZExpc3RlbmVyID0gZnVuY3Rpb24oZikge1xuICAgIHRoaXMubGlzdGVuZXJzLnB1c2goZik7XG59XG5cbk9ic2VydmVkLnByb3RvdHlwZS5hZGRMaXN0ZW5lckFuZEZpcmUgPSBmdW5jdGlvbihmKSB7XG4gICAgdGhpcy5saXN0ZW5lcnMucHVzaChmKTtcbiAgICBmKHRoaXMudmFsdWUpO1xufVxuXG5PYnNlcnZlZC5wcm90b3R5cGUucmVtb3ZlTGlzdGVuZXIgPSBmdW5jdGlvbihmKSB7XG4gICAgYXJyYXlSZW1vdmUodGhpcy5saXN0ZW5lcnMsIGYpO1xufVxuXG5PYnNlcnZlZC5wcm90b3R5cGUuZ2V0ID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIHRoaXMudmFsdWU7XG59XG5cbk9ic2VydmVkLnByb3RvdHlwZS5zZXQgPSBmdW5jdGlvbih4KSB7XG4gICAgdGhpcy52YWx1ZSA9IHg7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLmxpc3RlbmVycy5sZW5ndGg7ICsraSkge1xuICAgICAgICB0aGlzLmxpc3RlbmVyc1tpXSh4KTtcbiAgICB9XG59XG5cbmZ1bmN0aW9uIEF3YWl0ZWQoKSB7XG4gICAgdGhpcy5xdWV1ZSA9IFtdO1xufVxuXG5Bd2FpdGVkLnByb3RvdHlwZS5wcm92aWRlID0gZnVuY3Rpb24oeCkge1xuICAgIGlmICh0aGlzLnJlcyAhPT0gdW5kZWZpbmVkKSB7XG4gICAgICAgIHRocm93IFwiUmVzb3VyY2UgaGFzIGFscmVhZHkgYmVlbiBwcm92aWRlZC5cIjtcbiAgICB9XG5cbiAgICB0aGlzLnJlcyA9IHg7XG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCB0aGlzLnF1ZXVlLmxlbmd0aDsgKytpKSB7XG4gICAgICAgIHRoaXMucXVldWVbaV0oeCk7XG4gICAgfVxuICAgIHRoaXMucXVldWUgPSBudWxsOyAgIC8vIGF2b2lkIGxlYWtpbmcgY2xvc3VyZXMuXG59XG5cbkF3YWl0ZWQucHJvdG90eXBlLmF3YWl0ID0gZnVuY3Rpb24oZikge1xuICAgIGlmICh0aGlzLnJlcyAhPT0gdW5kZWZpbmVkKSB7XG4gICAgICAgIGYodGhpcy5yZXMpO1xuICAgICAgICByZXR1cm4gdGhpcy5yZXM7XG4gICAgfSBlbHNlIHtcbiAgICAgICAgdGhpcy5xdWV1ZS5wdXNoKGYpO1xuICAgIH1cbn1cblxudmFyIF9fZGFsbGlhbmNlX3NhbHRTZWVkID0gMDtcblxuZnVuY3Rpb24gc2FsdFVSTCh1cmwpIHtcbiAgICByZXR1cm4gdXJsICsgJz9zYWx0PScgKyBiNjRfc2hhMSgnJyArIERhdGUubm93KCkgKyAnLCcgKyAoKytfX2RhbGxpYW5jZV9zYWx0U2VlZCkpO1xufVxuXG5mdW5jdGlvbiB0ZXh0WEhSKHVybCwgY2FsbGJhY2ssIG9wdHMpIHtcbiAgICBpZiAob3B0cyAmJiBvcHRzLnNhbHQpIFxuICAgICAgICB1cmwgPSBzYWx0VVJMKHVybCk7XG5cbiAgICB0cnkge1xuICAgICAgICB2YXIgdGltZW91dDtcbiAgICAgICAgaWYgKG9wdHMudGltZW91dCkge1xuICAgICAgICAgICAgdGltZW91dCA9IHNldFRpbWVvdXQoXG4gICAgICAgICAgICAgICAgZnVuY3Rpb24oKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUubG9nKCd0aW1pbmcgb3V0ICcgKyB1cmwpO1xuICAgICAgICAgICAgICAgICAgICByZXEuYWJvcnQoKTtcbiAgICAgICAgICAgICAgICAgICAgcmV0dXJuIGNhbGxiYWNrKG51bGwsICdUaW1lb3V0Jyk7XG4gICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICBvcHRzLnRpbWVvdXRcbiAgICAgICAgICAgICk7XG4gICAgICAgIH1cblxuICAgICAgICB2YXIgcmVxID0gbmV3IFhNTEh0dHBSZXF1ZXN0KCk7XG4gICAgICAgIHJlcS5vbnJlYWR5c3RhdGVjaGFuZ2UgPSBmdW5jdGlvbigpIHtcbiAgICBcdCAgICBpZiAocmVxLnJlYWR5U3RhdGUgPT0gNCkge1xuICAgICAgICAgICAgICAgIGlmICh0aW1lb3V0KVxuICAgICAgICAgICAgICAgICAgICBjbGVhclRpbWVvdXQodGltZW91dCk7XG4gICAgXHQgICAgICAgIGlmIChyZXEuc3RhdHVzIDwgMjAwIHx8IHJlcS5zdGF0dXMgPj0gMzAwKSB7XG4gICAgXHRcdCAgICBjYWxsYmFjayhudWxsLCAnRXJyb3IgY29kZSAnICsgcmVxLnN0YXR1cyk7XG4gICAgXHQgICAgICAgIH0gZWxzZSB7XG4gICAgXHRcdCAgICBjYWxsYmFjayhyZXEucmVzcG9uc2VUZXh0KTtcbiAgICBcdCAgICAgICAgfVxuICAgIFx0ICAgIH1cbiAgICAgICAgfTtcbiAgICAgICAgXG4gICAgICAgIHJlcS5vcGVuKCdHRVQnLCB1cmwsIHRydWUpO1xuICAgICAgICByZXEucmVzcG9uc2VUeXBlID0gJ3RleHQnO1xuXG4gICAgICAgIGlmIChvcHRzICYmIG9wdHMuY3JlZGVudGlhbHMpIHtcbiAgICAgICAgICAgIHJlcS53aXRoQ3JlZGVudGlhbHMgPSB0cnVlO1xuICAgICAgICB9XG4gICAgICAgIHJlcS5zZW5kKCcnKTtcbiAgICB9IGNhdGNoIChlKSB7XG4gICAgICAgIGNhbGxiYWNrKG51bGwsICdFeGNlcHRpb24gJyArIGUpO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gcmVsYXRpdmVVUkwoYmFzZSwgcmVsKSB7XG4gICAgLy8gRklYTUUgcXVpdGUgbmFpdmUgLS0gZ29vZCBlbm91Z2ggZm9yIHRyYWNraHVicz9cblxuICAgIGlmIChyZWwuaW5kZXhPZignaHR0cDonKSA9PSAwIHx8IHJlbC5pbmRleE9mKCdodHRwczonKSA9PSAwKSB7XG4gICAgICAgIHJldHVybiByZWw7XG4gICAgfVxuXG4gICAgdmFyIGxpID0gYmFzZS5sYXN0SW5kZXhPZignLycpO1xuICAgIGlmIChsaSA+PSAwKSB7XG4gICAgICAgIHJldHVybiBiYXNlLnN1YnN0cigwLCBsaSArIDEpICsgcmVsO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHJldHVybiByZWw7XG4gICAgfVxufVxuXG52YXIgQU1JTk9fQUNJRF9UUkFOU0xBVElPTiA9IHtcbiAgICAnVFRUJzogJ0YnLFxuICAgICdUVEMnOiAnRicsXG4gICAgJ1RUQSc6ICdMJyxcbiAgICAnVFRHJzogJ0wnLFxuICAgICdDVFQnOiAnTCcsXG4gICAgJ0NUQyc6ICdMJyxcbiAgICAnQ1RBJzogJ0wnLFxuICAgICdDVEcnOiAnTCcsXG4gICAgJ0FUVCc6ICdJJyxcbiAgICAnQVRDJzogJ0knLFxuICAgICdBVEEnOiAnSScsXG4gICAgJ0FURyc6ICdNJyxcbiAgICAnR1RUJzogJ1YnLFxuICAgICdHVEMnOiAnVicsXG4gICAgJ0dUQSc6ICdWJyxcbiAgICAnR1RHJzogJ1YnLFxuICAgICdUQ1QnOiAnUycsXG4gICAgJ1RDQyc6ICdTJyxcbiAgICAnVENBJzogJ1MnLFxuICAgICdUQ0cnOiAnUycsXG4gICAgJ0NDVCc6ICdQJyxcbiAgICAnQ0NDJzogJ1AnLFxuICAgICdDQ0EnOiAnUCcsXG4gICAgJ0NDRyc6ICdQJyxcbiAgICAnQUNUJzogJ1QnLFxuICAgICdBQ0MnOiAnVCcsXG4gICAgJ0FDQSc6ICdUJyxcbiAgICAnQUNHJzogJ1QnLFxuICAgICdHQ1QnOiAnQScsXG4gICAgJ0dDQyc6ICdBJyxcbiAgICAnR0NBJzogJ0EnLFxuICAgICdHQ0cnOiAnQScsXG4gICAgJ1RBVCc6ICdZJyxcbiAgICAnVEFDJzogJ1knLFxuICAgICdUQUEnOiAnKicsICAvLyBzdG9wXG4gICAgJ1RBRyc6ICcqJywgIC8vIHN0b3BcbiAgICAnQ0FUJzogJ0gnLFxuICAgICdDQUMnOiAnSCcsXG4gICAgJ0NBQSc6ICdRJyxcbiAgICAnQ0FHJzogJ1EnLFxuICAgICdBQVQnOiAnTicsXG4gICAgJ0FBQyc6ICdOJyxcbiAgICAnQUFBJzogJ0snLFxuICAgICdBQUcnOiAnSycsXG4gICAgJ0dBVCc6ICdEJyxcbiAgICAnR0FDJzogJ0QnLFxuICAgICdHQUEnOiAnRScsXG4gICAgJ0dBRyc6ICdFJyxcbiAgICAnVEdUJzogJ0MnLFxuICAgICdUR0MnOiAnQycsXG4gICAgJ1RHQSc6ICcqJywgIC8vIHN0b3BcbiAgICAnVEdHJzogJ1cnLFxuICAgICdDR1QnOiAnUicsXG4gICAgJ0NHQyc6ICdSJyxcbiAgICAnQ0dBJzogJ1InLFxuICAgICdDR0cnOiAnUicsXG4gICAgJ0FHVCc6ICdTJyxcbiAgICAnQUdDJzogJ1MnLFxuICAgICdBR0EnOiAnUicsXG4gICAgJ0FHRyc6ICdSJyxcbiAgICAnR0dUJzogJ0cnLFxuICAgICdHR0MnOiAnRycsXG4gICAgJ0dHQSc6ICdHJyxcbiAgICAnR0dHJzogJ0cnXG59XG5cbmZ1bmN0aW9uIHJlc29sdmVVcmxUb1BhZ2UocmVsKSB7XG4gICAgcmV0dXJuIG1ha2VFbGVtZW50KCdhJywgbnVsbCwge2hyZWY6IHJlbH0pLmhyZWY7XG59XG5cbi8vXG4vLyBNaXNzaW5nIEFQSXNcbi8vIFxuXG5pZiAoISgndHJpbScgaW4gU3RyaW5nLnByb3RvdHlwZSkpIHtcbiAgICBTdHJpbmcucHJvdG90eXBlLnRyaW0gPSBmdW5jdGlvbigpIHtcbiAgICAgICAgcmV0dXJuIHRoaXMucmVwbGFjZSgvXlxccysvLCAnJykucmVwbGFjZSgvXFxzKyQvLCAnJyk7XG4gICAgfTtcbn1cblxuaWYgKHR5cGVvZihtb2R1bGUpICE9PSAndW5kZWZpbmVkJykge1xuICAgIG1vZHVsZS5leHBvcnRzID0ge1xuICAgICAgICB0ZXh0WEhSOiB0ZXh0WEhSLFxuICAgICAgICByZWxhdGl2ZVVSTDogcmVsYXRpdmVVUkwsXG4gICAgICAgIHJlc29sdmVVcmxUb1BhZ2U6IHJlc29sdmVVcmxUb1BhZ2UsXG4gICAgICAgIHNoYWxsb3dDb3B5OiBzaGFsbG93Q29weSxcbiAgICAgICAgcHVzaG86IHB1c2hvLFxuICAgICAgICBwdXNobmV3OiBwdXNobmV3LFxuICAgICAgICBwdXNobmV3bzogcHVzaG5ld28sXG4gICAgICAgIGFycmF5SW5kZXhPZjogYXJyYXlJbmRleE9mLFxuICAgICAgICBwaWNrOiBwaWNrLFxuXG4gICAgICAgIG1ha2VFbGVtZW50OiBtYWtlRWxlbWVudCxcbiAgICAgICAgbWFrZUVsZW1lbnROUzogbWFrZUVsZW1lbnROUyxcbiAgICAgICAgcmVtb3ZlQ2hpbGRyZW46IHJlbW92ZUNoaWxkcmVuLFxuXG4gICAgICAgIG1pbmlKU09OaWZ5OiBtaW5pSlNPTmlmeSxcblxuICAgICAgICBPYnNlcnZlZDogT2JzZXJ2ZWQsXG4gICAgICAgIEF3YWl0ZWQ6IEF3YWl0ZWQsXG5cbiAgICAgICAgQU1JTk9fQUNJRF9UUkFOU0xBVElPTjogQU1JTk9fQUNJRF9UUkFOU0xBVElPTlxuICAgIH1cbn1cbiIsIlwidXNlIHN0cmljdFwiO1xudmFyIFByb21pc2UgPSByZXF1aXJlKFwiLi9wcm9taXNlL3Byb21pc2VcIikuUHJvbWlzZTtcbnZhciBwb2x5ZmlsbCA9IHJlcXVpcmUoXCIuL3Byb21pc2UvcG9seWZpbGxcIikucG9seWZpbGw7XG5leHBvcnRzLlByb21pc2UgPSBQcm9taXNlO1xuZXhwb3J0cy5wb2x5ZmlsbCA9IHBvbHlmaWxsOyIsIlwidXNlIHN0cmljdFwiO1xuLyogZ2xvYmFsIHRvU3RyaW5nICovXG5cbnZhciBpc0FycmF5ID0gcmVxdWlyZShcIi4vdXRpbHNcIikuaXNBcnJheTtcbnZhciBpc0Z1bmN0aW9uID0gcmVxdWlyZShcIi4vdXRpbHNcIikuaXNGdW5jdGlvbjtcblxuLyoqXG4gIFJldHVybnMgYSBwcm9taXNlIHRoYXQgaXMgZnVsZmlsbGVkIHdoZW4gYWxsIHRoZSBnaXZlbiBwcm9taXNlcyBoYXZlIGJlZW5cbiAgZnVsZmlsbGVkLCBvciByZWplY3RlZCBpZiBhbnkgb2YgdGhlbSBiZWNvbWUgcmVqZWN0ZWQuIFRoZSByZXR1cm4gcHJvbWlzZVxuICBpcyBmdWxmaWxsZWQgd2l0aCBhbiBhcnJheSB0aGF0IGdpdmVzIGFsbCB0aGUgdmFsdWVzIGluIHRoZSBvcmRlciB0aGV5IHdlcmVcbiAgcGFzc2VkIGluIHRoZSBgcHJvbWlzZXNgIGFycmF5IGFyZ3VtZW50LlxuXG4gIEV4YW1wbGU6XG5cbiAgYGBgamF2YXNjcmlwdFxuICB2YXIgcHJvbWlzZTEgPSBSU1ZQLnJlc29sdmUoMSk7XG4gIHZhciBwcm9taXNlMiA9IFJTVlAucmVzb2x2ZSgyKTtcbiAgdmFyIHByb21pc2UzID0gUlNWUC5yZXNvbHZlKDMpO1xuICB2YXIgcHJvbWlzZXMgPSBbIHByb21pc2UxLCBwcm9taXNlMiwgcHJvbWlzZTMgXTtcblxuICBSU1ZQLmFsbChwcm9taXNlcykudGhlbihmdW5jdGlvbihhcnJheSl7XG4gICAgLy8gVGhlIGFycmF5IGhlcmUgd291bGQgYmUgWyAxLCAyLCAzIF07XG4gIH0pO1xuICBgYGBcblxuICBJZiBhbnkgb2YgdGhlIGBwcm9taXNlc2AgZ2l2ZW4gdG8gYFJTVlAuYWxsYCBhcmUgcmVqZWN0ZWQsIHRoZSBmaXJzdCBwcm9taXNlXG4gIHRoYXQgaXMgcmVqZWN0ZWQgd2lsbCBiZSBnaXZlbiBhcyBhbiBhcmd1bWVudCB0byB0aGUgcmV0dXJuZWQgcHJvbWlzZXMnc1xuICByZWplY3Rpb24gaGFuZGxlci4gRm9yIGV4YW1wbGU6XG5cbiAgRXhhbXBsZTpcblxuICBgYGBqYXZhc2NyaXB0XG4gIHZhciBwcm9taXNlMSA9IFJTVlAucmVzb2x2ZSgxKTtcbiAgdmFyIHByb21pc2UyID0gUlNWUC5yZWplY3QobmV3IEVycm9yKFwiMlwiKSk7XG4gIHZhciBwcm9taXNlMyA9IFJTVlAucmVqZWN0KG5ldyBFcnJvcihcIjNcIikpO1xuICB2YXIgcHJvbWlzZXMgPSBbIHByb21pc2UxLCBwcm9taXNlMiwgcHJvbWlzZTMgXTtcblxuICBSU1ZQLmFsbChwcm9taXNlcykudGhlbihmdW5jdGlvbihhcnJheSl7XG4gICAgLy8gQ29kZSBoZXJlIG5ldmVyIHJ1bnMgYmVjYXVzZSB0aGVyZSBhcmUgcmVqZWN0ZWQgcHJvbWlzZXMhXG4gIH0sIGZ1bmN0aW9uKGVycm9yKSB7XG4gICAgLy8gZXJyb3IubWVzc2FnZSA9PT0gXCIyXCJcbiAgfSk7XG4gIGBgYFxuXG4gIEBtZXRob2QgYWxsXG4gIEBmb3IgUlNWUFxuICBAcGFyYW0ge0FycmF5fSBwcm9taXNlc1xuICBAcGFyYW0ge1N0cmluZ30gbGFiZWxcbiAgQHJldHVybiB7UHJvbWlzZX0gcHJvbWlzZSB0aGF0IGlzIGZ1bGZpbGxlZCB3aGVuIGFsbCBgcHJvbWlzZXNgIGhhdmUgYmVlblxuICBmdWxmaWxsZWQsIG9yIHJlamVjdGVkIGlmIGFueSBvZiB0aGVtIGJlY29tZSByZWplY3RlZC5cbiovXG5mdW5jdGlvbiBhbGwocHJvbWlzZXMpIHtcbiAgLypqc2hpbnQgdmFsaWR0aGlzOnRydWUgKi9cbiAgdmFyIFByb21pc2UgPSB0aGlzO1xuXG4gIGlmICghaXNBcnJheShwcm9taXNlcykpIHtcbiAgICB0aHJvdyBuZXcgVHlwZUVycm9yKCdZb3UgbXVzdCBwYXNzIGFuIGFycmF5IHRvIGFsbC4nKTtcbiAgfVxuXG4gIHJldHVybiBuZXcgUHJvbWlzZShmdW5jdGlvbihyZXNvbHZlLCByZWplY3QpIHtcbiAgICB2YXIgcmVzdWx0cyA9IFtdLCByZW1haW5pbmcgPSBwcm9taXNlcy5sZW5ndGgsXG4gICAgcHJvbWlzZTtcblxuICAgIGlmIChyZW1haW5pbmcgPT09IDApIHtcbiAgICAgIHJlc29sdmUoW10pO1xuICAgIH1cblxuICAgIGZ1bmN0aW9uIHJlc29sdmVyKGluZGV4KSB7XG4gICAgICByZXR1cm4gZnVuY3Rpb24odmFsdWUpIHtcbiAgICAgICAgcmVzb2x2ZUFsbChpbmRleCwgdmFsdWUpO1xuICAgICAgfTtcbiAgICB9XG5cbiAgICBmdW5jdGlvbiByZXNvbHZlQWxsKGluZGV4LCB2YWx1ZSkge1xuICAgICAgcmVzdWx0c1tpbmRleF0gPSB2YWx1ZTtcbiAgICAgIGlmICgtLXJlbWFpbmluZyA9PT0gMCkge1xuICAgICAgICByZXNvbHZlKHJlc3VsdHMpO1xuICAgICAgfVxuICAgIH1cblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgcHJvbWlzZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgIHByb21pc2UgPSBwcm9taXNlc1tpXTtcblxuICAgICAgaWYgKHByb21pc2UgJiYgaXNGdW5jdGlvbihwcm9taXNlLnRoZW4pKSB7XG4gICAgICAgIHByb21pc2UudGhlbihyZXNvbHZlcihpKSwgcmVqZWN0KTtcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIHJlc29sdmVBbGwoaSwgcHJvbWlzZSk7XG4gICAgICB9XG4gICAgfVxuICB9KTtcbn1cblxuZXhwb3J0cy5hbGwgPSBhbGw7IiwiKGZ1bmN0aW9uIChwcm9jZXNzLGdsb2JhbCl7XG5cInVzZSBzdHJpY3RcIjtcbnZhciBicm93c2VyR2xvYmFsID0gKHR5cGVvZiB3aW5kb3cgIT09ICd1bmRlZmluZWQnKSA/IHdpbmRvdyA6IHt9O1xudmFyIEJyb3dzZXJNdXRhdGlvbk9ic2VydmVyID0gYnJvd3Nlckdsb2JhbC5NdXRhdGlvbk9ic2VydmVyIHx8IGJyb3dzZXJHbG9iYWwuV2ViS2l0TXV0YXRpb25PYnNlcnZlcjtcbnZhciBsb2NhbCA9ICh0eXBlb2YgZ2xvYmFsICE9PSAndW5kZWZpbmVkJykgPyBnbG9iYWwgOiAodGhpcyA9PT0gdW5kZWZpbmVkPyB3aW5kb3c6dGhpcyk7XG5cbi8vIG5vZGVcbmZ1bmN0aW9uIHVzZU5leHRUaWNrKCkge1xuICByZXR1cm4gZnVuY3Rpb24oKSB7XG4gICAgcHJvY2Vzcy5uZXh0VGljayhmbHVzaCk7XG4gIH07XG59XG5cbmZ1bmN0aW9uIHVzZU11dGF0aW9uT2JzZXJ2ZXIoKSB7XG4gIHZhciBpdGVyYXRpb25zID0gMDtcbiAgdmFyIG9ic2VydmVyID0gbmV3IEJyb3dzZXJNdXRhdGlvbk9ic2VydmVyKGZsdXNoKTtcbiAgdmFyIG5vZGUgPSBkb2N1bWVudC5jcmVhdGVUZXh0Tm9kZSgnJyk7XG4gIG9ic2VydmVyLm9ic2VydmUobm9kZSwgeyBjaGFyYWN0ZXJEYXRhOiB0cnVlIH0pO1xuXG4gIHJldHVybiBmdW5jdGlvbigpIHtcbiAgICBub2RlLmRhdGEgPSAoaXRlcmF0aW9ucyA9ICsraXRlcmF0aW9ucyAlIDIpO1xuICB9O1xufVxuXG5mdW5jdGlvbiB1c2VTZXRUaW1lb3V0KCkge1xuICByZXR1cm4gZnVuY3Rpb24oKSB7XG4gICAgbG9jYWwuc2V0VGltZW91dChmbHVzaCwgMSk7XG4gIH07XG59XG5cbnZhciBxdWV1ZSA9IFtdO1xuZnVuY3Rpb24gZmx1c2goKSB7XG4gIGZvciAodmFyIGkgPSAwOyBpIDwgcXVldWUubGVuZ3RoOyBpKyspIHtcbiAgICB2YXIgdHVwbGUgPSBxdWV1ZVtpXTtcbiAgICB2YXIgY2FsbGJhY2sgPSB0dXBsZVswXSwgYXJnID0gdHVwbGVbMV07XG4gICAgY2FsbGJhY2soYXJnKTtcbiAgfVxuICBxdWV1ZSA9IFtdO1xufVxuXG52YXIgc2NoZWR1bGVGbHVzaDtcblxuLy8gRGVjaWRlIHdoYXQgYXN5bmMgbWV0aG9kIHRvIHVzZSB0byB0cmlnZ2VyaW5nIHByb2Nlc3Npbmcgb2YgcXVldWVkIGNhbGxiYWNrczpcbmlmICh0eXBlb2YgcHJvY2VzcyAhPT0gJ3VuZGVmaW5lZCcgJiYge30udG9TdHJpbmcuY2FsbChwcm9jZXNzKSA9PT0gJ1tvYmplY3QgcHJvY2Vzc10nKSB7XG4gIHNjaGVkdWxlRmx1c2ggPSB1c2VOZXh0VGljaygpO1xufSBlbHNlIGlmIChCcm93c2VyTXV0YXRpb25PYnNlcnZlcikge1xuICBzY2hlZHVsZUZsdXNoID0gdXNlTXV0YXRpb25PYnNlcnZlcigpO1xufSBlbHNlIHtcbiAgc2NoZWR1bGVGbHVzaCA9IHVzZVNldFRpbWVvdXQoKTtcbn1cblxuZnVuY3Rpb24gYXNhcChjYWxsYmFjaywgYXJnKSB7XG4gIHZhciBsZW5ndGggPSBxdWV1ZS5wdXNoKFtjYWxsYmFjaywgYXJnXSk7XG4gIGlmIChsZW5ndGggPT09IDEpIHtcbiAgICAvLyBJZiBsZW5ndGggaXMgMSwgdGhhdCBtZWFucyB0aGF0IHdlIG5lZWQgdG8gc2NoZWR1bGUgYW4gYXN5bmMgZmx1c2guXG4gICAgLy8gSWYgYWRkaXRpb25hbCBjYWxsYmFja3MgYXJlIHF1ZXVlZCBiZWZvcmUgdGhlIHF1ZXVlIGlzIGZsdXNoZWQsIHRoZXlcbiAgICAvLyB3aWxsIGJlIHByb2Nlc3NlZCBieSB0aGlzIGZsdXNoIHRoYXQgd2UgYXJlIHNjaGVkdWxpbmcuXG4gICAgc2NoZWR1bGVGbHVzaCgpO1xuICB9XG59XG5cbmV4cG9ydHMuYXNhcCA9IGFzYXA7XG59KS5jYWxsKHRoaXMscmVxdWlyZShcIjFZaVo1U1wiKSx0eXBlb2Ygc2VsZiAhPT0gXCJ1bmRlZmluZWRcIiA/IHNlbGYgOiB0eXBlb2Ygd2luZG93ICE9PSBcInVuZGVmaW5lZFwiID8gd2luZG93IDoge30pIiwiXCJ1c2Ugc3RyaWN0XCI7XG4vKipcbiAgYFJTVlAuUHJvbWlzZS5jYXN0YCByZXR1cm5zIHRoZSBzYW1lIHByb21pc2UgaWYgdGhhdCBwcm9taXNlIHNoYXJlcyBhIGNvbnN0cnVjdG9yXG4gIHdpdGggdGhlIHByb21pc2UgYmVpbmcgY2FzdGVkLlxuXG4gIEV4YW1wbGU6XG5cbiAgYGBgamF2YXNjcmlwdFxuICB2YXIgcHJvbWlzZSA9IFJTVlAucmVzb2x2ZSgxKTtcbiAgdmFyIGNhc3RlZCA9IFJTVlAuUHJvbWlzZS5jYXN0KHByb21pc2UpO1xuXG4gIGNvbnNvbGUubG9nKHByb21pc2UgPT09IGNhc3RlZCk7IC8vIHRydWVcbiAgYGBgXG5cbiAgSW4gdGhlIGNhc2Ugb2YgYSBwcm9taXNlIHdob3NlIGNvbnN0cnVjdG9yIGRvZXMgbm90IG1hdGNoLCBpdCBpcyBhc3NpbWlsYXRlZC5cbiAgVGhlIHJlc3VsdGluZyBwcm9taXNlIHdpbGwgZnVsZmlsbCBvciByZWplY3QgYmFzZWQgb24gdGhlIG91dGNvbWUgb2YgdGhlXG4gIHByb21pc2UgYmVpbmcgY2FzdGVkLlxuXG4gIEluIHRoZSBjYXNlIG9mIGEgbm9uLXByb21pc2UsIGEgcHJvbWlzZSB3aGljaCB3aWxsIGZ1bGZpbGwgd2l0aCB0aGF0IHZhbHVlIGlzXG4gIHJldHVybmVkLlxuXG4gIEV4YW1wbGU6XG5cbiAgYGBgamF2YXNjcmlwdFxuICB2YXIgdmFsdWUgPSAxOyAvLyBjb3VsZCBiZSBhIG51bWJlciwgYm9vbGVhbiwgc3RyaW5nLCB1bmRlZmluZWQuLi5cbiAgdmFyIGNhc3RlZCA9IFJTVlAuUHJvbWlzZS5jYXN0KHZhbHVlKTtcblxuICBjb25zb2xlLmxvZyh2YWx1ZSA9PT0gY2FzdGVkKTsgLy8gZmFsc2VcbiAgY29uc29sZS5sb2coY2FzdGVkIGluc3RhbmNlb2YgUlNWUC5Qcm9taXNlKSAvLyB0cnVlXG5cbiAgY2FzdGVkLnRoZW4oZnVuY3Rpb24odmFsKSB7XG4gICAgdmFsID09PSB2YWx1ZSAvLyA9PiB0cnVlXG4gIH0pO1xuICBgYGBcblxuICBgUlNWUC5Qcm9taXNlLmNhc3RgIGlzIHNpbWlsYXIgdG8gYFJTVlAucmVzb2x2ZWAsIGJ1dCBgUlNWUC5Qcm9taXNlLmNhc3RgIGRpZmZlcnMgaW4gdGhlXG4gIGZvbGxvd2luZyB3YXlzOlxuICAqIGBSU1ZQLlByb21pc2UuY2FzdGAgc2VydmVzIGFzIGEgbWVtb3J5LWVmZmljaWVudCB3YXkgb2YgZ2V0dGluZyBhIHByb21pc2UsIHdoZW4geW91XG4gIGhhdmUgc29tZXRoaW5nIHRoYXQgY291bGQgZWl0aGVyIGJlIGEgcHJvbWlzZSBvciBhIHZhbHVlLiBSU1ZQLnJlc29sdmVcbiAgd2lsbCBoYXZlIHRoZSBzYW1lIGVmZmVjdCBidXQgd2lsbCBjcmVhdGUgYSBuZXcgcHJvbWlzZSB3cmFwcGVyIGlmIHRoZVxuICBhcmd1bWVudCBpcyBhIHByb21pc2UuXG4gICogYFJTVlAuUHJvbWlzZS5jYXN0YCBpcyBhIHdheSBvZiBjYXN0aW5nIGluY29taW5nIHRoZW5hYmxlcyBvciBwcm9taXNlIHN1YmNsYXNzZXMgdG9cbiAgcHJvbWlzZXMgb2YgdGhlIGV4YWN0IGNsYXNzIHNwZWNpZmllZCwgc28gdGhhdCB0aGUgcmVzdWx0aW5nIG9iamVjdCdzIGB0aGVuYCBpc1xuICBlbnN1cmVkIHRvIGhhdmUgdGhlIGJlaGF2aW9yIG9mIHRoZSBjb25zdHJ1Y3RvciB5b3UgYXJlIGNhbGxpbmcgY2FzdCBvbiAoaS5lLiwgUlNWUC5Qcm9taXNlKS5cblxuICBAbWV0aG9kIGNhc3RcbiAgQGZvciBSU1ZQXG4gIEBwYXJhbSB7T2JqZWN0fSBvYmplY3QgdG8gYmUgY2FzdGVkXG4gIEByZXR1cm4ge1Byb21pc2V9IHByb21pc2UgdGhhdCBpcyBmdWxmaWxsZWQgd2hlbiBhbGwgcHJvcGVydGllcyBvZiBgcHJvbWlzZXNgXG4gIGhhdmUgYmVlbiBmdWxmaWxsZWQsIG9yIHJlamVjdGVkIGlmIGFueSBvZiB0aGVtIGJlY29tZSByZWplY3RlZC5cbiovXG5cblxuZnVuY3Rpb24gY2FzdChvYmplY3QpIHtcbiAgLypqc2hpbnQgdmFsaWR0aGlzOnRydWUgKi9cbiAgaWYgKG9iamVjdCAmJiB0eXBlb2Ygb2JqZWN0ID09PSAnb2JqZWN0JyAmJiBvYmplY3QuY29uc3RydWN0b3IgPT09IHRoaXMpIHtcbiAgICByZXR1cm4gb2JqZWN0O1xuICB9XG5cbiAgdmFyIFByb21pc2UgPSB0aGlzO1xuXG4gIHJldHVybiBuZXcgUHJvbWlzZShmdW5jdGlvbihyZXNvbHZlKSB7XG4gICAgcmVzb2x2ZShvYmplY3QpO1xuICB9KTtcbn1cblxuZXhwb3J0cy5jYXN0ID0gY2FzdDsiLCJcInVzZSBzdHJpY3RcIjtcbnZhciBjb25maWcgPSB7XG4gIGluc3RydW1lbnQ6IGZhbHNlXG59O1xuXG5mdW5jdGlvbiBjb25maWd1cmUobmFtZSwgdmFsdWUpIHtcbiAgaWYgKGFyZ3VtZW50cy5sZW5ndGggPT09IDIpIHtcbiAgICBjb25maWdbbmFtZV0gPSB2YWx1ZTtcbiAgfSBlbHNlIHtcbiAgICByZXR1cm4gY29uZmlnW25hbWVdO1xuICB9XG59XG5cbmV4cG9ydHMuY29uZmlnID0gY29uZmlnO1xuZXhwb3J0cy5jb25maWd1cmUgPSBjb25maWd1cmU7IiwiKGZ1bmN0aW9uIChnbG9iYWwpe1xuXCJ1c2Ugc3RyaWN0XCI7XG4vKmdsb2JhbCBzZWxmKi9cbnZhciBSU1ZQUHJvbWlzZSA9IHJlcXVpcmUoXCIuL3Byb21pc2VcIikuUHJvbWlzZTtcbnZhciBpc0Z1bmN0aW9uID0gcmVxdWlyZShcIi4vdXRpbHNcIikuaXNGdW5jdGlvbjtcblxuZnVuY3Rpb24gcG9seWZpbGwoKSB7XG4gIHZhciBsb2NhbDtcblxuICBpZiAodHlwZW9mIGdsb2JhbCAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgICBsb2NhbCA9IGdsb2JhbDtcbiAgfSBlbHNlIGlmICh0eXBlb2Ygd2luZG93ICE9PSAndW5kZWZpbmVkJyAmJiB3aW5kb3cuZG9jdW1lbnQpIHtcbiAgICBsb2NhbCA9IHdpbmRvdztcbiAgfSBlbHNlIHtcbiAgICBsb2NhbCA9IHNlbGY7XG4gIH1cblxuICB2YXIgZXM2UHJvbWlzZVN1cHBvcnQgPSBcbiAgICBcIlByb21pc2VcIiBpbiBsb2NhbCAmJlxuICAgIC8vIFNvbWUgb2YgdGhlc2UgbWV0aG9kcyBhcmUgbWlzc2luZyBmcm9tXG4gICAgLy8gRmlyZWZveC9DaHJvbWUgZXhwZXJpbWVudGFsIGltcGxlbWVudGF0aW9uc1xuICAgIFwiY2FzdFwiIGluIGxvY2FsLlByb21pc2UgJiZcbiAgICBcInJlc29sdmVcIiBpbiBsb2NhbC5Qcm9taXNlICYmXG4gICAgXCJyZWplY3RcIiBpbiBsb2NhbC5Qcm9taXNlICYmXG4gICAgXCJhbGxcIiBpbiBsb2NhbC5Qcm9taXNlICYmXG4gICAgXCJyYWNlXCIgaW4gbG9jYWwuUHJvbWlzZSAmJlxuICAgIC8vIE9sZGVyIHZlcnNpb24gb2YgdGhlIHNwZWMgaGFkIGEgcmVzb2x2ZXIgb2JqZWN0XG4gICAgLy8gYXMgdGhlIGFyZyByYXRoZXIgdGhhbiBhIGZ1bmN0aW9uXG4gICAgKGZ1bmN0aW9uKCkge1xuICAgICAgdmFyIHJlc29sdmU7XG4gICAgICBuZXcgbG9jYWwuUHJvbWlzZShmdW5jdGlvbihyKSB7IHJlc29sdmUgPSByOyB9KTtcbiAgICAgIHJldHVybiBpc0Z1bmN0aW9uKHJlc29sdmUpO1xuICAgIH0oKSk7XG5cbiAgaWYgKCFlczZQcm9taXNlU3VwcG9ydCkge1xuICAgIGxvY2FsLlByb21pc2UgPSBSU1ZQUHJvbWlzZTtcbiAgfVxufVxuXG5leHBvcnRzLnBvbHlmaWxsID0gcG9seWZpbGw7XG59KS5jYWxsKHRoaXMsdHlwZW9mIHNlbGYgIT09IFwidW5kZWZpbmVkXCIgPyBzZWxmIDogdHlwZW9mIHdpbmRvdyAhPT0gXCJ1bmRlZmluZWRcIiA/IHdpbmRvdyA6IHt9KSIsIlwidXNlIHN0cmljdFwiO1xudmFyIGNvbmZpZyA9IHJlcXVpcmUoXCIuL2NvbmZpZ1wiKS5jb25maWc7XG52YXIgY29uZmlndXJlID0gcmVxdWlyZShcIi4vY29uZmlnXCIpLmNvbmZpZ3VyZTtcbnZhciBvYmplY3RPckZ1bmN0aW9uID0gcmVxdWlyZShcIi4vdXRpbHNcIikub2JqZWN0T3JGdW5jdGlvbjtcbnZhciBpc0Z1bmN0aW9uID0gcmVxdWlyZShcIi4vdXRpbHNcIikuaXNGdW5jdGlvbjtcbnZhciBub3cgPSByZXF1aXJlKFwiLi91dGlsc1wiKS5ub3c7XG52YXIgY2FzdCA9IHJlcXVpcmUoXCIuL2Nhc3RcIikuY2FzdDtcbnZhciBhbGwgPSByZXF1aXJlKFwiLi9hbGxcIikuYWxsO1xudmFyIHJhY2UgPSByZXF1aXJlKFwiLi9yYWNlXCIpLnJhY2U7XG52YXIgc3RhdGljUmVzb2x2ZSA9IHJlcXVpcmUoXCIuL3Jlc29sdmVcIikucmVzb2x2ZTtcbnZhciBzdGF0aWNSZWplY3QgPSByZXF1aXJlKFwiLi9yZWplY3RcIikucmVqZWN0O1xudmFyIGFzYXAgPSByZXF1aXJlKFwiLi9hc2FwXCIpLmFzYXA7XG5cbnZhciBjb3VudGVyID0gMDtcblxuY29uZmlnLmFzeW5jID0gYXNhcDsgLy8gZGVmYXVsdCBhc3luYyBpcyBhc2FwO1xuXG5mdW5jdGlvbiBQcm9taXNlKHJlc29sdmVyKSB7XG4gIGlmICghaXNGdW5jdGlvbihyZXNvbHZlcikpIHtcbiAgICB0aHJvdyBuZXcgVHlwZUVycm9yKCdZb3UgbXVzdCBwYXNzIGEgcmVzb2x2ZXIgZnVuY3Rpb24gYXMgdGhlIGZpcnN0IGFyZ3VtZW50IHRvIHRoZSBwcm9taXNlIGNvbnN0cnVjdG9yJyk7XG4gIH1cblxuICBpZiAoISh0aGlzIGluc3RhbmNlb2YgUHJvbWlzZSkpIHtcbiAgICB0aHJvdyBuZXcgVHlwZUVycm9yKFwiRmFpbGVkIHRvIGNvbnN0cnVjdCAnUHJvbWlzZSc6IFBsZWFzZSB1c2UgdGhlICduZXcnIG9wZXJhdG9yLCB0aGlzIG9iamVjdCBjb25zdHJ1Y3RvciBjYW5ub3QgYmUgY2FsbGVkIGFzIGEgZnVuY3Rpb24uXCIpO1xuICB9XG5cbiAgdGhpcy5fc3Vic2NyaWJlcnMgPSBbXTtcblxuICBpbnZva2VSZXNvbHZlcihyZXNvbHZlciwgdGhpcyk7XG59XG5cbmZ1bmN0aW9uIGludm9rZVJlc29sdmVyKHJlc29sdmVyLCBwcm9taXNlKSB7XG4gIGZ1bmN0aW9uIHJlc29sdmVQcm9taXNlKHZhbHVlKSB7XG4gICAgcmVzb2x2ZShwcm9taXNlLCB2YWx1ZSk7XG4gIH1cblxuICBmdW5jdGlvbiByZWplY3RQcm9taXNlKHJlYXNvbikge1xuICAgIHJlamVjdChwcm9taXNlLCByZWFzb24pO1xuICB9XG5cbiAgdHJ5IHtcbiAgICByZXNvbHZlcihyZXNvbHZlUHJvbWlzZSwgcmVqZWN0UHJvbWlzZSk7XG4gIH0gY2F0Y2goZSkge1xuICAgIHJlamVjdFByb21pc2UoZSk7XG4gIH1cbn1cblxuZnVuY3Rpb24gaW52b2tlQ2FsbGJhY2soc2V0dGxlZCwgcHJvbWlzZSwgY2FsbGJhY2ssIGRldGFpbCkge1xuICB2YXIgaGFzQ2FsbGJhY2sgPSBpc0Z1bmN0aW9uKGNhbGxiYWNrKSxcbiAgICAgIHZhbHVlLCBlcnJvciwgc3VjY2VlZGVkLCBmYWlsZWQ7XG5cbiAgaWYgKGhhc0NhbGxiYWNrKSB7XG4gICAgdHJ5IHtcbiAgICAgIHZhbHVlID0gY2FsbGJhY2soZGV0YWlsKTtcbiAgICAgIHN1Y2NlZWRlZCA9IHRydWU7XG4gICAgfSBjYXRjaChlKSB7XG4gICAgICBmYWlsZWQgPSB0cnVlO1xuICAgICAgZXJyb3IgPSBlO1xuICAgIH1cbiAgfSBlbHNlIHtcbiAgICB2YWx1ZSA9IGRldGFpbDtcbiAgICBzdWNjZWVkZWQgPSB0cnVlO1xuICB9XG5cbiAgaWYgKGhhbmRsZVRoZW5hYmxlKHByb21pc2UsIHZhbHVlKSkge1xuICAgIHJldHVybjtcbiAgfSBlbHNlIGlmIChoYXNDYWxsYmFjayAmJiBzdWNjZWVkZWQpIHtcbiAgICByZXNvbHZlKHByb21pc2UsIHZhbHVlKTtcbiAgfSBlbHNlIGlmIChmYWlsZWQpIHtcbiAgICByZWplY3QocHJvbWlzZSwgZXJyb3IpO1xuICB9IGVsc2UgaWYgKHNldHRsZWQgPT09IEZVTEZJTExFRCkge1xuICAgIHJlc29sdmUocHJvbWlzZSwgdmFsdWUpO1xuICB9IGVsc2UgaWYgKHNldHRsZWQgPT09IFJFSkVDVEVEKSB7XG4gICAgcmVqZWN0KHByb21pc2UsIHZhbHVlKTtcbiAgfVxufVxuXG52YXIgUEVORElORyAgID0gdm9pZCAwO1xudmFyIFNFQUxFRCAgICA9IDA7XG52YXIgRlVMRklMTEVEID0gMTtcbnZhciBSRUpFQ1RFRCAgPSAyO1xuXG5mdW5jdGlvbiBzdWJzY3JpYmUocGFyZW50LCBjaGlsZCwgb25GdWxmaWxsbWVudCwgb25SZWplY3Rpb24pIHtcbiAgdmFyIHN1YnNjcmliZXJzID0gcGFyZW50Ll9zdWJzY3JpYmVycztcbiAgdmFyIGxlbmd0aCA9IHN1YnNjcmliZXJzLmxlbmd0aDtcblxuICBzdWJzY3JpYmVyc1tsZW5ndGhdID0gY2hpbGQ7XG4gIHN1YnNjcmliZXJzW2xlbmd0aCArIEZVTEZJTExFRF0gPSBvbkZ1bGZpbGxtZW50O1xuICBzdWJzY3JpYmVyc1tsZW5ndGggKyBSRUpFQ1RFRF0gID0gb25SZWplY3Rpb247XG59XG5cbmZ1bmN0aW9uIHB1Ymxpc2gocHJvbWlzZSwgc2V0dGxlZCkge1xuICB2YXIgY2hpbGQsIGNhbGxiYWNrLCBzdWJzY3JpYmVycyA9IHByb21pc2UuX3N1YnNjcmliZXJzLCBkZXRhaWwgPSBwcm9taXNlLl9kZXRhaWw7XG5cbiAgZm9yICh2YXIgaSA9IDA7IGkgPCBzdWJzY3JpYmVycy5sZW5ndGg7IGkgKz0gMykge1xuICAgIGNoaWxkID0gc3Vic2NyaWJlcnNbaV07XG4gICAgY2FsbGJhY2sgPSBzdWJzY3JpYmVyc1tpICsgc2V0dGxlZF07XG5cbiAgICBpbnZva2VDYWxsYmFjayhzZXR0bGVkLCBjaGlsZCwgY2FsbGJhY2ssIGRldGFpbCk7XG4gIH1cblxuICBwcm9taXNlLl9zdWJzY3JpYmVycyA9IG51bGw7XG59XG5cblByb21pc2UucHJvdG90eXBlID0ge1xuICBjb25zdHJ1Y3RvcjogUHJvbWlzZSxcblxuICBfc3RhdGU6IHVuZGVmaW5lZCxcbiAgX2RldGFpbDogdW5kZWZpbmVkLFxuICBfc3Vic2NyaWJlcnM6IHVuZGVmaW5lZCxcblxuICB0aGVuOiBmdW5jdGlvbihvbkZ1bGZpbGxtZW50LCBvblJlamVjdGlvbikge1xuICAgIHZhciBwcm9taXNlID0gdGhpcztcblxuICAgIHZhciB0aGVuUHJvbWlzZSA9IG5ldyB0aGlzLmNvbnN0cnVjdG9yKGZ1bmN0aW9uKCkge30pO1xuXG4gICAgaWYgKHRoaXMuX3N0YXRlKSB7XG4gICAgICB2YXIgY2FsbGJhY2tzID0gYXJndW1lbnRzO1xuICAgICAgY29uZmlnLmFzeW5jKGZ1bmN0aW9uIGludm9rZVByb21pc2VDYWxsYmFjaygpIHtcbiAgICAgICAgaW52b2tlQ2FsbGJhY2socHJvbWlzZS5fc3RhdGUsIHRoZW5Qcm9taXNlLCBjYWxsYmFja3NbcHJvbWlzZS5fc3RhdGUgLSAxXSwgcHJvbWlzZS5fZGV0YWlsKTtcbiAgICAgIH0pO1xuICAgIH0gZWxzZSB7XG4gICAgICBzdWJzY3JpYmUodGhpcywgdGhlblByb21pc2UsIG9uRnVsZmlsbG1lbnQsIG9uUmVqZWN0aW9uKTtcbiAgICB9XG5cbiAgICByZXR1cm4gdGhlblByb21pc2U7XG4gIH0sXG5cbiAgJ2NhdGNoJzogZnVuY3Rpb24ob25SZWplY3Rpb24pIHtcbiAgICByZXR1cm4gdGhpcy50aGVuKG51bGwsIG9uUmVqZWN0aW9uKTtcbiAgfVxufTtcblxuUHJvbWlzZS5hbGwgPSBhbGw7XG5Qcm9taXNlLmNhc3QgPSBjYXN0O1xuUHJvbWlzZS5yYWNlID0gcmFjZTtcblByb21pc2UucmVzb2x2ZSA9IHN0YXRpY1Jlc29sdmU7XG5Qcm9taXNlLnJlamVjdCA9IHN0YXRpY1JlamVjdDtcblxuZnVuY3Rpb24gaGFuZGxlVGhlbmFibGUocHJvbWlzZSwgdmFsdWUpIHtcbiAgdmFyIHRoZW4gPSBudWxsLFxuICByZXNvbHZlZDtcblxuICB0cnkge1xuICAgIGlmIChwcm9taXNlID09PSB2YWx1ZSkge1xuICAgICAgdGhyb3cgbmV3IFR5cGVFcnJvcihcIkEgcHJvbWlzZXMgY2FsbGJhY2sgY2Fubm90IHJldHVybiB0aGF0IHNhbWUgcHJvbWlzZS5cIik7XG4gICAgfVxuXG4gICAgaWYgKG9iamVjdE9yRnVuY3Rpb24odmFsdWUpKSB7XG4gICAgICB0aGVuID0gdmFsdWUudGhlbjtcblxuICAgICAgaWYgKGlzRnVuY3Rpb24odGhlbikpIHtcbiAgICAgICAgdGhlbi5jYWxsKHZhbHVlLCBmdW5jdGlvbih2YWwpIHtcbiAgICAgICAgICBpZiAocmVzb2x2ZWQpIHsgcmV0dXJuIHRydWU7IH1cbiAgICAgICAgICByZXNvbHZlZCA9IHRydWU7XG5cbiAgICAgICAgICBpZiAodmFsdWUgIT09IHZhbCkge1xuICAgICAgICAgICAgcmVzb2x2ZShwcm9taXNlLCB2YWwpO1xuICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBmdWxmaWxsKHByb21pc2UsIHZhbCk7XG4gICAgICAgICAgfVxuICAgICAgICB9LCBmdW5jdGlvbih2YWwpIHtcbiAgICAgICAgICBpZiAocmVzb2x2ZWQpIHsgcmV0dXJuIHRydWU7IH1cbiAgICAgICAgICByZXNvbHZlZCA9IHRydWU7XG5cbiAgICAgICAgICByZWplY3QocHJvbWlzZSwgdmFsKTtcbiAgICAgICAgfSk7XG5cbiAgICAgICAgcmV0dXJuIHRydWU7XG4gICAgICB9XG4gICAgfVxuICB9IGNhdGNoIChlcnJvcikge1xuICAgIGlmIChyZXNvbHZlZCkgeyByZXR1cm4gdHJ1ZTsgfVxuICAgIHJlamVjdChwcm9taXNlLCBlcnJvcik7XG4gICAgcmV0dXJuIHRydWU7XG4gIH1cblxuICByZXR1cm4gZmFsc2U7XG59XG5cbmZ1bmN0aW9uIHJlc29sdmUocHJvbWlzZSwgdmFsdWUpIHtcbiAgaWYgKHByb21pc2UgPT09IHZhbHVlKSB7XG4gICAgZnVsZmlsbChwcm9taXNlLCB2YWx1ZSk7XG4gIH0gZWxzZSBpZiAoIWhhbmRsZVRoZW5hYmxlKHByb21pc2UsIHZhbHVlKSkge1xuICAgIGZ1bGZpbGwocHJvbWlzZSwgdmFsdWUpO1xuICB9XG59XG5cbmZ1bmN0aW9uIGZ1bGZpbGwocHJvbWlzZSwgdmFsdWUpIHtcbiAgaWYgKHByb21pc2UuX3N0YXRlICE9PSBQRU5ESU5HKSB7IHJldHVybjsgfVxuICBwcm9taXNlLl9zdGF0ZSA9IFNFQUxFRDtcbiAgcHJvbWlzZS5fZGV0YWlsID0gdmFsdWU7XG5cbiAgY29uZmlnLmFzeW5jKHB1Ymxpc2hGdWxmaWxsbWVudCwgcHJvbWlzZSk7XG59XG5cbmZ1bmN0aW9uIHJlamVjdChwcm9taXNlLCByZWFzb24pIHtcbiAgaWYgKHByb21pc2UuX3N0YXRlICE9PSBQRU5ESU5HKSB7IHJldHVybjsgfVxuICBwcm9taXNlLl9zdGF0ZSA9IFNFQUxFRDtcbiAgcHJvbWlzZS5fZGV0YWlsID0gcmVhc29uO1xuXG4gIGNvbmZpZy5hc3luYyhwdWJsaXNoUmVqZWN0aW9uLCBwcm9taXNlKTtcbn1cblxuZnVuY3Rpb24gcHVibGlzaEZ1bGZpbGxtZW50KHByb21pc2UpIHtcbiAgcHVibGlzaChwcm9taXNlLCBwcm9taXNlLl9zdGF0ZSA9IEZVTEZJTExFRCk7XG59XG5cbmZ1bmN0aW9uIHB1Ymxpc2hSZWplY3Rpb24ocHJvbWlzZSkge1xuICBwdWJsaXNoKHByb21pc2UsIHByb21pc2UuX3N0YXRlID0gUkVKRUNURUQpO1xufVxuXG5leHBvcnRzLlByb21pc2UgPSBQcm9taXNlOyIsIlwidXNlIHN0cmljdFwiO1xuLyogZ2xvYmFsIHRvU3RyaW5nICovXG52YXIgaXNBcnJheSA9IHJlcXVpcmUoXCIuL3V0aWxzXCIpLmlzQXJyYXk7XG5cbi8qKlxuICBgUlNWUC5yYWNlYCBhbGxvd3MgeW91IHRvIHdhdGNoIGEgc2VyaWVzIG9mIHByb21pc2VzIGFuZCBhY3QgYXMgc29vbiBhcyB0aGVcbiAgZmlyc3QgcHJvbWlzZSBnaXZlbiB0byB0aGUgYHByb21pc2VzYCBhcmd1bWVudCBmdWxmaWxscyBvciByZWplY3RzLlxuXG4gIEV4YW1wbGU6XG5cbiAgYGBgamF2YXNjcmlwdFxuICB2YXIgcHJvbWlzZTEgPSBuZXcgUlNWUC5Qcm9taXNlKGZ1bmN0aW9uKHJlc29sdmUsIHJlamVjdCl7XG4gICAgc2V0VGltZW91dChmdW5jdGlvbigpe1xuICAgICAgcmVzb2x2ZShcInByb21pc2UgMVwiKTtcbiAgICB9LCAyMDApO1xuICB9KTtcblxuICB2YXIgcHJvbWlzZTIgPSBuZXcgUlNWUC5Qcm9taXNlKGZ1bmN0aW9uKHJlc29sdmUsIHJlamVjdCl7XG4gICAgc2V0VGltZW91dChmdW5jdGlvbigpe1xuICAgICAgcmVzb2x2ZShcInByb21pc2UgMlwiKTtcbiAgICB9LCAxMDApO1xuICB9KTtcblxuICBSU1ZQLnJhY2UoW3Byb21pc2UxLCBwcm9taXNlMl0pLnRoZW4oZnVuY3Rpb24ocmVzdWx0KXtcbiAgICAvLyByZXN1bHQgPT09IFwicHJvbWlzZSAyXCIgYmVjYXVzZSBpdCB3YXMgcmVzb2x2ZWQgYmVmb3JlIHByb21pc2UxXG4gICAgLy8gd2FzIHJlc29sdmVkLlxuICB9KTtcbiAgYGBgXG5cbiAgYFJTVlAucmFjZWAgaXMgZGV0ZXJtaW5pc3RpYyBpbiB0aGF0IG9ubHkgdGhlIHN0YXRlIG9mIHRoZSBmaXJzdCBjb21wbGV0ZWRcbiAgcHJvbWlzZSBtYXR0ZXJzLiBGb3IgZXhhbXBsZSwgZXZlbiBpZiBvdGhlciBwcm9taXNlcyBnaXZlbiB0byB0aGUgYHByb21pc2VzYFxuICBhcnJheSBhcmd1bWVudCBhcmUgcmVzb2x2ZWQsIGJ1dCB0aGUgZmlyc3QgY29tcGxldGVkIHByb21pc2UgaGFzIGJlY29tZVxuICByZWplY3RlZCBiZWZvcmUgdGhlIG90aGVyIHByb21pc2VzIGJlY2FtZSBmdWxmaWxsZWQsIHRoZSByZXR1cm5lZCBwcm9taXNlXG4gIHdpbGwgYmVjb21lIHJlamVjdGVkOlxuXG4gIGBgYGphdmFzY3JpcHRcbiAgdmFyIHByb21pc2UxID0gbmV3IFJTVlAuUHJvbWlzZShmdW5jdGlvbihyZXNvbHZlLCByZWplY3Qpe1xuICAgIHNldFRpbWVvdXQoZnVuY3Rpb24oKXtcbiAgICAgIHJlc29sdmUoXCJwcm9taXNlIDFcIik7XG4gICAgfSwgMjAwKTtcbiAgfSk7XG5cbiAgdmFyIHByb21pc2UyID0gbmV3IFJTVlAuUHJvbWlzZShmdW5jdGlvbihyZXNvbHZlLCByZWplY3Qpe1xuICAgIHNldFRpbWVvdXQoZnVuY3Rpb24oKXtcbiAgICAgIHJlamVjdChuZXcgRXJyb3IoXCJwcm9taXNlIDJcIikpO1xuICAgIH0sIDEwMCk7XG4gIH0pO1xuXG4gIFJTVlAucmFjZShbcHJvbWlzZTEsIHByb21pc2UyXSkudGhlbihmdW5jdGlvbihyZXN1bHQpe1xuICAgIC8vIENvZGUgaGVyZSBuZXZlciBydW5zIGJlY2F1c2UgdGhlcmUgYXJlIHJlamVjdGVkIHByb21pc2VzIVxuICB9LCBmdW5jdGlvbihyZWFzb24pe1xuICAgIC8vIHJlYXNvbi5tZXNzYWdlID09PSBcInByb21pc2UyXCIgYmVjYXVzZSBwcm9taXNlIDIgYmVjYW1lIHJlamVjdGVkIGJlZm9yZVxuICAgIC8vIHByb21pc2UgMSBiZWNhbWUgZnVsZmlsbGVkXG4gIH0pO1xuICBgYGBcblxuICBAbWV0aG9kIHJhY2VcbiAgQGZvciBSU1ZQXG4gIEBwYXJhbSB7QXJyYXl9IHByb21pc2VzIGFycmF5IG9mIHByb21pc2VzIHRvIG9ic2VydmVcbiAgQHBhcmFtIHtTdHJpbmd9IGxhYmVsIG9wdGlvbmFsIHN0cmluZyBmb3IgZGVzY3JpYmluZyB0aGUgcHJvbWlzZSByZXR1cm5lZC5cbiAgVXNlZnVsIGZvciB0b29saW5nLlxuICBAcmV0dXJuIHtQcm9taXNlfSBhIHByb21pc2UgdGhhdCBiZWNvbWVzIGZ1bGZpbGxlZCB3aXRoIHRoZSB2YWx1ZSB0aGUgZmlyc3RcbiAgY29tcGxldGVkIHByb21pc2VzIGlzIHJlc29sdmVkIHdpdGggaWYgdGhlIGZpcnN0IGNvbXBsZXRlZCBwcm9taXNlIHdhc1xuICBmdWxmaWxsZWQsIG9yIHJlamVjdGVkIHdpdGggdGhlIHJlYXNvbiB0aGF0IHRoZSBmaXJzdCBjb21wbGV0ZWQgcHJvbWlzZVxuICB3YXMgcmVqZWN0ZWQgd2l0aC5cbiovXG5mdW5jdGlvbiByYWNlKHByb21pc2VzKSB7XG4gIC8qanNoaW50IHZhbGlkdGhpczp0cnVlICovXG4gIHZhciBQcm9taXNlID0gdGhpcztcblxuICBpZiAoIWlzQXJyYXkocHJvbWlzZXMpKSB7XG4gICAgdGhyb3cgbmV3IFR5cGVFcnJvcignWW91IG11c3QgcGFzcyBhbiBhcnJheSB0byByYWNlLicpO1xuICB9XG4gIHJldHVybiBuZXcgUHJvbWlzZShmdW5jdGlvbihyZXNvbHZlLCByZWplY3QpIHtcbiAgICB2YXIgcmVzdWx0cyA9IFtdLCBwcm9taXNlO1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBwcm9taXNlcy5sZW5ndGg7IGkrKykge1xuICAgICAgcHJvbWlzZSA9IHByb21pc2VzW2ldO1xuXG4gICAgICBpZiAocHJvbWlzZSAmJiB0eXBlb2YgcHJvbWlzZS50aGVuID09PSAnZnVuY3Rpb24nKSB7XG4gICAgICAgIHByb21pc2UudGhlbihyZXNvbHZlLCByZWplY3QpO1xuICAgICAgfSBlbHNlIHtcbiAgICAgICAgcmVzb2x2ZShwcm9taXNlKTtcbiAgICAgIH1cbiAgICB9XG4gIH0pO1xufVxuXG5leHBvcnRzLnJhY2UgPSByYWNlOyIsIlwidXNlIHN0cmljdFwiO1xuLyoqXG4gIGBSU1ZQLnJlamVjdGAgcmV0dXJucyBhIHByb21pc2UgdGhhdCB3aWxsIGJlY29tZSByZWplY3RlZCB3aXRoIHRoZSBwYXNzZWRcbiAgYHJlYXNvbmAuIGBSU1ZQLnJlamVjdGAgaXMgZXNzZW50aWFsbHkgc2hvcnRoYW5kIGZvciB0aGUgZm9sbG93aW5nOlxuXG4gIGBgYGphdmFzY3JpcHRcbiAgdmFyIHByb21pc2UgPSBuZXcgUlNWUC5Qcm9taXNlKGZ1bmN0aW9uKHJlc29sdmUsIHJlamVjdCl7XG4gICAgcmVqZWN0KG5ldyBFcnJvcignV0hPT1BTJykpO1xuICB9KTtcblxuICBwcm9taXNlLnRoZW4oZnVuY3Rpb24odmFsdWUpe1xuICAgIC8vIENvZGUgaGVyZSBkb2Vzbid0IHJ1biBiZWNhdXNlIHRoZSBwcm9taXNlIGlzIHJlamVjdGVkIVxuICB9LCBmdW5jdGlvbihyZWFzb24pe1xuICAgIC8vIHJlYXNvbi5tZXNzYWdlID09PSAnV0hPT1BTJ1xuICB9KTtcbiAgYGBgXG5cbiAgSW5zdGVhZCBvZiB3cml0aW5nIHRoZSBhYm92ZSwgeW91ciBjb2RlIG5vdyBzaW1wbHkgYmVjb21lcyB0aGUgZm9sbG93aW5nOlxuXG4gIGBgYGphdmFzY3JpcHRcbiAgdmFyIHByb21pc2UgPSBSU1ZQLnJlamVjdChuZXcgRXJyb3IoJ1dIT09QUycpKTtcblxuICBwcm9taXNlLnRoZW4oZnVuY3Rpb24odmFsdWUpe1xuICAgIC8vIENvZGUgaGVyZSBkb2Vzbid0IHJ1biBiZWNhdXNlIHRoZSBwcm9taXNlIGlzIHJlamVjdGVkIVxuICB9LCBmdW5jdGlvbihyZWFzb24pe1xuICAgIC8vIHJlYXNvbi5tZXNzYWdlID09PSAnV0hPT1BTJ1xuICB9KTtcbiAgYGBgXG5cbiAgQG1ldGhvZCByZWplY3RcbiAgQGZvciBSU1ZQXG4gIEBwYXJhbSB7QW55fSByZWFzb24gdmFsdWUgdGhhdCB0aGUgcmV0dXJuZWQgcHJvbWlzZSB3aWxsIGJlIHJlamVjdGVkIHdpdGguXG4gIEBwYXJhbSB7U3RyaW5nfSBsYWJlbCBvcHRpb25hbCBzdHJpbmcgZm9yIGlkZW50aWZ5aW5nIHRoZSByZXR1cm5lZCBwcm9taXNlLlxuICBVc2VmdWwgZm9yIHRvb2xpbmcuXG4gIEByZXR1cm4ge1Byb21pc2V9IGEgcHJvbWlzZSB0aGF0IHdpbGwgYmVjb21lIHJlamVjdGVkIHdpdGggdGhlIGdpdmVuXG4gIGByZWFzb25gLlxuKi9cbmZ1bmN0aW9uIHJlamVjdChyZWFzb24pIHtcbiAgLypqc2hpbnQgdmFsaWR0aGlzOnRydWUgKi9cbiAgdmFyIFByb21pc2UgPSB0aGlzO1xuXG4gIHJldHVybiBuZXcgUHJvbWlzZShmdW5jdGlvbiAocmVzb2x2ZSwgcmVqZWN0KSB7XG4gICAgcmVqZWN0KHJlYXNvbik7XG4gIH0pO1xufVxuXG5leHBvcnRzLnJlamVjdCA9IHJlamVjdDsiLCJcInVzZSBzdHJpY3RcIjtcbi8qKlxuICBgUlNWUC5yZXNvbHZlYCByZXR1cm5zIGEgcHJvbWlzZSB0aGF0IHdpbGwgYmVjb21lIGZ1bGZpbGxlZCB3aXRoIHRoZSBwYXNzZWRcbiAgYHZhbHVlYC4gYFJTVlAucmVzb2x2ZWAgaXMgZXNzZW50aWFsbHkgc2hvcnRoYW5kIGZvciB0aGUgZm9sbG93aW5nOlxuXG4gIGBgYGphdmFzY3JpcHRcbiAgdmFyIHByb21pc2UgPSBuZXcgUlNWUC5Qcm9taXNlKGZ1bmN0aW9uKHJlc29sdmUsIHJlamVjdCl7XG4gICAgcmVzb2x2ZSgxKTtcbiAgfSk7XG5cbiAgcHJvbWlzZS50aGVuKGZ1bmN0aW9uKHZhbHVlKXtcbiAgICAvLyB2YWx1ZSA9PT0gMVxuICB9KTtcbiAgYGBgXG5cbiAgSW5zdGVhZCBvZiB3cml0aW5nIHRoZSBhYm92ZSwgeW91ciBjb2RlIG5vdyBzaW1wbHkgYmVjb21lcyB0aGUgZm9sbG93aW5nOlxuXG4gIGBgYGphdmFzY3JpcHRcbiAgdmFyIHByb21pc2UgPSBSU1ZQLnJlc29sdmUoMSk7XG5cbiAgcHJvbWlzZS50aGVuKGZ1bmN0aW9uKHZhbHVlKXtcbiAgICAvLyB2YWx1ZSA9PT0gMVxuICB9KTtcbiAgYGBgXG5cbiAgQG1ldGhvZCByZXNvbHZlXG4gIEBmb3IgUlNWUFxuICBAcGFyYW0ge0FueX0gdmFsdWUgdmFsdWUgdGhhdCB0aGUgcmV0dXJuZWQgcHJvbWlzZSB3aWxsIGJlIHJlc29sdmVkIHdpdGhcbiAgQHBhcmFtIHtTdHJpbmd9IGxhYmVsIG9wdGlvbmFsIHN0cmluZyBmb3IgaWRlbnRpZnlpbmcgdGhlIHJldHVybmVkIHByb21pc2UuXG4gIFVzZWZ1bCBmb3IgdG9vbGluZy5cbiAgQHJldHVybiB7UHJvbWlzZX0gYSBwcm9taXNlIHRoYXQgd2lsbCBiZWNvbWUgZnVsZmlsbGVkIHdpdGggdGhlIGdpdmVuXG4gIGB2YWx1ZWBcbiovXG5mdW5jdGlvbiByZXNvbHZlKHZhbHVlKSB7XG4gIC8qanNoaW50IHZhbGlkdGhpczp0cnVlICovXG4gIHZhciBQcm9taXNlID0gdGhpcztcbiAgcmV0dXJuIG5ldyBQcm9taXNlKGZ1bmN0aW9uKHJlc29sdmUsIHJlamVjdCkge1xuICAgIHJlc29sdmUodmFsdWUpO1xuICB9KTtcbn1cblxuZXhwb3J0cy5yZXNvbHZlID0gcmVzb2x2ZTsiLCJcInVzZSBzdHJpY3RcIjtcbmZ1bmN0aW9uIG9iamVjdE9yRnVuY3Rpb24oeCkge1xuICByZXR1cm4gaXNGdW5jdGlvbih4KSB8fCAodHlwZW9mIHggPT09IFwib2JqZWN0XCIgJiYgeCAhPT0gbnVsbCk7XG59XG5cbmZ1bmN0aW9uIGlzRnVuY3Rpb24oeCkge1xuICByZXR1cm4gdHlwZW9mIHggPT09IFwiZnVuY3Rpb25cIjtcbn1cblxuZnVuY3Rpb24gaXNBcnJheSh4KSB7XG4gIHJldHVybiBPYmplY3QucHJvdG90eXBlLnRvU3RyaW5nLmNhbGwoeCkgPT09IFwiW29iamVjdCBBcnJheV1cIjtcbn1cblxuLy8gRGF0ZS5ub3cgaXMgbm90IGF2YWlsYWJsZSBpbiBicm93c2VycyA8IElFOVxuLy8gaHR0cHM6Ly9kZXZlbG9wZXIubW96aWxsYS5vcmcvZW4tVVMvZG9jcy9XZWIvSmF2YVNjcmlwdC9SZWZlcmVuY2UvR2xvYmFsX09iamVjdHMvRGF0ZS9ub3cjQ29tcGF0aWJpbGl0eVxudmFyIG5vdyA9IERhdGUubm93IHx8IGZ1bmN0aW9uKCkgeyByZXR1cm4gbmV3IERhdGUoKS5nZXRUaW1lKCk7IH07XG5cblxuZXhwb3J0cy5vYmplY3RPckZ1bmN0aW9uID0gb2JqZWN0T3JGdW5jdGlvbjtcbmV4cG9ydHMuaXNGdW5jdGlvbiA9IGlzRnVuY3Rpb247XG5leHBvcnRzLmlzQXJyYXkgPSBpc0FycmF5O1xuZXhwb3J0cy5ub3cgPSBub3c7IiwiLy8gc2hpbSBmb3IgdXNpbmcgcHJvY2VzcyBpbiBicm93c2VyXG5cbnZhciBwcm9jZXNzID0gbW9kdWxlLmV4cG9ydHMgPSB7fTtcblxucHJvY2Vzcy5uZXh0VGljayA9IChmdW5jdGlvbiAoKSB7XG4gICAgdmFyIGNhblNldEltbWVkaWF0ZSA9IHR5cGVvZiB3aW5kb3cgIT09ICd1bmRlZmluZWQnXG4gICAgJiYgd2luZG93LnNldEltbWVkaWF0ZTtcbiAgICB2YXIgY2FuUG9zdCA9IHR5cGVvZiB3aW5kb3cgIT09ICd1bmRlZmluZWQnXG4gICAgJiYgd2luZG93LnBvc3RNZXNzYWdlICYmIHdpbmRvdy5hZGRFdmVudExpc3RlbmVyXG4gICAgO1xuXG4gICAgaWYgKGNhblNldEltbWVkaWF0ZSkge1xuICAgICAgICByZXR1cm4gZnVuY3Rpb24gKGYpIHsgcmV0dXJuIHdpbmRvdy5zZXRJbW1lZGlhdGUoZikgfTtcbiAgICB9XG5cbiAgICBpZiAoY2FuUG9zdCkge1xuICAgICAgICB2YXIgcXVldWUgPSBbXTtcbiAgICAgICAgd2luZG93LmFkZEV2ZW50TGlzdGVuZXIoJ21lc3NhZ2UnLCBmdW5jdGlvbiAoZXYpIHtcbiAgICAgICAgICAgIHZhciBzb3VyY2UgPSBldi5zb3VyY2U7XG4gICAgICAgICAgICBpZiAoKHNvdXJjZSA9PT0gd2luZG93IHx8IHNvdXJjZSA9PT0gbnVsbCkgJiYgZXYuZGF0YSA9PT0gJ3Byb2Nlc3MtdGljaycpIHtcbiAgICAgICAgICAgICAgICBldi5zdG9wUHJvcGFnYXRpb24oKTtcbiAgICAgICAgICAgICAgICBpZiAocXVldWUubGVuZ3RoID4gMCkge1xuICAgICAgICAgICAgICAgICAgICB2YXIgZm4gPSBxdWV1ZS5zaGlmdCgpO1xuICAgICAgICAgICAgICAgICAgICBmbigpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH1cbiAgICAgICAgfSwgdHJ1ZSk7XG5cbiAgICAgICAgcmV0dXJuIGZ1bmN0aW9uIG5leHRUaWNrKGZuKSB7XG4gICAgICAgICAgICBxdWV1ZS5wdXNoKGZuKTtcbiAgICAgICAgICAgIHdpbmRvdy5wb3N0TWVzc2FnZSgncHJvY2Vzcy10aWNrJywgJyonKTtcbiAgICAgICAgfTtcbiAgICB9XG5cbiAgICByZXR1cm4gZnVuY3Rpb24gbmV4dFRpY2soZm4pIHtcbiAgICAgICAgc2V0VGltZW91dChmbiwgMCk7XG4gICAgfTtcbn0pKCk7XG5cbnByb2Nlc3MudGl0bGUgPSAnYnJvd3Nlcic7XG5wcm9jZXNzLmJyb3dzZXIgPSB0cnVlO1xucHJvY2Vzcy5lbnYgPSB7fTtcbnByb2Nlc3MuYXJndiA9IFtdO1xuXG5mdW5jdGlvbiBub29wKCkge31cblxucHJvY2Vzcy5vbiA9IG5vb3A7XG5wcm9jZXNzLmFkZExpc3RlbmVyID0gbm9vcDtcbnByb2Nlc3Mub25jZSA9IG5vb3A7XG5wcm9jZXNzLm9mZiA9IG5vb3A7XG5wcm9jZXNzLnJlbW92ZUxpc3RlbmVyID0gbm9vcDtcbnByb2Nlc3MucmVtb3ZlQWxsTGlzdGVuZXJzID0gbm9vcDtcbnByb2Nlc3MuZW1pdCA9IG5vb3A7XG5cbnByb2Nlc3MuYmluZGluZyA9IGZ1bmN0aW9uIChuYW1lKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKCdwcm9jZXNzLmJpbmRpbmcgaXMgbm90IHN1cHBvcnRlZCcpO1xufVxuXG4vLyBUT0RPKHNodHlsbWFuKVxucHJvY2Vzcy5jd2QgPSBmdW5jdGlvbiAoKSB7IHJldHVybiAnLycgfTtcbnByb2Nlc3MuY2hkaXIgPSBmdW5jdGlvbiAoZGlyKSB7XG4gICAgdGhyb3cgbmV3IEVycm9yKCdwcm9jZXNzLmNoZGlyIGlzIG5vdCBzdXBwb3J0ZWQnKTtcbn07XG4iLCIvKiAtKi0gbW9kZTogamF2YXNjcmlwdDsgYy1iYXNpYy1vZmZzZXQ6IDQ7IGluZGVudC10YWJzLW1vZGU6IG5pbCAtKi0gKi9cblxuLy8gXG4vLyBKYXZhc2NyaXB0IFpMaWJcbi8vIEJ5IFRob21hcyBEb3duIDIwMTAtMjAxMVxuLy9cbi8vIEJhc2VkIHZlcnkgaGVhdmlseSBvbiBwb3J0aW9ucyBvZiBqemxpYiAoYnkgeW1ua0BqY3JhZnQuY29tKSwgd2hvIGluXG4vLyB0dXJuIGNyZWRpdHMgSmVhbi1sb3VwIEdhaWxseSBhbmQgTWFyayBBZGxlciBmb3IgdGhlIG9yaWdpbmFsIHpsaWIgY29kZS5cbi8vXG4vLyBpbmZsYXRlLmpzOiBaTGliIGluZmxhdGUgY29kZVxuLy9cblxuLy9cbi8vIFNoYXJlZCBjb25zdGFudHNcbi8vXG5cbnZhciBNQVhfV0JJVFM9MTU7IC8vIDMySyBMWjc3IHdpbmRvd1xudmFyIERFRl9XQklUUz1NQVhfV0JJVFM7XG52YXIgTUFYX01FTV9MRVZFTD05O1xudmFyIE1BTlk9MTQ0MDtcbnZhciBCTUFYID0gMTU7XG5cbi8vIHByZXNldCBkaWN0aW9uYXJ5IGZsYWcgaW4gemxpYiBoZWFkZXJcbnZhciBQUkVTRVRfRElDVD0weDIwO1xuXG52YXIgWl9OT19GTFVTSD0wO1xudmFyIFpfUEFSVElBTF9GTFVTSD0xO1xudmFyIFpfU1lOQ19GTFVTSD0yO1xudmFyIFpfRlVMTF9GTFVTSD0zO1xudmFyIFpfRklOSVNIPTQ7XG5cbnZhciBaX0RFRkxBVEVEPTg7XG5cbnZhciBaX09LPTA7XG52YXIgWl9TVFJFQU1fRU5EPTE7XG52YXIgWl9ORUVEX0RJQ1Q9MjtcbnZhciBaX0VSUk5PPS0xO1xudmFyIFpfU1RSRUFNX0VSUk9SPS0yO1xudmFyIFpfREFUQV9FUlJPUj0tMztcbnZhciBaX01FTV9FUlJPUj0tNDtcbnZhciBaX0JVRl9FUlJPUj0tNTtcbnZhciBaX1ZFUlNJT05fRVJST1I9LTY7XG5cbnZhciBNRVRIT0Q9MDsgICAvLyB3YWl0aW5nIGZvciBtZXRob2QgYnl0ZVxudmFyIEZMQUc9MTsgICAgIC8vIHdhaXRpbmcgZm9yIGZsYWcgYnl0ZVxudmFyIERJQ1Q0PTI7ICAgIC8vIGZvdXIgZGljdGlvbmFyeSBjaGVjayBieXRlcyB0byBnb1xudmFyIERJQ1QzPTM7ICAgIC8vIHRocmVlIGRpY3Rpb25hcnkgY2hlY2sgYnl0ZXMgdG8gZ29cbnZhciBESUNUMj00OyAgICAvLyB0d28gZGljdGlvbmFyeSBjaGVjayBieXRlcyB0byBnb1xudmFyIERJQ1QxPTU7ICAgIC8vIG9uZSBkaWN0aW9uYXJ5IGNoZWNrIGJ5dGUgdG8gZ29cbnZhciBESUNUMD02OyAgICAvLyB3YWl0aW5nIGZvciBpbmZsYXRlU2V0RGljdGlvbmFyeVxudmFyIEJMT0NLUz03OyAgIC8vIGRlY29tcHJlc3NpbmcgYmxvY2tzXG52YXIgQ0hFQ0s0PTg7ICAgLy8gZm91ciBjaGVjayBieXRlcyB0byBnb1xudmFyIENIRUNLMz05OyAgIC8vIHRocmVlIGNoZWNrIGJ5dGVzIHRvIGdvXG52YXIgQ0hFQ0syPTEwOyAgLy8gdHdvIGNoZWNrIGJ5dGVzIHRvIGdvXG52YXIgQ0hFQ0sxPTExOyAgLy8gb25lIGNoZWNrIGJ5dGUgdG8gZ29cbnZhciBET05FPTEyOyAgICAvLyBmaW5pc2hlZCBjaGVjaywgZG9uZVxudmFyIEJBRD0xMzsgICAgIC8vIGdvdCBhbiBlcnJvci0tc3RheSBoZXJlXG5cbnZhciBpbmZsYXRlX21hc2sgPSBbMHgwMDAwMDAwMCwgMHgwMDAwMDAwMSwgMHgwMDAwMDAwMywgMHgwMDAwMDAwNywgMHgwMDAwMDAwZiwgMHgwMDAwMDAxZiwgMHgwMDAwMDAzZiwgMHgwMDAwMDA3ZiwgMHgwMDAwMDBmZiwgMHgwMDAwMDFmZiwgMHgwMDAwMDNmZiwgMHgwMDAwMDdmZiwgMHgwMDAwMGZmZiwgMHgwMDAwMWZmZiwgMHgwMDAwM2ZmZiwgMHgwMDAwN2ZmZiwgMHgwMDAwZmZmZl07XG5cbnZhciBJQl9UWVBFPTA7ICAvLyBnZXQgdHlwZSBiaXRzICgzLCBpbmNsdWRpbmcgZW5kIGJpdClcbnZhciBJQl9MRU5TPTE7ICAvLyBnZXQgbGVuZ3RocyBmb3Igc3RvcmVkXG52YXIgSUJfU1RPUkVEPTI7Ly8gcHJvY2Vzc2luZyBzdG9yZWQgYmxvY2tcbnZhciBJQl9UQUJMRT0zOyAvLyBnZXQgdGFibGUgbGVuZ3Roc1xudmFyIElCX0JUUkVFPTQ7IC8vIGdldCBiaXQgbGVuZ3RocyB0cmVlIGZvciBhIGR5bmFtaWMgYmxvY2tcbnZhciBJQl9EVFJFRT01OyAvLyBnZXQgbGVuZ3RoLCBkaXN0YW5jZSB0cmVlcyBmb3IgYSBkeW5hbWljIGJsb2NrXG52YXIgSUJfQ09ERVM9NjsgLy8gcHJvY2Vzc2luZyBmaXhlZCBvciBkeW5hbWljIGJsb2NrXG52YXIgSUJfRFJZPTc7ICAgLy8gb3V0cHV0IHJlbWFpbmluZyB3aW5kb3cgYnl0ZXNcbnZhciBJQl9ET05FPTg7ICAvLyBmaW5pc2hlZCBsYXN0IGJsb2NrLCBkb25lXG52YXIgSUJfQkFEPTk7ICAgLy8gb3QgYSBkYXRhIGVycm9yLS1zdHVjayBoZXJlXG5cbnZhciBmaXhlZF9ibCA9IDk7XG52YXIgZml4ZWRfYmQgPSA1O1xuXG52YXIgZml4ZWRfdGwgPSBbXG4gICAgOTYsNywyNTYsIDAsOCw4MCwgMCw4LDE2LCA4NCw4LDExNSxcbiAgICA4Miw3LDMxLCAwLDgsMTEyLCAwLDgsNDgsIDAsOSwxOTIsXG4gICAgODAsNywxMCwgMCw4LDk2LCAwLDgsMzIsIDAsOSwxNjAsXG4gICAgMCw4LDAsIDAsOCwxMjgsIDAsOCw2NCwgMCw5LDIyNCxcbiAgICA4MCw3LDYsIDAsOCw4OCwgMCw4LDI0LCAwLDksMTQ0LFxuICAgIDgzLDcsNTksIDAsOCwxMjAsIDAsOCw1NiwgMCw5LDIwOCxcbiAgICA4MSw3LDE3LCAwLDgsMTA0LCAwLDgsNDAsIDAsOSwxNzYsXG4gICAgMCw4LDgsIDAsOCwxMzYsIDAsOCw3MiwgMCw5LDI0MCxcbiAgICA4MCw3LDQsIDAsOCw4NCwgMCw4LDIwLCA4NSw4LDIyNyxcbiAgICA4Myw3LDQzLCAwLDgsMTE2LCAwLDgsNTIsIDAsOSwyMDAsXG4gICAgODEsNywxMywgMCw4LDEwMCwgMCw4LDM2LCAwLDksMTY4LFxuICAgIDAsOCw0LCAwLDgsMTMyLCAwLDgsNjgsIDAsOSwyMzIsXG4gICAgODAsNyw4LCAwLDgsOTIsIDAsOCwyOCwgMCw5LDE1MixcbiAgICA4NCw3LDgzLCAwLDgsMTI0LCAwLDgsNjAsIDAsOSwyMTYsXG4gICAgODIsNywyMywgMCw4LDEwOCwgMCw4LDQ0LCAwLDksMTg0LFxuICAgIDAsOCwxMiwgMCw4LDE0MCwgMCw4LDc2LCAwLDksMjQ4LFxuICAgIDgwLDcsMywgMCw4LDgyLCAwLDgsMTgsIDg1LDgsMTYzLFxuICAgIDgzLDcsMzUsIDAsOCwxMTQsIDAsOCw1MCwgMCw5LDE5NixcbiAgICA4MSw3LDExLCAwLDgsOTgsIDAsOCwzNCwgMCw5LDE2NCxcbiAgICAwLDgsMiwgMCw4LDEzMCwgMCw4LDY2LCAwLDksMjI4LFxuICAgIDgwLDcsNywgMCw4LDkwLCAwLDgsMjYsIDAsOSwxNDgsXG4gICAgODQsNyw2NywgMCw4LDEyMiwgMCw4LDU4LCAwLDksMjEyLFxuICAgIDgyLDcsMTksIDAsOCwxMDYsIDAsOCw0MiwgMCw5LDE4MCxcbiAgICAwLDgsMTAsIDAsOCwxMzgsIDAsOCw3NCwgMCw5LDI0NCxcbiAgICA4MCw3LDUsIDAsOCw4NiwgMCw4LDIyLCAxOTIsOCwwLFxuICAgIDgzLDcsNTEsIDAsOCwxMTgsIDAsOCw1NCwgMCw5LDIwNCxcbiAgICA4MSw3LDE1LCAwLDgsMTAyLCAwLDgsMzgsIDAsOSwxNzIsXG4gICAgMCw4LDYsIDAsOCwxMzQsIDAsOCw3MCwgMCw5LDIzNixcbiAgICA4MCw3LDksIDAsOCw5NCwgMCw4LDMwLCAwLDksMTU2LFxuICAgIDg0LDcsOTksIDAsOCwxMjYsIDAsOCw2MiwgMCw5LDIyMCxcbiAgICA4Miw3LDI3LCAwLDgsMTEwLCAwLDgsNDYsIDAsOSwxODgsXG4gICAgMCw4LDE0LCAwLDgsMTQyLCAwLDgsNzgsIDAsOSwyNTIsXG4gICAgOTYsNywyNTYsIDAsOCw4MSwgMCw4LDE3LCA4NSw4LDEzMSxcbiAgICA4Miw3LDMxLCAwLDgsMTEzLCAwLDgsNDksIDAsOSwxOTQsXG4gICAgODAsNywxMCwgMCw4LDk3LCAwLDgsMzMsIDAsOSwxNjIsXG4gICAgMCw4LDEsIDAsOCwxMjksIDAsOCw2NSwgMCw5LDIyNixcbiAgICA4MCw3LDYsIDAsOCw4OSwgMCw4LDI1LCAwLDksMTQ2LFxuICAgIDgzLDcsNTksIDAsOCwxMjEsIDAsOCw1NywgMCw5LDIxMCxcbiAgICA4MSw3LDE3LCAwLDgsMTA1LCAwLDgsNDEsIDAsOSwxNzgsXG4gICAgMCw4LDksIDAsOCwxMzcsIDAsOCw3MywgMCw5LDI0MixcbiAgICA4MCw3LDQsIDAsOCw4NSwgMCw4LDIxLCA4MCw4LDI1OCxcbiAgICA4Myw3LDQzLCAwLDgsMTE3LCAwLDgsNTMsIDAsOSwyMDIsXG4gICAgODEsNywxMywgMCw4LDEwMSwgMCw4LDM3LCAwLDksMTcwLFxuICAgIDAsOCw1LCAwLDgsMTMzLCAwLDgsNjksIDAsOSwyMzQsXG4gICAgODAsNyw4LCAwLDgsOTMsIDAsOCwyOSwgMCw5LDE1NCxcbiAgICA4NCw3LDgzLCAwLDgsMTI1LCAwLDgsNjEsIDAsOSwyMTgsXG4gICAgODIsNywyMywgMCw4LDEwOSwgMCw4LDQ1LCAwLDksMTg2LFxuICAgIDAsOCwxMywgMCw4LDE0MSwgMCw4LDc3LCAwLDksMjUwLFxuICAgIDgwLDcsMywgMCw4LDgzLCAwLDgsMTksIDg1LDgsMTk1LFxuICAgIDgzLDcsMzUsIDAsOCwxMTUsIDAsOCw1MSwgMCw5LDE5OCxcbiAgICA4MSw3LDExLCAwLDgsOTksIDAsOCwzNSwgMCw5LDE2NixcbiAgICAwLDgsMywgMCw4LDEzMSwgMCw4LDY3LCAwLDksMjMwLFxuICAgIDgwLDcsNywgMCw4LDkxLCAwLDgsMjcsIDAsOSwxNTAsXG4gICAgODQsNyw2NywgMCw4LDEyMywgMCw4LDU5LCAwLDksMjE0LFxuICAgIDgyLDcsMTksIDAsOCwxMDcsIDAsOCw0MywgMCw5LDE4MixcbiAgICAwLDgsMTEsIDAsOCwxMzksIDAsOCw3NSwgMCw5LDI0NixcbiAgICA4MCw3LDUsIDAsOCw4NywgMCw4LDIzLCAxOTIsOCwwLFxuICAgIDgzLDcsNTEsIDAsOCwxMTksIDAsOCw1NSwgMCw5LDIwNixcbiAgICA4MSw3LDE1LCAwLDgsMTAzLCAwLDgsMzksIDAsOSwxNzQsXG4gICAgMCw4LDcsIDAsOCwxMzUsIDAsOCw3MSwgMCw5LDIzOCxcbiAgICA4MCw3LDksIDAsOCw5NSwgMCw4LDMxLCAwLDksMTU4LFxuICAgIDg0LDcsOTksIDAsOCwxMjcsIDAsOCw2MywgMCw5LDIyMixcbiAgICA4Miw3LDI3LCAwLDgsMTExLCAwLDgsNDcsIDAsOSwxOTAsXG4gICAgMCw4LDE1LCAwLDgsMTQzLCAwLDgsNzksIDAsOSwyNTQsXG4gICAgOTYsNywyNTYsIDAsOCw4MCwgMCw4LDE2LCA4NCw4LDExNSxcbiAgICA4Miw3LDMxLCAwLDgsMTEyLCAwLDgsNDgsIDAsOSwxOTMsXG5cbiAgICA4MCw3LDEwLCAwLDgsOTYsIDAsOCwzMiwgMCw5LDE2MSxcbiAgICAwLDgsMCwgMCw4LDEyOCwgMCw4LDY0LCAwLDksMjI1LFxuICAgIDgwLDcsNiwgMCw4LDg4LCAwLDgsMjQsIDAsOSwxNDUsXG4gICAgODMsNyw1OSwgMCw4LDEyMCwgMCw4LDU2LCAwLDksMjA5LFxuICAgIDgxLDcsMTcsIDAsOCwxMDQsIDAsOCw0MCwgMCw5LDE3NyxcbiAgICAwLDgsOCwgMCw4LDEzNiwgMCw4LDcyLCAwLDksMjQxLFxuICAgIDgwLDcsNCwgMCw4LDg0LCAwLDgsMjAsIDg1LDgsMjI3LFxuICAgIDgzLDcsNDMsIDAsOCwxMTYsIDAsOCw1MiwgMCw5LDIwMSxcbiAgICA4MSw3LDEzLCAwLDgsMTAwLCAwLDgsMzYsIDAsOSwxNjksXG4gICAgMCw4LDQsIDAsOCwxMzIsIDAsOCw2OCwgMCw5LDIzMyxcbiAgICA4MCw3LDgsIDAsOCw5MiwgMCw4LDI4LCAwLDksMTUzLFxuICAgIDg0LDcsODMsIDAsOCwxMjQsIDAsOCw2MCwgMCw5LDIxNyxcbiAgICA4Miw3LDIzLCAwLDgsMTA4LCAwLDgsNDQsIDAsOSwxODUsXG4gICAgMCw4LDEyLCAwLDgsMTQwLCAwLDgsNzYsIDAsOSwyNDksXG4gICAgODAsNywzLCAwLDgsODIsIDAsOCwxOCwgODUsOCwxNjMsXG4gICAgODMsNywzNSwgMCw4LDExNCwgMCw4LDUwLCAwLDksMTk3LFxuICAgIDgxLDcsMTEsIDAsOCw5OCwgMCw4LDM0LCAwLDksMTY1LFxuICAgIDAsOCwyLCAwLDgsMTMwLCAwLDgsNjYsIDAsOSwyMjksXG4gICAgODAsNyw3LCAwLDgsOTAsIDAsOCwyNiwgMCw5LDE0OSxcbiAgICA4NCw3LDY3LCAwLDgsMTIyLCAwLDgsNTgsIDAsOSwyMTMsXG4gICAgODIsNywxOSwgMCw4LDEwNiwgMCw4LDQyLCAwLDksMTgxLFxuICAgIDAsOCwxMCwgMCw4LDEzOCwgMCw4LDc0LCAwLDksMjQ1LFxuICAgIDgwLDcsNSwgMCw4LDg2LCAwLDgsMjIsIDE5Miw4LDAsXG4gICAgODMsNyw1MSwgMCw4LDExOCwgMCw4LDU0LCAwLDksMjA1LFxuICAgIDgxLDcsMTUsIDAsOCwxMDIsIDAsOCwzOCwgMCw5LDE3MyxcbiAgICAwLDgsNiwgMCw4LDEzNCwgMCw4LDcwLCAwLDksMjM3LFxuICAgIDgwLDcsOSwgMCw4LDk0LCAwLDgsMzAsIDAsOSwxNTcsXG4gICAgODQsNyw5OSwgMCw4LDEyNiwgMCw4LDYyLCAwLDksMjIxLFxuICAgIDgyLDcsMjcsIDAsOCwxMTAsIDAsOCw0NiwgMCw5LDE4OSxcbiAgICAwLDgsMTQsIDAsOCwxNDIsIDAsOCw3OCwgMCw5LDI1MyxcbiAgICA5Niw3LDI1NiwgMCw4LDgxLCAwLDgsMTcsIDg1LDgsMTMxLFxuICAgIDgyLDcsMzEsIDAsOCwxMTMsIDAsOCw0OSwgMCw5LDE5NSxcbiAgICA4MCw3LDEwLCAwLDgsOTcsIDAsOCwzMywgMCw5LDE2MyxcbiAgICAwLDgsMSwgMCw4LDEyOSwgMCw4LDY1LCAwLDksMjI3LFxuICAgIDgwLDcsNiwgMCw4LDg5LCAwLDgsMjUsIDAsOSwxNDcsXG4gICAgODMsNyw1OSwgMCw4LDEyMSwgMCw4LDU3LCAwLDksMjExLFxuICAgIDgxLDcsMTcsIDAsOCwxMDUsIDAsOCw0MSwgMCw5LDE3OSxcbiAgICAwLDgsOSwgMCw4LDEzNywgMCw4LDczLCAwLDksMjQzLFxuICAgIDgwLDcsNCwgMCw4LDg1LCAwLDgsMjEsIDgwLDgsMjU4LFxuICAgIDgzLDcsNDMsIDAsOCwxMTcsIDAsOCw1MywgMCw5LDIwMyxcbiAgICA4MSw3LDEzLCAwLDgsMTAxLCAwLDgsMzcsIDAsOSwxNzEsXG4gICAgMCw4LDUsIDAsOCwxMzMsIDAsOCw2OSwgMCw5LDIzNSxcbiAgICA4MCw3LDgsIDAsOCw5MywgMCw4LDI5LCAwLDksMTU1LFxuICAgIDg0LDcsODMsIDAsOCwxMjUsIDAsOCw2MSwgMCw5LDIxOSxcbiAgICA4Miw3LDIzLCAwLDgsMTA5LCAwLDgsNDUsIDAsOSwxODcsXG4gICAgMCw4LDEzLCAwLDgsMTQxLCAwLDgsNzcsIDAsOSwyNTEsXG4gICAgODAsNywzLCAwLDgsODMsIDAsOCwxOSwgODUsOCwxOTUsXG4gICAgODMsNywzNSwgMCw4LDExNSwgMCw4LDUxLCAwLDksMTk5LFxuICAgIDgxLDcsMTEsIDAsOCw5OSwgMCw4LDM1LCAwLDksMTY3LFxuICAgIDAsOCwzLCAwLDgsMTMxLCAwLDgsNjcsIDAsOSwyMzEsXG4gICAgODAsNyw3LCAwLDgsOTEsIDAsOCwyNywgMCw5LDE1MSxcbiAgICA4NCw3LDY3LCAwLDgsMTIzLCAwLDgsNTksIDAsOSwyMTUsXG4gICAgODIsNywxOSwgMCw4LDEwNywgMCw4LDQzLCAwLDksMTgzLFxuICAgIDAsOCwxMSwgMCw4LDEzOSwgMCw4LDc1LCAwLDksMjQ3LFxuICAgIDgwLDcsNSwgMCw4LDg3LCAwLDgsMjMsIDE5Miw4LDAsXG4gICAgODMsNyw1MSwgMCw4LDExOSwgMCw4LDU1LCAwLDksMjA3LFxuICAgIDgxLDcsMTUsIDAsOCwxMDMsIDAsOCwzOSwgMCw5LDE3NSxcbiAgICAwLDgsNywgMCw4LDEzNSwgMCw4LDcxLCAwLDksMjM5LFxuICAgIDgwLDcsOSwgMCw4LDk1LCAwLDgsMzEsIDAsOSwxNTksXG4gICAgODQsNyw5OSwgMCw4LDEyNywgMCw4LDYzLCAwLDksMjIzLFxuICAgIDgyLDcsMjcsIDAsOCwxMTEsIDAsOCw0NywgMCw5LDE5MSxcbiAgICAwLDgsMTUsIDAsOCwxNDMsIDAsOCw3OSwgMCw5LDI1NVxuXTtcbnZhciBmaXhlZF90ZCA9IFtcbiAgICA4MCw1LDEsIDg3LDUsMjU3LCA4Myw1LDE3LCA5MSw1LDQwOTcsXG4gICAgODEsNSw1LCA4OSw1LDEwMjUsIDg1LDUsNjUsIDkzLDUsMTYzODUsXG4gICAgODAsNSwzLCA4OCw1LDUxMywgODQsNSwzMywgOTIsNSw4MTkzLFxuICAgIDgyLDUsOSwgOTAsNSwyMDQ5LCA4Niw1LDEyOSwgMTkyLDUsMjQ1NzcsXG4gICAgODAsNSwyLCA4Nyw1LDM4NSwgODMsNSwyNSwgOTEsNSw2MTQ1LFxuICAgIDgxLDUsNywgODksNSwxNTM3LCA4NSw1LDk3LCA5Myw1LDI0NTc3LFxuICAgIDgwLDUsNCwgODgsNSw3NjksIDg0LDUsNDksIDkyLDUsMTIyODksXG4gICAgODIsNSwxMywgOTAsNSwzMDczLCA4Niw1LDE5MywgMTkyLDUsMjQ1Nzdcbl07XG5cbiAgLy8gVGFibGVzIGZvciBkZWZsYXRlIGZyb20gUEtaSVAncyBhcHBub3RlLnR4dC5cbiAgdmFyIGNwbGVucyA9IFsgLy8gQ29weSBsZW5ndGhzIGZvciBsaXRlcmFsIGNvZGVzIDI1Ny4uMjg1XG4gICAgICAgIDMsIDQsIDUsIDYsIDcsIDgsIDksIDEwLCAxMSwgMTMsIDE1LCAxNywgMTksIDIzLCAyNywgMzEsXG4gICAgICAgIDM1LCA0MywgNTEsIDU5LCA2NywgODMsIDk5LCAxMTUsIDEzMSwgMTYzLCAxOTUsIDIyNywgMjU4LCAwLCAwXG4gIF07XG5cbiAgLy8gc2VlIG5vdGUgIzEzIGFib3ZlIGFib3V0IDI1OFxuICB2YXIgY3BsZXh0ID0gWyAvLyBFeHRyYSBiaXRzIGZvciBsaXRlcmFsIGNvZGVzIDI1Ny4uMjg1XG4gICAgICAgIDAsIDAsIDAsIDAsIDAsIDAsIDAsIDAsIDEsIDEsIDEsIDEsIDIsIDIsIDIsIDIsXG4gICAgICAgIDMsIDMsIDMsIDMsIDQsIDQsIDQsIDQsIDUsIDUsIDUsIDUsIDAsIDExMiwgMTEyICAvLyAxMTI9PWludmFsaWRcbiAgXTtcblxuIHZhciBjcGRpc3QgPSBbIC8vIENvcHkgb2Zmc2V0cyBmb3IgZGlzdGFuY2UgY29kZXMgMC4uMjlcbiAgICAgICAgMSwgMiwgMywgNCwgNSwgNywgOSwgMTMsIDE3LCAyNSwgMzMsIDQ5LCA2NSwgOTcsIDEyOSwgMTkzLFxuICAgICAgICAyNTcsIDM4NSwgNTEzLCA3NjksIDEwMjUsIDE1MzcsIDIwNDksIDMwNzMsIDQwOTcsIDYxNDUsXG4gICAgICAgIDgxOTMsIDEyMjg5LCAxNjM4NSwgMjQ1NzdcbiAgXTtcblxuICB2YXIgY3BkZXh0ID0gWyAvLyBFeHRyYSBiaXRzIGZvciBkaXN0YW5jZSBjb2Rlc1xuICAgICAgICAwLCAwLCAwLCAwLCAxLCAxLCAyLCAyLCAzLCAzLCA0LCA0LCA1LCA1LCA2LCA2LFxuICAgICAgICA3LCA3LCA4LCA4LCA5LCA5LCAxMCwgMTAsIDExLCAxMSxcbiAgICAgICAgMTIsIDEyLCAxMywgMTNdO1xuXG4vL1xuLy8gWlN0cmVhbS5qYXZhXG4vL1xuXG5mdW5jdGlvbiBaU3RyZWFtKCkge1xufVxuXG5cblpTdHJlYW0ucHJvdG90eXBlLmluZmxhdGVJbml0ID0gZnVuY3Rpb24odywgbm93cmFwKSB7XG4gICAgaWYgKCF3KSB7XG5cdHcgPSBERUZfV0JJVFM7XG4gICAgfVxuICAgIGlmIChub3dyYXApIHtcblx0bm93cmFwID0gZmFsc2U7XG4gICAgfVxuICAgIHRoaXMuaXN0YXRlID0gbmV3IEluZmxhdGUoKTtcbiAgICByZXR1cm4gdGhpcy5pc3RhdGUuaW5mbGF0ZUluaXQodGhpcywgbm93cmFwPy13OncpO1xufVxuXG5aU3RyZWFtLnByb3RvdHlwZS5pbmZsYXRlID0gZnVuY3Rpb24oZikge1xuICAgIGlmKHRoaXMuaXN0YXRlPT1udWxsKSByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgcmV0dXJuIHRoaXMuaXN0YXRlLmluZmxhdGUodGhpcywgZik7XG59XG5cblpTdHJlYW0ucHJvdG90eXBlLmluZmxhdGVFbmQgPSBmdW5jdGlvbigpe1xuICAgIGlmKHRoaXMuaXN0YXRlPT1udWxsKSByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgdmFyIHJldD1pc3RhdGUuaW5mbGF0ZUVuZCh0aGlzKTtcbiAgICB0aGlzLmlzdGF0ZSA9IG51bGw7XG4gICAgcmV0dXJuIHJldDtcbn1cblpTdHJlYW0ucHJvdG90eXBlLmluZmxhdGVTeW5jID0gZnVuY3Rpb24oKXtcbiAgICAvLyBpZihpc3RhdGUgPT0gbnVsbCkgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIHJldHVybiBpc3RhdGUuaW5mbGF0ZVN5bmModGhpcyk7XG59XG5aU3RyZWFtLnByb3RvdHlwZS5pbmZsYXRlU2V0RGljdGlvbmFyeSA9IGZ1bmN0aW9uKGRpY3Rpb25hcnksIGRpY3RMZW5ndGgpe1xuICAgIC8vIGlmKGlzdGF0ZSA9PSBudWxsKSByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgcmV0dXJuIGlzdGF0ZS5pbmZsYXRlU2V0RGljdGlvbmFyeSh0aGlzLCBkaWN0aW9uYXJ5LCBkaWN0TGVuZ3RoKTtcbn1cblxuLypcblxuICBwdWJsaWMgaW50IGRlZmxhdGVJbml0KGludCBsZXZlbCl7XG4gICAgcmV0dXJuIGRlZmxhdGVJbml0KGxldmVsLCBNQVhfV0JJVFMpO1xuICB9XG4gIHB1YmxpYyBpbnQgZGVmbGF0ZUluaXQoaW50IGxldmVsLCBib29sZWFuIG5vd3JhcCl7XG4gICAgcmV0dXJuIGRlZmxhdGVJbml0KGxldmVsLCBNQVhfV0JJVFMsIG5vd3JhcCk7XG4gIH1cbiAgcHVibGljIGludCBkZWZsYXRlSW5pdChpbnQgbGV2ZWwsIGludCBiaXRzKXtcbiAgICByZXR1cm4gZGVmbGF0ZUluaXQobGV2ZWwsIGJpdHMsIGZhbHNlKTtcbiAgfVxuICBwdWJsaWMgaW50IGRlZmxhdGVJbml0KGludCBsZXZlbCwgaW50IGJpdHMsIGJvb2xlYW4gbm93cmFwKXtcbiAgICBkc3RhdGU9bmV3IERlZmxhdGUoKTtcbiAgICByZXR1cm4gZHN0YXRlLmRlZmxhdGVJbml0KHRoaXMsIGxldmVsLCBub3dyYXA/LWJpdHM6Yml0cyk7XG4gIH1cbiAgcHVibGljIGludCBkZWZsYXRlKGludCBmbHVzaCl7XG4gICAgaWYoZHN0YXRlPT1udWxsKXtcbiAgICAgIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICB9XG4gICAgcmV0dXJuIGRzdGF0ZS5kZWZsYXRlKHRoaXMsIGZsdXNoKTtcbiAgfVxuICBwdWJsaWMgaW50IGRlZmxhdGVFbmQoKXtcbiAgICBpZihkc3RhdGU9PW51bGwpIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICBpbnQgcmV0PWRzdGF0ZS5kZWZsYXRlRW5kKCk7XG4gICAgZHN0YXRlPW51bGw7XG4gICAgcmV0dXJuIHJldDtcbiAgfVxuICBwdWJsaWMgaW50IGRlZmxhdGVQYXJhbXMoaW50IGxldmVsLCBpbnQgc3RyYXRlZ3kpe1xuICAgIGlmKGRzdGF0ZT09bnVsbCkgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIHJldHVybiBkc3RhdGUuZGVmbGF0ZVBhcmFtcyh0aGlzLCBsZXZlbCwgc3RyYXRlZ3kpO1xuICB9XG4gIHB1YmxpYyBpbnQgZGVmbGF0ZVNldERpY3Rpb25hcnkgKGJ5dGVbXSBkaWN0aW9uYXJ5LCBpbnQgZGljdExlbmd0aCl7XG4gICAgaWYoZHN0YXRlID09IG51bGwpXG4gICAgICByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgcmV0dXJuIGRzdGF0ZS5kZWZsYXRlU2V0RGljdGlvbmFyeSh0aGlzLCBkaWN0aW9uYXJ5LCBkaWN0TGVuZ3RoKTtcbiAgfVxuXG4qL1xuXG4vKlxuICAvLyBGbHVzaCBhcyBtdWNoIHBlbmRpbmcgb3V0cHV0IGFzIHBvc3NpYmxlLiBBbGwgZGVmbGF0ZSgpIG91dHB1dCBnb2VzXG4gIC8vIHRocm91Z2ggdGhpcyBmdW5jdGlvbiBzbyBzb21lIGFwcGxpY2F0aW9ucyBtYXkgd2lzaCB0byBtb2RpZnkgaXRcbiAgLy8gdG8gYXZvaWQgYWxsb2NhdGluZyBhIGxhcmdlIHN0cm0tPm5leHRfb3V0IGJ1ZmZlciBhbmQgY29weWluZyBpbnRvIGl0LlxuICAvLyAoU2VlIGFsc28gcmVhZF9idWYoKSkuXG4gIHZvaWQgZmx1c2hfcGVuZGluZygpe1xuICAgIGludCBsZW49ZHN0YXRlLnBlbmRpbmc7XG5cbiAgICBpZihsZW4+YXZhaWxfb3V0KSBsZW49YXZhaWxfb3V0O1xuICAgIGlmKGxlbj09MCkgcmV0dXJuO1xuXG4gICAgaWYoZHN0YXRlLnBlbmRpbmdfYnVmLmxlbmd0aDw9ZHN0YXRlLnBlbmRpbmdfb3V0IHx8XG4gICAgICAgbmV4dF9vdXQubGVuZ3RoPD1uZXh0X291dF9pbmRleCB8fFxuICAgICAgIGRzdGF0ZS5wZW5kaW5nX2J1Zi5sZW5ndGg8KGRzdGF0ZS5wZW5kaW5nX291dCtsZW4pIHx8XG4gICAgICAgbmV4dF9vdXQubGVuZ3RoPChuZXh0X291dF9pbmRleCtsZW4pKXtcbiAgICAgIFN5c3RlbS5vdXQucHJpbnRsbihkc3RhdGUucGVuZGluZ19idWYubGVuZ3RoK1wiLCBcIitkc3RhdGUucGVuZGluZ19vdXQrXG5cdFx0XHQgXCIsIFwiK25leHRfb3V0Lmxlbmd0aCtcIiwgXCIrbmV4dF9vdXRfaW5kZXgrXCIsIFwiK2xlbik7XG4gICAgICBTeXN0ZW0ub3V0LnByaW50bG4oXCJhdmFpbF9vdXQ9XCIrYXZhaWxfb3V0KTtcbiAgICB9XG5cbiAgICBTeXN0ZW0uYXJyYXljb3B5KGRzdGF0ZS5wZW5kaW5nX2J1ZiwgZHN0YXRlLnBlbmRpbmdfb3V0LFxuXHRcdCAgICAgbmV4dF9vdXQsIG5leHRfb3V0X2luZGV4LCBsZW4pO1xuXG4gICAgbmV4dF9vdXRfaW5kZXgrPWxlbjtcbiAgICBkc3RhdGUucGVuZGluZ19vdXQrPWxlbjtcbiAgICB0b3RhbF9vdXQrPWxlbjtcbiAgICBhdmFpbF9vdXQtPWxlbjtcbiAgICBkc3RhdGUucGVuZGluZy09bGVuO1xuICAgIGlmKGRzdGF0ZS5wZW5kaW5nPT0wKXtcbiAgICAgIGRzdGF0ZS5wZW5kaW5nX291dD0wO1xuICAgIH1cbiAgfVxuXG4gIC8vIFJlYWQgYSBuZXcgYnVmZmVyIGZyb20gdGhlIGN1cnJlbnQgaW5wdXQgc3RyZWFtLCB1cGRhdGUgdGhlIGFkbGVyMzJcbiAgLy8gYW5kIHRvdGFsIG51bWJlciBvZiBieXRlcyByZWFkLiAgQWxsIGRlZmxhdGUoKSBpbnB1dCBnb2VzIHRocm91Z2hcbiAgLy8gdGhpcyBmdW5jdGlvbiBzbyBzb21lIGFwcGxpY2F0aW9ucyBtYXkgd2lzaCB0byBtb2RpZnkgaXQgdG8gYXZvaWRcbiAgLy8gYWxsb2NhdGluZyBhIGxhcmdlIHN0cm0tPm5leHRfaW4gYnVmZmVyIGFuZCBjb3B5aW5nIGZyb20gaXQuXG4gIC8vIChTZWUgYWxzbyBmbHVzaF9wZW5kaW5nKCkpLlxuICBpbnQgcmVhZF9idWYoYnl0ZVtdIGJ1ZiwgaW50IHN0YXJ0LCBpbnQgc2l6ZSkge1xuICAgIGludCBsZW49YXZhaWxfaW47XG5cbiAgICBpZihsZW4+c2l6ZSkgbGVuPXNpemU7XG4gICAgaWYobGVuPT0wKSByZXR1cm4gMDtcblxuICAgIGF2YWlsX2luLT1sZW47XG5cbiAgICBpZihkc3RhdGUubm9oZWFkZXI9PTApIHtcbiAgICAgIGFkbGVyPV9hZGxlci5hZGxlcjMyKGFkbGVyLCBuZXh0X2luLCBuZXh0X2luX2luZGV4LCBsZW4pO1xuICAgIH1cbiAgICBTeXN0ZW0uYXJyYXljb3B5KG5leHRfaW4sIG5leHRfaW5faW5kZXgsIGJ1Ziwgc3RhcnQsIGxlbik7XG4gICAgbmV4dF9pbl9pbmRleCAgKz0gbGVuO1xuICAgIHRvdGFsX2luICs9IGxlbjtcbiAgICByZXR1cm4gbGVuO1xuICB9XG5cbiAgcHVibGljIHZvaWQgZnJlZSgpe1xuICAgIG5leHRfaW49bnVsbDtcbiAgICBuZXh0X291dD1udWxsO1xuICAgIG1zZz1udWxsO1xuICAgIF9hZGxlcj1udWxsO1xuICB9XG59XG4qL1xuXG5cbi8vXG4vLyBJbmZsYXRlLmphdmFcbi8vXG5cbmZ1bmN0aW9uIEluZmxhdGUoKSB7XG4gICAgdGhpcy53YXMgPSBbMF07XG59XG5cbkluZmxhdGUucHJvdG90eXBlLmluZmxhdGVSZXNldCA9IGZ1bmN0aW9uKHopIHtcbiAgICBpZih6ID09IG51bGwgfHwgei5pc3RhdGUgPT0gbnVsbCkgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIFxuICAgIHoudG90YWxfaW4gPSB6LnRvdGFsX291dCA9IDA7XG4gICAgei5tc2cgPSBudWxsO1xuICAgIHouaXN0YXRlLm1vZGUgPSB6LmlzdGF0ZS5ub3dyYXAhPTAgPyBCTE9DS1MgOiBNRVRIT0Q7XG4gICAgei5pc3RhdGUuYmxvY2tzLnJlc2V0KHosIG51bGwpO1xuICAgIHJldHVybiBaX09LO1xufVxuXG5JbmZsYXRlLnByb3RvdHlwZS5pbmZsYXRlRW5kID0gZnVuY3Rpb24oeil7XG4gICAgaWYodGhpcy5ibG9ja3MgIT0gbnVsbClcbiAgICAgIHRoaXMuYmxvY2tzLmZyZWUoeik7XG4gICAgdGhpcy5ibG9ja3M9bnVsbDtcbiAgICByZXR1cm4gWl9PSztcbn1cblxuSW5mbGF0ZS5wcm90b3R5cGUuaW5mbGF0ZUluaXQgPSBmdW5jdGlvbih6LCB3KXtcbiAgICB6Lm1zZyA9IG51bGw7XG4gICAgdGhpcy5ibG9ja3MgPSBudWxsO1xuXG4gICAgLy8gaGFuZGxlIHVuZG9jdW1lbnRlZCBub3dyYXAgb3B0aW9uIChubyB6bGliIGhlYWRlciBvciBjaGVjaylcbiAgICBub3dyYXAgPSAwO1xuICAgIGlmKHcgPCAwKXtcbiAgICAgIHcgPSAtIHc7XG4gICAgICBub3dyYXAgPSAxO1xuICAgIH1cblxuICAgIC8vIHNldCB3aW5kb3cgc2l6ZVxuICAgIGlmKHc8OCB8fHc+MTUpe1xuICAgICAgdGhpcy5pbmZsYXRlRW5kKHopO1xuICAgICAgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIH1cbiAgICB0aGlzLndiaXRzPXc7XG5cbiAgICB6LmlzdGF0ZS5ibG9ja3M9bmV3IEluZkJsb2Nrcyh6LCBcblx0XHRcdFx0ICB6LmlzdGF0ZS5ub3dyYXAhPTAgPyBudWxsIDogdGhpcyxcblx0XHRcdFx0ICAxPDx3KTtcblxuICAgIC8vIHJlc2V0IHN0YXRlXG4gICAgdGhpcy5pbmZsYXRlUmVzZXQoeik7XG4gICAgcmV0dXJuIFpfT0s7XG4gIH1cblxuSW5mbGF0ZS5wcm90b3R5cGUuaW5mbGF0ZSA9IGZ1bmN0aW9uKHosIGYpe1xuICAgIHZhciByLCBiO1xuXG4gICAgaWYoeiA9PSBudWxsIHx8IHouaXN0YXRlID09IG51bGwgfHwgei5uZXh0X2luID09IG51bGwpXG4gICAgICByZXR1cm4gWl9TVFJFQU1fRVJST1I7XG4gICAgZiA9IGYgPT0gWl9GSU5JU0ggPyBaX0JVRl9FUlJPUiA6IFpfT0s7XG4gICAgciA9IFpfQlVGX0VSUk9SO1xuICAgIHdoaWxlICh0cnVlKXtcbiAgICAgIHN3aXRjaCAoei5pc3RhdGUubW9kZSl7XG4gICAgICBjYXNlIE1FVEhPRDpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgaWYoKCh6LmlzdGF0ZS5tZXRob2QgPSB6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdKSYweGYpIT1aX0RFRkxBVEVEKXtcbiAgICAgICAgICB6LmlzdGF0ZS5tb2RlID0gQkFEO1xuICAgICAgICAgIHoubXNnPVwidW5rbm93biBjb21wcmVzc2lvbiBtZXRob2RcIjtcbiAgICAgICAgICB6LmlzdGF0ZS5tYXJrZXIgPSA1OyAgICAgICAvLyBjYW4ndCB0cnkgaW5mbGF0ZVN5bmNcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgfVxuICAgICAgICBpZigoei5pc3RhdGUubWV0aG9kPj40KSs4PnouaXN0YXRlLndiaXRzKXtcbiAgICAgICAgICB6LmlzdGF0ZS5tb2RlID0gQkFEO1xuICAgICAgICAgIHoubXNnPVwiaW52YWxpZCB3aW5kb3cgc2l6ZVwiO1xuICAgICAgICAgIHouaXN0YXRlLm1hcmtlciA9IDU7ICAgICAgIC8vIGNhbid0IHRyeSBpbmZsYXRlU3luY1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG4gICAgICAgIHouaXN0YXRlLm1vZGU9RkxBRztcbiAgICAgIGNhc2UgRkxBRzpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgYiA9ICh6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdKSYweGZmO1xuXG4gICAgICAgIGlmKCgoKHouaXN0YXRlLm1ldGhvZCA8PCA4KStiKSAlIDMxKSE9MCl7XG4gICAgICAgICAgei5pc3RhdGUubW9kZSA9IEJBRDtcbiAgICAgICAgICB6Lm1zZyA9IFwiaW5jb3JyZWN0IGhlYWRlciBjaGVja1wiO1xuICAgICAgICAgIHouaXN0YXRlLm1hcmtlciA9IDU7ICAgICAgIC8vIGNhbid0IHRyeSBpbmZsYXRlU3luY1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG5cbiAgICAgICAgaWYoKGImUFJFU0VUX0RJQ1QpPT0wKXtcbiAgICAgICAgICB6LmlzdGF0ZS5tb2RlID0gQkxPQ0tTO1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG4gICAgICAgIHouaXN0YXRlLm1vZGUgPSBESUNUNDtcbiAgICAgIGNhc2UgRElDVDQ6XG5cbiAgICAgICAgaWYoei5hdmFpbF9pbj09MClyZXR1cm4gcjtyPWY7XG5cbiAgICAgICAgei5hdmFpbF9pbi0tOyB6LnRvdGFsX2luKys7XG4gICAgICAgIHouaXN0YXRlLm5lZWQ9KCh6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdJjB4ZmYpPDwyNCkmMHhmZjAwMDAwMDtcbiAgICAgICAgei5pc3RhdGUubW9kZT1ESUNUMztcbiAgICAgIGNhc2UgRElDVDM6XG5cbiAgICAgICAgaWYoei5hdmFpbF9pbj09MClyZXR1cm4gcjtyPWY7XG5cbiAgICAgICAgei5hdmFpbF9pbi0tOyB6LnRvdGFsX2luKys7XG4gICAgICAgIHouaXN0YXRlLm5lZWQrPSgoei5uZXh0X2luW3oubmV4dF9pbl9pbmRleCsrXSYweGZmKTw8MTYpJjB4ZmYwMDAwO1xuICAgICAgICB6LmlzdGF0ZS5tb2RlPURJQ1QyO1xuICAgICAgY2FzZSBESUNUMjpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgei5pc3RhdGUubmVlZCs9KCh6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdJjB4ZmYpPDw4KSYweGZmMDA7XG4gICAgICAgIHouaXN0YXRlLm1vZGU9RElDVDE7XG4gICAgICBjYXNlIERJQ1QxOlxuXG4gICAgICAgIGlmKHouYXZhaWxfaW49PTApcmV0dXJuIHI7cj1mO1xuXG4gICAgICAgIHouYXZhaWxfaW4tLTsgei50b3RhbF9pbisrO1xuICAgICAgICB6LmlzdGF0ZS5uZWVkICs9ICh6Lm5leHRfaW5bei5uZXh0X2luX2luZGV4KytdJjB4ZmYpO1xuICAgICAgICB6LmFkbGVyID0gei5pc3RhdGUubmVlZDtcbiAgICAgICAgei5pc3RhdGUubW9kZSA9IERJQ1QwO1xuICAgICAgICByZXR1cm4gWl9ORUVEX0RJQ1Q7XG4gICAgICBjYXNlIERJQ1QwOlxuICAgICAgICB6LmlzdGF0ZS5tb2RlID0gQkFEO1xuICAgICAgICB6Lm1zZyA9IFwibmVlZCBkaWN0aW9uYXJ5XCI7XG4gICAgICAgIHouaXN0YXRlLm1hcmtlciA9IDA7ICAgICAgIC8vIGNhbiB0cnkgaW5mbGF0ZVN5bmNcbiAgICAgICAgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgICAgY2FzZSBCTE9DS1M6XG5cbiAgICAgICAgciA9IHouaXN0YXRlLmJsb2Nrcy5wcm9jKHosIHIpO1xuICAgICAgICBpZihyID09IFpfREFUQV9FUlJPUil7XG4gICAgICAgICAgei5pc3RhdGUubW9kZSA9IEJBRDtcbiAgICAgICAgICB6LmlzdGF0ZS5tYXJrZXIgPSAwOyAgICAgLy8gY2FuIHRyeSBpbmZsYXRlU3luY1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG4gICAgICAgIGlmKHIgPT0gWl9PSyl7XG4gICAgICAgICAgciA9IGY7XG4gICAgICAgIH1cbiAgICAgICAgaWYociAhPSBaX1NUUkVBTV9FTkQpe1xuICAgICAgICAgIHJldHVybiByO1xuICAgICAgICB9XG4gICAgICAgIHIgPSBmO1xuICAgICAgICB6LmlzdGF0ZS5ibG9ja3MucmVzZXQoeiwgei5pc3RhdGUud2FzKTtcbiAgICAgICAgaWYoei5pc3RhdGUubm93cmFwIT0wKXtcbiAgICAgICAgICB6LmlzdGF0ZS5tb2RlPURPTkU7XG4gICAgICAgICAgYnJlYWs7XG4gICAgICAgIH1cbiAgICAgICAgei5pc3RhdGUubW9kZT1DSEVDSzQ7XG4gICAgICBjYXNlIENIRUNLNDpcblxuICAgICAgICBpZih6LmF2YWlsX2luPT0wKXJldHVybiByO3I9ZjtcblxuICAgICAgICB6LmF2YWlsX2luLS07IHoudG90YWxfaW4rKztcbiAgICAgICAgei5pc3RhdGUubmVlZD0oKHoubmV4dF9pblt6Lm5leHRfaW5faW5kZXgrK10mMHhmZik8PDI0KSYweGZmMDAwMDAwO1xuICAgICAgICB6LmlzdGF0ZS5tb2RlPUNIRUNLMztcbiAgICAgIGNhc2UgQ0hFQ0szOlxuXG4gICAgICAgIGlmKHouYXZhaWxfaW49PTApcmV0dXJuIHI7cj1mO1xuXG4gICAgICAgIHouYXZhaWxfaW4tLTsgei50b3RhbF9pbisrO1xuICAgICAgICB6LmlzdGF0ZS5uZWVkKz0oKHoubmV4dF9pblt6Lm5leHRfaW5faW5kZXgrK10mMHhmZik8PDE2KSYweGZmMDAwMDtcbiAgICAgICAgei5pc3RhdGUubW9kZSA9IENIRUNLMjtcbiAgICAgIGNhc2UgQ0hFQ0syOlxuXG4gICAgICAgIGlmKHouYXZhaWxfaW49PTApcmV0dXJuIHI7cj1mO1xuXG4gICAgICAgIHouYXZhaWxfaW4tLTsgei50b3RhbF9pbisrO1xuICAgICAgICB6LmlzdGF0ZS5uZWVkKz0oKHoubmV4dF9pblt6Lm5leHRfaW5faW5kZXgrK10mMHhmZik8PDgpJjB4ZmYwMDtcbiAgICAgICAgei5pc3RhdGUubW9kZSA9IENIRUNLMTtcbiAgICAgIGNhc2UgQ0hFQ0sxOlxuXG4gICAgICAgIGlmKHouYXZhaWxfaW49PTApcmV0dXJuIHI7cj1mO1xuXG4gICAgICAgIHouYXZhaWxfaW4tLTsgei50b3RhbF9pbisrO1xuICAgICAgICB6LmlzdGF0ZS5uZWVkKz0oei5uZXh0X2luW3oubmV4dF9pbl9pbmRleCsrXSYweGZmKTtcblxuICAgICAgICBpZigoKHouaXN0YXRlLndhc1swXSkpICE9ICgoei5pc3RhdGUubmVlZCkpKXtcbiAgICAgICAgICB6LmlzdGF0ZS5tb2RlID0gQkFEO1xuICAgICAgICAgIHoubXNnID0gXCJpbmNvcnJlY3QgZGF0YSBjaGVja1wiO1xuICAgICAgICAgIHouaXN0YXRlLm1hcmtlciA9IDU7ICAgICAgIC8vIGNhbid0IHRyeSBpbmZsYXRlU3luY1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICB9XG5cbiAgICAgICAgei5pc3RhdGUubW9kZSA9IERPTkU7XG4gICAgICBjYXNlIERPTkU6XG4gICAgICAgIHJldHVybiBaX1NUUkVBTV9FTkQ7XG4gICAgICBjYXNlIEJBRDpcbiAgICAgICAgcmV0dXJuIFpfREFUQV9FUlJPUjtcbiAgICAgIGRlZmF1bHQ6XG4gICAgICAgIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcbiAgICAgIH1cbiAgICB9XG4gIH1cblxuXG5JbmZsYXRlLnByb3RvdHlwZS5pbmZsYXRlU2V0RGljdGlvbmFyeSA9IGZ1bmN0aW9uKHosICBkaWN0aW9uYXJ5LCBkaWN0TGVuZ3RoKSB7XG4gICAgdmFyIGluZGV4PTA7XG4gICAgdmFyIGxlbmd0aCA9IGRpY3RMZW5ndGg7XG4gICAgaWYoej09bnVsbCB8fCB6LmlzdGF0ZSA9PSBudWxsfHwgei5pc3RhdGUubW9kZSAhPSBESUNUMClcbiAgICAgIHJldHVybiBaX1NUUkVBTV9FUlJPUjtcblxuICAgIGlmKHouX2FkbGVyLmFkbGVyMzIoMSwgZGljdGlvbmFyeSwgMCwgZGljdExlbmd0aCkhPXouYWRsZXIpe1xuICAgICAgcmV0dXJuIFpfREFUQV9FUlJPUjtcbiAgICB9XG5cbiAgICB6LmFkbGVyID0gei5fYWRsZXIuYWRsZXIzMigwLCBudWxsLCAwLCAwKTtcblxuICAgIGlmKGxlbmd0aCA+PSAoMTw8ei5pc3RhdGUud2JpdHMpKXtcbiAgICAgIGxlbmd0aCA9ICgxPDx6LmlzdGF0ZS53Yml0cyktMTtcbiAgICAgIGluZGV4PWRpY3RMZW5ndGggLSBsZW5ndGg7XG4gICAgfVxuICAgIHouaXN0YXRlLmJsb2Nrcy5zZXRfZGljdGlvbmFyeShkaWN0aW9uYXJ5LCBpbmRleCwgbGVuZ3RoKTtcbiAgICB6LmlzdGF0ZS5tb2RlID0gQkxPQ0tTO1xuICAgIHJldHVybiBaX09LO1xuICB9XG5cbi8vICBzdGF0aWMgcHJpdmF0ZSBieXRlW10gbWFyayA9IHsoYnl0ZSkwLCAoYnl0ZSkwLCAoYnl0ZSkweGZmLCAoYnl0ZSkweGZmfTtcbnZhciBtYXJrID0gWzAsIDAsIDI1NSwgMjU1XVxuXG5JbmZsYXRlLnByb3RvdHlwZS5pbmZsYXRlU3luYyA9IGZ1bmN0aW9uKHope1xuICAgIHZhciBuOyAgICAgICAvLyBudW1iZXIgb2YgYnl0ZXMgdG8gbG9vayBhdFxuICAgIHZhciBwOyAgICAgICAvLyBwb2ludGVyIHRvIGJ5dGVzXG4gICAgdmFyIG07ICAgICAgIC8vIG51bWJlciBvZiBtYXJrZXIgYnl0ZXMgZm91bmQgaW4gYSByb3dcbiAgICB2YXIgciwgdzsgICAvLyB0ZW1wb3JhcmllcyB0byBzYXZlIHRvdGFsX2luIGFuZCB0b3RhbF9vdXRcblxuICAgIC8vIHNldCB1cFxuICAgIGlmKHogPT0gbnVsbCB8fCB6LmlzdGF0ZSA9PSBudWxsKVxuICAgICAgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIGlmKHouaXN0YXRlLm1vZGUgIT0gQkFEKXtcbiAgICAgIHouaXN0YXRlLm1vZGUgPSBCQUQ7XG4gICAgICB6LmlzdGF0ZS5tYXJrZXIgPSAwO1xuICAgIH1cbiAgICBpZigobj16LmF2YWlsX2luKT09MClcbiAgICAgIHJldHVybiBaX0JVRl9FUlJPUjtcbiAgICBwPXoubmV4dF9pbl9pbmRleDtcbiAgICBtPXouaXN0YXRlLm1hcmtlcjtcblxuICAgIC8vIHNlYXJjaFxuICAgIHdoaWxlIChuIT0wICYmIG0gPCA0KXtcbiAgICAgIGlmKHoubmV4dF9pbltwXSA9PSBtYXJrW21dKXtcbiAgICAgICAgbSsrO1xuICAgICAgfVxuICAgICAgZWxzZSBpZih6Lm5leHRfaW5bcF0hPTApe1xuICAgICAgICBtID0gMDtcbiAgICAgIH1cbiAgICAgIGVsc2V7XG4gICAgICAgIG0gPSA0IC0gbTtcbiAgICAgIH1cbiAgICAgIHArKzsgbi0tO1xuICAgIH1cblxuICAgIC8vIHJlc3RvcmVcbiAgICB6LnRvdGFsX2luICs9IHAtei5uZXh0X2luX2luZGV4O1xuICAgIHoubmV4dF9pbl9pbmRleCA9IHA7XG4gICAgei5hdmFpbF9pbiA9IG47XG4gICAgei5pc3RhdGUubWFya2VyID0gbTtcblxuICAgIC8vIHJldHVybiBubyBqb3kgb3Igc2V0IHVwIHRvIHJlc3RhcnQgb24gYSBuZXcgYmxvY2tcbiAgICBpZihtICE9IDQpe1xuICAgICAgcmV0dXJuIFpfREFUQV9FUlJPUjtcbiAgICB9XG4gICAgcj16LnRvdGFsX2luOyAgdz16LnRvdGFsX291dDtcbiAgICB0aGlzLmluZmxhdGVSZXNldCh6KTtcbiAgICB6LnRvdGFsX2luPXI7ICB6LnRvdGFsX291dCA9IHc7XG4gICAgei5pc3RhdGUubW9kZSA9IEJMT0NLUztcbiAgICByZXR1cm4gWl9PSztcbn1cblxuICAvLyBSZXR1cm5zIHRydWUgaWYgaW5mbGF0ZSBpcyBjdXJyZW50bHkgYXQgdGhlIGVuZCBvZiBhIGJsb2NrIGdlbmVyYXRlZFxuICAvLyBieSBaX1NZTkNfRkxVU0ggb3IgWl9GVUxMX0ZMVVNILiBUaGlzIGZ1bmN0aW9uIGlzIHVzZWQgYnkgb25lIFBQUFxuICAvLyBpbXBsZW1lbnRhdGlvbiB0byBwcm92aWRlIGFuIGFkZGl0aW9uYWwgc2FmZXR5IGNoZWNrLiBQUFAgdXNlcyBaX1NZTkNfRkxVU0hcbiAgLy8gYnV0IHJlbW92ZXMgdGhlIGxlbmd0aCBieXRlcyBvZiB0aGUgcmVzdWx0aW5nIGVtcHR5IHN0b3JlZCBibG9jay4gV2hlblxuICAvLyBkZWNvbXByZXNzaW5nLCBQUFAgY2hlY2tzIHRoYXQgYXQgdGhlIGVuZCBvZiBpbnB1dCBwYWNrZXQsIGluZmxhdGUgaXNcbiAgLy8gd2FpdGluZyBmb3IgdGhlc2UgbGVuZ3RoIGJ5dGVzLlxuSW5mbGF0ZS5wcm90b3R5cGUuaW5mbGF0ZVN5bmNQb2ludCA9IGZ1bmN0aW9uKHope1xuICAgIGlmKHogPT0gbnVsbCB8fCB6LmlzdGF0ZSA9PSBudWxsIHx8IHouaXN0YXRlLmJsb2NrcyA9PSBudWxsKVxuICAgICAgcmV0dXJuIFpfU1RSRUFNX0VSUk9SO1xuICAgIHJldHVybiB6LmlzdGF0ZS5ibG9ja3Muc3luY19wb2ludCgpO1xufVxuXG5cbi8vXG4vLyBJbmZCbG9ja3MuamF2YVxuLy9cblxudmFyIElORkJMT0NLU19CT1JERVIgPSBbMTYsIDE3LCAxOCwgMCwgOCwgNywgOSwgNiwgMTAsIDUsIDExLCA0LCAxMiwgMywgMTMsIDIsIDE0LCAxLCAxNV07XG5cbmZ1bmN0aW9uIEluZkJsb2Nrcyh6LCBjaGVja2ZuLCB3KSB7XG4gICAgdGhpcy5odWZ0cz1uZXcgSW50MzJBcnJheShNQU5ZKjMpO1xuICAgIHRoaXMud2luZG93PW5ldyBVaW50OEFycmF5KHcpO1xuICAgIHRoaXMuZW5kPXc7XG4gICAgdGhpcy5jaGVja2ZuID0gY2hlY2tmbjtcbiAgICB0aGlzLm1vZGUgPSBJQl9UWVBFO1xuICAgIHRoaXMucmVzZXQoeiwgbnVsbCk7XG5cbiAgICB0aGlzLmxlZnQgPSAwOyAgICAgICAgICAgIC8vIGlmIFNUT1JFRCwgYnl0ZXMgbGVmdCB0byBjb3B5IFxuXG4gICAgdGhpcy50YWJsZSA9IDA7ICAgICAgICAgICAvLyB0YWJsZSBsZW5ndGhzICgxNCBiaXRzKSBcbiAgICB0aGlzLmluZGV4ID0gMDsgICAgICAgICAgIC8vIGluZGV4IGludG8gYmxlbnMgKG9yIGJvcmRlcikgXG4gICAgdGhpcy5ibGVucyA9IG51bGw7ICAgICAgICAgLy8gYml0IGxlbmd0aHMgb2YgY29kZXMgXG4gICAgdGhpcy5iYj1uZXcgSW50MzJBcnJheSgxKTsgLy8gYml0IGxlbmd0aCB0cmVlIGRlcHRoIFxuICAgIHRoaXMudGI9bmV3IEludDMyQXJyYXkoMSk7IC8vIGJpdCBsZW5ndGggZGVjb2RpbmcgdHJlZSBcblxuICAgIHRoaXMuY29kZXMgPSBuZXcgSW5mQ29kZXMoKTtcblxuICAgIHRoaXMubGFzdCA9IDA7ICAgICAgICAgICAgLy8gdHJ1ZSBpZiB0aGlzIGJsb2NrIGlzIHRoZSBsYXN0IGJsb2NrIFxuXG4gIC8vIG1vZGUgaW5kZXBlbmRlbnQgaW5mb3JtYXRpb24gXG4gICAgdGhpcy5iaXRrID0gMDsgICAgICAgICAgICAvLyBiaXRzIGluIGJpdCBidWZmZXIgXG4gICAgdGhpcy5iaXRiID0gMDsgICAgICAgICAgICAvLyBiaXQgYnVmZmVyIFxuICAgIHRoaXMucmVhZCA9IDA7ICAgICAgICAgICAgLy8gd2luZG93IHJlYWQgcG9pbnRlciBcbiAgICB0aGlzLndyaXRlID0gMDsgICAgICAgICAgIC8vIHdpbmRvdyB3cml0ZSBwb2ludGVyIFxuICAgIHRoaXMuY2hlY2sgPSAwOyAgICAgICAgICAvLyBjaGVjayBvbiBvdXRwdXQgXG5cbiAgICB0aGlzLmluZnRyZWU9bmV3IEluZlRyZWUoKTtcbn1cblxuXG5cblxuSW5mQmxvY2tzLnByb3RvdHlwZS5yZXNldCA9IGZ1bmN0aW9uKHosIGMpe1xuICAgIGlmKGMpIGNbMF09dGhpcy5jaGVjaztcbiAgICBpZih0aGlzLm1vZGU9PUlCX0NPREVTKXtcbiAgICAgIHRoaXMuY29kZXMuZnJlZSh6KTtcbiAgICB9XG4gICAgdGhpcy5tb2RlPUlCX1RZUEU7XG4gICAgdGhpcy5iaXRrPTA7XG4gICAgdGhpcy5iaXRiPTA7XG4gICAgdGhpcy5yZWFkPXRoaXMud3JpdGU9MDtcblxuICAgIGlmKHRoaXMuY2hlY2tmbilcbiAgICAgIHouYWRsZXI9dGhpcy5jaGVjaz16Ll9hZGxlci5hZGxlcjMyKDAsIG51bGwsIDAsIDApO1xuICB9XG5cbiBJbmZCbG9ja3MucHJvdG90eXBlLnByb2MgPSBmdW5jdGlvbih6LCByKXtcbiAgICB2YXIgdDsgICAgICAgICAgICAgIC8vIHRlbXBvcmFyeSBzdG9yYWdlXG4gICAgdmFyIGI7ICAgICAgICAgICAgICAvLyBiaXQgYnVmZmVyXG4gICAgdmFyIGs7ICAgICAgICAgICAgICAvLyBiaXRzIGluIGJpdCBidWZmZXJcbiAgICB2YXIgcDsgICAgICAgICAgICAgIC8vIGlucHV0IGRhdGEgcG9pbnRlclxuICAgIHZhciBuOyAgICAgICAgICAgICAgLy8gYnl0ZXMgYXZhaWxhYmxlIHRoZXJlXG4gICAgdmFyIHE7ICAgICAgICAgICAgICAvLyBvdXRwdXQgd2luZG93IHdyaXRlIHBvaW50ZXJcbiAgICB2YXIgbTsgICAgICAgICAgICAgIC8vIGJ5dGVzIHRvIGVuZCBvZiB3aW5kb3cgb3IgcmVhZCBwb2ludGVyXG5cbiAgICAvLyBjb3B5IGlucHV0L291dHB1dCBpbmZvcm1hdGlvbiB0byBsb2NhbHMgKFVQREFURSBtYWNybyByZXN0b3JlcylcbiAgICB7cD16Lm5leHRfaW5faW5kZXg7bj16LmF2YWlsX2luO2I9dGhpcy5iaXRiO2s9dGhpcy5iaXRrO31cbiAgICB7cT10aGlzLndyaXRlO209KHE8dGhpcy5yZWFkID8gdGhpcy5yZWFkLXEtMSA6IHRoaXMuZW5kLXEpO31cblxuICAgIC8vIHByb2Nlc3MgaW5wdXQgYmFzZWQgb24gY3VycmVudCBzdGF0ZVxuICAgIHdoaWxlKHRydWUpe1xuICAgICAgc3dpdGNoICh0aGlzLm1vZGUpe1xuICAgICAgY2FzZSBJQl9UWVBFOlxuXG5cdHdoaWxlKGs8KDMpKXtcblx0ICBpZihuIT0wKXtcblx0ICAgIHI9Wl9PSztcblx0ICB9XG5cdCAgZWxzZXtcblx0ICAgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgICB6LmF2YWlsX2luPW47XG5cdCAgICB6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgIHRoaXMud3JpdGU9cTtcblx0ICAgIHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICB9O1xuXHQgIG4tLTtcblx0ICBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7XG5cdCAgays9ODtcblx0fVxuXHR0ID0gKGIgJiA3KTtcblx0dGhpcy5sYXN0ID0gdCAmIDE7XG5cblx0c3dpdGNoICh0ID4+PiAxKXtcbiAgICAgICAgY2FzZSAwOiAgICAgICAgICAgICAgICAgICAgICAgICAvLyBzdG9yZWQgXG4gICAgICAgICAge2I+Pj49KDMpO2stPSgzKTt9XG4gICAgICAgICAgdCA9IGsgJiA3OyAgICAgICAgICAgICAgICAgICAgLy8gZ28gdG8gYnl0ZSBib3VuZGFyeVxuXG4gICAgICAgICAge2I+Pj49KHQpO2stPSh0KTt9XG4gICAgICAgICAgdGhpcy5tb2RlID0gSUJfTEVOUzsgICAgICAgICAgICAgICAgICAvLyBnZXQgbGVuZ3RoIG9mIHN0b3JlZCBibG9ja1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICBjYXNlIDE6ICAgICAgICAgICAgICAgICAgICAgICAgIC8vIGZpeGVkXG4gICAgICAgICAge1xuICAgICAgICAgICAgICB2YXIgYmw9bmV3IEludDMyQXJyYXkoMSk7XG5cdCAgICAgIHZhciBiZD1uZXcgSW50MzJBcnJheSgxKTtcbiAgICAgICAgICAgICAgdmFyIHRsPVtdO1xuXHQgICAgICB2YXIgdGQ9W107XG5cblx0ICAgICAgaW5mbGF0ZV90cmVlc19maXhlZChibCwgYmQsIHRsLCB0ZCwgeik7XG4gICAgICAgICAgICAgIHRoaXMuY29kZXMuaW5pdChibFswXSwgYmRbMF0sIHRsWzBdLCAwLCB0ZFswXSwgMCwgeik7XG4gICAgICAgICAgfVxuXG4gICAgICAgICAge2I+Pj49KDMpO2stPSgzKTt9XG5cbiAgICAgICAgICB0aGlzLm1vZGUgPSBJQl9DT0RFUztcbiAgICAgICAgICBicmVhaztcbiAgICAgICAgY2FzZSAyOiAgICAgICAgICAgICAgICAgICAgICAgICAvLyBkeW5hbWljXG5cbiAgICAgICAgICB7Yj4+Pj0oMyk7ay09KDMpO31cblxuICAgICAgICAgIHRoaXMubW9kZSA9IElCX1RBQkxFO1xuICAgICAgICAgIGJyZWFrO1xuICAgICAgICBjYXNlIDM6ICAgICAgICAgICAgICAgICAgICAgICAgIC8vIGlsbGVnYWxcblxuICAgICAgICAgIHtiPj4+PSgzKTtrLT0oMyk7fVxuICAgICAgICAgIHRoaXMubW9kZSA9IEJBRDtcbiAgICAgICAgICB6Lm1zZyA9IFwiaW52YWxpZCBibG9jayB0eXBlXCI7XG4gICAgICAgICAgciA9IFpfREFUQV9FUlJPUjtcblxuXHQgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgIHRoaXMud3JpdGU9cTtcblx0ICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdH1cblx0YnJlYWs7XG4gICAgICBjYXNlIElCX0xFTlM6XG5cdHdoaWxlKGs8KDMyKSl7XG5cdCAgaWYobiE9MCl7XG5cdCAgICByPVpfT0s7XG5cdCAgfVxuXHQgIGVsc2V7XG5cdCAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgei5hdmFpbF9pbj1uO1xuXHQgICAgei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICB0aGlzLndyaXRlPXE7XG5cdCAgICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgfTtcblx0ICBuLS07XG5cdCAgYnw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO1xuXHQgIGsrPTg7XG5cdH1cblxuXHRpZiAoKCgofmIpID4+PiAxNikgJiAweGZmZmYpICE9IChiICYgMHhmZmZmKSl7XG5cdCAgdGhpcy5tb2RlID0gQkFEO1xuXHQgIHoubXNnID0gXCJpbnZhbGlkIHN0b3JlZCBibG9jayBsZW5ndGhzXCI7XG5cdCAgciA9IFpfREFUQV9FUlJPUjtcblxuXHQgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgIHRoaXMud3JpdGU9cTtcblx0ICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdH1cblx0dGhpcy5sZWZ0ID0gKGIgJiAweGZmZmYpO1xuXHRiID0gayA9IDA7ICAgICAgICAgICAgICAgICAgICAgICAvLyBkdW1wIGJpdHNcblx0dGhpcy5tb2RlID0gdGhpcy5sZWZ0IT0wID8gSUJfU1RPUkVEIDogKHRoaXMubGFzdCE9MCA/IElCX0RSWSA6IElCX1RZUEUpO1xuXHRicmVhaztcbiAgICAgIGNhc2UgSUJfU1RPUkVEOlxuXHRpZiAobiA9PSAwKXtcblx0ICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICB3cml0ZT1xO1xuXHQgIHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeixyKTtcblx0fVxuXG5cdGlmKG09PTApe1xuXHQgIGlmKHE9PWVuZCYmcmVhZCE9MCl7XG5cdCAgICBxPTA7IG09KHE8dGhpcy5yZWFkID8gdGhpcy5yZWFkLXEtMSA6IHRoaXMuZW5kLXEpO1xuXHQgIH1cblx0ICBpZihtPT0wKXtcblx0ICAgIHRoaXMud3JpdGU9cTsgXG5cdCAgICByPXRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgcT10aGlzLndyaXRlOyBtID0gKHEgPCB0aGlzLnJlYWQgPyB0aGlzLnJlYWQtcS0xIDogdGhpcy5lbmQtcSk7XG5cdCAgICBpZihxPT10aGlzLmVuZCAmJiB0aGlzLnJlYWQgIT0gMCl7XG5cdCAgICAgIHE9MDsgbSA9IChxIDwgdGhpcy5yZWFkID8gdGhpcy5yZWFkLXEtMSA6IHRoaXMuZW5kLXEpO1xuXHQgICAgfVxuXHQgICAgaWYobT09MCl7XG5cdCAgICAgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgICAgdGhpcy53cml0ZT1xO1xuXHQgICAgICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgICB9XG5cdCAgfVxuXHR9XG5cdHI9Wl9PSztcblxuXHR0ID0gdGhpcy5sZWZ0O1xuXHRpZih0Pm4pIHQgPSBuO1xuXHRpZih0Pm0pIHQgPSBtO1xuXHRhcnJheUNvcHkoei5uZXh0X2luLCBwLCB0aGlzLndpbmRvdywgcSwgdCk7XG5cdHAgKz0gdDsgIG4gLT0gdDtcblx0cSArPSB0OyAgbSAtPSB0O1xuXHRpZiAoKHRoaXMubGVmdCAtPSB0KSAhPSAwKVxuXHQgIGJyZWFrO1xuXHR0aGlzLm1vZGUgPSAodGhpcy5sYXN0ICE9IDAgPyBJQl9EUlkgOiBJQl9UWVBFKTtcblx0YnJlYWs7XG4gICAgICBjYXNlIElCX1RBQkxFOlxuXG5cdHdoaWxlKGs8KDE0KSl7XG5cdCAgaWYobiE9MCl7XG5cdCAgICByPVpfT0s7XG5cdCAgfVxuXHQgIGVsc2V7XG5cdCAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgei5hdmFpbF9pbj1uO1xuXHQgICAgei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICB0aGlzLndyaXRlPXE7XG5cdCAgICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgfTtcblx0ICBuLS07XG5cdCAgYnw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO1xuXHQgIGsrPTg7XG5cdH1cblxuXHR0aGlzLnRhYmxlID0gdCA9IChiICYgMHgzZmZmKTtcblx0aWYgKCh0ICYgMHgxZikgPiAyOSB8fCAoKHQgPj4gNSkgJiAweDFmKSA+IDI5KVxuXHQgIHtcblx0ICAgIHRoaXMubW9kZSA9IElCX0JBRDtcblx0ICAgIHoubXNnID0gXCJ0b28gbWFueSBsZW5ndGggb3IgZGlzdGFuY2Ugc3ltYm9sc1wiO1xuXHQgICAgciA9IFpfREFUQV9FUlJPUjtcblxuXHQgICAgdGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgIHRoaXMud3JpdGU9cTtcblx0ICAgIHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICB9XG5cdHQgPSAyNTggKyAodCAmIDB4MWYpICsgKCh0ID4+IDUpICYgMHgxZik7XG5cdGlmKHRoaXMuYmxlbnM9PW51bGwgfHwgdGhpcy5ibGVucy5sZW5ndGg8dCl7XG5cdCAgICB0aGlzLmJsZW5zPW5ldyBJbnQzMkFycmF5KHQpO1xuXHR9XG5cdGVsc2V7XG5cdCAgZm9yKHZhciBpPTA7IGk8dDsgaSsrKXtcbiAgICAgICAgICAgICAgdGhpcy5ibGVuc1tpXT0wO1xuICAgICAgICAgIH1cblx0fVxuXG5cdHtiPj4+PSgxNCk7ay09KDE0KTt9XG5cblx0dGhpcy5pbmRleCA9IDA7XG5cdG1vZGUgPSBJQl9CVFJFRTtcbiAgICAgIGNhc2UgSUJfQlRSRUU6XG5cdHdoaWxlICh0aGlzLmluZGV4IDwgNCArICh0aGlzLnRhYmxlID4+PiAxMCkpe1xuXHQgIHdoaWxlKGs8KDMpKXtcblx0ICAgIGlmKG4hPTApe1xuXHQgICAgICByPVpfT0s7XG5cdCAgICB9XG5cdCAgICBlbHNle1xuXHQgICAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgICB6LmF2YWlsX2luPW47XG5cdCAgICAgIHoudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgICB0aGlzLndyaXRlPXE7XG5cdCAgICAgIHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICAgIH07XG5cdCAgICBuLS07XG5cdCAgICBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7XG5cdCAgICBrKz04O1xuXHQgIH1cblxuXHQgIHRoaXMuYmxlbnNbSU5GQkxPQ0tTX0JPUkRFUlt0aGlzLmluZGV4KytdXSA9IGImNztcblxuXHQgIHtiPj4+PSgzKTtrLT0oMyk7fVxuXHR9XG5cblx0d2hpbGUodGhpcy5pbmRleCA8IDE5KXtcblx0ICB0aGlzLmJsZW5zW0lORkJMT0NLU19CT1JERVJbdGhpcy5pbmRleCsrXV0gPSAwO1xuXHR9XG5cblx0dGhpcy5iYlswXSA9IDc7XG5cdHQgPSB0aGlzLmluZnRyZWUuaW5mbGF0ZV90cmVlc19iaXRzKHRoaXMuYmxlbnMsIHRoaXMuYmIsIHRoaXMudGIsIHRoaXMuaHVmdHMsIHopO1xuXHRpZiAodCAhPSBaX09LKXtcblx0ICByID0gdDtcblx0ICBpZiAociA9PSBaX0RBVEFfRVJST1Ipe1xuXHQgICAgdGhpcy5ibGVucz1udWxsO1xuXHQgICAgdGhpcy5tb2RlID0gSUJfQkFEO1xuXHQgIH1cblxuXHQgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgIHdyaXRlPXE7XG5cdCAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHR9XG5cblx0dGhpcy5pbmRleCA9IDA7XG5cdHRoaXMubW9kZSA9IElCX0RUUkVFO1xuICAgICAgY2FzZSBJQl9EVFJFRTpcblx0d2hpbGUgKHRydWUpe1xuXHQgIHQgPSB0aGlzLnRhYmxlO1xuXHQgIGlmKCEodGhpcy5pbmRleCA8IDI1OCArICh0ICYgMHgxZikgKyAoKHQgPj4gNSkgJiAweDFmKSkpe1xuXHQgICAgYnJlYWs7XG5cdCAgfVxuXG5cdCAgdmFyIGg7IC8vaW50W11cblx0ICB2YXIgaSwgaiwgYztcblxuXHQgIHQgPSB0aGlzLmJiWzBdO1xuXG5cdCAgd2hpbGUoazwodCkpe1xuXHQgICAgaWYobiE9MCl7XG5cdCAgICAgIHI9Wl9PSztcblx0ICAgIH1cblx0ICAgIGVsc2V7XG5cdCAgICAgIHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdCAgICAgIHouYXZhaWxfaW49bjtcblx0ICAgICAgei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICAgIHRoaXMud3JpdGU9cTtcblx0ICAgICAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgfTtcblx0ICAgIG4tLTtcblx0ICAgIGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztcblx0ICAgIGsrPTg7XG5cdCAgfVxuXG4vL1x0ICBpZiAodGhpcy50YlswXT09LTEpe1xuLy8gICAgICAgICAgICBkbG9nKFwibnVsbC4uLlwiKTtcbi8vXHQgIH1cblxuXHQgIHQ9dGhpcy5odWZ0c1sodGhpcy50YlswXSsoYiAmIGluZmxhdGVfbWFza1t0XSkpKjMrMV07XG5cdCAgYz10aGlzLmh1ZnRzWyh0aGlzLnRiWzBdKyhiICYgaW5mbGF0ZV9tYXNrW3RdKSkqMysyXTtcblxuXHQgIGlmIChjIDwgMTYpe1xuXHQgICAgYj4+Pj0odCk7ay09KHQpO1xuXHQgICAgdGhpcy5ibGVuc1t0aGlzLmluZGV4KytdID0gYztcblx0ICB9XG5cdCAgZWxzZSB7IC8vIGMgPT0gMTYuLjE4XG5cdCAgICBpID0gYyA9PSAxOCA/IDcgOiBjIC0gMTQ7XG5cdCAgICBqID0gYyA9PSAxOCA/IDExIDogMztcblxuXHQgICAgd2hpbGUoazwodCtpKSl7XG5cdCAgICAgIGlmKG4hPTApe1xuXHRcdHI9Wl9PSztcblx0ICAgICAgfVxuXHQgICAgICBlbHNle1xuXHRcdHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdFx0ei5hdmFpbF9pbj1uO1xuXHRcdHoudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHRcdHRoaXMud3JpdGU9cTtcblx0XHRyZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgICAgIH07XG5cdCAgICAgIG4tLTtcblx0ICAgICAgYnw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO1xuXHQgICAgICBrKz04O1xuXHQgICAgfVxuXG5cdCAgICBiPj4+PSh0KTtrLT0odCk7XG5cblx0ICAgIGogKz0gKGIgJiBpbmZsYXRlX21hc2tbaV0pO1xuXG5cdCAgICBiPj4+PShpKTtrLT0oaSk7XG5cblx0ICAgIGkgPSB0aGlzLmluZGV4O1xuXHQgICAgdCA9IHRoaXMudGFibGU7XG5cdCAgICBpZiAoaSArIGogPiAyNTggKyAodCAmIDB4MWYpICsgKCh0ID4+IDUpICYgMHgxZikgfHxcblx0XHQoYyA9PSAxNiAmJiBpIDwgMSkpe1xuXHQgICAgICB0aGlzLmJsZW5zPW51bGw7XG5cdCAgICAgIHRoaXMubW9kZSA9IElCX0JBRDtcblx0ICAgICAgei5tc2cgPSBcImludmFsaWQgYml0IGxlbmd0aCByZXBlYXRcIjtcblx0ICAgICAgciA9IFpfREFUQV9FUlJPUjtcblxuXHQgICAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICAgIHRoaXMud3JpdGU9cTtcblx0ICAgICAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgfVxuXG5cdCAgICBjID0gYyA9PSAxNiA/IHRoaXMuYmxlbnNbaS0xXSA6IDA7XG5cdCAgICBkb3tcblx0ICAgICAgdGhpcy5ibGVuc1tpKytdID0gYztcblx0ICAgIH1cblx0ICAgIHdoaWxlICgtLWohPTApO1xuXHQgICAgdGhpcy5pbmRleCA9IGk7XG5cdCAgfVxuXHR9XG5cblx0dGhpcy50YlswXT0tMTtcblx0e1xuXHQgICAgdmFyIGJsPW5ldyBJbnQzMkFycmF5KDEpO1xuXHQgICAgdmFyIGJkPW5ldyBJbnQzMkFycmF5KDEpO1xuXHQgICAgdmFyIHRsPW5ldyBJbnQzMkFycmF5KDEpO1xuXHQgICAgdmFyIHRkPW5ldyBJbnQzMkFycmF5KDEpO1xuXHQgICAgYmxbMF0gPSA5OyAgICAgICAgIC8vIG11c3QgYmUgPD0gOSBmb3IgbG9va2FoZWFkIGFzc3VtcHRpb25zXG5cdCAgICBiZFswXSA9IDY7ICAgICAgICAgLy8gbXVzdCBiZSA8PSA5IGZvciBsb29rYWhlYWQgYXNzdW1wdGlvbnNcblxuXHQgICAgdCA9IHRoaXMudGFibGU7XG5cdCAgICB0ID0gdGhpcy5pbmZ0cmVlLmluZmxhdGVfdHJlZXNfZHluYW1pYygyNTcgKyAodCAmIDB4MWYpLCBcblx0XHRcdFx0XHQgICAgICAxICsgKCh0ID4+IDUpICYgMHgxZiksXG5cdFx0XHRcdFx0ICAgICAgdGhpcy5ibGVucywgYmwsIGJkLCB0bCwgdGQsIHRoaXMuaHVmdHMsIHopO1xuXG5cdCAgICBpZiAodCAhPSBaX09LKXtcblx0ICAgICAgICBpZiAodCA9PSBaX0RBVEFfRVJST1Ipe1xuXHQgICAgICAgICAgICB0aGlzLmJsZW5zPW51bGw7XG5cdCAgICAgICAgICAgIHRoaXMubW9kZSA9IEJBRDtcblx0ICAgICAgICB9XG5cdCAgICAgICAgciA9IHQ7XG5cblx0ICAgICAgICB0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHQgICAgICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgICAgICB0aGlzLndyaXRlPXE7XG5cdCAgICAgICAgcmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgfVxuXHQgICAgdGhpcy5jb2Rlcy5pbml0KGJsWzBdLCBiZFswXSwgdGhpcy5odWZ0cywgdGxbMF0sIHRoaXMuaHVmdHMsIHRkWzBdLCB6KTtcblx0fVxuXHR0aGlzLm1vZGUgPSBJQl9DT0RFUztcbiAgICAgIGNhc2UgSUJfQ09ERVM6XG5cdHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9aztcblx0ei5hdmFpbF9pbj1uOyB6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0dGhpcy53cml0ZT1xO1xuXG5cdGlmICgociA9IHRoaXMuY29kZXMucHJvYyh0aGlzLCB6LCByKSkgIT0gWl9TVFJFQU1fRU5EKXtcblx0ICByZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHosIHIpO1xuXHR9XG5cdHIgPSBaX09LO1xuXHR0aGlzLmNvZGVzLmZyZWUoeik7XG5cblx0cD16Lm5leHRfaW5faW5kZXg7IG49ei5hdmFpbF9pbjtiPXRoaXMuYml0YjtrPXRoaXMuYml0aztcblx0cT10aGlzLndyaXRlO20gPSAocSA8IHRoaXMucmVhZCA/IHRoaXMucmVhZC1xLTEgOiB0aGlzLmVuZC1xKTtcblxuXHRpZiAodGhpcy5sYXN0PT0wKXtcblx0ICB0aGlzLm1vZGUgPSBJQl9UWVBFO1xuXHQgIGJyZWFrO1xuXHR9XG5cdHRoaXMubW9kZSA9IElCX0RSWTtcbiAgICAgIGNhc2UgSUJfRFJZOlxuXHR0aGlzLndyaXRlPXE7IFxuXHRyID0gdGhpcy5pbmZsYXRlX2ZsdXNoKHosIHIpOyBcblx0cT10aGlzLndyaXRlOyBtID0gKHEgPCB0aGlzLnJlYWQgPyB0aGlzLnJlYWQtcS0xIDogdGhpcy5lbmQtcSk7XG5cdGlmICh0aGlzLnJlYWQgIT0gdGhpcy53cml0ZSl7XG5cdCAgdGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgdGhpcy53cml0ZT1xO1xuXHQgIHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeiwgcik7XG5cdH1cblx0bW9kZSA9IERPTkU7XG4gICAgICBjYXNlIElCX0RPTkU6XG5cdHIgPSBaX1NUUkVBTV9FTkQ7XG5cblx0dGhpcy5iaXRiPWI7IHRoaXMuYml0az1rOyBcblx0ei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHR0aGlzLndyaXRlPXE7XG5cdHJldHVybiB0aGlzLmluZmxhdGVfZmx1c2goeiwgcik7XG4gICAgICBjYXNlIElCX0JBRDpcblx0ciA9IFpfREFUQV9FUlJPUjtcblxuXHR0aGlzLmJpdGI9YjsgdGhpcy5iaXRrPWs7IFxuXHR6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdHRoaXMud3JpdGU9cTtcblx0cmV0dXJuIHRoaXMuaW5mbGF0ZV9mbHVzaCh6LCByKTtcblxuICAgICAgZGVmYXVsdDpcblx0ciA9IFpfU1RSRUFNX0VSUk9SO1xuXG5cdHRoaXMuYml0Yj1iOyB0aGlzLmJpdGs9azsgXG5cdHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0dGhpcy53cml0ZT1xO1xuXHRyZXR1cm4gdGhpcy5pbmZsYXRlX2ZsdXNoKHosIHIpO1xuICAgICAgfVxuICAgIH1cbiAgfVxuXG5JbmZCbG9ja3MucHJvdG90eXBlLmZyZWUgPSBmdW5jdGlvbih6KXtcbiAgICB0aGlzLnJlc2V0KHosIG51bGwpO1xuICAgIHRoaXMud2luZG93PW51bGw7XG4gICAgdGhpcy5odWZ0cz1udWxsO1xufVxuXG5JbmZCbG9ja3MucHJvdG90eXBlLnNldF9kaWN0aW9uYXJ5ID0gZnVuY3Rpb24oZCwgc3RhcnQsIG4pe1xuICAgIGFycmF5Q29weShkLCBzdGFydCwgd2luZG93LCAwLCBuKTtcbiAgICB0aGlzLnJlYWQgPSB0aGlzLndyaXRlID0gbjtcbn1cblxuICAvLyBSZXR1cm5zIHRydWUgaWYgaW5mbGF0ZSBpcyBjdXJyZW50bHkgYXQgdGhlIGVuZCBvZiBhIGJsb2NrIGdlbmVyYXRlZFxuICAvLyBieSBaX1NZTkNfRkxVU0ggb3IgWl9GVUxMX0ZMVVNILiBcbkluZkJsb2Nrcy5wcm90b3R5cGUuc3luY19wb2ludCA9IGZ1bmN0aW9uKCl7XG4gICAgcmV0dXJuIHRoaXMubW9kZSA9PSBJQl9MRU5TO1xufVxuXG4gIC8vIGNvcHkgYXMgbXVjaCBhcyBwb3NzaWJsZSBmcm9tIHRoZSBzbGlkaW5nIHdpbmRvdyB0byB0aGUgb3V0cHV0IGFyZWFcbkluZkJsb2Nrcy5wcm90b3R5cGUuaW5mbGF0ZV9mbHVzaCA9IGZ1bmN0aW9uKHosIHIpe1xuICAgIHZhciBuO1xuICAgIHZhciBwO1xuICAgIHZhciBxO1xuXG4gICAgLy8gbG9jYWwgY29waWVzIG9mIHNvdXJjZSBhbmQgZGVzdGluYXRpb24gcG9pbnRlcnNcbiAgICBwID0gei5uZXh0X291dF9pbmRleDtcbiAgICBxID0gdGhpcy5yZWFkO1xuXG4gICAgLy8gY29tcHV0ZSBudW1iZXIgb2YgYnl0ZXMgdG8gY29weSBhcyBmYXIgYXMgZW5kIG9mIHdpbmRvd1xuICAgIG4gPSAoKHEgPD0gdGhpcy53cml0ZSA/IHRoaXMud3JpdGUgOiB0aGlzLmVuZCkgLSBxKTtcbiAgICBpZiAobiA+IHouYXZhaWxfb3V0KSBuID0gei5hdmFpbF9vdXQ7XG4gICAgaWYgKG4hPTAgJiYgciA9PSBaX0JVRl9FUlJPUikgciA9IFpfT0s7XG5cbiAgICAvLyB1cGRhdGUgY291bnRlcnNcbiAgICB6LmF2YWlsX291dCAtPSBuO1xuICAgIHoudG90YWxfb3V0ICs9IG47XG5cbiAgICAvLyB1cGRhdGUgY2hlY2sgaW5mb3JtYXRpb25cbiAgICBpZih0aGlzLmNoZWNrZm4gIT0gbnVsbClcbiAgICAgIHouYWRsZXI9dGhpcy5jaGVjaz16Ll9hZGxlci5hZGxlcjMyKHRoaXMuY2hlY2ssIHRoaXMud2luZG93LCBxLCBuKTtcblxuICAgIC8vIGNvcHkgYXMgZmFyIGFzIGVuZCBvZiB3aW5kb3dcbiAgICBhcnJheUNvcHkodGhpcy53aW5kb3csIHEsIHoubmV4dF9vdXQsIHAsIG4pO1xuICAgIHAgKz0gbjtcbiAgICBxICs9IG47XG5cbiAgICAvLyBzZWUgaWYgbW9yZSB0byBjb3B5IGF0IGJlZ2lubmluZyBvZiB3aW5kb3dcbiAgICBpZiAocSA9PSB0aGlzLmVuZCl7XG4gICAgICAvLyB3cmFwIHBvaW50ZXJzXG4gICAgICBxID0gMDtcbiAgICAgIGlmICh0aGlzLndyaXRlID09IHRoaXMuZW5kKVxuICAgICAgICB0aGlzLndyaXRlID0gMDtcblxuICAgICAgLy8gY29tcHV0ZSBieXRlcyB0byBjb3B5XG4gICAgICBuID0gdGhpcy53cml0ZSAtIHE7XG4gICAgICBpZiAobiA+IHouYXZhaWxfb3V0KSBuID0gei5hdmFpbF9vdXQ7XG4gICAgICBpZiAobiE9MCAmJiByID09IFpfQlVGX0VSUk9SKSByID0gWl9PSztcblxuICAgICAgLy8gdXBkYXRlIGNvdW50ZXJzXG4gICAgICB6LmF2YWlsX291dCAtPSBuO1xuICAgICAgei50b3RhbF9vdXQgKz0gbjtcblxuICAgICAgLy8gdXBkYXRlIGNoZWNrIGluZm9ybWF0aW9uXG4gICAgICBpZih0aGlzLmNoZWNrZm4gIT0gbnVsbClcblx0ei5hZGxlcj10aGlzLmNoZWNrPXouX2FkbGVyLmFkbGVyMzIodGhpcy5jaGVjaywgdGhpcy53aW5kb3csIHEsIG4pO1xuXG4gICAgICAvLyBjb3B5XG4gICAgICBhcnJheUNvcHkodGhpcy53aW5kb3csIHEsIHoubmV4dF9vdXQsIHAsIG4pO1xuICAgICAgcCArPSBuO1xuICAgICAgcSArPSBuO1xuICAgIH1cblxuICAgIC8vIHVwZGF0ZSBwb2ludGVyc1xuICAgIHoubmV4dF9vdXRfaW5kZXggPSBwO1xuICAgIHRoaXMucmVhZCA9IHE7XG5cbiAgICAvLyBkb25lXG4gICAgcmV0dXJuIHI7XG4gIH1cblxuLy9cbi8vIEluZkNvZGVzLmphdmFcbi8vXG5cbnZhciBJQ19TVEFSVD0wOyAgLy8geDogc2V0IHVwIGZvciBMRU5cbnZhciBJQ19MRU49MTsgICAgLy8gaTogZ2V0IGxlbmd0aC9saXRlcmFsL2VvYiBuZXh0XG52YXIgSUNfTEVORVhUPTI7IC8vIGk6IGdldHRpbmcgbGVuZ3RoIGV4dHJhIChoYXZlIGJhc2UpXG52YXIgSUNfRElTVD0zOyAgIC8vIGk6IGdldCBkaXN0YW5jZSBuZXh0XG52YXIgSUNfRElTVEVYVD00Oy8vIGk6IGdldHRpbmcgZGlzdGFuY2UgZXh0cmFcbnZhciBJQ19DT1BZPTU7ICAgLy8gbzogY29weWluZyBieXRlcyBpbiB3aW5kb3csIHdhaXRpbmcgZm9yIHNwYWNlXG52YXIgSUNfTElUPTY7ICAgIC8vIG86IGdvdCBsaXRlcmFsLCB3YWl0aW5nIGZvciBvdXRwdXQgc3BhY2VcbnZhciBJQ19XQVNIPTc7ICAgLy8gbzogZ290IGVvYiwgcG9zc2libHkgc3RpbGwgb3V0cHV0IHdhaXRpbmdcbnZhciBJQ19FTkQ9ODsgICAgLy8geDogZ290IGVvYiBhbmQgYWxsIGRhdGEgZmx1c2hlZFxudmFyIElDX0JBRENPREU9OTsvLyB4OiBnb3QgZXJyb3JcblxuZnVuY3Rpb24gSW5mQ29kZXMoKSB7XG59XG5cbkluZkNvZGVzLnByb3RvdHlwZS5pbml0ID0gZnVuY3Rpb24oYmwsIGJkLCB0bCwgdGxfaW5kZXgsIHRkLCB0ZF9pbmRleCwgeikge1xuICAgIHRoaXMubW9kZT1JQ19TVEFSVDtcbiAgICB0aGlzLmxiaXRzPWJsO1xuICAgIHRoaXMuZGJpdHM9YmQ7XG4gICAgdGhpcy5sdHJlZT10bDtcbiAgICB0aGlzLmx0cmVlX2luZGV4PXRsX2luZGV4O1xuICAgIHRoaXMuZHRyZWUgPSB0ZDtcbiAgICB0aGlzLmR0cmVlX2luZGV4PXRkX2luZGV4O1xuICAgIHRoaXMudHJlZT1udWxsO1xufVxuXG5JbmZDb2Rlcy5wcm90b3R5cGUucHJvYyA9IGZ1bmN0aW9uKHMsIHosIHIpeyBcbiAgICB2YXIgajsgICAgICAgICAgICAgIC8vIHRlbXBvcmFyeSBzdG9yYWdlXG4gICAgdmFyIHQ7ICAgICAgICAgICAgICAvLyB0ZW1wb3JhcnkgcG9pbnRlciAoaW50W10pXG4gICAgdmFyIHRpbmRleDsgICAgICAgICAvLyB0ZW1wb3JhcnkgcG9pbnRlclxuICAgIHZhciBlOyAgICAgICAgICAgICAgLy8gZXh0cmEgYml0cyBvciBvcGVyYXRpb25cbiAgICB2YXIgYj0wOyAgICAgICAgICAgIC8vIGJpdCBidWZmZXJcbiAgICB2YXIgaz0wOyAgICAgICAgICAgIC8vIGJpdHMgaW4gYml0IGJ1ZmZlclxuICAgIHZhciBwPTA7ICAgICAgICAgICAgLy8gaW5wdXQgZGF0YSBwb2ludGVyXG4gICAgdmFyIG47ICAgICAgICAgICAgICAvLyBieXRlcyBhdmFpbGFibGUgdGhlcmVcbiAgICB2YXIgcTsgICAgICAgICAgICAgIC8vIG91dHB1dCB3aW5kb3cgd3JpdGUgcG9pbnRlclxuICAgIHZhciBtOyAgICAgICAgICAgICAgLy8gYnl0ZXMgdG8gZW5kIG9mIHdpbmRvdyBvciByZWFkIHBvaW50ZXJcbiAgICB2YXIgZjsgICAgICAgICAgICAgIC8vIHBvaW50ZXIgdG8gY29weSBzdHJpbmdzIGZyb21cblxuICAgIC8vIGNvcHkgaW5wdXQvb3V0cHV0IGluZm9ybWF0aW9uIHRvIGxvY2FscyAoVVBEQVRFIG1hY3JvIHJlc3RvcmVzKVxuICAgIHA9ei5uZXh0X2luX2luZGV4O249ei5hdmFpbF9pbjtiPXMuYml0YjtrPXMuYml0aztcbiAgICBxPXMud3JpdGU7bT1xPHMucmVhZD9zLnJlYWQtcS0xOnMuZW5kLXE7XG5cbiAgICAvLyBwcm9jZXNzIGlucHV0IGFuZCBvdXRwdXQgYmFzZWQgb24gY3VycmVudCBzdGF0ZVxuICAgIHdoaWxlICh0cnVlKXtcbiAgICAgIHN3aXRjaCAodGhpcy5tb2RlKXtcblx0Ly8gd2FpdGluZyBmb3IgXCJpOlwiPWlucHV0LCBcIm86XCI9b3V0cHV0LCBcIng6XCI9bm90aGluZ1xuICAgICAgY2FzZSBJQ19TVEFSVDogICAgICAgICAvLyB4OiBzZXQgdXAgZm9yIExFTlxuXHRpZiAobSA+PSAyNTggJiYgbiA+PSAxMCl7XG5cblx0ICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgcy53cml0ZT1xO1xuXHQgIHIgPSB0aGlzLmluZmxhdGVfZmFzdCh0aGlzLmxiaXRzLCB0aGlzLmRiaXRzLCBcblx0XHRcdCAgIHRoaXMubHRyZWUsIHRoaXMubHRyZWVfaW5kZXgsIFxuXHRcdFx0ICAgdGhpcy5kdHJlZSwgdGhpcy5kdHJlZV9pbmRleCxcblx0XHRcdCAgIHMsIHopO1xuXG5cdCAgcD16Lm5leHRfaW5faW5kZXg7bj16LmF2YWlsX2luO2I9cy5iaXRiO2s9cy5iaXRrO1xuXHQgIHE9cy53cml0ZTttPXE8cy5yZWFkP3MucmVhZC1xLTE6cy5lbmQtcTtcblxuXHQgIGlmIChyICE9IFpfT0spe1xuXHQgICAgdGhpcy5tb2RlID0gciA9PSBaX1NUUkVBTV9FTkQgPyBJQ19XQVNIIDogSUNfQkFEQ09ERTtcblx0ICAgIGJyZWFrO1xuXHQgIH1cblx0fVxuXHR0aGlzLm5lZWQgPSB0aGlzLmxiaXRzO1xuXHR0aGlzLnRyZWUgPSB0aGlzLmx0cmVlO1xuXHR0aGlzLnRyZWVfaW5kZXg9dGhpcy5sdHJlZV9pbmRleDtcblxuXHR0aGlzLm1vZGUgPSBJQ19MRU47XG4gICAgICBjYXNlIElDX0xFTjogICAgICAgICAgIC8vIGk6IGdldCBsZW5ndGgvbGl0ZXJhbC9lb2IgbmV4dFxuXHRqID0gdGhpcy5uZWVkO1xuXG5cdHdoaWxlKGs8KGopKXtcblx0ICBpZihuIT0wKXI9Wl9PSztcblx0ICBlbHNle1xuXG5cdCAgICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgIHMud3JpdGU9cTtcblx0ICAgIHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICB9XG5cdCAgbi0tO1xuXHQgIGJ8PSh6Lm5leHRfaW5bcCsrXSYweGZmKTw8aztcblx0ICBrKz04O1xuXHR9XG5cblx0dGluZGV4PSh0aGlzLnRyZWVfaW5kZXgrKGImaW5mbGF0ZV9tYXNrW2pdKSkqMztcblxuXHRiPj4+PSh0aGlzLnRyZWVbdGluZGV4KzFdKTtcblx0ay09KHRoaXMudHJlZVt0aW5kZXgrMV0pO1xuXG5cdGU9dGhpcy50cmVlW3RpbmRleF07XG5cblx0aWYoZSA9PSAwKXsgICAgICAgICAgICAgICAvLyBsaXRlcmFsXG5cdCAgdGhpcy5saXQgPSB0aGlzLnRyZWVbdGluZGV4KzJdO1xuXHQgIHRoaXMubW9kZSA9IElDX0xJVDtcblx0ICBicmVhaztcblx0fVxuXHRpZigoZSAmIDE2KSE9MCApeyAgICAgICAgICAvLyBsZW5ndGhcblx0ICB0aGlzLmdldCA9IGUgJiAxNTtcblx0ICB0aGlzLmxlbiA9IHRoaXMudHJlZVt0aW5kZXgrMl07XG5cdCAgdGhpcy5tb2RlID0gSUNfTEVORVhUO1xuXHQgIGJyZWFrO1xuXHR9XG5cdGlmICgoZSAmIDY0KSA9PSAwKXsgICAgICAgIC8vIG5leHQgdGFibGVcblx0ICB0aGlzLm5lZWQgPSBlO1xuXHQgIHRoaXMudHJlZV9pbmRleCA9IHRpbmRleC8zICsgdGhpcy50cmVlW3RpbmRleCsyXTtcblx0ICBicmVhaztcblx0fVxuXHRpZiAoKGUgJiAzMikhPTApeyAgICAgICAgICAgICAgIC8vIGVuZCBvZiBibG9ja1xuXHQgIHRoaXMubW9kZSA9IElDX1dBU0g7XG5cdCAgYnJlYWs7XG5cdH1cblx0dGhpcy5tb2RlID0gSUNfQkFEQ09ERTsgICAgICAgIC8vIGludmFsaWQgY29kZVxuXHR6Lm1zZyA9IFwiaW52YWxpZCBsaXRlcmFsL2xlbmd0aCBjb2RlXCI7XG5cdHIgPSBaX0RBVEFfRVJST1I7XG5cblx0cy5iaXRiPWI7cy5iaXRrPWs7XG5cdHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0cy53cml0ZT1xO1xuXHRyZXR1cm4gcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cbiAgICAgIGNhc2UgSUNfTEVORVhUOiAgICAgICAgLy8gaTogZ2V0dGluZyBsZW5ndGggZXh0cmEgKGhhdmUgYmFzZSlcblx0aiA9IHRoaXMuZ2V0O1xuXG5cdHdoaWxlKGs8KGopKXtcblx0ICBpZihuIT0wKXI9Wl9PSztcblx0ICBlbHNle1xuXG5cdCAgICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgIHMud3JpdGU9cTtcblx0ICAgIHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICB9XG5cdCAgbi0tOyBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7XG5cdCAgays9ODtcblx0fVxuXG5cdHRoaXMubGVuICs9IChiICYgaW5mbGF0ZV9tYXNrW2pdKTtcblxuXHRiPj49ajtcblx0ay09ajtcblxuXHR0aGlzLm5lZWQgPSB0aGlzLmRiaXRzO1xuXHR0aGlzLnRyZWUgPSB0aGlzLmR0cmVlO1xuXHR0aGlzLnRyZWVfaW5kZXggPSB0aGlzLmR0cmVlX2luZGV4O1xuXHR0aGlzLm1vZGUgPSBJQ19ESVNUO1xuICAgICAgY2FzZSBJQ19ESVNUOiAgICAgICAgICAvLyBpOiBnZXQgZGlzdGFuY2UgbmV4dFxuXHRqID0gdGhpcy5uZWVkO1xuXG5cdHdoaWxlKGs8KGopKXtcblx0ICBpZihuIT0wKXI9Wl9PSztcblx0ICBlbHNle1xuXG5cdCAgICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICAgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICAgIHMud3JpdGU9cTtcblx0ICAgIHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICB9XG5cdCAgbi0tOyBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7XG5cdCAgays9ODtcblx0fVxuXG5cdHRpbmRleD0odGhpcy50cmVlX2luZGV4KyhiICYgaW5mbGF0ZV9tYXNrW2pdKSkqMztcblxuXHRiPj49dGhpcy50cmVlW3RpbmRleCsxXTtcblx0ay09dGhpcy50cmVlW3RpbmRleCsxXTtcblxuXHRlID0gKHRoaXMudHJlZVt0aW5kZXhdKTtcblx0aWYoKGUgJiAxNikhPTApeyAgICAgICAgICAgICAgIC8vIGRpc3RhbmNlXG5cdCAgdGhpcy5nZXQgPSBlICYgMTU7XG5cdCAgdGhpcy5kaXN0ID0gdGhpcy50cmVlW3RpbmRleCsyXTtcblx0ICB0aGlzLm1vZGUgPSBJQ19ESVNURVhUO1xuXHQgIGJyZWFrO1xuXHR9XG5cdGlmICgoZSAmIDY0KSA9PSAwKXsgICAgICAgIC8vIG5leHQgdGFibGVcblx0ICB0aGlzLm5lZWQgPSBlO1xuXHQgIHRoaXMudHJlZV9pbmRleCA9IHRpbmRleC8zICsgdGhpcy50cmVlW3RpbmRleCsyXTtcblx0ICBicmVhaztcblx0fVxuXHR0aGlzLm1vZGUgPSBJQ19CQURDT0RFOyAgICAgICAgLy8gaW52YWxpZCBjb2RlXG5cdHoubXNnID0gXCJpbnZhbGlkIGRpc3RhbmNlIGNvZGVcIjtcblx0ciA9IFpfREFUQV9FUlJPUjtcblxuXHRzLmJpdGI9YjtzLmJpdGs9aztcblx0ei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHRzLndyaXRlPXE7XG5cdHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblxuICAgICAgY2FzZSBJQ19ESVNURVhUOiAgICAgICAvLyBpOiBnZXR0aW5nIGRpc3RhbmNlIGV4dHJhXG5cdGogPSB0aGlzLmdldDtcblxuXHR3aGlsZShrPChqKSl7XG5cdCAgaWYobiE9MClyPVpfT0s7XG5cdCAgZWxzZXtcblxuXHQgICAgcy5iaXRiPWI7cy5iaXRrPWs7XG5cdCAgICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICBzLndyaXRlPXE7XG5cdCAgICByZXR1cm4gcy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgfVxuXHQgIG4tLTsgYnw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO1xuXHQgIGsrPTg7XG5cdH1cblxuXHR0aGlzLmRpc3QgKz0gKGIgJiBpbmZsYXRlX21hc2tbal0pO1xuXG5cdGI+Pj1qO1xuXHRrLT1qO1xuXG5cdHRoaXMubW9kZSA9IElDX0NPUFk7XG4gICAgICBjYXNlIElDX0NPUFk6ICAgICAgICAgIC8vIG86IGNvcHlpbmcgYnl0ZXMgaW4gd2luZG93LCB3YWl0aW5nIGZvciBzcGFjZVxuICAgICAgICBmID0gcSAtIHRoaXMuZGlzdDtcbiAgICAgICAgd2hpbGUoZiA8IDApeyAgICAgLy8gbW9kdWxvIHdpbmRvdyBzaXplLVwid2hpbGVcIiBpbnN0ZWFkXG4gICAgICAgICAgZiArPSBzLmVuZDsgICAgIC8vIG9mIFwiaWZcIiBoYW5kbGVzIGludmFsaWQgZGlzdGFuY2VzXG5cdH1cblx0d2hpbGUgKHRoaXMubGVuIT0wKXtcblxuXHQgIGlmKG09PTApe1xuXHQgICAgaWYocT09cy5lbmQmJnMucmVhZCE9MCl7cT0wO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO31cblx0ICAgIGlmKG09PTApe1xuXHQgICAgICBzLndyaXRlPXE7IHI9cy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdCAgICAgIHE9cy53cml0ZTttPXE8cy5yZWFkP3MucmVhZC1xLTE6cy5lbmQtcTtcblxuXHQgICAgICBpZihxPT1zLmVuZCYmcy5yZWFkIT0wKXtxPTA7bT1xPHMucmVhZD9zLnJlYWQtcS0xOnMuZW5kLXE7fVxuXG5cdCAgICAgIGlmKG09PTApe1xuXHRcdHMuYml0Yj1iO3MuYml0az1rO1xuXHRcdHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0XHRzLndyaXRlPXE7XG5cdFx0cmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHQgICAgICB9ICBcblx0ICAgIH1cblx0ICB9XG5cblx0ICBzLndpbmRvd1txKytdPXMud2luZG93W2YrK107IG0tLTtcblxuXHQgIGlmIChmID09IHMuZW5kKVxuICAgICAgICAgICAgZiA9IDA7XG5cdCAgdGhpcy5sZW4tLTtcblx0fVxuXHR0aGlzLm1vZGUgPSBJQ19TVEFSVDtcblx0YnJlYWs7XG4gICAgICBjYXNlIElDX0xJVDogICAgICAgICAgIC8vIG86IGdvdCBsaXRlcmFsLCB3YWl0aW5nIGZvciBvdXRwdXQgc3BhY2Vcblx0aWYobT09MCl7XG5cdCAgaWYocT09cy5lbmQmJnMucmVhZCE9MCl7cT0wO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO31cblx0ICBpZihtPT0wKXtcblx0ICAgIHMud3JpdGU9cTsgcj1zLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICAgIHE9cy53cml0ZTttPXE8cy5yZWFkP3MucmVhZC1xLTE6cy5lbmQtcTtcblxuXHQgICAgaWYocT09cy5lbmQmJnMucmVhZCE9MCl7cT0wO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO31cblx0ICAgIGlmKG09PTApe1xuXHQgICAgICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICAgICAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgICAgICBzLndyaXRlPXE7XG5cdCAgICAgIHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblx0ICAgIH1cblx0ICB9XG5cdH1cblx0cj1aX09LO1xuXG5cdHMud2luZG93W3ErK109dGhpcy5saXQ7IG0tLTtcblxuXHR0aGlzLm1vZGUgPSBJQ19TVEFSVDtcblx0YnJlYWs7XG4gICAgICBjYXNlIElDX1dBU0g6ICAgICAgICAgICAvLyBvOiBnb3QgZW9iLCBwb3NzaWJseSBtb3JlIG91dHB1dFxuXHRpZiAoayA+IDcpeyAgICAgICAgLy8gcmV0dXJuIHVudXNlZCBieXRlLCBpZiBhbnlcblx0ICBrIC09IDg7XG5cdCAgbisrO1xuXHQgIHAtLTsgICAgICAgICAgICAgLy8gY2FuIGFsd2F5cyByZXR1cm4gb25lXG5cdH1cblxuXHRzLndyaXRlPXE7IHI9cy5pbmZsYXRlX2ZsdXNoKHoscik7XG5cdHE9cy53cml0ZTttPXE8cy5yZWFkP3MucmVhZC1xLTE6cy5lbmQtcTtcblxuXHRpZiAocy5yZWFkICE9IHMud3JpdGUpe1xuXHQgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgIHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0ICBzLndyaXRlPXE7XG5cdCAgcmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXHR9XG5cdHRoaXMubW9kZSA9IElDX0VORDtcbiAgICAgIGNhc2UgSUNfRU5EOlxuXHRyID0gWl9TVFJFQU1fRU5EO1xuXHRzLmJpdGI9YjtzLmJpdGs9aztcblx0ei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHRzLndyaXRlPXE7XG5cdHJldHVybiBzLmluZmxhdGVfZmx1c2goeixyKTtcblxuICAgICAgY2FzZSBJQ19CQURDT0RFOiAgICAgICAvLyB4OiBnb3QgZXJyb3JcblxuXHRyID0gWl9EQVRBX0VSUk9SO1xuXG5cdHMuYml0Yj1iO3MuYml0az1rO1xuXHR6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdHMud3JpdGU9cTtcblx0cmV0dXJuIHMuaW5mbGF0ZV9mbHVzaCh6LHIpO1xuXG4gICAgICBkZWZhdWx0OlxuXHRyID0gWl9TVFJFQU1fRVJST1I7XG5cblx0cy5iaXRiPWI7cy5iaXRrPWs7XG5cdHouYXZhaWxfaW49bjt6LnRvdGFsX2luKz1wLXoubmV4dF9pbl9pbmRleDt6Lm5leHRfaW5faW5kZXg9cDtcblx0cy53cml0ZT1xO1xuXHRyZXR1cm4gcy5pbmZsYXRlX2ZsdXNoKHoscik7XG4gICAgICB9XG4gICAgfVxuICB9XG5cbkluZkNvZGVzLnByb3RvdHlwZS5mcmVlID0gZnVuY3Rpb24oeil7XG4gICAgLy8gIFpGUkVFKHosIGMpO1xufVxuXG4gIC8vIENhbGxlZCB3aXRoIG51bWJlciBvZiBieXRlcyBsZWZ0IHRvIHdyaXRlIGluIHdpbmRvdyBhdCBsZWFzdCAyNThcbiAgLy8gKHRoZSBtYXhpbXVtIHN0cmluZyBsZW5ndGgpIGFuZCBudW1iZXIgb2YgaW5wdXQgYnl0ZXMgYXZhaWxhYmxlXG4gIC8vIGF0IGxlYXN0IHRlbi4gIFRoZSB0ZW4gYnl0ZXMgYXJlIHNpeCBieXRlcyBmb3IgdGhlIGxvbmdlc3QgbGVuZ3RoL1xuICAvLyBkaXN0YW5jZSBwYWlyIHBsdXMgZm91ciBieXRlcyBmb3Igb3ZlcmxvYWRpbmcgdGhlIGJpdCBidWZmZXIuXG5cbkluZkNvZGVzLnByb3RvdHlwZS5pbmZsYXRlX2Zhc3QgPSBmdW5jdGlvbihibCwgYmQsIHRsLCB0bF9pbmRleCwgdGQsIHRkX2luZGV4LCBzLCB6KSB7XG4gICAgdmFyIHQ7ICAgICAgICAgICAgICAgIC8vIHRlbXBvcmFyeSBwb2ludGVyXG4gICAgdmFyICAgdHA7ICAgICAgICAgICAgIC8vIHRlbXBvcmFyeSBwb2ludGVyIChpbnRbXSlcbiAgICB2YXIgdHBfaW5kZXg7ICAgICAgICAgLy8gdGVtcG9yYXJ5IHBvaW50ZXJcbiAgICB2YXIgZTsgICAgICAgICAgICAgICAgLy8gZXh0cmEgYml0cyBvciBvcGVyYXRpb25cbiAgICB2YXIgYjsgICAgICAgICAgICAgICAgLy8gYml0IGJ1ZmZlclxuICAgIHZhciBrOyAgICAgICAgICAgICAgICAvLyBiaXRzIGluIGJpdCBidWZmZXJcbiAgICB2YXIgcDsgICAgICAgICAgICAgICAgLy8gaW5wdXQgZGF0YSBwb2ludGVyXG4gICAgdmFyIG47ICAgICAgICAgICAgICAgIC8vIGJ5dGVzIGF2YWlsYWJsZSB0aGVyZVxuICAgIHZhciBxOyAgICAgICAgICAgICAgICAvLyBvdXRwdXQgd2luZG93IHdyaXRlIHBvaW50ZXJcbiAgICB2YXIgbTsgICAgICAgICAgICAgICAgLy8gYnl0ZXMgdG8gZW5kIG9mIHdpbmRvdyBvciByZWFkIHBvaW50ZXJcbiAgICB2YXIgbWw7ICAgICAgICAgICAgICAgLy8gbWFzayBmb3IgbGl0ZXJhbC9sZW5ndGggdHJlZVxuICAgIHZhciBtZDsgICAgICAgICAgICAgICAvLyBtYXNrIGZvciBkaXN0YW5jZSB0cmVlXG4gICAgdmFyIGM7ICAgICAgICAgICAgICAgIC8vIGJ5dGVzIHRvIGNvcHlcbiAgICB2YXIgZDsgICAgICAgICAgICAgICAgLy8gZGlzdGFuY2UgYmFjayB0byBjb3B5IGZyb21cbiAgICB2YXIgcjsgICAgICAgICAgICAgICAgLy8gY29weSBzb3VyY2UgcG9pbnRlclxuXG4gICAgdmFyIHRwX2luZGV4X3RfMzsgICAgIC8vICh0cF9pbmRleCt0KSozXG5cbiAgICAvLyBsb2FkIGlucHV0LCBvdXRwdXQsIGJpdCB2YWx1ZXNcbiAgICBwPXoubmV4dF9pbl9pbmRleDtuPXouYXZhaWxfaW47Yj1zLmJpdGI7az1zLmJpdGs7XG4gICAgcT1zLndyaXRlO209cTxzLnJlYWQ/cy5yZWFkLXEtMTpzLmVuZC1xO1xuXG4gICAgLy8gaW5pdGlhbGl6ZSBtYXNrc1xuICAgIG1sID0gaW5mbGF0ZV9tYXNrW2JsXTtcbiAgICBtZCA9IGluZmxhdGVfbWFza1tiZF07XG5cbiAgICAvLyBkbyB1bnRpbCBub3QgZW5vdWdoIGlucHV0IG9yIG91dHB1dCBzcGFjZSBmb3IgZmFzdCBsb29wXG4gICAgZG8geyAgICAgICAgICAgICAgICAgICAgICAgICAgLy8gYXNzdW1lIGNhbGxlZCB3aXRoIG0gPj0gMjU4ICYmIG4gPj0gMTBcbiAgICAgIC8vIGdldCBsaXRlcmFsL2xlbmd0aCBjb2RlXG4gICAgICB3aGlsZShrPCgyMCkpeyAgICAgICAgICAgICAgLy8gbWF4IGJpdHMgZm9yIGxpdGVyYWwvbGVuZ3RoIGNvZGVcblx0bi0tO1xuXHRifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7ays9ODtcbiAgICAgIH1cblxuICAgICAgdD0gYiZtbDtcbiAgICAgIHRwPXRsOyBcbiAgICAgIHRwX2luZGV4PXRsX2luZGV4O1xuICAgICAgdHBfaW5kZXhfdF8zPSh0cF9pbmRleCt0KSozO1xuICAgICAgaWYgKChlID0gdHBbdHBfaW5kZXhfdF8zXSkgPT0gMCl7XG5cdGI+Pj0odHBbdHBfaW5kZXhfdF8zKzFdKTsgay09KHRwW3RwX2luZGV4X3RfMysxXSk7XG5cblx0cy53aW5kb3dbcSsrXSA9IHRwW3RwX2luZGV4X3RfMysyXTtcblx0bS0tO1xuXHRjb250aW51ZTtcbiAgICAgIH1cbiAgICAgIGRvIHtcblxuXHRiPj49KHRwW3RwX2luZGV4X3RfMysxXSk7IGstPSh0cFt0cF9pbmRleF90XzMrMV0pO1xuXG5cdGlmKChlJjE2KSE9MCl7XG5cdCAgZSAmPSAxNTtcblx0ICBjID0gdHBbdHBfaW5kZXhfdF8zKzJdICsgKGIgJiBpbmZsYXRlX21hc2tbZV0pO1xuXG5cdCAgYj4+PWU7IGstPWU7XG5cblx0ICAvLyBkZWNvZGUgZGlzdGFuY2UgYmFzZSBvZiBibG9jayB0byBjb3B5XG5cdCAgd2hpbGUoazwoMTUpKXsgICAgICAgICAgIC8vIG1heCBiaXRzIGZvciBkaXN0YW5jZSBjb2RlXG5cdCAgICBuLS07XG5cdCAgICBifD0oei5uZXh0X2luW3ArK10mMHhmZik8PGs7ays9ODtcblx0ICB9XG5cblx0ICB0PSBiJm1kO1xuXHQgIHRwPXRkO1xuXHQgIHRwX2luZGV4PXRkX2luZGV4O1xuICAgICAgICAgIHRwX2luZGV4X3RfMz0odHBfaW5kZXgrdCkqMztcblx0ICBlID0gdHBbdHBfaW5kZXhfdF8zXTtcblxuXHQgIGRvIHtcblxuXHQgICAgYj4+PSh0cFt0cF9pbmRleF90XzMrMV0pOyBrLT0odHBbdHBfaW5kZXhfdF8zKzFdKTtcblxuXHQgICAgaWYoKGUmMTYpIT0wKXtcblx0ICAgICAgLy8gZ2V0IGV4dHJhIGJpdHMgdG8gYWRkIHRvIGRpc3RhbmNlIGJhc2Vcblx0ICAgICAgZSAmPSAxNTtcblx0ICAgICAgd2hpbGUoazwoZSkpeyAgICAgICAgIC8vIGdldCBleHRyYSBiaXRzICh1cCB0byAxMylcblx0XHRuLS07XG5cdFx0Ynw9KHoubmV4dF9pbltwKytdJjB4ZmYpPDxrO2srPTg7XG5cdCAgICAgIH1cblxuXHQgICAgICBkID0gdHBbdHBfaW5kZXhfdF8zKzJdICsgKGImaW5mbGF0ZV9tYXNrW2VdKTtcblxuXHQgICAgICBiPj49KGUpOyBrLT0oZSk7XG5cblx0ICAgICAgLy8gZG8gdGhlIGNvcHlcblx0ICAgICAgbSAtPSBjO1xuXHQgICAgICBpZiAocSA+PSBkKXsgICAgICAgICAgICAgICAgLy8gb2Zmc2V0IGJlZm9yZSBkZXN0XG5cdFx0Ly8gIGp1c3QgY29weVxuXHRcdHI9cS1kO1xuXHRcdGlmKHEtcj4wICYmIDI+KHEtcikpeyAgICAgICAgICAgXG5cdFx0ICBzLndpbmRvd1txKytdPXMud2luZG93W3IrK107IC8vIG1pbmltdW0gY291bnQgaXMgdGhyZWUsXG5cdFx0ICBzLndpbmRvd1txKytdPXMud2luZG93W3IrK107IC8vIHNvIHVucm9sbCBsb29wIGEgbGl0dGxlXG5cdFx0ICBjLT0yO1xuXHRcdH1cblx0XHRlbHNle1xuXHRcdCAgcy53aW5kb3dbcSsrXT1zLndpbmRvd1tyKytdOyAvLyBtaW5pbXVtIGNvdW50IGlzIHRocmVlLFxuXHRcdCAgcy53aW5kb3dbcSsrXT1zLndpbmRvd1tyKytdOyAvLyBzbyB1bnJvbGwgbG9vcCBhIGxpdHRsZVxuXHRcdCAgYy09Mjtcblx0XHR9XG5cdCAgICAgIH1cblx0ICAgICAgZWxzZXsgICAgICAgICAgICAgICAgICAvLyBlbHNlIG9mZnNldCBhZnRlciBkZXN0aW5hdGlvblxuICAgICAgICAgICAgICAgIHI9cS1kO1xuICAgICAgICAgICAgICAgIGRve1xuICAgICAgICAgICAgICAgICAgcis9cy5lbmQ7ICAgICAgICAgIC8vIGZvcmNlIHBvaW50ZXIgaW4gd2luZG93XG4gICAgICAgICAgICAgICAgfXdoaWxlKHI8MCk7ICAgICAgICAgLy8gY292ZXJzIGludmFsaWQgZGlzdGFuY2VzXG5cdFx0ZT1zLmVuZC1yO1xuXHRcdGlmKGM+ZSl7ICAgICAgICAgICAgIC8vIGlmIHNvdXJjZSBjcm9zc2VzLFxuXHRcdCAgYy09ZTsgICAgICAgICAgICAgIC8vIHdyYXBwZWQgY29weVxuXHRcdCAgaWYocS1yPjAgJiYgZT4ocS1yKSl7ICAgICAgICAgICBcblx0XHQgICAgZG97cy53aW5kb3dbcSsrXSA9IHMud2luZG93W3IrK107fVxuXHRcdCAgICB3aGlsZSgtLWUhPTApO1xuXHRcdCAgfVxuXHRcdCAgZWxzZXtcblx0XHQgICAgYXJyYXlDb3B5KHMud2luZG93LCByLCBzLndpbmRvdywgcSwgZSk7XG5cdFx0ICAgIHErPWU7IHIrPWU7IGU9MDtcblx0XHQgIH1cblx0XHQgIHIgPSAwOyAgICAgICAgICAgICAgICAgIC8vIGNvcHkgcmVzdCBmcm9tIHN0YXJ0IG9mIHdpbmRvd1xuXHRcdH1cblxuXHQgICAgICB9XG5cblx0ICAgICAgLy8gY29weSBhbGwgb3Igd2hhdCdzIGxlZnRcbiAgICAgICAgICAgICAgZG97cy53aW5kb3dbcSsrXSA9IHMud2luZG93W3IrK107fVxuXHRcdHdoaWxlKC0tYyE9MCk7XG5cdCAgICAgIGJyZWFrO1xuXHQgICAgfVxuXHQgICAgZWxzZSBpZigoZSY2NCk9PTApe1xuXHQgICAgICB0Kz10cFt0cF9pbmRleF90XzMrMl07XG5cdCAgICAgIHQrPShiJmluZmxhdGVfbWFza1tlXSk7XG5cdCAgICAgIHRwX2luZGV4X3RfMz0odHBfaW5kZXgrdCkqMztcblx0ICAgICAgZT10cFt0cF9pbmRleF90XzNdO1xuXHQgICAgfVxuXHQgICAgZWxzZXtcblx0ICAgICAgei5tc2cgPSBcImludmFsaWQgZGlzdGFuY2UgY29kZVwiO1xuXG5cdCAgICAgIGM9ei5hdmFpbF9pbi1uO2M9KGs+PjMpPGM/az4+MzpjO24rPWM7cC09YztrLT1jPDwzO1xuXG5cdCAgICAgIHMuYml0Yj1iO3MuYml0az1rO1xuXHQgICAgICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgICAgIHMud3JpdGU9cTtcblxuXHQgICAgICByZXR1cm4gWl9EQVRBX0VSUk9SO1xuXHQgICAgfVxuXHQgIH1cblx0ICB3aGlsZSh0cnVlKTtcblx0ICBicmVhaztcblx0fVxuXG5cdGlmKChlJjY0KT09MCl7XG5cdCAgdCs9dHBbdHBfaW5kZXhfdF8zKzJdO1xuXHQgIHQrPShiJmluZmxhdGVfbWFza1tlXSk7XG5cdCAgdHBfaW5kZXhfdF8zPSh0cF9pbmRleCt0KSozO1xuXHQgIGlmKChlPXRwW3RwX2luZGV4X3RfM10pPT0wKXtcblxuXHQgICAgYj4+PSh0cFt0cF9pbmRleF90XzMrMV0pOyBrLT0odHBbdHBfaW5kZXhfdF8zKzFdKTtcblxuXHQgICAgcy53aW5kb3dbcSsrXT10cFt0cF9pbmRleF90XzMrMl07XG5cdCAgICBtLS07XG5cdCAgICBicmVhaztcblx0ICB9XG5cdH1cblx0ZWxzZSBpZigoZSYzMikhPTApe1xuXG5cdCAgYz16LmF2YWlsX2luLW47Yz0oaz4+Myk8Yz9rPj4zOmM7bis9YztwLT1jO2stPWM8PDM7XG4gXG5cdCAgcy5iaXRiPWI7cy5iaXRrPWs7XG5cdCAgei5hdmFpbF9pbj1uO3oudG90YWxfaW4rPXAtei5uZXh0X2luX2luZGV4O3oubmV4dF9pbl9pbmRleD1wO1xuXHQgIHMud3JpdGU9cTtcblxuXHQgIHJldHVybiBaX1NUUkVBTV9FTkQ7XG5cdH1cblx0ZWxzZXtcblx0ICB6Lm1zZz1cImludmFsaWQgbGl0ZXJhbC9sZW5ndGggY29kZVwiO1xuXG5cdCAgYz16LmF2YWlsX2luLW47Yz0oaz4+Myk8Yz9rPj4zOmM7bis9YztwLT1jO2stPWM8PDM7XG5cblx0ICBzLmJpdGI9YjtzLmJpdGs9aztcblx0ICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG5cdCAgcy53cml0ZT1xO1xuXG5cdCAgcmV0dXJuIFpfREFUQV9FUlJPUjtcblx0fVxuICAgICAgfSBcbiAgICAgIHdoaWxlKHRydWUpO1xuICAgIH0gXG4gICAgd2hpbGUobT49MjU4ICYmIG4+PSAxMCk7XG5cbiAgICAvLyBub3QgZW5vdWdoIGlucHV0IG9yIG91dHB1dC0tcmVzdG9yZSBwb2ludGVycyBhbmQgcmV0dXJuXG4gICAgYz16LmF2YWlsX2luLW47Yz0oaz4+Myk8Yz9rPj4zOmM7bis9YztwLT1jO2stPWM8PDM7XG5cbiAgICBzLmJpdGI9YjtzLmJpdGs9aztcbiAgICB6LmF2YWlsX2luPW47ei50b3RhbF9pbis9cC16Lm5leHRfaW5faW5kZXg7ei5uZXh0X2luX2luZGV4PXA7XG4gICAgcy53cml0ZT1xO1xuXG4gICAgcmV0dXJuIFpfT0s7XG59XG5cbi8vXG4vLyBJbmZUcmVlLmphdmFcbi8vXG5cbmZ1bmN0aW9uIEluZlRyZWUoKSB7XG59XG5cbkluZlRyZWUucHJvdG90eXBlLmh1ZnRfYnVpbGQgPSBmdW5jdGlvbihiLCBiaW5kZXgsIG4sIHMsIGQsIGUsIHQsIG0sIGhwLCBobiwgdikge1xuXG4gICAgLy8gR2l2ZW4gYSBsaXN0IG9mIGNvZGUgbGVuZ3RocyBhbmQgYSBtYXhpbXVtIHRhYmxlIHNpemUsIG1ha2UgYSBzZXQgb2ZcbiAgICAvLyB0YWJsZXMgdG8gZGVjb2RlIHRoYXQgc2V0IG9mIGNvZGVzLiAgUmV0dXJuIFpfT0sgb24gc3VjY2VzcywgWl9CVUZfRVJST1JcbiAgICAvLyBpZiB0aGUgZ2l2ZW4gY29kZSBzZXQgaXMgaW5jb21wbGV0ZSAodGhlIHRhYmxlcyBhcmUgc3RpbGwgYnVpbHQgaW4gdGhpc1xuICAgIC8vIGNhc2UpLCBaX0RBVEFfRVJST1IgaWYgdGhlIGlucHV0IGlzIGludmFsaWQgKGFuIG92ZXItc3Vic2NyaWJlZCBzZXQgb2ZcbiAgICAvLyBsZW5ndGhzKSwgb3IgWl9NRU1fRVJST1IgaWYgbm90IGVub3VnaCBtZW1vcnkuXG5cbiAgICB2YXIgYTsgICAgICAgICAgICAgICAgICAgICAgIC8vIGNvdW50ZXIgZm9yIGNvZGVzIG9mIGxlbmd0aCBrXG4gICAgdmFyIGY7ICAgICAgICAgICAgICAgICAgICAgICAvLyBpIHJlcGVhdHMgaW4gdGFibGUgZXZlcnkgZiBlbnRyaWVzXG4gICAgdmFyIGc7ICAgICAgICAgICAgICAgICAgICAgICAvLyBtYXhpbXVtIGNvZGUgbGVuZ3RoXG4gICAgdmFyIGg7ICAgICAgICAgICAgICAgICAgICAgICAvLyB0YWJsZSBsZXZlbFxuICAgIHZhciBpOyAgICAgICAgICAgICAgICAgICAgICAgLy8gY291bnRlciwgY3VycmVudCBjb2RlXG4gICAgdmFyIGo7ICAgICAgICAgICAgICAgICAgICAgICAvLyBjb3VudGVyXG4gICAgdmFyIGs7ICAgICAgICAgICAgICAgICAgICAgICAvLyBudW1iZXIgb2YgYml0cyBpbiBjdXJyZW50IGNvZGVcbiAgICB2YXIgbDsgICAgICAgICAgICAgICAgICAgICAgIC8vIGJpdHMgcGVyIHRhYmxlIChyZXR1cm5lZCBpbiBtKVxuICAgIHZhciBtYXNrOyAgICAgICAgICAgICAgICAgICAgLy8gKDEgPDwgdykgLSAxLCB0byBhdm9pZCBjYyAtTyBidWcgb24gSFBcbiAgICB2YXIgcDsgICAgICAgICAgICAgICAgICAgICAgIC8vIHBvaW50ZXIgaW50byBjW10sIGJbXSwgb3IgdltdXG4gICAgdmFyIHE7ICAgICAgICAgICAgICAgICAgICAgICAvLyBwb2ludHMgdG8gY3VycmVudCB0YWJsZVxuICAgIHZhciB3OyAgICAgICAgICAgICAgICAgICAgICAgLy8gYml0cyBiZWZvcmUgdGhpcyB0YWJsZSA9PSAobCAqIGgpXG4gICAgdmFyIHhwOyAgICAgICAgICAgICAgICAgICAgICAvLyBwb2ludGVyIGludG8geFxuICAgIHZhciB5OyAgICAgICAgICAgICAgICAgICAgICAgLy8gbnVtYmVyIG9mIGR1bW15IGNvZGVzIGFkZGVkXG4gICAgdmFyIHo7ICAgICAgICAgICAgICAgICAgICAgICAvLyBudW1iZXIgb2YgZW50cmllcyBpbiBjdXJyZW50IHRhYmxlXG5cbiAgICAvLyBHZW5lcmF0ZSBjb3VudHMgZm9yIGVhY2ggYml0IGxlbmd0aFxuXG4gICAgcCA9IDA7IGkgPSBuO1xuICAgIGRvIHtcbiAgICAgIHRoaXMuY1tiW2JpbmRleCtwXV0rKzsgcCsrOyBpLS07ICAgLy8gYXNzdW1lIGFsbCBlbnRyaWVzIDw9IEJNQVhcbiAgICB9d2hpbGUoaSE9MCk7XG5cbiAgICBpZih0aGlzLmNbMF0gPT0gbil7ICAgICAgICAgICAgICAgIC8vIG51bGwgaW5wdXQtLWFsbCB6ZXJvIGxlbmd0aCBjb2Rlc1xuICAgICAgdFswXSA9IC0xO1xuICAgICAgbVswXSA9IDA7XG4gICAgICByZXR1cm4gWl9PSztcbiAgICB9XG5cbiAgICAvLyBGaW5kIG1pbmltdW0gYW5kIG1heGltdW0gbGVuZ3RoLCBib3VuZCAqbSBieSB0aG9zZVxuICAgIGwgPSBtWzBdO1xuICAgIGZvciAoaiA9IDE7IGogPD0gQk1BWDsgaisrKVxuICAgICAgaWYodGhpcy5jW2pdIT0wKSBicmVhaztcbiAgICBrID0gajsgICAgICAgICAgICAgICAgICAgICAgICAvLyBtaW5pbXVtIGNvZGUgbGVuZ3RoXG4gICAgaWYobCA8IGope1xuICAgICAgbCA9IGo7XG4gICAgfVxuICAgIGZvciAoaSA9IEJNQVg7IGkhPTA7IGktLSl7XG4gICAgICBpZih0aGlzLmNbaV0hPTApIGJyZWFrO1xuICAgIH1cbiAgICBnID0gaTsgICAgICAgICAgICAgICAgICAgICAgICAvLyBtYXhpbXVtIGNvZGUgbGVuZ3RoXG4gICAgaWYobCA+IGkpe1xuICAgICAgbCA9IGk7XG4gICAgfVxuICAgIG1bMF0gPSBsO1xuXG4gICAgLy8gQWRqdXN0IGxhc3QgbGVuZ3RoIGNvdW50IHRvIGZpbGwgb3V0IGNvZGVzLCBpZiBuZWVkZWRcbiAgICBmb3IgKHkgPSAxIDw8IGo7IGogPCBpOyBqKyssIHkgPDw9IDEpe1xuICAgICAgaWYgKCh5IC09IHRoaXMuY1tqXSkgPCAwKXtcbiAgICAgICAgcmV0dXJuIFpfREFUQV9FUlJPUjtcbiAgICAgIH1cbiAgICB9XG4gICAgaWYgKCh5IC09IHRoaXMuY1tpXSkgPCAwKXtcbiAgICAgIHJldHVybiBaX0RBVEFfRVJST1I7XG4gICAgfVxuICAgIHRoaXMuY1tpXSArPSB5O1xuXG4gICAgLy8gR2VuZXJhdGUgc3RhcnRpbmcgb2Zmc2V0cyBpbnRvIHRoZSB2YWx1ZSB0YWJsZSBmb3IgZWFjaCBsZW5ndGhcbiAgICB0aGlzLnhbMV0gPSBqID0gMDtcbiAgICBwID0gMTsgIHhwID0gMjtcbiAgICB3aGlsZSAoLS1pIT0wKSB7ICAgICAgICAgICAgICAgICAvLyBub3RlIHRoYXQgaSA9PSBnIGZyb20gYWJvdmVcbiAgICAgIHRoaXMueFt4cF0gPSAoaiArPSB0aGlzLmNbcF0pO1xuICAgICAgeHArKztcbiAgICAgIHArKztcbiAgICB9XG5cbiAgICAvLyBNYWtlIGEgdGFibGUgb2YgdmFsdWVzIGluIG9yZGVyIG9mIGJpdCBsZW5ndGhzXG4gICAgaSA9IDA7IHAgPSAwO1xuICAgIGRvIHtcbiAgICAgIGlmICgoaiA9IGJbYmluZGV4K3BdKSAhPSAwKXtcbiAgICAgICAgdGhpcy52W3RoaXMueFtqXSsrXSA9IGk7XG4gICAgICB9XG4gICAgICBwKys7XG4gICAgfVxuICAgIHdoaWxlICgrK2kgPCBuKTtcbiAgICBuID0gdGhpcy54W2ddOyAgICAgICAgICAgICAgICAgICAgIC8vIHNldCBuIHRvIGxlbmd0aCBvZiB2XG5cbiAgICAvLyBHZW5lcmF0ZSB0aGUgSHVmZm1hbiBjb2RlcyBhbmQgZm9yIGVhY2gsIG1ha2UgdGhlIHRhYmxlIGVudHJpZXNcbiAgICB0aGlzLnhbMF0gPSBpID0gMDsgICAgICAgICAgICAgICAgIC8vIGZpcnN0IEh1ZmZtYW4gY29kZSBpcyB6ZXJvXG4gICAgcCA9IDA7ICAgICAgICAgICAgICAgICAgICAgICAgLy8gZ3JhYiB2YWx1ZXMgaW4gYml0IG9yZGVyXG4gICAgaCA9IC0xOyAgICAgICAgICAgICAgICAgICAgICAgLy8gbm8gdGFibGVzIHlldC0tbGV2ZWwgLTFcbiAgICB3ID0gLWw7ICAgICAgICAgICAgICAgICAgICAgICAvLyBiaXRzIGRlY29kZWQgPT0gKGwgKiBoKVxuICAgIHRoaXMudVswXSA9IDA7ICAgICAgICAgICAgICAgICAgICAgLy8ganVzdCB0byBrZWVwIGNvbXBpbGVycyBoYXBweVxuICAgIHEgPSAwOyAgICAgICAgICAgICAgICAgICAgICAgIC8vIGRpdHRvXG4gICAgeiA9IDA7ICAgICAgICAgICAgICAgICAgICAgICAgLy8gZGl0dG9cblxuICAgIC8vIGdvIHRocm91Z2ggdGhlIGJpdCBsZW5ndGhzIChrIGFscmVhZHkgaXMgYml0cyBpbiBzaG9ydGVzdCBjb2RlKVxuICAgIGZvciAoOyBrIDw9IGc7IGsrKyl7XG4gICAgICBhID0gdGhpcy5jW2tdO1xuICAgICAgd2hpbGUgKGEtLSE9MCl7XG5cdC8vIGhlcmUgaSBpcyB0aGUgSHVmZm1hbiBjb2RlIG9mIGxlbmd0aCBrIGJpdHMgZm9yIHZhbHVlICpwXG5cdC8vIG1ha2UgdGFibGVzIHVwIHRvIHJlcXVpcmVkIGxldmVsXG4gICAgICAgIHdoaWxlIChrID4gdyArIGwpe1xuICAgICAgICAgIGgrKztcbiAgICAgICAgICB3ICs9IGw7ICAgICAgICAgICAgICAgICAvLyBwcmV2aW91cyB0YWJsZSBhbHdheXMgbCBiaXRzXG5cdCAgLy8gY29tcHV0ZSBtaW5pbXVtIHNpemUgdGFibGUgbGVzcyB0aGFuIG9yIGVxdWFsIHRvIGwgYml0c1xuICAgICAgICAgIHogPSBnIC0gdztcbiAgICAgICAgICB6ID0gKHogPiBsKSA/IGwgOiB6OyAgICAgICAgLy8gdGFibGUgc2l6ZSB1cHBlciBsaW1pdFxuICAgICAgICAgIGlmKChmPTE8PChqPWstdykpPmErMSl7ICAgICAvLyB0cnkgYSBrLXcgYml0IHRhYmxlXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIC8vIHRvbyBmZXcgY29kZXMgZm9yIGstdyBiaXQgdGFibGVcbiAgICAgICAgICAgIGYgLT0gYSArIDE7ICAgICAgICAgICAgICAgLy8gZGVkdWN0IGNvZGVzIGZyb20gcGF0dGVybnMgbGVmdFxuICAgICAgICAgICAgeHAgPSBrO1xuICAgICAgICAgICAgaWYoaiA8IHope1xuICAgICAgICAgICAgICB3aGlsZSAoKytqIDwgeil7ICAgICAgICAvLyB0cnkgc21hbGxlciB0YWJsZXMgdXAgdG8geiBiaXRzXG4gICAgICAgICAgICAgICAgaWYoKGYgPDw9IDEpIDw9IHRoaXMuY1srK3hwXSlcbiAgICAgICAgICAgICAgICAgIGJyZWFrOyAgICAgICAgICAgICAgLy8gZW5vdWdoIGNvZGVzIHRvIHVzZSB1cCBqIGJpdHNcbiAgICAgICAgICAgICAgICBmIC09IHRoaXMuY1t4cF07ICAgICAgICAgICAvLyBlbHNlIGRlZHVjdCBjb2RlcyBmcm9tIHBhdHRlcm5zXG4gICAgICAgICAgICAgIH1cblx0ICAgIH1cbiAgICAgICAgICB9XG4gICAgICAgICAgeiA9IDEgPDwgajsgICAgICAgICAgICAgICAgIC8vIHRhYmxlIGVudHJpZXMgZm9yIGotYml0IHRhYmxlXG5cblx0ICAvLyBhbGxvY2F0ZSBuZXcgdGFibGVcbiAgICAgICAgICBpZiAodGhpcy5oblswXSArIHogPiBNQU5ZKXsgICAgICAgLy8gKG5vdGU6IGRvZXNuJ3QgbWF0dGVyIGZvciBmaXhlZClcbiAgICAgICAgICAgIHJldHVybiBaX0RBVEFfRVJST1I7ICAgICAgIC8vIG92ZXJmbG93IG9mIE1BTllcbiAgICAgICAgICB9XG4gICAgICAgICAgdGhpcy51W2hdID0gcSA9IC8qaHArKi8gdGhpcy5oblswXTsgICAvLyBERUJVR1xuICAgICAgICAgIHRoaXMuaG5bMF0gKz0gejtcbiBcblx0ICAvLyBjb25uZWN0IHRvIGxhc3QgdGFibGUsIGlmIHRoZXJlIGlzIG9uZVxuXHQgIGlmKGghPTApe1xuICAgICAgICAgICAgdGhpcy54W2hdPWk7ICAgICAgICAgICAvLyBzYXZlIHBhdHRlcm4gZm9yIGJhY2tpbmcgdXBcbiAgICAgICAgICAgIHRoaXMuclswXT1qOyAgICAgLy8gYml0cyBpbiB0aGlzIHRhYmxlXG4gICAgICAgICAgICB0aGlzLnJbMV09bDsgICAgIC8vIGJpdHMgdG8gZHVtcCBiZWZvcmUgdGhpcyB0YWJsZVxuICAgICAgICAgICAgaj1pPj4+KHcgLSBsKTtcbiAgICAgICAgICAgIHRoaXMuclsyXSA9IChxIC0gdGhpcy51W2gtMV0gLSBqKTsgICAgICAgICAgICAgICAvLyBvZmZzZXQgdG8gdGhpcyB0YWJsZVxuICAgICAgICAgICAgYXJyYXlDb3B5KHRoaXMuciwgMCwgaHAsICh0aGlzLnVbaC0xXStqKSozLCAzKTsgLy8gY29ubmVjdCB0byBsYXN0IHRhYmxlXG4gICAgICAgICAgfVxuICAgICAgICAgIGVsc2V7XG4gICAgICAgICAgICB0WzBdID0gcTsgICAgICAgICAgICAgICAvLyBmaXJzdCB0YWJsZSBpcyByZXR1cm5lZCByZXN1bHRcblx0ICB9XG4gICAgICAgIH1cblxuXHQvLyBzZXQgdXAgdGFibGUgZW50cnkgaW4gclxuICAgICAgICB0aGlzLnJbMV0gPSAoayAtIHcpO1xuICAgICAgICBpZiAocCA+PSBuKXtcbiAgICAgICAgICB0aGlzLnJbMF0gPSAxMjggKyA2NDsgICAgICAvLyBvdXQgb2YgdmFsdWVzLS1pbnZhbGlkIGNvZGVcblx0fVxuICAgICAgICBlbHNlIGlmICh2W3BdIDwgcyl7XG4gICAgICAgICAgdGhpcy5yWzBdID0gKHRoaXMudltwXSA8IDI1NiA/IDAgOiAzMiArIDY0KTsgIC8vIDI1NiBpcyBlbmQtb2YtYmxvY2tcbiAgICAgICAgICB0aGlzLnJbMl0gPSB0aGlzLnZbcCsrXTsgICAgICAgICAgLy8gc2ltcGxlIGNvZGUgaXMganVzdCB0aGUgdmFsdWVcbiAgICAgICAgfVxuICAgICAgICBlbHNle1xuICAgICAgICAgIHRoaXMuclswXT0oZVt0aGlzLnZbcF0tc10rMTYrNjQpOyAvLyBub24tc2ltcGxlLS1sb29rIHVwIGluIGxpc3RzXG4gICAgICAgICAgdGhpcy5yWzJdPWRbdGhpcy52W3ArK10gLSBzXTtcbiAgICAgICAgfVxuXG4gICAgICAgIC8vIGZpbGwgY29kZS1saWtlIGVudHJpZXMgd2l0aCByXG4gICAgICAgIGY9MTw8KGstdyk7XG4gICAgICAgIGZvciAoaj1pPj4+dztqPHo7ais9Zil7XG4gICAgICAgICAgYXJyYXlDb3B5KHRoaXMuciwgMCwgaHAsIChxK2opKjMsIDMpO1xuXHR9XG5cblx0Ly8gYmFja3dhcmRzIGluY3JlbWVudCB0aGUgay1iaXQgY29kZSBpXG4gICAgICAgIGZvciAoaiA9IDEgPDwgKGsgLSAxKTsgKGkgJiBqKSE9MDsgaiA+Pj49IDEpe1xuICAgICAgICAgIGkgXj0gajtcblx0fVxuICAgICAgICBpIF49IGo7XG5cblx0Ly8gYmFja3VwIG92ZXIgZmluaXNoZWQgdGFibGVzXG4gICAgICAgIG1hc2sgPSAoMSA8PCB3KSAtIDE7ICAgICAgLy8gbmVlZGVkIG9uIEhQLCBjYyAtTyBidWdcbiAgICAgICAgd2hpbGUgKChpICYgbWFzaykgIT0gdGhpcy54W2hdKXtcbiAgICAgICAgICBoLS07ICAgICAgICAgICAgICAgICAgICAvLyBkb24ndCBuZWVkIHRvIHVwZGF0ZSBxXG4gICAgICAgICAgdyAtPSBsO1xuICAgICAgICAgIG1hc2sgPSAoMSA8PCB3KSAtIDE7XG4gICAgICAgIH1cbiAgICAgIH1cbiAgICB9XG4gICAgLy8gUmV0dXJuIFpfQlVGX0VSUk9SIGlmIHdlIHdlcmUgZ2l2ZW4gYW4gaW5jb21wbGV0ZSB0YWJsZVxuICAgIHJldHVybiB5ICE9IDAgJiYgZyAhPSAxID8gWl9CVUZfRVJST1IgOiBaX09LO1xufVxuXG5JbmZUcmVlLnByb3RvdHlwZS5pbmZsYXRlX3RyZWVzX2JpdHMgPSBmdW5jdGlvbihjLCBiYiwgdGIsIGhwLCB6KSB7XG4gICAgdmFyIHJlc3VsdDtcbiAgICB0aGlzLmluaXRXb3JrQXJlYSgxOSk7XG4gICAgdGhpcy5oblswXT0wO1xuICAgIHJlc3VsdCA9IHRoaXMuaHVmdF9idWlsZChjLCAwLCAxOSwgMTksIG51bGwsIG51bGwsIHRiLCBiYiwgaHAsIHRoaXMuaG4sIHRoaXMudik7XG5cbiAgICBpZihyZXN1bHQgPT0gWl9EQVRBX0VSUk9SKXtcbiAgICAgIHoubXNnID0gXCJvdmVyc3Vic2NyaWJlZCBkeW5hbWljIGJpdCBsZW5ndGhzIHRyZWVcIjtcbiAgICB9XG4gICAgZWxzZSBpZihyZXN1bHQgPT0gWl9CVUZfRVJST1IgfHwgYmJbMF0gPT0gMCl7XG4gICAgICB6Lm1zZyA9IFwiaW5jb21wbGV0ZSBkeW5hbWljIGJpdCBsZW5ndGhzIHRyZWVcIjtcbiAgICAgIHJlc3VsdCA9IFpfREFUQV9FUlJPUjtcbiAgICB9XG4gICAgcmV0dXJuIHJlc3VsdDtcbn1cblxuSW5mVHJlZS5wcm90b3R5cGUuaW5mbGF0ZV90cmVlc19keW5hbWljID0gZnVuY3Rpb24obmwsIG5kLCBjLCBibCwgYmQsIHRsLCB0ZCwgaHAsIHopIHtcbiAgICB2YXIgcmVzdWx0O1xuXG4gICAgLy8gYnVpbGQgbGl0ZXJhbC9sZW5ndGggdHJlZVxuICAgIHRoaXMuaW5pdFdvcmtBcmVhKDI4OCk7XG4gICAgdGhpcy5oblswXT0wO1xuICAgIHJlc3VsdCA9IHRoaXMuaHVmdF9idWlsZChjLCAwLCBubCwgMjU3LCBjcGxlbnMsIGNwbGV4dCwgdGwsIGJsLCBocCwgdGhpcy5obiwgdGhpcy52KTtcbiAgICBpZiAocmVzdWx0ICE9IFpfT0sgfHwgYmxbMF0gPT0gMCl7XG4gICAgICBpZihyZXN1bHQgPT0gWl9EQVRBX0VSUk9SKXtcbiAgICAgICAgei5tc2cgPSBcIm92ZXJzdWJzY3JpYmVkIGxpdGVyYWwvbGVuZ3RoIHRyZWVcIjtcbiAgICAgIH1cbiAgICAgIGVsc2UgaWYgKHJlc3VsdCAhPSBaX01FTV9FUlJPUil7XG4gICAgICAgIHoubXNnID0gXCJpbmNvbXBsZXRlIGxpdGVyYWwvbGVuZ3RoIHRyZWVcIjtcbiAgICAgICAgcmVzdWx0ID0gWl9EQVRBX0VSUk9SO1xuICAgICAgfVxuICAgICAgcmV0dXJuIHJlc3VsdDtcbiAgICB9XG5cbiAgICAvLyBidWlsZCBkaXN0YW5jZSB0cmVlXG4gICAgdGhpcy5pbml0V29ya0FyZWEoMjg4KTtcbiAgICByZXN1bHQgPSB0aGlzLmh1ZnRfYnVpbGQoYywgbmwsIG5kLCAwLCBjcGRpc3QsIGNwZGV4dCwgdGQsIGJkLCBocCwgdGhpcy5obiwgdGhpcy52KTtcblxuICAgIGlmIChyZXN1bHQgIT0gWl9PSyB8fCAoYmRbMF0gPT0gMCAmJiBubCA+IDI1Nykpe1xuICAgICAgaWYgKHJlc3VsdCA9PSBaX0RBVEFfRVJST1Ipe1xuICAgICAgICB6Lm1zZyA9IFwib3ZlcnN1YnNjcmliZWQgZGlzdGFuY2UgdHJlZVwiO1xuICAgICAgfVxuICAgICAgZWxzZSBpZiAocmVzdWx0ID09IFpfQlVGX0VSUk9SKSB7XG4gICAgICAgIHoubXNnID0gXCJpbmNvbXBsZXRlIGRpc3RhbmNlIHRyZWVcIjtcbiAgICAgICAgcmVzdWx0ID0gWl9EQVRBX0VSUk9SO1xuICAgICAgfVxuICAgICAgZWxzZSBpZiAocmVzdWx0ICE9IFpfTUVNX0VSUk9SKXtcbiAgICAgICAgei5tc2cgPSBcImVtcHR5IGRpc3RhbmNlIHRyZWUgd2l0aCBsZW5ndGhzXCI7XG4gICAgICAgIHJlc3VsdCA9IFpfREFUQV9FUlJPUjtcbiAgICAgIH1cbiAgICAgIHJldHVybiByZXN1bHQ7XG4gICAgfVxuXG4gICAgcmV0dXJuIFpfT0s7XG59XG4vKlxuICBzdGF0aWMgaW50IGluZmxhdGVfdHJlZXNfZml4ZWQoaW50W10gYmwsICAvL2xpdGVyYWwgZGVzaXJlZC9hY3R1YWwgYml0IGRlcHRoXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbnRbXSBiZCwgIC8vZGlzdGFuY2UgZGVzaXJlZC9hY3R1YWwgYml0IGRlcHRoXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBpbnRbXVtdIHRsLC8vbGl0ZXJhbC9sZW5ndGggdHJlZSByZXN1bHRcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGludFtdW10gdGQsLy9kaXN0YW5jZSB0cmVlIHJlc3VsdCBcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFpTdHJlYW0geiAgLy9mb3IgbWVtb3J5IGFsbG9jYXRpb25cblx0XHRcdFx0ICl7XG5cbiovXG5cbmZ1bmN0aW9uIGluZmxhdGVfdHJlZXNfZml4ZWQoYmwsIGJkLCB0bCwgdGQsIHopIHtcbiAgICBibFswXT1maXhlZF9ibDtcbiAgICBiZFswXT1maXhlZF9iZDtcbiAgICB0bFswXT1maXhlZF90bDtcbiAgICB0ZFswXT1maXhlZF90ZDtcbiAgICByZXR1cm4gWl9PSztcbn1cblxuSW5mVHJlZS5wcm90b3R5cGUuaW5pdFdvcmtBcmVhID0gZnVuY3Rpb24odnNpemUpe1xuICAgIGlmKHRoaXMuaG49PW51bGwpe1xuICAgICAgICB0aGlzLmhuPW5ldyBJbnQzMkFycmF5KDEpO1xuICAgICAgICB0aGlzLnY9bmV3IEludDMyQXJyYXkodnNpemUpO1xuICAgICAgICB0aGlzLmM9bmV3IEludDMyQXJyYXkoQk1BWCsxKTtcbiAgICAgICAgdGhpcy5yPW5ldyBJbnQzMkFycmF5KDMpO1xuICAgICAgICB0aGlzLnU9bmV3IEludDMyQXJyYXkoQk1BWCk7XG4gICAgICAgIHRoaXMueD1uZXcgSW50MzJBcnJheShCTUFYKzEpO1xuICAgIH1cbiAgICBpZih0aGlzLnYubGVuZ3RoPHZzaXplKXsgXG4gICAgICAgIHRoaXMudj1uZXcgSW50MzJBcnJheSh2c2l6ZSk7IFxuICAgIH1cbiAgICBmb3IodmFyIGk9MDsgaTx2c2l6ZTsgaSsrKXt0aGlzLnZbaV09MDt9XG4gICAgZm9yKHZhciBpPTA7IGk8Qk1BWCsxOyBpKyspe3RoaXMuY1tpXT0wO31cbiAgICBmb3IodmFyIGk9MDsgaTwzOyBpKyspe3RoaXMucltpXT0wO31cbi8vICBmb3IoaW50IGk9MDsgaTxCTUFYOyBpKyspe3VbaV09MDt9XG4gICAgYXJyYXlDb3B5KHRoaXMuYywgMCwgdGhpcy51LCAwLCBCTUFYKTtcbi8vICBmb3IoaW50IGk9MDsgaTxCTUFYKzE7IGkrKyl7eFtpXT0wO31cbiAgICBhcnJheUNvcHkodGhpcy5jLCAwLCB0aGlzLngsIDAsIEJNQVgrMSk7XG59XG5cbnZhciB0ZXN0QXJyYXkgPSBuZXcgVWludDhBcnJheSgxKTtcbnZhciBoYXNTdWJhcnJheSA9ICh0eXBlb2YgdGVzdEFycmF5LnN1YmFycmF5ID09PSAnZnVuY3Rpb24nKTtcbnZhciBoYXNTbGljZSA9IGZhbHNlOyAvKiAodHlwZW9mIHRlc3RBcnJheS5zbGljZSA9PT0gJ2Z1bmN0aW9uJyk7ICovIC8vIENocm9tZSBzbGljZSBwZXJmb3JtYW5jZSBpcyBzbyBkaXJlIHRoYXQgd2UncmUgY3VycmVudGx5IG5vdCB1c2luZyBpdC4uLlxuXG5mdW5jdGlvbiBhcnJheUNvcHkoc3JjLCBzcmNPZmZzZXQsIGRlc3QsIGRlc3RPZmZzZXQsIGNvdW50KSB7XG4gICAgaWYgKGNvdW50ID09IDApIHtcbiAgICAgICAgcmV0dXJuO1xuICAgIH0gXG4gICAgaWYgKCFzcmMpIHtcbiAgICAgICAgdGhyb3cgXCJVbmRlZiBzcmNcIjtcbiAgICB9IGVsc2UgaWYgKCFkZXN0KSB7XG4gICAgICAgIHRocm93IFwiVW5kZWYgZGVzdFwiO1xuICAgIH1cblxuICAgIGlmIChzcmNPZmZzZXQgPT0gMCAmJiBjb3VudCA9PSBzcmMubGVuZ3RoKSB7XG4gICAgICAgIGFycmF5Q29weV9mYXN0KHNyYywgZGVzdCwgZGVzdE9mZnNldCk7XG4gICAgfSBlbHNlIGlmIChoYXNTdWJhcnJheSkge1xuICAgICAgICBhcnJheUNvcHlfZmFzdChzcmMuc3ViYXJyYXkoc3JjT2Zmc2V0LCBzcmNPZmZzZXQgKyBjb3VudCksIGRlc3QsIGRlc3RPZmZzZXQpOyBcbiAgICB9IGVsc2UgaWYgKHNyYy5CWVRFU19QRVJfRUxFTUVOVCA9PSAxICYmIGNvdW50ID4gMTAwKSB7XG4gICAgICAgIGFycmF5Q29weV9mYXN0KG5ldyBVaW50OEFycmF5KHNyYy5idWZmZXIsIHNyYy5ieXRlT2Zmc2V0ICsgc3JjT2Zmc2V0LCBjb3VudCksIGRlc3QsIGRlc3RPZmZzZXQpO1xuICAgIH0gZWxzZSB7IFxuICAgICAgICBhcnJheUNvcHlfc2xvdyhzcmMsIHNyY09mZnNldCwgZGVzdCwgZGVzdE9mZnNldCwgY291bnQpO1xuICAgIH1cblxufVxuXG5mdW5jdGlvbiBhcnJheUNvcHlfc2xvdyhzcmMsIHNyY09mZnNldCwgZGVzdCwgZGVzdE9mZnNldCwgY291bnQpIHtcblxuICAgIC8vIGRsb2coJ19zbG93IGNhbGw6IHNyY09mZnNldD0nICsgc3JjT2Zmc2V0ICsgJzsgZGVzdE9mZnNldD0nICsgZGVzdE9mZnNldCArICc7IGNvdW50PScgKyBjb3VudCk7XG5cbiAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBjb3VudDsgKytpKSB7XG4gICAgICAgIGRlc3RbZGVzdE9mZnNldCArIGldID0gc3JjW3NyY09mZnNldCArIGldO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gYXJyYXlDb3B5X2Zhc3Qoc3JjLCBkZXN0LCBkZXN0T2Zmc2V0KSB7XG4gICAgZGVzdC5zZXQoc3JjLCBkZXN0T2Zmc2V0KTtcbn1cblxuXG4gIC8vIGxhcmdlc3QgcHJpbWUgc21hbGxlciB0aGFuIDY1NTM2XG52YXIgQURMRVJfQkFTRT02NTUyMTsgXG4gIC8vIE5NQVggaXMgdGhlIGxhcmdlc3QgbiBzdWNoIHRoYXQgMjU1bihuKzEpLzIgKyAobisxKShCQVNFLTEpIDw9IDJeMzItMVxudmFyIEFETEVSX05NQVg9NTU1MjtcblxuZnVuY3Rpb24gYWRsZXIzMihhZGxlciwgLyogYnl0ZVtdICovIGJ1ZiwgIGluZGV4LCBsZW4pe1xuICAgIGlmKGJ1ZiA9PSBudWxsKXsgcmV0dXJuIDE7IH1cblxuICAgIHZhciBzMT1hZGxlciYweGZmZmY7XG4gICAgdmFyIHMyPShhZGxlcj4+MTYpJjB4ZmZmZjtcbiAgICB2YXIgaztcblxuICAgIHdoaWxlKGxlbiA+IDApIHtcbiAgICAgIGs9bGVuPEFETEVSX05NQVg/bGVuOkFETEVSX05NQVg7XG4gICAgICBsZW4tPWs7XG4gICAgICB3aGlsZShrPj0xNil7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgczErPWJ1ZltpbmRleCsrXSYweGZmOyBzMis9czE7XG4gICAgICAgIHMxKz1idWZbaW5kZXgrK10mMHhmZjsgczIrPXMxO1xuICAgICAgICBrLT0xNjtcbiAgICAgIH1cbiAgICAgIGlmKGshPTApe1xuICAgICAgICBkb3tcbiAgICAgICAgICBzMSs9YnVmW2luZGV4KytdJjB4ZmY7IHMyKz1zMTtcbiAgICAgICAgfVxuICAgICAgICB3aGlsZSgtLWshPTApO1xuICAgICAgfVxuICAgICAgczElPUFETEVSX0JBU0U7XG4gICAgICBzMiU9QURMRVJfQkFTRTtcbiAgICB9XG4gICAgcmV0dXJuIChzMjw8MTYpfHMxO1xufVxuXG5cblxuZnVuY3Rpb24ganN6bGliX2luZmxhdGVfYnVmZmVyKGJ1ZmZlciwgc3RhcnQsIGxlbmd0aCwgYWZ0ZXJVbmNPZmZzZXQpIHtcbiAgICBpZiAoIXN0YXJ0KSB7XG4gICAgICAgIGJ1ZmZlciA9IG5ldyBVaW50OEFycmF5KGJ1ZmZlcik7XG4gICAgfSBlbHNlIGlmICghbGVuZ3RoKSB7XG4gICAgICAgIGJ1ZmZlciA9IG5ldyBVaW50OEFycmF5KGJ1ZmZlciwgc3RhcnQsIGJ1ZmZlci5ieXRlTGVuZ3RoIC0gc3RhcnQpO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIGJ1ZmZlciA9IG5ldyBVaW50OEFycmF5KGJ1ZmZlciwgc3RhcnQsIGxlbmd0aCk7XG4gICAgfVxuXG4gICAgdmFyIHogPSBuZXcgWlN0cmVhbSgpO1xuICAgIHouaW5mbGF0ZUluaXQoREVGX1dCSVRTLCB0cnVlKTtcbiAgICB6Lm5leHRfaW4gPSBidWZmZXI7XG4gICAgei5uZXh0X2luX2luZGV4ID0gMDtcbiAgICB6LmF2YWlsX2luID0gYnVmZmVyLmxlbmd0aDtcblxuICAgIHZhciBvQmxvY2tMaXN0ID0gW107XG4gICAgdmFyIHRvdGFsU2l6ZSA9IDA7XG4gICAgd2hpbGUgKHRydWUpIHtcbiAgICAgICAgdmFyIG9idWYgPSBuZXcgVWludDhBcnJheSgzMjAwMCk7XG4gICAgICAgIHoubmV4dF9vdXQgPSBvYnVmO1xuICAgICAgICB6Lm5leHRfb3V0X2luZGV4ID0gMDtcbiAgICAgICAgei5hdmFpbF9vdXQgPSBvYnVmLmxlbmd0aDtcbiAgICAgICAgdmFyIHN0YXR1cyA9IHouaW5mbGF0ZShaX05PX0ZMVVNIKTtcbiAgICAgICAgaWYgKHN0YXR1cyAhPSBaX09LICYmIHN0YXR1cyAhPSBaX1NUUkVBTV9FTkQgJiYgc3RhdHVzICE9IFpfQlVGX0VSUk9SKSB7XG4gICAgICAgICAgICB0aHJvdyB6Lm1zZztcbiAgICAgICAgfVxuICAgICAgICBpZiAoei5hdmFpbF9vdXQgIT0gMCkge1xuICAgICAgICAgICAgdmFyIG5ld29iID0gbmV3IFVpbnQ4QXJyYXkob2J1Zi5sZW5ndGggLSB6LmF2YWlsX291dCk7XG4gICAgICAgICAgICBhcnJheUNvcHkob2J1ZiwgMCwgbmV3b2IsIDAsIChvYnVmLmxlbmd0aCAtIHouYXZhaWxfb3V0KSk7XG4gICAgICAgICAgICBvYnVmID0gbmV3b2I7XG4gICAgICAgIH1cbiAgICAgICAgb0Jsb2NrTGlzdC5wdXNoKG9idWYpO1xuICAgICAgICB0b3RhbFNpemUgKz0gb2J1Zi5sZW5ndGg7XG4gICAgICAgIGlmIChzdGF0dXMgPT0gWl9TVFJFQU1fRU5EIHx8IHN0YXR1cyA9PSBaX0JVRl9FUlJPUikge1xuICAgICAgICAgICAgYnJlYWs7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBpZiAoYWZ0ZXJVbmNPZmZzZXQpIHtcbiAgICAgICAgYWZ0ZXJVbmNPZmZzZXRbMF0gPSAoc3RhcnQgfHwgMCkgKyB6Lm5leHRfaW5faW5kZXg7XG4gICAgfVxuXG4gICAgaWYgKG9CbG9ja0xpc3QubGVuZ3RoID09IDEpIHtcbiAgICAgICAgcmV0dXJuIG9CbG9ja0xpc3RbMF0uYnVmZmVyO1xuICAgIH0gZWxzZSB7XG4gICAgICAgIHZhciBvdXQgPSBuZXcgVWludDhBcnJheSh0b3RhbFNpemUpO1xuICAgICAgICB2YXIgY3Vyc29yID0gMDtcbiAgICAgICAgZm9yICh2YXIgaSA9IDA7IGkgPCBvQmxvY2tMaXN0Lmxlbmd0aDsgKytpKSB7XG4gICAgICAgICAgICB2YXIgYiA9IG9CbG9ja0xpc3RbaV07XG4gICAgICAgICAgICBhcnJheUNvcHkoYiwgMCwgb3V0LCBjdXJzb3IsIGIubGVuZ3RoKTtcbiAgICAgICAgICAgIGN1cnNvciArPSBiLmxlbmd0aDtcbiAgICAgICAgfVxuICAgICAgICByZXR1cm4gb3V0LmJ1ZmZlcjtcbiAgICB9XG59XG5cbmlmICh0eXBlb2YobW9kdWxlKSAhPT0gJ3VuZGVmaW5lZCcpIHtcbiAgbW9kdWxlLmV4cG9ydHMgPSB7XG4gICAgaW5mbGF0ZUJ1ZmZlcjoganN6bGliX2luZmxhdGVfYnVmZmVyLFxuICAgIGFycmF5Q29weTogYXJyYXlDb3B5XG4gIH07XG59XG4iXX0=
