class CramData {
	constructor(fil, option){
		this.file_ = fil;
		this.option_ = (option !== undefined)? option: {};
		this.cramReader_ = new CramReader(fil[0], fil[1], option.localFlg, fil[2], fil[3], "js/cram-reader-worker.min.js");
	}

	readWaitReader(chr, start, end, callback, reject, option) {
		if(option === undefined) option = {};
		if(option.timeout === undefined) option.timeout = 300;

		var func = (reads)=>{
			var result = [];
			for (var read of reads) {
				var arr = [
					read.position,
					read.mappingQuality,
					read.bf_,
					read.readName,
					read.cigar,
					read.positionEnd,
					read.seq,
					read.qualityScore,
					read.mateRefId,
					read.matePos - 1,
					read.templateSize,
					read.samString,
					read.cigarLn,
					read.cigarOp
				];
				result.push(arr);
			}
			callback(result);
		};

		try{
			var timeoutID = setTimeout(() => {throw "timeout cram.js"});
			this.cramReader_.getRecords(chr, start, end, func);
			clearTimeout(timeoutID);
		}catch (e){
			reject;
		}
	}
}

class CramReader{constructor(r,e,t,i=undefined,R=undefined,s="cram-reader-worker.min.js"){if(!window.Worker){throw"Web Workers API is needed"}if(!r||!e){throw"Files are Falsy"}this.t=new Map;this.i=this.R(s);this.i.cram_reader=this;this.i.onmessage=function(r){this.cram_reader.h(r)};this.o("init",[r,e,t,i,R],undefined)}getRecords(r,e,t,i){this.o("read",[r,e,t],i)}setOnerror(r){this.i.onerror=r}R(r){try{var e=window.location.href.replace(/\/[^\/]*$/g,"/");var t=e+r;var i=new Blob(["importScripts('"+t+"');"],{type:"text/javascript"});var R=URL.createObjectURL(i);return new Worker(R)}catch{return new Worker(r)}}o(r,e,t){var i;while(true){var R=this.v();if(!this.t.has(R)){i=R;break}}this.t.set(i,t);this.i.postMessage([i,r,e])}h(r){var e=this.t.get(r.data[0]);this.t.delete(r.data[0]);if(e){e(r.data[1])}}v(){var r="RRRRRRRR-RRRR-4RRR-rRRR-RRRRRRRRRRRR";for(var e=0;e<r.length;e++){switch(r.charAt(e)){case"R":r=r.replace("R",Math.floor(Math.random()*16).toString(16));break;case"r":r=r.replace("r",(Math.floor(Math.random()*4)+8).toString(16));break}}return r}}
