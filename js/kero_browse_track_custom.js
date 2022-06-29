//nameは自分のID(オブジェクトごとユニークにする)
var WgCram = function(name, dispName, bamUrl, option) {
	this.option = (option !== undefined)? option: {};
	if(this.option.uriDirFlg === undefined) this.option.uriDirFlg = true;
	
	this.name = name;
	this.dispName = dispName;
	
	this.charPx = 10;
	
	//this.bamUrl = bamUrl;
	
	this.minHeight = 30;
	this.y;
	
	//配列データのデータ取得先
	this.getSeqData = (this.option.seqUrl === undefined)? {
		"url": "https://dbtss.hgc.jp/cgi-bin/dbtss_db_json.cgi/9606",
		"param": "SEE=1&UID=2"
	}: this.option.seqUrl;
	
	this.powMax = 1;
	this.margin = 1;
	//this.eachStep = 15;
	this.eachStep = 10;
	//readの表示はここまでですの表示を入れるためのスペース(mexHeightに含まれる)
	this.limitShowHeight = 5;
	this.maxHeight = this.eachStep * 30 + this.limitShowHeight;
	//readが積み重なりすぎていて全部表示できない領域
	this.overReg = {};
	this.overRegList = [];
	this.bam = new CramData(bamUrl, {localFlg: this.option.localFlg});
};
WgCram.prototype = new WgRoot();
//描画する(Y位置, 裏画面の幅, chr, 裏画面のstart, 裏画面のend)
//実際には裏画面を横に3つに区切ったうちの真ん中が表示される
WgCram.prototype.paint = function(y, width, chr, start, end, strand) {
	//複数回カウントを避ける
	var existId = {};
	var status = [];
	
	var pow = Math.floor(Math.LOG10E * Math.log((end - start + 1) / width));
	if(pow < 0) pow = 0;
	var reg = Math.pow(10, pow) * POW_REG;
	
	var binStart = Math.floor((start - 1) / reg);
	var binEnd = Math.floor((end - 1) / reg);
	
	this.imgObj.font = this.charPx + "px 'Helvetica'";
	
	if(pow > this.powMax) {
		var y1 = y;
		var y2 = y1 + this.height - 1;
		var y3 = (y1 + y2) / 2 + 3;

		this.imgObj.fillStyle = "#EEEEEE";
		this.imgObj.fillRect(0, y1, width, y2 - y1 + 1);
		
		var tooLargeStr = "Too large region for showing alignment";
		var strWidth = this.imgObj.measureText(tooLargeStr).width;
		this.imgObj.fillStyle = "#888888";
		x3 = (width - strWidth) / 2;
		this.imgObj.fillText(tooLargeStr, x3, y3);
		
		return [];
	}
	
	let pileup = {};
	
	for(var i = binStart; i <= binEnd; i ++) {
		if(
			this.ojson[pow] && this.ojson[pow][chr + "|" + i] &&
			this.ojson[pow][chr + "|" + i][this.name]
		) {
			var each = this.ojson[pow][chr + "|" + i][this.name];
			var read = each.reads;
			
			if(!read) continue;
			
			for(var j = 0; j < read.length; j ++) {
				var readS = read[j].pos;
				var readE = read[j].posEnd;
				var flag = read[j].flag;
				var stepNo = read[j].step;
				var cigar = read[j].cigar;
				var cigarOp = read[j].cigarOp;
				var cigarLn = read[j].cigarLn;
				var seq = read[j].seq;
				var seqPoi = 0;
				
				//stepNoがundefinedもしくは0未満のデータは表示しない
				if(start <= readE && readS <= end) {
					if(existId[read[j].id]) {
						
						continue;
					} else {
						existId[read[j].id] = true;
					}
					var x1 = (readS - start) * (width - 1) / (end - start + 1);
					var x2 = (readE - start + 1) * (width - 1) / (end - start + 1);
					if(strand == "-") {var tmp = width - 1 - x1; x1 = width - 1 - x2; x2 = tmp;}
					
					var y1 = y + this.eachStep * stepNo;
					var y2 = y1 + this.eachStep - 3;
					if(stepNo >= 0 && false) {
					//簡易表示
						this.imgObj.fillStyle = (flag & 16)? "#BBBBFF": "#FFBBBB";
						this.imgObj.fillRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
					} else {
						let inPosList = [];
						let nowPos = readS;
						//let ops = cigar.split(/\d+/);
						//let lns = cigar.split(/\D/);
						let ops = cigarOp;
						let lns = cigarLn;
						for(let k = 1; k < ops.length; k ++) {
							let op = ops[k];
							let ln = parseInt(lns[k - 1]);
							if(stepNo >= 0 && op == "I") {
								if(pow == 0) {
									let x3 = (nowPos - start) * (width - 1) / (end - start + 1);
									inPosList.push(x3);
								}
								seqPoi += ln;
							} else if(stepNo >= 0 && op == "N") {
								let x3 = (nowPos - start) * (width - 1) / (end - start + 1);
								let x4 = (nowPos + ln - start) * (width - 1) / (end - start + 1);
								if(strand == "-") {var tmp = width - 1 - x3; x3 = width - 1 - x4; x4 = tmp;}
								let y3 = (y1 + y2) / 2 - 1;
								this.imgObj.fillStyle = (flag & 16)? "#BBBBFF": "#FFBBBB";
								this.imgObj.fillRect(x3, y3, x4 - x3 + 1, 2);
								nowPos += ln;
							} else if(stepNo >= 0 && op == "H") {
							} else if(stepNo >= 0 && op == "S") {
								seqPoi += ln;
							} else if(stepNo >= 0 && op == "P") {
								alert("not supported cigar P");
							} else if(op == "M" || op == "X" || op == "=" || op == "D") {
								if(stepNo >= 0) {
									let x3 = (nowPos - start) * (width - 1) / (end - start + 1);
									let x4 = (nowPos + ln - start) * (width - 1) / (end - start + 1);
									if(strand == "-") {var tmp = width - 1 - x3; x3 = width - 1 - x4; x4 = tmp;}
									this.imgObj.fillStyle = (flag & 16)? "#BBBBFF": "#FFBBBB";
									if(op == "D") {
										this.imgObj.fillRect(x3, (y1 + y2) / 2 - 1, x4 - x3 + 1, 2);
										nowPos += ln;
										continue;
									} else {
										this.imgObj.fillRect(x3, y1, x4 - x3 + 1, y2 - y1 + 1);
									}
								}
								
								//mutationの表示
								var sSeq = "";
								var binAlnS = Math.floor((nowPos - 1) / reg);
								var strStart = (nowPos - 1) % reg;
								var nowLn = ln;
								while(sSeq.length < ln) {
									var rest = reg - strStart;
									var kn = (rest < nowLn)? rest: nowLn;
									if(
										this.ojson[pow][chr + "|" + binAlnS] && this.ojson[pow][chr + "|" + binAlnS][this.name] &&
										this.ojson[pow][chr + "|" + binAlnS][this.name].sequence
									) {
										var gSeq = this.ojson[pow][chr + "|" + binAlnS][this.name].sequence;
										sSeq += gSeq.substr(strStart, kn).toUpperCase();
									} else {
										sSeq += "*".repeat(kn);
									}
									binAlnS ++;
									strStart = 0;
									nowLn -= kn;
								}
								
								var qSeq = seq.substr(seqPoi, ln);
								
								var cpPos = nowPos;
								for(let l = 0; l < sSeq.length; l ++) {
									if(cpPos < start || end < cpPos) {
										cpPos ++;
										continue;
									}
									let sBase = sSeq.charAt(l);
									let qBase = qSeq.charAt(l);
									let char = qBase;
									if(strand == "-") {
										var tmp = char;
										if(char == "A") tmp = "T"; if(char == "a") tmp = "t";
										if(char == "C") tmp = "G"; if(char == "c") tmp = "g";
										if(char == "G") tmp = "C"; if(char == "g") tmp = "c";
										if(char == "T") tmp = "A"; if(char == "t") tmp = "a";
										char = tmp;
									}
									
									//if(pileup[cpPos] === undefined) pileup[cpPos] = {ref:sBase, pile:{}};
									//pileup[cpPos][char] ++;
									
									if(stepNo >= 0 && sBase != qBase && sBase != "*" && qBase != "-") {
										let x5 = (cpPos - start) * (width - 1) / (end - start + 1);
										let x6 = (cpPos - start + 1) * (width - 1) / (end - start + 1);
										if(strand == "-") {var tmp = width - 1 - x5; x5 = width - 1 - x6; x6 = tmp;}
										this.imgObj.fillStyle =
											(char == "A" || char == "a")? "#88FF88":
											(char == "C" || char == "c")? "#8888FF":
											(char == "G" || char == "g")? "#FF8800":
											(char == "T" || char == "t")? "#FF4488": "#AAAAAA";
										this.imgObj.fillRect(x5, y1, x6 - x5 + 1, y2 - y1 + 1);
										if(x6 - x5 > this.charPx && this.eachStep >= 10) {
											this.imgObj.fillStyle = "#000000";
											this.imgObj.fillText(char, (x5 + x6) / 2 - 2, (y1 + y2) / 2 + 4);
										}
									}
									
									cpPos ++;
								}
								//mutationの表示（ここまで）
								nowPos += ln;
								seqPoi += ln;
							}
						}
						this.imgObj.fillStyle = "#333333";
						for(let k = 0; k < inPosList.length; k ++) {
							var x3 = inPosList[k];
							if(strand == "-") {var x3 = width - 1 - x3;}
							this.imgObj.fillRect(x3 - 1, y1, 2, y2 - y1 + 1);
						}
					}
				}
			}
		} else {
			if(i >= 0) status.push(i);
			this.paintLoading(y, width, start, end, strand, i);
		}
	}
	
	//積み重なりが大きくて表示できない領域
	this.imgObj.fillStyle = "#AAAAAA";
	for(let regS in this.overReg) {
		var regE = this.overReg[regS];
		var x1 = (regS - start) * (width - 1) / (end - start + 1);
		var x2 = (regE - start + 1) * (width - 1) / (end - start + 1);
		if(strand == "-") {var tmp = width - 1 - x1; x1 = width - 1 - x2; x2 = tmp;}
		var y1 = y + this.maxHeight - this.limitShowHeight;
		var y2 = y1 + 3;
		this.imgObj.fillRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
	}
	
	return status;
};
//step情報が毎回付加更新されます！
WgCram.prototype.setHeight = function(width, chr, start, end) {
	//複数回表示を避ける
	var existId = {};
	
	var pow = Math.floor(Math.LOG10E * Math.log((end - start + 1) / width));
	if(pow < 0) pow = 0;
	var reg = Math.pow(10, pow) * POW_REG;
	
	if(pow > this.powMax) {
		this.height = this.minHeight;
		return;
	}
	
	var binStart = Math.floor((start - 1) / reg);
	var binEnd = Math.floor((end - 1) / reg);
	
	var step = [];
	var overReg = this.overReg;
	var overRegList = this.overRegList;
	var stepMax = 0;
	for(var i = binStart; i <= binEnd; i ++) {
		if(
			this.ojson[pow] && this.ojson[pow][chr + "|" + i] &&
			this.ojson[pow][chr + "|" + i][this.name] && this.ojson[pow][chr + "|" + i][this.name].reads
		) {
			var read = this.ojson[pow][chr + "|" + i][this.name].reads;
			for(var j = 0; j < read.length; j ++) {
				let stepNo;
				var readS = read[j].pos;
				var readE = read[j].posEnd;
				var flag = read[j].flag;
				var cigar = read[j].cigar;
				//cigarないデータはスキップ
				if(cigar == "") continue;
				
				if(start <= readE && readS <= end) {
					if(existId[read[j].id]) {
						read[j].step = undefined;	//stepNo = undefined
						continue;
					} else {
						existId[read[j].id] = true;
					}
					
					if(read[j].step !== undefined) {
					//file:///E:/Documents/work/tokyo_univ/2017/js_test/local_file_access/index.html#chr1:24822488-24887847:- "hg38.illumina.adrenal.1.bam"
					//高速化のため
						var k = read[j].step;
						if(k != -1) {
							if(step[k] === undefined) step[k] = [];
							if(stepMax < k) stepMax = k;
							
							var ovFlg = false;
							for(var l = 0; l < step[k].length; l ++) {
								var checkS = step[k][l][0];
								var checkE = step[k][l][1];
								if(readS - this.margin <= checkE && checkS <= readE + this.margin) {
									ovFlg = true;
									break;
								}
							}
							
							if(!ovFlg) {
								step[k].push([readS, readE]);
								stepNo = k;
							} else {
								for(k = 0; k < step.length; k ++) {
									ovFlg = false;
									if(step[k] === undefined) step[k] = [];
									for(var l = 0; l < step[k].length; l ++) {
										var checkS = step[k][l][0];
										var checkE = step[k][l][1];
										if(readS - this.margin <= checkE && checkS <= readE + this.margin) {
											ovFlg = true;
											break;
										}
									}
									if(!ovFlg) {
										step[k].push([readS, readE]);
										stepNo = k;
										read[j].step = stepNo;
										break;
									}
								}
							}
						}
					} else {
						for(var k = 0; k < step.length; k ++) {
							if(step[k] === undefined) step[k] = [];
							var ovFlg = false;
							for(var l = 0; l < step[k].length; l ++) {
								var checkS = step[k][l][0];
								var checkE = step[k][l][1];
								if(readS - this.margin <= checkE && checkS <= readE + this.margin) {
									ovFlg = true;
									break;
								}
							}
							if(!ovFlg) {
								step[k].push([readS, readE]);
								stepNo = k;
								break;
							}
						}
						read[j].step = stepNo;
					}
					
					if(stepNo === undefined || stepNo == -1) {
						if(
							(step.length + 1) * this.eachStep <= this.maxHeight - this.limitShowHeight &&
							stepNo != -1
						) {
							stepNo = step.length;
							step[stepNo] = [[readS, readE]];
							stepMax = stepNo;
						} else {
						//stepが上限を超えている場合
							if(overRegList.length == 0) {
								overRegList[0] = readS;
								overReg[readS] = readE;
							} else {
								var findFlg = false;
								var targetK = overRegList.length - 1;
								for(var k = 0; k < overRegList.length; k ++) {
									var checkS = overRegList[k];
									var checkE = overReg[checkS];
									
									if(checkS <= readS && readE <= checkE) {
									//領域に含まれている
										findFlg = true;
										break;
									} else if(readS <= checkE && checkS <= readE) {
									//重なっている(かつリードが領域を飛び出ている)
										if(checkS <= readE && readE <= checkE) {
										//リードが左に出てて、右には出てない
											overRegList.splice(k, 1, readS);
											delete overReg[checkS];
											overReg[readS] = readE;
										} else {
										//リードが右に出てる(左は出てる出てない両方)
											var newS = (checkS < readS)? checkS: readS;
											var newE = readE;
											var delCnt = 1;
											delete overReg[checkS];
											for(var l = k + 1; l < overRegList.length; l ++) {
												var checkS2 = overRegList[l];
												var checkE2 = overReg[checkS2];
												if(readE < checkS2) break;
												delCnt ++;
												delete overReg[checkS2];
												if(readE < checkE2) {
													newE = checkE2;
													break;
												}
											}
											overRegList.splice(k, delCnt, newS);
											overReg[newS] = newE;
										}
										findFlg = true;
										break;
									} else if(readE < checkS) {
									//リードは領域間に含まれる領域である
										overRegList.splice(k, 0, readS);
										overReg[readS] = readE;
										findFlg = true;
										break;
									}
								}
								
								if(!findFlg) {
									overRegList.push(readS);
									overReg[readS] = readE;
								}
							}
							
							read[j].step = -1; //stepNo = -1
							continue;
						}
						read[j].step = stepNo;
					}
				}
			}
		} else {
			//データのロードが終わってない
		}
	}
	
	var height = this.eachStep * (stepMax + 1);
	if(this.overRegList.length) height += this.limitShowHeight;
	this.height = (height > this.minHeight)? height: this.minHeight;
	
};
//ポップアップ用データを返す
WgCram.prototype.getPopupData = function() {
	return;
};
//自分のID
WgCram.prototype.getName = function() {
	return this.name;
};
WgCram.prototype.getItemDispName = function() {
	var dispName = (this.dispName)? this.dispName: this.name;
	return dispName;
};
WgCram.prototype.getPannelBgcolor = function() {
	return (this.option.pannelBgcolor)? this.option.pannelBgcolor: "#EEEEEE";
};
WgCram.prototype.accessObj = function(chr, binStart, binEnd, powP, accDefault) {
	var m = this;
	var bpPerPixel = Math.pow(10, powP);
	
	var bpStart = binStart * bpPerPixel * POW_REG + 1;
	var bpEnd = (parseInt(binStart) + 1) * bpPerPixel * POW_REG;
	
	if(powP > this.powMax) {
		return;
	}
	
	var magic = [];
	var data = {};
	if(this.option.uriDirFlg) {
		data[this.name] = {};
		data[this.name].reads = magic;
	} else {
		data[powP] = {};
		data[powP][chr + "|" + binStart] = {};
		data[powP][chr + "|" + binStart][this.name] = {};
		data[powP][chr + "|" + binStart][this.name].reads = magic;
	}
	
	if(binStart == binEnd && this.getSeqData) {
	//配列を取得する場合
		var jsonUrl = this.getSeqData.url + '/' + powP + '/' + chr + '/' + binStart + '/sequence/?' + this.getSeqData.param;
		//var jsonUrl = "https://dbtss.hgc.jp/cgi-bin/dbtss_db_json.cgi";
		//jsonUrl += '/9606/' + powP + '/' + chr + '/' + binStart + '/sequence/?SEE=1&UID=2';
		//var jsonUrl = "http://fullmal.hgc.jp/cgi-bin/fullmal_db_json.cgi";
		//jsonUrl += '?bin_region=' + chr + ':' + binStart + '-' + binStart + '&pow=' + powP + '&kinds=sequence&SPID=0301&UID=7&SEE=1';
		
		
		var ajaxParam = {
			success: function(seqData) {
				if(seqData.sequence !== undefined) data[m.name].sequence = seqData.sequence;
				//if(seqData[powP][chr + "|" + binStart] !== undefined)
				//	data[powP][chr + "|" + binStart][m.name].sequence = seqData[powP][chr + "|" + binStart].sequence;
				m.accessBamFile(chr, bpStart, bpEnd, accDefault, [data, magic]);
			},
			error: function(XMLHttpRequest, textStatus, errorThrown) {
				m.accessBamFile(chr, bpStart, bpEnd, accDefault, [data, magic]);
			},
			complete: function() {
			}
		};
		ajaxParam.url = jsonUrl;
		ajaxParam.xhrFields = {withCredentials: true};
		ajaxParam.dataType = 'json';
		$.ajax(ajaxParam);
	} else {
		this.accessBamFile(chr, bpStart, bpEnd, accDefault, [data, magic]);
	}
	
	
};
WgCram.prototype.accessBamFile = function(chr, bpStart, bpEnd, accDefault, dataMagic) {
	var m = this;
	this.bam.readWaitReader(chr, bpStart, bpEnd, function(result) {
		var counter = 0;
		var reads = dataMagic[1];
		for(var alnEach of result) {
			reads.push({
				qname: alnEach[3], flag: alnEach[2], pos: alnEach[0], mapq: alnEach[1],
				cigar: alnEach[4], seq: alnEach[6], posEnd: alnEach[5], id: alnEach[11],
				cigarLn: alnEach[12], cigarOp: alnEach[13]
			});
			counter ++;
			if(counter % 100000 == 0) {
				alert("Too many read (>= 100000) and some reads truncated. " + chr + ":" + bpStart + "-" + bpEnd + " will be incomplete display.");
				break;
			}
		}
		accDefault.success(dataMagic[0]);
		accDefault.complete();
	}, function(err) {
		accDefault.error(err, err, err);
		accDefault.complete();
	});
};

var WgFasta = function(fa, fai, option) {
	this.imgObj, this.ojson;
	
	this.option = (option !== undefined)? option: {};
	this.option.frameFlg = 
		(this.option.frameFlg === undefined)? true: this.option.frameFlg;
	this.option.inColorFlg = 
		(this.option.inColorFlg === undefined)? true: this.option.inColorFlg;
	this.option.buttonOnOffFlg = 
		(this.option.buttonOnOffFlg === undefined)? false: option.buttonOnOffFlg;
	
	this.charPx = (this.option.charPx === undefined)? 10: this.option.charPx;
	this.height = 12;
	this.y;
	
	this.fasta = new Fasta(fa, fai);
};
WgFasta.prototype = new WgRoot();
WgFasta.prototype.paint = function(y, width, chr, start, end, strand) {
	//描写状況:存在しないデータのbin情報を入れる
	var status = [];
	
	var statusFor = {};
	//1px未満になったら表示しない
	if(end - start + 1 < width) {
		var reg = 1 * POW_REG;
		var y1 = y;
		var y2 = y1 + 10;
		this.imgObj.font = this.charPx + "px 'Helvetica'";
		for(var i = parseInt(start); i <= parseInt(end + 1); i ++) {
			var bin = Math.floor((i - 1) / reg);
			//以下はもっと効率よくできるかも
			if(this.ojson["0"] && this.ojson["0"][chr + "|" + bin] && this.ojson["0"][chr + "|" + bin]["fasta"]) {
				var char = this.ojson["0"][chr + "|" + bin]["fasta"].charAt(i - 1 - bin * reg);
				if(strand == "-") {
					var tmp = char;
					if(char == "A") tmp = "T"; if(char == "a") tmp = "t";
					if(char == "C") tmp = "G"; if(char == "c") tmp = "g";
					if(char == "G") tmp = "C"; if(char == "g") tmp = "c";
					if(char == "T") tmp = "A"; if(char == "t") tmp = "a";
					char = tmp;
				}
				var x1 = (i - start) * (width - 1) / (end - start + 1);
				var x2 = (i - start + 1) * (width - 1) / (end - start + 1);
				if(strand == "-") {var tmp = width - 1 - x1; x1 = width - 1 - x2; x2 = tmp;}
				if(this.option.inColorFlg) {
					this.imgObj.fillStyle = 
						(char == "A" || char == "a")? "#88FF88":
						(char == "C" || char == "c")? "#8888FF":
						(char == "G" || char == "g")? "#FF8800":
						(char == "T" || char == "t")? "#FF4488": "#AAAAAA";
					this.imgObj.fillRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
				}
				if(x2 - x1 > this.charPx) {
					this.imgObj.fillStyle = "#000000";
					this.imgObj.fillText(char, (x1 + x2) / 2 - 2, y2 - 1);
					if(this.option.frameFlg) 
						this.imgObj.strokeRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
				}
			} else {
				statusFor[bin] = 1;
			}
		}
	}
	for(var bin in statusFor) if(bin >= 0) status.push(bin);
	
	return status;
};
WgFasta.prototype.getMenuPopup = function() {
	var htmlStr = "";
	return htmlStr;
};
WgFasta.prototype.getName = function() {
	return "fasta";
};
WgFasta.prototype.getItemDispName = function() {
	return "Fasta";
};
WgFasta.prototype.getButtonInfo = function() {
	var buttonColor = (this.option.buttonColor === undefined)? 
		{"color": ["eeeeee", "ffcccc"], "onOff": this.option.buttonOnOffFlg}: 
		{"color": this.option.buttonColor, "onOff": this.option.buttonOnOffFlg};
	
	return buttonColor;
};
WgFasta.prototype.accessObj = function(chr, binStart, binEnd, powP, accDefault) {
	var m = this;
	var bpPerPixel = Math.pow(10, powP);
	
	var bpStart = binStart * bpPerPixel * POW_REG + 1;
	var bpEnd = (parseInt(binStart) + 1) * bpPerPixel * POW_REG;
	
	this.fasta.readWaitReader(chr, bpStart, bpEnd, function(seq) {
		var data = {};
		data[m.getName()] = seq;
		accDefault.success(data);
		accDefault.complete();
	}, function(err) {
		throw new Error("data not found.");
	});
};

var WgFastaAmino = function(fa, fai, option) {
	this.imgObj, this.ojson;
	
	this.option = (option !== undefined)? option: {};
	this.option.frameFlg = 
		(this.option.frameFlg === undefined)? true: this.option.frameFlg;
	this.option.inColorFlg = 
		(this.option.inColorFlg === undefined)? true: this.option.inColorFlg;
	this.option.buttonOnOffFlg = 
		(this.option.buttonOnOffFlg === undefined)? false: option.buttonOnOffFlg;
	
	this.charPx = (this.option.charPx === undefined)? 10: this.option.charPx;
	this.height = 82;
	this.y;
	this.showType = (this.option.showType)? this.option.showType: "forward";
	this.translationTable = (this.option.translationTable) ? this.option.translationTable: "1";

	this.fasta = new Fasta(fa, fai);
};
WgFastaAmino.prototype = new WgRoot();
WgFastaAmino.prototype.paint = function(y, width, chr, start, end, strand) {
	//描写状況:存在しないデータのbin情報を入れる
	var status = [];
	
	var statusFor = {};
	//1px未満になったら表示しない
	if(end - start + 1 < width) {
		var reg = 1 * POW_REG;
		var y1 = y;
		var y2 = y1 + 10;
		this.imgObj.font = this.charPx + "px 'Helvetica'";
		for(var i = parseInt(start); i <= parseInt(end + 1); i ++) {
			var bin = Math.floor((i - 1) / reg);
			var rev_comp = function(char){
				var changeBase = {
					A: "T", T: "A", G: "C", C: "G",
					R: "Y", Y: "R", M: "K", K: "M",
					B: "V", V: "B", D: "H", H: "D",
					a: "t", t: "a", g: "c", c: "g",
					r: "y", y: "r", m: "k", k: "m",
					b: "v", v: "b", d: "h", h: "d",
				};
				return changeBase[char];
			}
			var genetic_code = {
				"1": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Ter", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"2": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Met", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Ter", AGG: "Ter",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"3": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Thr", CUC: "Thr", CUA: "Thr", CUG: "Thr",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Met", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"4": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"5": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Met", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Ser", AGG: "Ser",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"6": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Gln", UAG: "Gln",
					UGU: "Cys", UGC: "Cys", UGA: "Ter", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"9": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Asn", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Ser", AGG: "Ser",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"10": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Cys", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"11": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Ter", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"12": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Ter", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Ser",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"13": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Met", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Gly", AGG: "Gly",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"14": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Tyr", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Asn", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Ser", AGG: "Ser",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"16": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Leu",
					UGU: "Cys", UGC: "Cys", UGA: "Ter", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"21": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Trp", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Met", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Asn", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Ser", AGG: "Ser",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"22": {
					UUU: "Phe", UUC: "Phe", UUA: "Leu", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ter", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Leu",
					UGU: "Cys", UGC: "Cys", UGA: "Ter", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
				"23": {
					UUU: "Phe", UUC: "Phe", UUA: "Ter", UUG: "Leu",
					UCU: "Ser", UCC: "Ser", UCA: "Ser", UCG: "Ser",
					UAU: "Tyr", UAC: "Tyr", UAA: "Ter", UAG: "Ter",
					UGU: "Cys", UGC: "Cys", UGA: "Ter", UGG: "Trp",
					CUU: "Leu", CUC: "Leu", CUA: "Leu", CUG: "Leu",
					CCU: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
					CAU: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
					CGU: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
					AUU: "Ile", AUC: "Ile", AUA: "Ile", AUG: "Met",
					ACU: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
					AAU: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
					AGU: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
					GUU: "Val", GUC: "Val", GUA: "Val", GUG: "Val",
					GCU: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
					GAU: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
					GGU: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
				},
			};
			var amino_code = {
				"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
				"Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
				"Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
				"Tyr": "Y", "Val": "V", "Ter": "*",
			};

			if(this.ojson["0"]
					&& this.ojson["0"][chr + "|" + bin]
					&& this.ojson["0"][chr + "|" + bin]["fastaamino"]) {
				var char0 = this.ojson["0"][chr + "|" + bin]["fastaamino"].charAt(i - 2 - bin * reg);
				var char1 = this.ojson["0"][chr + "|" + bin]["fastaamino"].charAt(i - 1 - bin * reg);
				var char2 = this.ojson["0"][chr + "|" + bin]["fastaamino"].charAt(i - 0 - bin * reg);
				if (i - 2 - bin * reg < 0) {
					char0 = this.ojson["0"][chr + "|" + (bin - 1)]
						&& this.ojson["0"][chr + "|" + (bin - 1)]["fastaamino"]
						? this.ojson["0"][chr + "|" + (bin - 1)]["fastaamino"].charAt(
							this.ojson["0"][chr + "|" + (bin - 1)]["fastaamino"].length - 1)
						: undefined;
				} else if (i - 0 - bin * reg >= this.ojson["0"][chr + "|" + bin]["fastaamino"].length) {
					char2 = this.ojson["0"][chr + "|" + (bin + 1)]
						&& this.ojson["0"][chr + "|" + (bin + 1)]["fastaamino"]
						? this.ojson["0"][chr + "|" + (bin + 1)]["fastaamino"].charAt(0)
						: undefined;
				}
				if(strand == "-") {
					var tmp = char0;
					char0 = rev_comp(char2);
					char1 = rev_comp(char1);
					char2 = rev_comp(tmp);
				}

				var x1 = (i - start) * (width - 1) / (end - start + 1);
				var x2 = (i - start + 1) * (width - 1) / (end - start + 1);
				if(strand == "-") {
					var tmp = width - 1 - x1;
					x1 = width - 1 - x2;
					x2 = tmp;
				}
				var x3 = x1 - (width - 1) / (end - start + 1);
				var x4 = x2 + (width - 1) / (end - start + 1);
				var y3 = y1 + 12 + 10 * (i % 3);
				var y4 = y3 + 8;
				var y5 = y3 + 30
				var y6 = y4 + 30;

				// paint sequence
				if(this.option.inColorFlg) {
					this.imgObj.fillStyle =
						(char1 == "A" || char1 == "a")? "#88FF88":
						(char1 == "C" || char1 == "c")? "#8888FF":
						(char1 == "G" || char1 == "g")? "#FF8800":
						(char1 == "T" || char1 == "t")? "#FF4488": "#AAAAAA";
					this.imgObj.fillRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
				}
				if(x2 - x1 > this.charPx) {
					this.imgObj.fillStyle = "#000000";
					this.imgObj.fillText(char1, (x1 + x2) / 2 - 2, y2 - 1);
				}
				if(x2 - x1 > this.charPx && this.option.frameFlg) {
					this.imgObj.strokeStyle = "#000000";
					this.imgObj.strokeRect(x1, y1, x2 - x1 + 1, y2 - y1 + 1);
				}
				if (this.showType == "onlyseq" || !(char0 && char1 && char2)) {
					continue;
				}
				// paint amino acid (+strand)
				var codon = "";
				codon += char0.toUpperCase() == "T" ? "U" : char0.toUpperCase();
				codon += char1.toUpperCase() == "T" ? "U" : char1.toUpperCase();
				codon += char2.toUpperCase() == "T" ? "U" : char2.toUpperCase();
				var codon_rev = "";
				codon_rev += rev_comp(char2).toUpperCase() == "T" ? "U" : rev_comp(char2).toUpperCase();
				codon_rev += rev_comp(char1).toUpperCase() == "T" ? "U" : rev_comp(char1).toUpperCase();
				codon_rev += rev_comp(char0).toUpperCase() == "T" ? "U" : rev_comp(char0).toUpperCase();
				var amino = genetic_code[this.translationTable][codon];
				var amino_rev = genetic_code[this.translationTable][codon_rev];

				if(this.option.inColorFlg) {
					this.imgObj.fillStyle =
						amino == "Met" ? "#88FF88" :
						amino == "Ter" ? "#FF4488" :
						i % 6 < 3 ? "#AAAAAA" : "#FFFFFF";
					this.imgObj.fillRect(x3, y3, x4 - x3 + 1, y4 - y3 + 1);
				}
				if(x4 - x3 > this.charPx) {
					this.imgObj.fillStyle = "#000000";
					this.imgObj.fillText(
							x4 - x3 > this.charPx * 3 ? amino : amino_code[amino],
							(x3 + x4) / 2 - this.charPx * 0.3 - this.charPx * 0.6 * (x4 - x3 > this.charPx * 3),
							y4);
				}
				if(x4 - x3 > this.charPx && this.option.frameFlg) {
					this.imgObj.strokeStyle = "#000000";
					this.imgObj.strokeRect(x3, y3, x4 - x3 + 1, y4 - y3 + 1);
				}
				if (this.showType == "forward") {
					continue;
				}
				// paint reverse strand
				if(this.option.inColorFlg) {
					this.imgObj.fillStyle = 
						(rev_comp(char1) == "A" || rev_comp(char1) == "a")? "#88FF88":
						(rev_comp(char1) == "C" || rev_comp(char1) == "c")? "#8888FF":
						(rev_comp(char1) == "G" || rev_comp(char1) == "g")? "#FF8800":
						(rev_comp(char1) == "T" || rev_comp(char1) == "t")? "#FF4488": "#AAAAAA";
					this.imgObj.fillRect(x1, y1 + 72, x2 - x1 + 1, 10);
				}
				if(x2 - x1 > this.charPx) {
					this.imgObj.fillStyle = "#000000";
					this.imgObj.fillText(rev_comp(char1), (x1 + x2) / 2 - 2, y2 + 71);
				}
				if(x2 - x1 > this.charPx && this.option.frameFlg) {
					this.imgObj.strokeStyle = "#000000";
					this.imgObj.strokeRect(x1, y1 + 72, x2 - x1 + 1, y2 - y1 + 1);
				}
				if(this.option.inColorFlg) {
					this.imgObj.fillStyle =
						amino_rev == "Met" ? "#88FF88" :
						amino_rev == "Ter" ? "#FF4488" :
						i % 6 < 3 ? "#FFFFFF" : "#AAAAAA";
					this.imgObj.fillRect(x3, y5, x4 - x3 + 1, y6 - y5 + 1);
				}
				if(x4 - x3 > this.charPx) {
					this.imgObj.fillStyle = "#000000";
					this.imgObj.fillText(
							x4 - x3 > this.charPx * 3 ? amino_rev : amino_code[amino_rev],
							(x3 + x4) / 2 - this.charPx * 0.3 - this.charPx * 0.6 * (x4 - x3 > this.charPx * 3),
							y6);
				}
				if(x4 - x3 > this.charPx && this.option.frameFlg) {
					this.imgObj.strokeStyle = "#000000";
					this.imgObj.strokeRect(x3, y5, x4 - x3 + 1, y6 - y5 + 1);
				}
			} else {
				statusFor[bin] = 1;
			}
		}
	}
	for(var bin in statusFor) if(bin >= 0) status.push(bin);
	
	return status;
};
WgFastaAmino.prototype.setHeight = function(dt, width, chr, start, end) {
	if (this.showType == "onlyseq") {
		this.height = 12;
	} else if (this.showType == "forward") {
		this.height = 42;
	} else {
		this.height = 82;
	}
};
WgFastaAmino.prototype.getMenuPopup = function() {
	var checked = new Array("", "", "");
	if(this.showType == "onlyseq") checked[0] = "checked=\"checked\"";
	if(this.showType == "withrev") checked[1] = "checked=\"checked\"";
	if(this.showType == "forward") checked[2] = "checked=\"checked\"";
	var htmlStr = "";
	htmlStr += "<div style=\"border:1px solid\"><table border=\"0\" width=\"100%\"><tr><th align=\"left\" bgcolor=\"#aaaaaa\" colspan=\"2\">Setting:</th></tr><tr>";
	htmlStr += "<td bgcolor=\"#aaaaaa\">&nbsp;</td><td><form>";
	htmlStr += "<div><input type=\"radio\" name=\"show_type\" value=\"onlyseq\" id=\"onlyseq\" ";
	htmlStr += checked[0] + " /><label for=\"onlyseq\">Only Sequence</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"show_type\" value=\"withrev\" id=\"withrev\" ";
	htmlStr += checked[1] + " /><label for=\"withrev\">With reverse strand</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"show_type\" value=\"forward\" id=\"forward\" ";
	htmlStr += checked[2] + " /><label for=\"forward\">Default</label></div>";
	htmlStr += "</form>";
	htmlStr += "<hr /><a href=\"#\" class=\"det\">Change translation table</a>";
	htmlStr += "</td></tr></table></div>";
	// htmlStr += "<hr /><a href=\"#\" class=\"help\">Help</a>";
	
	return htmlStr;
};
WgFastaAmino.prototype.getMenuDetail = function() {
	var checked = new Array("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "");
	if(this.translationTable == "1") checked[0] = "checked=\"checked\"";
	if(this.translationTable == "2") checked[1] = "checked=\"checked\"";
	if(this.translationTable == "3") checked[2] = "checked=\"checked\"";
	if(this.translationTable == "4") checked[3] = "checked=\"checked\"";
	if(this.translationTable == "5") checked[4] = "checked=\"checked\"";
	if(this.translationTable == "6") checked[5] = "checked=\"checked\"";
	if(this.translationTable == "9") checked[6] = "checked=\"checked\"";
	if(this.translationTable == "10") checked[7] = "checked=\"checked\"";
	if(this.translationTable == "11") checked[8] = "checked=\"checked\"";
	if(this.translationTable == "12") checked[9] = "checked=\"checked\"";
	if(this.translationTable == "13") checked[10] = "checked=\"checked\"";
	if(this.translationTable == "14") checked[11] = "checked=\"checked\"";
	if(this.translationTable == "16") checked[12] = "checked=\"checked\"";
	if(this.translationTable == "21") checked[13] = "checked=\"checked\"";
	if(this.translationTable == "23") checked[14] = "checked=\"checked\"";
	var htmlStr = "";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"1\" id=\"1\" ";
	htmlStr += checked[0] + " /><label for=\"1\">Standard</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"2\" id=\"2\" ";
	htmlStr += checked[1] + " /><label for=\"2\">Vertebrate Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"3\" id=\"3\" ";
	htmlStr += checked[2] + " /><label for=\"3\">Yeast Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"4\" id=\"4\" ";
	htmlStr += checked[3] + " /><label for=\"4\">Mold, Protozoan, and Coelenterate Mitochondrial and Mycoplasma/Spiroplasma</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"5\" id=\"5\" ";
	htmlStr += checked[4] + " /><label for=\"5\">Invertebrate Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"6\" id=\"6\" ";
	htmlStr += checked[5] + " /><label for=\"6\">Ciliate, Dasycladacean and Hexamita Nuclear</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"9\" id=\"9\" ";
	htmlStr += checked[6] + " /><label for=\"9\">Echinoderm and Flatworm Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"10\" id=\"10\" ";
	htmlStr += checked[7] + " /><label for=\"10\">Euplotid Nuclear</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"11\" id=\"11\" ";
	htmlStr += checked[8] + " /><label for=\"11\">Bacterial, Archaeal and Plant Plastid</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"12\" id=\"12\" ";
	htmlStr += checked[9] + " /><label for=\"12\">Alternative Yeast Nuclear</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"13\" id=\"13\" ";
	htmlStr += checked[10] + " /><label for=\"13\">Ascidian Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"14\" id=\"14\" ";
	htmlStr += checked[11] + " /><label for=\"14\">Alternative Flatworm Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"16\" id=\"16\" ";
	htmlStr += checked[12] + " /><label for=\"16\">Chlorophycean Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"21\" id=\"21\" ";
	htmlStr += checked[13] + " /><label for=\"21\">Trematode Mitochondrial</label></div>";
	htmlStr += "<div><input type=\"radio\" name=\"translationTable\" value=\"23\" id=\"23\" ";
	htmlStr += checked[14] + " /><label for=\"23\">Thraustochytrium Mitochondrial</label></div>";
	htmlStr += "<div class=\"modal_btn\"><input type=\"button\" id=\"apply_button\" value=\"Apply\" /><input type=\"button\" id=\"cancel_button\" value=\"Cancel\" /></div>";
	return htmlStr;
};
WgFastaAmino.prototype.setViewDetail = function(divId) {
	this.translationTable =
		$("div" + divId + " div#modal div.container #1")[0].checked ? "1" :
		$("div" + divId + " div#modal div.container #2")[0].checked ? "2" :
		$("div" + divId + " div#modal div.container #3")[0].checked ? "3" :
		$("div" + divId + " div#modal div.container #4")[0].checked ? "4" :
		$("div" + divId + " div#modal div.container #5")[0].checked ? "5" :
		$("div" + divId + " div#modal div.container #6")[0].checked ? "6" :
		$("div" + divId + " div#modal div.container #9")[0].checked ? "9" :
		$("div" + divId + " div#modal div.container #10")[0].checked ? "10" :
		$("div" + divId + " div#modal div.container #11")[0].checked ? "11" :
		$("div" + divId + " div#modal div.container #12")[0].checked ? "12" :
		$("div" + divId + " div#modal div.container #13")[0].checked ? "13" :
		$("div" + divId + " div#modal div.container #14")[0].checked ? "14" :
		$("div" + divId + " div#modal div.container #16")[0].checked ? "16" :
		$("div" + divId + " div#modal div.container #21")[0].checked ? "21" :
		$("div" + divId + " div#modal div.container #23")[0].checked ? "23" :
		"1";
};
WgFastaAmino.prototype.getName = function() {
	return "fastaamino";
};
WgFastaAmino.prototype.getItemDispName = function() {
	return "FastaAmino";
};
WgFastaAmino.prototype.getButtonInfo = function() {
	var buttonColor = (this.option.buttonColor === undefined)? 
		{"color": ["eeeeee", "ffcccc"], "onOff": this.option.buttonOnOffFlg}: 
		{"color": this.option.buttonColor, "onOff": this.option.buttonOnOffFlg};
	
	return buttonColor;
};
WgFastaAmino.prototype.accessObj = function(chr, binStart, binEnd, powP, accDefault) {
	var m = this;
	var bpPerPixel = Math.pow(10, powP);
	
	var bpStart = binStart * bpPerPixel * POW_REG + 1;
	var bpEnd = (parseInt(binStart) + 1) * bpPerPixel * POW_REG;
	
	this.fasta.readWaitReader(chr, bpStart, bpEnd, function(seq) {
		var data = {};
		data[m.getName()] = seq;
		accDefault.success(data);
		accDefault.complete();
	}, function(err) {
		throw new Error("data not found.");
	});
};

class Fasta {
	constructor(fa, fai, cache_size = 0) {
		// fa and fai are FileHandler
		if (!fa || !fai) {
			throw "Files are Falsy";
		}
		this.fa_ = new FileHandler(fa, false);
		this.fai_ = new FileHandler(fai, false);
		this.faindex_ = undefined;
		this.changeBase_ = {
			A: "T",
			T: "A",
			G: "C",
			C: "G",
			R: "Y",
			Y: "R",
			M: "K",
			K: "M",
			B: "V",
			V: "B",
			D: "H",
			H: "D",
			a: "t",
			t: "a",
			g: "c",
			c: "g",
			r: "y",
			y: "r",
			m: "k",
			k: "m",
			b: "v",
			v: "b",
			d: "h",
			h: "d",
		};
		this.cache_size_ = cache_size;
		if (this.cache_size_ > 0) {
			this.cache_ = [
				[
					new Date(),
					"",
					1,
					cache_size,
					new Promise((resolve) => {
						resolve("N".repeat(cache_size));
					}),
				],
			];
		}
	}

	readWaitReader(chr, start, end, callback, reject, option) {
		try{
			this.loadSequence_(chr, start, end).then((seq)=>{
				callback(seq);
			});
		} catch(e) {
			reject(e);
		}
	}

	loadSequence_(chr, start, end, strand = "+") {
		return this.loadSeq_(chr, start, end).then((seq) => {
			if (strand !== "-" && strand !== -1) {
				return seq;
			} else {
				return this.reverseComplement_(seq);
			}
		});
	}

	async loadSeq_(chr, start, end) {
		if (!this.cache_) {
			return await this.loadFromSource_(chr, start, end);
		}
		var fragments = [];
		this.cache_.sort((a, b) => {
			return a[2] - b[2];
		});
		for (var i = 0; i < this.cache_.length; i++) {
			var elem = this.cache_[i];
			if (elem[1] == chr && elem[2] <= start && elem[3] >= end) {
				// entire data is in this element of cache
				this.cache_[i][0] = new Date();
				var start_index = start - elem[2];
				var end_index = start_index + (end - start) + 1;
				var cache_seq = await elem[4];
				return cache_seq.slice(start_index, end_index);
			} else if (elem[1] == chr && elem[2] <= end && elem[3] >= start) {
				// part of data is in this element of cache
				fragments.push([elem, i]);
			}
		}

		if (fragments.length == 0) {
			// no data is in this element of cache
			var seq = this.loadFromSource_(chr, start, end);
			this.cache_.push([new Date(), chr, start, end, seq]);
			this.shrinkCache_(end - start + 1);
			return await seq;
		}

		fragments[0][0][4].then((seq) => {});
		var load_length = 0;
		if (fragments[0][0][2] > start) {
			var p_seq = this.loadFromSource_(chr, start, fragments[0][2] - 1);
			load_length += fragments[0][2] - start;
			this.cache_[fragments[0][1]][0] = new Date();
			this.cache_[fragments[0][1]][2] = start;
			this.cache_[fragments[0][1]][4] = Promise.all([
				p_seq,
				this.cache_[fragments[0][1]][4],
			]).then((values) => {
				return values[0].concat(values[1]);
			});
		}
		if (fragments.length >= 2) {
			for (var i = 1; i < fragments.length; i++) {
				var p_loaded_seq = this.loadFromSource_(
					chr,
					this.cache_[fragments[i - 1][1]][3] + 1,
					fragments[i][0][2] - 1
				);
				load_length +=
					fragments[i][0][2] - this.cache_[fragments[i - 1][1]][3] - 1;
				this.cache_[fragments[0][1]][0] = new Date();
				this.cache_[fragments[0][1]][3] = fragments[i][0][3];
				this.cache_[fragments[0][1]][4] = p_loaded_seq.then((loaded_seq) => {
					return this.cache_[fragments[0][1]][4].concat(loaded_seq);
				});
				this.cache_[fragments[0][1]][4] = Promise.all([
					this.cache_[fragments[0][1]][4],
					fragments[i][0][4],
				]).then((values) => {
					return values[0].concat(values[1]);
				});
			}
		}
		if (fragments[fragments.length - 1][0][3] < end) {
			var p_seq = this.loadFromSource_(
				chr,
				fragments[fragments.length - 1][0][3] + 1,
				end
			);
			load_length += end - fragments[fragments.length - 1][0][3];
			this.cache_[fragments[0][1]][0] = new Date();
			this.cache_[fragments[0][1]][3] = end;
			this.cache_[fragments[0][1]][4] = Promise.all([
				this.cache_[fragments[0][1]][4],
				p_seq,
			]).then((values) => {
				return values[0].concat(values[1]);
			});
		}

		if (fragments.length >= 2) {
			for (var i = 1; i < fragments.length; i++) {
				var index = this.cache_.indexOf(fragments[i][0]);
				this.cache_.splice(index, 1);
			}
		}
		this.shrinkCache_(load_length);
		return this.loadSeq_(chr, start, end);
	}

	getBytePos_(pos, offset, linebases, linewidth) {
		// pos is 1-start coordinate
		var m = (pos - 1) % linebases;
		var n = (pos - 1 - m) / linebases;
		return offset + linewidth * n + m;
	}

	loadFromSource_(chr, start, end) {
		if (typeof this.faindex_ === "undefined") {
			this.faindex_ = this.loadFai_();
		}
		return this.faindex_
			.then((faindex) => {
				const index = faindex.find((elem) => {
					return elem[0] == chr;
				});
				if (typeof index === "undefined") {
					throw "Chromosome (" + chr + ") is not found in faindex.";
				}
				if (start < 1 || end > index[1]) {
					throw "out of bounds";
				}
				var start_byte = this.getBytePos_(start, index[2], index[3], index[4]);
				var end_byte = this.getBytePos_(end, index[2], index[3], index[4]);
				return this.fa_.load(start_byte, end_byte - start_byte + 1);
			})
			.then((arraybuffer) => {
				var str = String.fromCharCode.apply("", new Uint8Array(arraybuffer));
				return str.replace(/\r?\n/g, "");
			});
	}

	reverseComplement_(seq) {
		var retval = new String("");
		for (var i = 0; i < seq.length; i++) {
			retval = this.changeBase_[seq.charAt(i)].concat(retval);
		}
		return retval;
	}

	shrinkCache_(length) {
		this.cache_.sort((a, b) => {
			return a[0] - b[0];
		});
		while (length >= this.cache_[0][3] - this.cache_[0][2] + 1) {
			length -= this.cache_[0][3] - this.cache_[0][2] + 1;
			this.cache_.splice(0, 1);
		}
		this.cache_[0][2] += length;
		this.cache_[0][4] = this.cache_[0][4].then((seq) => {
			return seq.slice(length);
		});
	}

	loadFai_() {
		return this.fai_.load().then((arraybuffer) => {
			var faindex = [];
			var fai = String.fromCharCode.apply("", new Uint8Array(arraybuffer));
			var lines = fai.split("\n");
			lines.forEach((line) => {
				var l = line.split("\t");
				if (l.length == 5) {
					faindex.push([
						l[0],
						parseInt(l[1], 10),
						parseInt(l[2], 10),
						parseInt(l[3], 10),
						parseInt(l[4], 10),
					]);
				}
			});
			return faindex;
		});
	}
}
class FileHandler {
  constructor(
    file /* instance of File object or String of URL */,
    local_flag = true
  ) {
    this.file_ = file;
    this.local_flag_ = local_flag;
  }

  load(pos = -1, length = -1) {
    return new Promise((resolve, reject) => {
      if (this.local_flag_) {
        if (pos >= 0 && length > 0) {
          var sliced = this.file_.slice(pos, pos + length);
          resolve(sliced.arrayBuffer());
        } else {
          resolve(this.file_.arrayBuffer());
        }
      } else {
        var oReq = new XMLHttpRequest();
        oReq.open("GET", this.file_);
        if (pos >= 0 && length > 0) {
          oReq.setRequestHeader(
            "Range",
            "bytes=" + pos + "-" + (pos + length - 1)
          );
        }
        oReq.responseType = "arraybuffer";
        oReq.onload = function (oEvent) {
          const ab = oReq.response;
          if (ab) {
            resolve(ab);
          } else {
            reject(oReq.statusText);
          }
        };
        oReq.onerror = function () {
          reject("An error occurred during HTTP access");
        };
        oReq.onabort = function () {
          reject("HTTP access is aborted.");
        };
        oReq.timeout = function () {
          reject("HTTP access timed out.");
        };
        oReq.send();
      }
    });
  }
}
