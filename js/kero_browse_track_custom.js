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
	
	this.wg.readWaitReader(chr, bpStart, bpEnd, function(fetcher) {
		var seq = "";
		for(var seqEach of fetcher()) {
			seq += seqEach;
		}
		var data = {};
		data[m.getName()] = seq;
		accDefault.success(data);
		accDefault.complete();
	}, function(err) {
		throw new Error("data not found.");
	});
};
