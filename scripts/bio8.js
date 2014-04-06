/**********************************
 *** Common routines from k8.js ***
 **********************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

function intv_ovlp(intv, bits)
{
	if (typeof bits == "undefined") bits = 13;
	intv.sort(function(a,b) {return a[0]-b[0];});
	// merge overlapping regions
	var j = 0;
	for (var i = 1; i < intv.length; ++i) {
		if (intv[j][1] >= intv[i][0])
			intv[j][1] = intv[j][1] > intv[i][1]? intv[j][1] : intv[i][1];
		else intv[++j] = [intv[i][0], intv[i][1]];
	}
	intv.length = j + 1;
	// create the index
	var idx = [], max = 0;
	for (var i = 0; i < intv.length; ++i) {
		var b = intv[i][0]>>bits;
		var e = (intv[i][1]-1)>>bits;
		if (b != e) {
			for (var j = b; j <= e; ++j)
				if (idx[j] == null) idx[j] = i;
		} else if (idx[b] == null) idx[b] = i;
		max = max > e? max : e;
	}
	return function(_b, _e) { // closure
		var x = _b >> bits;
		if (x > max) return false;
		var off = idx[x];
		if (off == null) {
			var i;
			for (i = ((_e - 1) >> bits) - 1; i >= 0; --i)
				if (idx[i] != null) break;
			off = i < 0? 0 : idx[i];
		}
		for (var i = off; i < intv.length && intv[i][0] < _e; ++i)
			if (intv[i][1] > _b) return true;
		return false;
	}
}

/*****************************
 *** Compare two BED files ***
 *****************************/

function read_bed(fn)
{
	var f = fn == '-'? new File() : new File(fn);
	var b = new Bytes();
	var reg = {}, idx = {};
	while (f.readline(b) >= 0) {
		var t = b.toString().split("\t");
		if (reg[t[0]] == null) reg[t[0]] = [];
		reg[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
	}
	var n_rec = 0;
	for (var chr in reg) {
		idx[chr] = intv_ovlp(reg[chr]);
		n_rec += reg[chr].length;
	}
	b.destroy();
	f.close();
	return [idx, n_rec];
}

function b8_bedovlp(args)
{
	if (args.length < 2) {
		print("Usage: k8 bio8.js bedovlp <read.bed> <inRAM.bed>");
		exit(1);
	}
	var n_rec;
	var f1 = new File(args[0]);
	var tmp = read_bed(args[1]);
	var bed = tmp[0], n_rec = tmp[1];
	var b1 = new Bytes();
	var c = [0, 0, 0];
	while (f1.readline(b1) >= 0) {
		var t = b1.toString().split("\t", 3);
		t[1] = parseInt(t[1]); t[2] = parseInt(t[2]);
		if (bed[t[0]] == null) ++c[0];
		else if (bed[t[0]](t[1], t[2])) ++c[1];
		else ++c[2];
	}
	b1.destroy();
	f1.close();
	print(args[0], args[1], c[0], c[1], c[2], n_rec);
}

/*********************
 *********************/

function b8_bedcmpm(args)
{
	function cmpfunc(a, b) {
		if (a[0] < b[0]) return -1;
		if (a[0] > b[0]) return 1;
		return a[1] - b[1];
	}

	var c, win_size = 0, print_joint = false;
	while ((c = getopt(args, 'pw:')) != null)
		if (c == 'w') win_size = parseInt(getopt.arg);
		else if (c == 'p') print_joint = true;

	var bed = [];
	var buf = new Bytes();
	var n = args.length - getopt.ind;
	for (var i = getopt.ind; i < args.length; ++i) {
		var file = new File(args[i]);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			t[1] = parseInt(t[1]);
			if (t.length > 2) t[2] = parseInt(t[2]);
			else t[2] = t[1], --t[1];
			t[1] = t[1] < win_size? 0 : t[1] - win_size;
			t[2] += win_size;
			bed.push([t[0], t[1], t[2], 1 << (i-getopt.ind)]);
		}
		file.close();
	}
	bed.sort(cmpfunc);
	var labels = [];
	for (var i = 0; i < 1<<n; ++i) labels[i] = 0;
	var last = null, start = 0, end = 0, label = 0;
	for (var j = 0; j < bed.length; ++j) {
		if (bed[j][0] != last || bed[j][1] >= end) { // no overlap
			if (last != null) {
				if (print_joint) print(last, start, end, label);
				++labels[label];
			}
			last = bed[j][0], start = bed[j][1], end = bed[j][2], label = bed[j][3];
		} else end = end > bed[j][2]? end : bed[j][2], label |= bed[j][3];
	}
	if (print_joint) print(last, start, end, label);
	++labels[label];
	for (var i = 0; i < 1<<n; ++i)
		print(i, labels[i]);
}

/******************************************
 *** Circular Binary Segmentation (CBS) ***
 ******************************************/

function cnv_cbs1(a0, start, end, winsize) // find the best cut for a single segment
{
	var a = a0.slice(start, end);
	a.sort(function(x,y){return x[1]-y[1]}); // sort by value
	var same = 1, T = 0;
	for (var i = 1; i <= a.length; ++i) { // compute the rank
		if (i == a.length || a[i-1][1] != a[i][1]) {
			var rank = i - .5 * (same - 1);
			for (var j = i - same; j <= i - 1; ++j)
				a[j][2] = rank;
			if (same > 1) T += (same - 1.) * same * (same + 1.), same = 1; // T to correct for ties
		} else ++same;
	}
	var z = ((a.length + 1) - T / a.length / (a.length - 1)) / 12.;
	a.sort(function(x,y){return x[0]-y[0]});
	var s = [], sum = 0;
	s[0] = 0;
	for (var i = 0; i < a.length; ++i) s[i+1] = s[i] + a[i][2];
	var max = [0, -1, -1, 0];
	for (var j = 1; j <= a.length; ++j) {
		for (var i = j > winsize? j - winsize : 0; i < j; ++i) { // compute the standardized Mann-Whitney U
			var n = j - i, m = a.length - n;
			var x = (s[j] - s[i] - .5 * n * (a.length + 1)) / Math.sqrt(m * n * z);
			var y = x > 0? x : -x;
			if (y > max[0]) max = [y, i, j, x]; // this is the best cut
		}
	}
	return max;
}

function cnv_cbs(chr, a, thres, winsize) // recursively cut
{ // FIXME: it does not consider the edge effect; should not be a big issue
	function print_intv(b, e) { // print an interval
		if (b == e) return;
		var t = a.slice(b, e);
		t.sort(function(x,y){return x[1]-y[1]}); // sort to compute the median
		print(chr, a[b][0] - 1, a[e-1][0], t.length, t[Math.floor(t.length/2)][1]);
	}

	warn("Processing chromsome " + chr + "...");
	a.sort(function(x,y){return x[0]-y[0]});
	var queue = [[0, a.length, -1]];
	while (queue.length) {
		var e = queue.shift();
		var start = e[0], end = e[1];
		var cut = cnv_cbs1(a, start, end, winsize);
		if (cut[0] < -thres || cut[0] > thres) { // a significant cut; then cut further
			if (cut[1] >= 4) queue.push([start, start + cut[1]]);
			else print_intv(start, start + cut[1]);
			if (cut[2] - cut[1] >= 4) queue.push([start + cut[1], start + cut[2]]);
			else print_intv(start + cut[1], start + cut[2]);
			if (end - start - cut[2] >= 4) queue.push([start + cut[2], end]);
			else print_intv(start + cut[2], end);
		} else print_intv(start, end); // no significant cut can be found
	}
}

function b8_cbs(args)
{
	var c, winsize = 200, thres = 3.891; // 3.891 is equivalent to P-value 1e-4
	while ((c = getopt(args, 't:w:')) != null) {
		if (c == 't') thres = parseFloat(getopt.arg);
		else if (c == 'w') winsize = parseInt(getopt.arg);
	}

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var last_chr = null, a = [];

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0] != last_chr) {
			if (a.length) cnv_cbs(last_chr, a, thres, winsize);
			last_chr = t[0];
			a = [];
		}
		a.push([parseInt(t[1]), parseFloat(t[2])]);
	}
	if (a.length) cnv_cbs(last_chr, a, thres, winsize);

	buf.destroy();
	file.close();
}

/****************************
 *** Parse one-sample VCF ***
 ****************************/

function b8_parse_vcf1(t) // t = vcf_line.split("\t")
{
	var a = [];
	if (t.length < 6) return null;
	t[1] = parseInt(t[1]) - 1; t[3] = t[3].toUpperCase(); t[4] = t[4].toUpperCase(); t[5] = parseFloat(t[5]);
	var s = t[4].split(","); // list of ALT alleles
	// find the genotype
	var gt, match = /^(\d+)(\/|\|)(\d+)/.exec(t[9]);
	if (match == null) { // special-casing samtools, as for single-sample, it effectively gives the genotype at FQ
		var m2 = /FQ=([^;\t]+)/.exec(t[7]);
		if (m2 == null) gt = /^\.\/\./.test(t[7])? 0 : -1; // special casing CG VCF
		else gt = parseFloat(m2[1]) > 0? 1 : 0;
	} else gt = parseInt(match[1]) != parseInt(match[3])? 1 : 0;
	// get CIGAR for freebayes
	var m3 = /CIGAR=([^;\t]+)/.exec(t[7]);
	var cigar = m3 != null? m3[1].split(",") : [];
	if (cigar.length && cigar.length != s.length) throw Error("Inconsistent ALT and CIGAR");
	// loop through each ALT allele
	for (var i = 0; i < s.length; ++i) {
		var type; // 0=ts, 1=tv, 2=mnp, 3=ins, 4=del
		if (t[3].length == 1 && s[i].length == 1) { // SNP
			type = ((t[3] == 'A' && s[i] == 'G') || (t[3] == 'G' && s[i] == 'A') || (t[3] == 'C' && s[i] == 'T') || (t[3] == 'T' && s[i] == 'C'))? 0 : 1;
			a.push([t[5], type, gt, t[1], 0]);
		} else if (cigar.length) { // MNP or INDEL from freebayes
			var x = 0, y = 0;
			var m4, re = /(\d+)([MXID])/g;
			while ((m4 = re.exec(cigar[i])) != null) {
				var l = parseInt(m4[1]);
				if (m4[2] == 'X') {
					for (var j = 0; j < l; ++j) {
						u = t[3].substr(x+j, 1), v = s[i].substr(y+j, 1);
						type = ((u == 'A' && v == 'G') || (u == 'G' && v == 'A') || (u == 'C' && v == 'T') || (u == 'T' && v == 'C'))? 0 : 1;
						a.push([t[5], type, gt, t[1] + x, 0]);
					}
					x += l, y += l;
				} else if (m4[2] == 'I') {
					a.push([t[5], 3, gt, t[1] + x, l]);
					y += l;
				} else if (m4[2] == 'D') {
					a.push([t[5], 4, gt, t[1] + x, -l]);
					x += l;
				} else if (m4[2] == 'M') x += l, y += l;
			}
		} else { // MNP or INDEL from Platypus and others
			var l = t[3].length < s[i].length? t[3].length : s[i].length;
			for (var j = 0; j < l; ++j) { // decompose long variants
				var u = t[3].substr(j, 1), v = s[i].substr(j, 1);
				if (u != v) {
					type = ((u == 'A' && v == 'G') || (u == 'G' && v == 'A') || (u == 'C' && v == 'T') || (u == 'T' && v == 'C'))? 0 : 1;
					a.push([t[5], type, gt, t[1] + j, 0]);
				}
			}
			if (t[3].length != s[i].length) { // INDEL
				type = t[3].length < s[i].length? 3 : 4;
				a.push([t[5], type, gt, t[1], s[i].length - t[3].length]);
			}
		}
	}
	return a; // [qual, ts/tv/ins/del, gt, pos, indelLen]
}

/*****************
 *** VCF stats ***
 *****************/

function b8_qst1(args)
{
	var c, b = 0.02, show_hdr = false, min_q = 0, het_num = null, tstv = null, indel_bed = false, snp_pos = false;
	while ((c = getopt(args, "SGhb:q:H:T:")) != null)
		if (c == 'b') b = parseFloat(getopt.arg);
		else if (c == 'q') min_q = parseFloat(getopt.arg);
		else if (c == 'h') show_hdr = true;
		else if (c == 'H') het_num = parseInt(getopt.arg);
		else if (c == 'T') tstv = parseFloat(getopt.arg);
		else if (c == 'G') indel_bed = true;
		else if (c == 'S') snp_pos = true;

	if (het_num != null && tstv != null) throw Error("-H and -T cannot be specified at the same time");

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var a = [];
	while (file.readline(buf) >= 0) {
		if (buf[0] == 35) continue;
		var t = buf.toString().split("\t");
		var u = b8_parse_vcf1(t);
		if (u == null || u.length == 0 || u[0][0] < min_q) continue;
		if (snp_pos) {
			for (var i = 0; i < u.length; ++i)
				if (u[i][1] == 0 || u[i][1] == 1)
					print(t[0], u[i][3], u[i][2]);
		} else if (indel_bed) {
			var indel_len = [];
			t[1] = parseInt(t[1]);
			for (var i = 0; i < u.length; ++i)
				if (u[i][1] == 3 || u[i][1] == 4) indel_len.push(u[i][4]);
			if (indel_len.length)
				print(t[0], t[1] - 1, t[1] - 1 + t[3].length, indel_len.sort().join(","), u[0][2]);
			continue;
		} else {
			if (u[0][2] < 0) continue;
			for (var i = 0; i < u.length; ++i)
				a.push(u[i].slice(0, 3));
		}
	}
	buf.destroy(buf);
	file.close();

	a.sort(function(x,y) {return y[0]-x[0]});
	if (het_num != null) {
		var i, s = 0;
		for (i = 0; i < a.length; ++i) {
			if (a[i][1] == 2) continue; // ignore MNPs as they are decomposed to SNPs
			if (a[i][2]) ++s;
			if (het_num != null && s > het_num) break;
		}
		print(a[i-1][0]);
	} else if (tstv != null) {
		var b = [], max = 10000;
		for (var q = 0; q <= max; ++q) b[q] = [0, 0];
		for (var i = 0; i < a.length; ++i) {
			var q = Math.floor(a[i][0] + .499);
			q = q < max? q : max;
			if (a[i][1] == 0) ++b[q][0];
			else if (a[i][1] == 1) ++b[q][1];
		}
		var score = 0, m = 0, mq = 0;
		for (var q = 0; q < max; ++q) {
			score += tstv * b[q][1] - b[q][0];
			if (score > m) m = score, mq = q;
		}
		print(mq);
	} else if (!indel_bed && !snp_pos) {
		if (show_hdr) print("Q", "#", "#tsHom", "#tsHet", "#tvHom", "#tvHet", "#mnpHom", "#mnpHet", "#insHom", "#insHet", "#delHom", "#delHet");
		var size = Math.floor(a.length * b + 1.);
		var lastq = -1;
		var c = [], ac = [];
		for (var j = 0; j < 11; ++j) c[j] = ac[j] = 0;
		for (var i = 0; i <= a.length; ++i) {
			if (i == a.length || (a[i][0] != lastq && c[0] > a.length * b)) {
				print(lastq, ac.join("\t"), c.join("\t"));
				if (i == a.length) break;
				for (var j = 0; j < 11; ++j) c[j] = 0;
			}
			++c[0]; ++ac[0];
			++c[a[i][1] * 2 + a[i][2] + 1];
			++ac[a[i][1] * 2 + a[i][2] + 1];
			lastq = a[i][0];
		}
	}
}

/************************************************
 *** recompute GQ and GT in single-sample VCF ***
 ************************************************/

function b8_upd1gt(args)
{
	var c, prior = 30;
	while ((c = getopt(args, "p:")) != null)
		if (c == 'p') prior = parseFloat(getopt.arg);

	// each value in GL corresponds to which genotype? Assuming diploid
	var max_alleles = 16, het_arr = [], gt_arr = [];
	for (var i = 0; i < max_alleles; ++i) {
		for (var j = 0; j <= i; ++j) {
			het_arr.push(i != j);
			gt_arr.push(j + '/' + i);
		}
	}

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var c1 = '#'.charCodeAt(0);
	var re = /^(\d+)(\/|\|)(\d+)/;
	var lineno = 0;

	while (file.readline(buf) >= 0) {
		++lineno;
		if (buf[0] == c1) {
			print(buf);
			continue;
		}
		var t = buf.toString().split("\t");
		var s1 = t[8].split(":");
		var s2 = t[9].split(":");
		var i;
		for (i = 0; i < s1.length; ++i)
			if (s1[i] == 'GL' || s1[i] == 'PL') break;
		if (i == s1.length) {
			warn("WARNING: no GL/PL at line " + lineno + ":\n" + buf.toString());
			print(buf);
			continue;
		}
		var u = s2[i].split(",");
		var n_alleles = t[4].split(",").length + 1;
		if (u.length != n_alleles * (n_alleles + 1) / 2) {
			warn("WARNING: inconsistent GL/PL at line " + lineno + ":\n" + buf.toString());
			print(buf);
			continue;
		}
		var is_GL = (s1[i] == 'GL');
		for (var j = 0; j < u.length; ++j) {
			u[j] = parseFloat(u[j]);
			if (is_GL) u[j] *= -10.;
			u[j] += (het_arr[j]? prior : 0);
		}
		var min = 1e37, min_j = -1, min2 = 1e37;
		for (var j = 0; j < u.length; ++j)
			if (min > u[j]) min2 = min, min = u[j], min_j = j;
			else if (min2 > u[j]) min2 = u[j]; 
		var GT = gt_arr[min_j], GQ = Math.floor(min2 - min + .499);
		if (s1[0] != 'GT') {
			s1.unshift('GT'); s2.unshift(GT);
		} else s2[0] = GT;
		var k;
		for (k = 0; k < s1.length; ++k)
			if (s1[k] == 'GQ') break;
		if (k == s1.length) {
			s1.push('GQ'); s2.push(GQ);
		} else s2[k] = GQ;
		t[8] = s1.join(":");
		t[9] = s2.join(":");
		print(t.join("\t"));
	}
	buf.destroy(); file.close();
}

/******************
 *** De-overlap ***
 ******************/

function b8_deovlp(args)
{
	var c;
	while ((c = getopt(args, "")) != null);

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var c1 = '#'.charCodeAt(0);
	var a = [];
	var n = 0;
	while (file.readline(buf) >= 0) {
		if (buf[0] == c1) {
			print(buf);
			continue;
		}
		var s = buf.toString();
		var t = s.split("\t", 6);
		t[1] = parseInt(t[1]) - 1; t[5] = parseFloat(t[5]);
		if (a.length) {
			var i;
			for (i = 0; i < a.length && (a[i][0] != t[0] || a[i][1] <= t[1]); ++i)
				if (a[i][3]) print(a[i][4]);
				else ++n;
			while (i--) a.shift();
		}
		var to_print = true;
		if (a.length) {
			for (var i = 0; i < a.length; ++i) {
				if (a[i][1] <= t[1] || !a[i][3]) continue;
				if (a[i][2] < t[5]) a[i][3] = false;
				else to_print = false;
			}
		}
		a.push([t[0], t[1] + t[3].length, t[5], to_print, s]);
	}
	for (var i = 0; i < a.length; ++i)
		if (a[i][3]) print(a[i][4]);
		else ++n;
	buf.destroy();
	file.close();
	warn(n + " variants have been dropped");
}

/************************************
 *** Convert pileup output to var ***
 ************************************/

function b8_plp2var(args)
{
	var file = args.length? new File(args[0]) : new File();
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (!/[ACGTacgt]/.test(t[4])) continue;
		t[2] = t[2].toUpperCase();
		var i = 0, j = 0;
		var alleles = {}, cnt = [], sum = 0;
		while (i < t[4].length && j < t[5].length) {
			var b = t[4].charAt(i);
			if (b == '$') { ++i; continue; }
			if (b == '^') { i += 2; continue; }
			if (b == '*') { ++i, ++j; continue; }

			// determine the allele sequence
			var a, q, forward;
			var match = /[.,A-Za-z]([+-](\d+)[A-Za-z])?/.exec(t[4].substr(i));
			if (b == '.') b = t[2].toUpperCase();
			else if (b == ',') b = t[2].toLowerCase();
			forward = (b.charCodeAt(0) < 97);
			q = t[5].charCodeAt(j) - 33;
			var l_int = 0, l = 0;
			if (match[1] != null) {
				l_int = match[2].length + 1; // including +/-
				l = parseInt(match[2]);
				a = (b + t[4].substr(i + 1, l_int +l)).toUpperCase();
			} else a = b.toUpperCase();
			i += 1 + l_int + l;
			++j;

			// count
			var ci;
			if (alleles[a] == null) alleles[a] = cnt.length, cnt.push([a, 0, 0, 0, 0]);
			ci = alleles[a];
			++cnt[ci][forward? 1 : 2];
			cnt[ci][forward? 3 : 4] += q;
			sum += q;
		}

		var out = [t[0], t[1], t[2], sum, cnt.length];
		for (var i = 0; i < cnt.length; ++i)
			for (var j = 0; j < 5; ++j)
				out.push(cnt[i][j]);
		print(out.join("\t"));
	}

	buf.destroy();
	file.close();
}

/*************************************
 *** Convert plp2var output to VCF ***
 *************************************/

function b8_var2vcf(args)
{
	var c, qdp = false;
	while ((c = getopt(args, 'q')) != null)
		if (c == 'q') qdp = true;
	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var max = 0, match, max_del = '';
		var t = buf.toString().split("\t");
		t[3] = parseInt(t[3]); t[4] = parseInt(t[4]);
		for (var i = 0; i < t[4]; ++i) {
			var match = /^[A-Z]-(\d+)([A-Z]+)/.exec(t[5*(i+1)]);
			if (match != null && max < parseInt(match[1]))
				max = parseInt(match[1]), max_del = match[2];
		}
		var alt = [], dp4 = [0, 0, 0, 0], q = [], qs4 = [0, 0, 0, 0];
		for (var i = 0; i < t[4]; ++i) {
			var a = t[5*(i+1)], match;
			if (a == t[2]) {
				dp4[0] += parseInt(t[5*(i+1) + 1]);
				dp4[1] += parseInt(t[5*(i+1) + 2]);
				qs4[0] += parseInt(t[5*(i+1) + 3]);
				qs4[1] += parseInt(t[5*(i+1) + 4]);
				continue; // identical to the reference
			} else {
				dp4[2] += parseInt(t[5*(i+1) + 1]);
				dp4[3] += parseInt(t[5*(i+1) + 2]);
				qs4[2] += parseInt(t[5*(i+1) + 3]);
				qs4[3] += parseInt(t[5*(i+1) + 4]);
			}
			if ((match = /^[A-Z]\+(\d+)([A-Z]+)/.exec(a)) != null) { // insertion
				alt.push(t[2] + match[2] + max_del);
			} else if ((match = /^[A-Z]-(\d+)([A-Z]+)/.exec(a)) != null) { // deletion
				alt.push(t[2] + max_del.substr(parseInt(match[1])));
			} else { // SNP
				alt.push(a);
			}
			q.push(qs4[2] + qs4[3]);
		}
		if (alt.length == 0) continue; // not a variant
		var alt_sum = 0;
		for (var i = 0; i < q.length; ++i) alt_sum += q[i];
		q.unshift(qs4[0] + qs4[1]);
		var gt;
		if (alt.length == 1 && qs4[0] + qs4[1] == 0) {
			gt = "1/1";
		} else {
			var max = -1, max2 = -1, max_i = -1, max2_i = -1;
			for (var i = 0; i < q.length; ++i) {
				if (max < q[i]) max2 = max, max2_i = max_i, max = q[i], max_i = i;
				else if (max2 < q[i]) max2 = q[i], max2_i = i;
			}
			if (max_i > max2_i) max_i ^= max2_i, max2_i ^= max_i, max_i ^= max2_i;
			gt = max_i + "/" + max2_i;
		}
		var info = qdp? 'DP4='+qs4.join(",")+';RD4='+dp4.join(",") : 'DP4='+dp4.join(",")+';QS4='+qs4.join(",");
		var out = [t[0], t[1], '.', t[2] + max_del, alt.join(","), alt_sum, '.', info, 'GT', gt];
		print(out.join("\t"));
	}

	buf.destroy();
	file.close();
}

function b8_cg2vcf(args)
{
	var file = args.length > 0? new File(args[0]) : new File();
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		var a = [];
		if (t.length < 5 || buf[0] < 48 || buf[0] > 57 || t[5] == 'no-call' || t[6] == 'ref') continue;
		if (t[7].length != t[8].length || t[7].length != t[9].length)
			for (var i = 7; i <= 9; ++i)
				if (t[i] != '?') t[i] = 'N' + t[i];
		var dp = t[24] == ''? 0 : parseInt(t[24]);
		var call = [], dr = 0, da = 0, gq = 1000000, qual = 0;
		for (var i = 0; i < 2; ++i) {
			var q = parseInt(t[10+i]);
			gq = gq < q? gq : q;
			if (t[i+8] == '?') continue;
			call.push(t[i+8]);
			var d = t[21+i] == ''? 0 : parseInt(t[21+i]);
			if (t[i+8] != t[7]) {
				da += i == 0 || t[9] != t[8]? d : 0;
				qual = qual > q? qual : q;
			} else dr += d;
		}
		var gt = './.', alt = '';
		if (call.length == 2) {
			if (call[0] == call[1]) {
				gt = call[0] == t[7]? '0/0' : '1/1';
				alt = call[0] == t[7]? '.' : call[0];
			} else if (call[0] == t[7] || call[1] == t[7]) { // biallelic
				gt = '0/1';
				alt = call[0] != t[7]? call[0] : call[1];
			} else { // triallelic
				gt = '1/2';
				alt = call.join(",");
			}
		} else alt = call[0];
		var o = [t[2], parseInt(t[3]) + 1, '.', t[7], alt, qual, '.', "DP="+dp+";"+"CGDP2="+dr+","+da, "GT:GQ", gt+":"+gq];
		print(o.join("\t"));
	}
	buf.destroy();
	file.close();
}

/***********************
 *** Main() function ***
 ***********************/

function main(args)
{
	if (args.length == 0) {
		print("\nUsage:    k8 bio8.js <command> [arguments]\n");
		print("Commands: qst1     vcf stats stratified by QUAL, one sample only");
		print("          upd1gt   update genotypes in single-sample VCFs");
		print("          deovlp   remove overlaps between variants");
		print("          plp2var  extract alleles from pileup output");
		print("          var2vcf  convert plp2var output to VCF");
		print("          cg2vcf   convert CG's masterVarBeta to VCF");
		print("          bedovlp  counts lines overlapping in a second BED");
		print("          bedcmpm  intersections between multiple sorted BED files");
		print("          cbs      circular binary segmentation");
		print("");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'upd1gt') b8_upd1gt(args);
	else if (cmd == 'qst1') b8_qst1(args);
	else if (cmd == 'deovlp') b8_deovlp(args);
	else if (cmd == 'plp2var') b8_plp2var(args);
	else if (cmd == 'var2vcf') b8_var2vcf(args);
	else if (cmd == 'cg2vcf') b8_cg2vcf(args);
	else if (cmd == 'bedovlp') b8_bedovlp(args);
	else if (cmd == 'bedcmpm') b8_bedcmpm(args);
	else if (cmd == 'cbs') b8_cbs(args);
	else warn("Unrecognized command");
}

main(arguments);
