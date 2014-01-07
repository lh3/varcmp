/***********************************************************
 *** WARNING: this script is NOT a generic VCF filter!!! ***
 ***********************************************************/

/******************
 *** From k8.js ***
 ******************/

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

/*********************
 *** Misc routines ***
 *********************/

function read_bed(fn)
{
	var f = new File(fn);
	var b = new Bytes();
	var reg = {}, idx = {};
	while (f.readline(b) >= 0) {
		var t = b.toString().split("\t");
		if (reg[t[0]] == null) reg[t[0]] = [];
		reg[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
	}
	for (var chr in reg)
		idx[chr] = intv_ovlp(reg[chr]);
	b.destroy();
	f.close();
	return idx;
}

/******************
 *** Main entry ***
 ******************/

var c, min_dp = 10, max_dp = -1, min_da = 5, min_r = .2, min_q = 0., min_sda = 1, str = null, mapa = null, dup = null;
while ((c = getopt(arguments, 'd:a:r:q:D:s:Q:t:m:c:')) != null) {
	if (c == 'd') min_dp = parseInt(getopt.arg);
	else if (c == 'a') min_da = parseInt(getopt.arg);
	else if (c == 'r') min_r = parseFloat(getopt.arg);
	else if (c == 'q') min_q = parseFloat(getopt.arg);
	else if (c == 's') min_sda = parseInt(getopt.arg);
	else if (c == 'D') max_dp = parseInt(getopt.arg);
	else if (c == 't') str = read_bed(getopt.arg);
	else if (c == 'm') mapa = read_bed(getopt.arg);
	else if (c == 'c') dup = read_bed(getopt.arg);
}

var file = arguments.length > getopt.ind? new File(arguments[getopt.ind]) : new File();
var buf = new Bytes();

var flag = { lowQual:0x01, lowDP:0x02, highDP:0x04, lowAB:0x08, lowDV:0x10, lowSDV:0x20, STR:0x40, mapa:0x80, dup:0x100 };

var c1 = '1'.charCodeAt(0);
var c9 = '9'.charCodeAt(0);
var lineno = 0;
while (file.readline(buf) >= 0) {
	++lineno;
	if (buf[0] < c1 || buf[0] > c9) continue; // skip header lines and non-autosomes
	var t = buf.toString().split("\t");
	t[1] = parseInt(t[1]);
	var qual = parseFloat(t[5]);
	var f = 0, nf = 0;
	if (qual < min_q) f |= flag.lowQual, ++nf;
	// extract depth information
	var m = /DP=(\d+)/.exec(t[7]);
	var depth = m != null? parseInt(m[1]) : -1; // get read depth
	var dr = -1, da = -1, daf = -1, dar = -1;
	if ((m = /DP4=(\d+),(\d+),(\d+),(\d+)/.exec(t[7])) != null) { // samtools
		dr = parseInt(m[1]) + parseInt(m[2]);
		daf = parseInt(m[3]);
		dar = parseInt(m[4]);
		da = daf + dar;
	} else if (/^GT:DP:RO:QR:AO/.test(t[8])) { // freebayes; no daf/dar
		m = /^\d\/\d:\d+:(\d+):\d+:(\d+)/.exec(t[9]);
		dr = parseInt(t[1]);
		da = parseInt(t[2]);
	} else if (/GT:AD/.test(t[8])) { // GATK; no daf/dar
		m = /^\d\/\d:(\d+),(\d+)/.exec(t[9]);
		dr = parseInt(t[1]);
		da = parseInt(t[2]);
	} else if (/^GT:GL:GOF:GQ:NR:NV/.test(t[8])) { // Platypus; depth inaccurate
		m = /NF=([\d,]+).*NR=([\d,]+)/.exec(t[7]);
		var m1 = m[1].split(",");
		var m2 = m[2].split(",");
		var s = t[9].split(":");
		var w = s[4].split(","); // spliting NR
		var max = -1, max_j = -1;
		if (m1.length != m2.length || m1.length != w.length)
			throw Error("Platypus: INFO:NF/NR and Sample:NV are inconsistent");
		for (var j = 0; j < m1.length && j < m2.length && j < w.length; ++j) {
			m1[j] = parseInt(m1[j]); m2[j] = parseInt(m2[j]); w[j] = parseInt(w[j]);
			var x = m1[j] + m2[j] + w[j];
			if (max < x) max = x, max_j = j;
		}
		daf = m1[max_j]; dar = m2[max_j];
		da = daf + dar;
		dr = w[max_j] - da;
	} else warn("Unrecognized format at line " + lineno + ":\n" + buf.toString());
	// set flag
	if (da >= 0 && dr >= 0) {
		if (depth < da + dr) depth = da + dr;
		if (da < min_da) f |= flag.lowDV, ++nf;
		if (da < depth * min_r) f |= flag.lowAB, ++nf;
		if (min_sda > 0 && daf >= 0 && dar >= 0 && (daf < min_sda || dar < min_sda)) f |= flag.lowSDV, ++nf;
	}
	if (depth < 0) throw Error("Cannot get the depth information at line " + lineno + ":\n" + buf.toString());
	if (depth < min_dp) f |= flag.lowDP, ++nf;
	if (max_dp > 0 && depth > max_dp) f |= flag.highDP, ++nf;
	// bed filters
	if (str != null && str[t[0]] != null) {
		if (str[t[0]](t[1] - 1, t[1] - 1 + t[3].length)) f |= flag.STR, ++nf;
	}
	if (mapa != null && mapa[t[0]] != null) {
		if (mapa[t[0]](t[1] - 1, t[1] - 1 + t[3].length)) f |= flag.mapa, ++nf;
	}
	if (dup != null && dup[t[0]] != null) {
		if (dup[t[0]](t[1] - 1, t[1] - 1 + t[3].length)) f |= flag.dup, ++nf;
	}
	// write FILTER
	var flt = "0x" + f.toString(16) + ';' + nf;
	if (t[6] == "PASS" || t[6] == ".") t[6] = flt;
	else t[6] = flt + ";" + t[6];
	print(t.join("\t"));
}

buf.destroy();
file.close();
