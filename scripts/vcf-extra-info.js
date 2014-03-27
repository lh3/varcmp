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
		if (intv[j][1] > intv[i][0])
			intv[j][1] = intv[j][1] > intv[i][1]? intv[j][1] : intv[i][1];
		else intv[++j] = intv[i].slice(0);
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
			if (intv[i][1] > _b) return intv[i];
		return null;
	}
}

Math.lgamma = function(z) {
	var x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return Math.log(x) - 5.58106146679532777 - z + (z-0.5) * Math.log(z+6.5);
}

Math.fisher_exact = function(n11, n12, n21, n22)
{
	function lbinom(n, k) {
		if (k == 0 || n == k) return 0;
		return Math.lgamma(n+1) - Math.lgamma(k+1) - Math.lgamma(n-k+1);
	}

	function hypergeo(n11, n1_, n_1, n) {
		return Math.exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
	}

	function hypergeo_acc(n11, n1_, n_1, n, aux) {
		if (n1_ || n_1 || n) {
			aux.n11 = n11; aux.n1_ = n1_; aux.n_1 = n_1; aux.n = n;
		} else { // then only n11 changed; the rest fixed
			if (n11%11 && n11 + aux.n - aux.n1_ - aux.n_1) {
				if (n11 == aux.n11 + 1) { // incremental
					aux.p *= (aux.n1_ - aux.n11) / n11
						* (aux.n_1 - aux.n11) / (n11 + aux.n - aux.n1_ - aux.n_1);
					aux.n11 = n11;
					return aux.p;
				}
				if (n11 == aux.n11 - 1) { // incremental
					aux.p *= aux.n11 / (aux.n1_ - n11)
						* (aux.n11 + aux.n - aux.n1_ - aux.n_1) / (aux.n_1 - n11);
					aux.n11 = n11;
					return aux.p;
				}
			}
			aux.n11 = n11;
		}
		aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
		return aux.p;
	}

	var i, j, max, min;
	var p, q, left, right, two;
	var _aux = { n11:0, n1_:0, n_1:0, n:0, p:0. };
	var n1_, n_1, n;

	n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
	max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
	min = n1_ + n_1 - n;
	if (min < 0) min = 0; // min n11, for left tail
	if (min == max) return [1., 1., 1.]; // no need to do test
	q = hypergeo_acc(n11, n1_, n_1, n, _aux); // the probability of the current table
	// left tail
	p = hypergeo_acc(min, 0, 0, 0, _aux);
	for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) // loop until underflow
		left += p, p = hypergeo_acc(i, 0, 0, 0, _aux);
	--i;
	if (p < 1.00000001 * q) left += p;
	else --i;
	// right tail
	p = hypergeo_acc(max, 0, 0, 0, _aux);
	for (right = 0., j = max - 1; p < 0.99999999 * q; --j) // loop until underflow
		right += p, p = hypergeo_acc(j, 0, 0, 0, _aux);
	++j;
	if (p < 1.00000001 * q) right += p;
	else ++j;
	// two-tail
	two = left + right;
	if (two > 1.) two = 1.;
	// adjust left and right
	if (Math.abs(i - n11) < Math.abs(j - n11)) right = 1. - left + q;
	else left = 1.0 - right + q;
	return [two, left, right];
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
		t[1] = parseInt(t[1]); t[2] = parseInt(t[2]);
		var chr = t.shift();
		reg[chr].push(t);
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

var c, bed_lc = null, bed_cbs = null, bed_hwe = null, depth_only = false, hard_Q = null;
while ((c = getopt(arguments, 'dl:c:Q:h:')) != null) {
	if (c == 'l') bed_lc = read_bed(getopt.arg);
	else if (c == 'c') bed_cbs = read_bed(getopt.arg);
	else if (c == 'h') bed_hwe = read_bed(getopt.arg);
	else if (c == 'd') depth_only = true;
	else if (c == 'Q') hard_Q = parseInt(getopt.arg);
}

var file = arguments.length > getopt.ind? new File(arguments[getopt.ind]) : new File();
var buf = new Bytes();

var csharp = '#'.charCodeAt(0);
var lineno = 0;
while (file.readline(buf) >= 0) {
	++lineno;
	if (buf[0] == csharp) continue; // skip header lines
	var t = buf.toString().split("\t");
	if (hard_Q != null && parseFloat(t[5]) < hard_Q) continue;
	t[1] = parseInt(t[1]);
	// extract depth information
	var m = /DP=(\d+)/.exec(t[7]);
	var depth = m != null? parseInt(m[1]) : -1; // get read depth
	var dp4 = [], dp_ref = null, dp_alt = null, dp_alt_for = null, dp_alt_rev = null, FS = null;
	var m4 = [];
	if (/RD4=\d+,\d+,\d+,\d+/.test(t[7]) || /QS4=\d+,\d+,\d+,\d+/.test(t[7])) { // fermi2 + mpileup + plp2var + var2vcf-q
		m = /DP4=(\d+),(\d+),(\d+),(\d+)/.exec(t[7]);
		for (var j = 1; j <= 4; ++j) dp4[j-1] = parseInt(m[j]);
		dp_ref = dp4[0] + dp4[1];
		dp_alt = dp4[2] + dp4[3];
		depth = dp_ref + dp_alt;
		dp4 = []; // unset dp4[] as we should not filter based on this false DP4
	} else if ((m = /DP4=(\d+),(\d+),(\d+),(\d+)/.exec(t[7])) != null) { // samtools
		for (var j = 1; j <= 4; ++j) dp4[j-1] = parseInt(m[j]);
	} else if ((m = /CGDP2=(\d+),(\d+)/.exec(t[7])) != null) { // CG vcf converted by bio8.js cg2vcf
		dp_ref = parseInt(m[1]);
		dp_alt = parseInt(m[2]);
	} else if ((m4[0] = /SRF=(\d+)/.exec(t[7])) != null && (m4[1] = /SRR=(\d+)/.exec(t[7])) != null
			&& (m4[2] = /SAF=(\d+)/.exec(t[7])) != null && (m4[3] = /SAR=(\d+)/.exec(t[7])) != null) { // freebayes; in four separate tags
		for (var j = 0; j < 4; ++j) dp4[j] = parseInt(m4[j][1]);
	} else if (/GT:AD/.test(t[8])) { // GATK; no strand information
		m = /^\d\/\d:(\d+),(\d+)(,(\d+))?/.exec(t[9]);
		dp_ref = parseInt(m[1]);
		dp_alt = parseInt(m[2]) + (m[4] != null? parseInt(m[4]) : 0);
		if ((m = /[;\t]FS=([^\t;]+)/.exec(t[7])) != null)
			FS = Math.pow(10, -.1 * parseFloat(m[1]));
	} else if ((m = /NF=([\d,]+).*NR=([\d,]+).*TCF=(\d+).*TCR=(\d+)/.exec(t[7])) != null) { // Platypus; in four tags, but I don't really know what they mean...
		var m1 = m[1].split(",");
		var m2 = m[2].split(",");
		var max = -1, max_j = -1;
		if (m1.length != m2.length)
			throw Error("Platypus: INFO:NF and INFO:NR are inconsistent");
		for (var j = 0; j < m1.length; ++j) {
			m1[j] = parseInt(m1[j]); m2[j] = parseInt(m2[j]);
			var x = m1[j] + m2[j]
			if (max < x) max = x, max_j = j;
		}
		dp4[2] = m1[max_j]; dp4[3] = m2[max_j];
		dp4[0] = parseInt(m[3]) - dp4[2];
		dp4[1] = parseInt(m[4]) - dp4[3];
		if (dp4[0] < 0 || dp4[1] < 0) {
			warn("Platypus: negative reference coverage at line " + lineno);
			if (dp4[0] < 0) dp4[0] = 0;
			if (dp4[1] < 0) dp4[1] = 0;
		}
	} else warn("Unrecognized format at line " + lineno + ":\n" + buf.toString());

	if (dp4.length) dp_ref = dp4[0] + dp4[1], dp_alt_for = dp4[2], dp_alt_rev = dp4[3], dp_alt = dp_alt_for + dp_alt_rev;
	if (dp_alt != null && dp_ref != null) depth = dp_alt + dp_ref;
	if (depth == null) throw Error("Cannot get the depth information at line " + lineno + ":\n" + buf.toString());
	if (FS == null && dp4.length) FS = Math.fisher_exact(dp4[0], dp4[1], dp4[2], dp4[3])[0];

	var labels = ['qLC', 'qRD', 'qDP', 'qDS', 'qFS', 'qAB', 'qHW'];
	//             0      1      2      3      4      5      6
	var values = [null,   null,  null,  null,  null,  null,  0];
	if (FS != null) // fisher strand bias
		values[4] = Math.floor((FS < 1e-10? 0 : 4.343 * Math.log(FS)) + .499);
	if (dp_alt_for != null && dp_alt_rev != null) // double-strand support
		values[3] = dp_alt_for < dp_alt_rev? dp_alt_for : dp_alt_rev;
	if (dp_alt != null) // allele balance
		values[5] = depth == 0? 0 : Math.floor(100 * dp_alt / depth + .499);
	values[2] = -depth; // total depth
	if (bed_lc != null && bed_lc[t[0]] != null) { // low-complexity regions
		var x = bed_lc[t[0]](t[1] - 1, t[1] - 1 + t[3].length);
		values[0] = x == null? 0 : -(x[1] - x[0]);
	}
	if (bed_cbs != null && bed_cbs[t[0]] != null) { // regional depth
		var x = bed_cbs[t[0]](t[1] - 1, t[1] - 1 + t[3].length);
		values[1] = x != null? -x[3] : -depth;
	}
	if (bed_hwe != null && bed_hwe[t[0]] != null) { // HWE
		var x = bed_hwe[t[0]](t[1] - 1, t[1] - 1 + t[3].length);
		values[6] = x == null? 0 : -(x[1] - x[0]);
	}
	// write INFO
	var extra_info = [];
	for (var i = 0; i < values.length; ++i)
		if (values[i] != null)
			extra_info.push(labels[i] + '=' + values[i]);
	t[7] = extra_info.join(';') + ';' + t[7];
	if (depth_only) print(t[0], t[1], depth);
	else print(t.join("\t"));
}

buf.destroy();
file.close();
