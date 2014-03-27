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

{
	var c, frac = 0.01, qthres = 30, tag = null, keep_LC = false, flt_DP = false, gap1_only = false;

	while ((c = getopt(arguments, 'Lt:b:q:D1')) != null) {
		if (c == 't') tag = getopt.arg;
		else if (c == 'b') frac = parseFloat(getopt.arg);
		else if (c == 'q') qthres = parseFloat(getopt.arg);
		else if (c == 'L') keep_LC = true;
		else if (c == 'D') flt_DP = true;
		else if (c == '1') gap1_only = true;
	}

	if (getopt.ind + 2 > arguments.length) {
		print("Usage: k8 qst1p.js <NA12878.vcf> <CHM1.vcf>");
		exit(1);
	}

	var a = [];

	function read_vcf(fn, which) {
		var file = new File(fn);
		var buf = new Bytes();
		var pat = tag != null? new RegExp(tag + '=(-?\\d+)') : null;
		var sum_dp = 0, n_dp = 0;
		var aa = [];
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			var c = b8_parse_vcf1(t);
			var m, q = null, depth = null;
			if (!keep_LC && tag != 'qLC') {
				m = /qLC=(-?\d+)/.exec(t[7]);
				if (m != null && m[1] < 0) continue;
			}
			m = /qDP=(-?\d+)/.exec(t[7]);
			if (m != null) depth = -parseInt(m[1]);
			if (pat != null) {
				var match = pat.exec(t[7]);
				if (match == null) continue;
				q = parseFloat(match[1]);
			}
			if (c[0][0] < qthres) continue;
			sum_dp += depth; ++n_dp;
			for (var i = 0; i < c.length; ++i) {
				var qual = q == null? c[i][0] : q;
				if (c[i][0] < qthres) continue;
				if (gap1_only) {
					if (c[i][1] != 3 && c[i][1] != 4) continue;
					if (c[i][4] != 1 && c[i][4] != -1) continue;
				}
				if (flt_DP) aa.push([qual, c[i][1], c[i][2], which, depth]);
				else a.push([qual, c[i][1], c[i][2], which, depth]);
			}
		}
		buf.destroy();
		file.close();
		if (flt_DP) {
			var avg_dp = sum_dp / n_dp;
			var max_dp = Math.sqrt(avg_dp) * 3. + avg_dp;
			warn(n_dp, avg_dp, max_dp);
			for (var i = 0; i < aa.length; ++i)
				if (aa[i][4] <= max_dp)
					a.push(aa[i]);
			aa = [];
		}
	}

	// read NA12878 calls
	warn('Reading the NA12878 calls... ');
	read_vcf(arguments[getopt.ind], 0);
	// read CHM1 calls
	warn('Reading the CHM1 calls...');
	read_vcf(arguments[getopt.ind + 1], 1);
	// stratify by quality
	warn('Stratifying by QUAL...');
	a.sort(function(x,y) {return y[0]-x[0]});
	var size = Math.floor(a.length * frac + 1.);
	var lastq = -1;
	var c = [[], []], ac = [[], []];
	for (var j = 0; j < 11; ++j) c[0][j] = ac[0][j] = c[1][j] = ac[1][j] = 0;
	for (var i = 0; i <= a.length; ++i) {
		if (i == a.length || (a[i][0] != lastq && c[0][0] + c[1][0] > a.length * frac)) {
			print('QR', lastq, ac[0].join('\t'), ac[1].join('\t'));
			if (i == a.length) break;
			for (var j = 0; j < 11; ++j) c[0][j] = c[1][j] = 0;
		}
		++c[a[i][3]][0];
		++ac[a[i][3]][0];
		++c[a[i][3]][a[i][1] * 2 + a[i][2] + 1];
		++ac[a[i][3]][a[i][1] * 2 + a[i][2] + 1];
		lastq = a[i][0];
	}
}
