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

/*****************
 *** VCF stats ***
 *****************/

function b8_qst1(args)
{
	var c, b = 0.02, show_hdr = false, min_q = 0, use_GQ = false, het_num = null, tstv = null;
	while ((c = getopt(args, "ghb:q:H:T:")) != null)
		if (c == 'b') b = parseFloat(getopt.arg);
		else if (c == 'q') min_q = parseFloat(getopt.arg);
		else if (c == 'h') show_hdr = true;
		else if (c == 'g') use_GQ = true;
		else if (c == 'H') het_num = parseInt(getopt.arg);
		else if (c == 'T') tstv = parseFloat(getopt.arg);

	if (het_num != null && tstv != null) throw Error("-H and -T cannot be specified at the same time");

	var file = args.length > getopt.ind? new File(args[getopt.ind]) : new File();
	var buf = new Bytes();
	var a = [];
	while (file.readline(buf) >= 0) {
		if (buf[0] == 35) continue;
		var t = buf.toString().split("\t");
		if (t.length < 6) continue;
		t[5] = parseFloat(t[5]);
		if (t[5] < min_q) continue;
		t[3] = t[3].toUpperCase(); t[4] = t[4].toUpperCase();
		var s = t[4].split(","); // list of ALT alleles
		// find the genotype
		var is_het;
		var match = /^(\d+)(\/|\|)(\d+)/.exec(t[9]);
		if (match == null) { // special-casing samtools, as for single-sample, it effectively gives the genotype at FQ
			var m2 = /FQ=([^;\t]+)/.exec(t[7]);
			if (m2 == null) continue;
			is_het = parseFloat(m2[1]) > 0? 1 : 0;
		} else is_het = (parseInt(match[1]) != parseInt(match[3]))? 1 : 0;
		// find GQ
		if (use_GQ) {
			var i, u = t[8].split(":"), v = t[9].split(":");
			for (i = 0; i < u.length; ++i)
				if (u[i] == 'GQ') break;
			if (i == u.length) throw Error("No GQ");
			t[5] = parseFloat(v[i]);
		}
		// get CIGAR for freebayes
		var m3 = /CIGAR=([^;\t]+)/.exec(t[7]);
		var cigar = m3 != null? m3[1].split(",") : [];
		if (cigar.length && cigar.length != s.length) throw Error("Inconsistent ALT and CIGAR");
		// loop through each ALT allele
		for (var i = 0; i < s.length; ++i) {
			var type; // 1=ts, 2=tv, 3=mnp, 4=ins, 5=del
			if (t[3].length == 1 && s[i].length == 1) { // SNP
				type = ((t[3] == 'A' && s[i] == 'G') || (t[3] == 'G' && s[i] == 'A') || (t[3] == 'C' && s[i] == 'T') || (t[3] == 'T' && s[i] == 'C'))? 0 : 1;
			} else if (cigar.length) { // MNP or INDEL from freebayes
				var x = 0, y = 0;
				var m4, re = /(\d+)([MXID])/g;
				while ((m4 = re.exec(cigar[i])) != null) {
					var l = parseInt(m4[1]);
					if (m4[2] == 'X') {
						for (var j = 0; j < l; ++j) {
							u = t[3].substr(x+j, 1), v = s[i].substr(y+j, 1);
							type = ((u == 'A' && v == 'G') || (u == 'G' && v == 'A') || (u == 'C' && v == 'T') || (u == 'T' && v == 'C'))? 0 : 1;
							a.push([parseFloat(t[5]), type, is_het]);
						}
						x += l, y += l;
					} else if (m4[2] == 'M') x += l, y += l;
					else if (m4[2] == 'I') y += l;
					else if (m4[2] == 'D') x += l;
				}
				if (t[3].length == s[i].length) type = 2; // MNP
				else if (t[3].length < s[i].length) type = 3; // insertion
				else type = 4; // deletion
			} else { // MNP or INDEL from Platypus
				var l = t[3].length < s[i].length? t[3].length : s[i].length;
				for (var j = 0; j < l; ++j) { // decompose long variants
					var u = t[3].substr(j, 1), v = s[i].substr(j, 1);
					if (u != v) {
						type = ((u == 'A' && v == 'G') || (u == 'G' && v == 'A') || (u == 'C' && v == 'T') || (u == 'T' && v == 'C'))? 0 : 1;
						a.push([parseFloat(t[5]), type, is_het]);
					}
				}
				if (t[3].length == s[i].length) type = 2; // MNP
				else if (t[3].length < s[i].length) type = 3; // insertion
				else type = 4; // deletion
			}
			a.push([t[5], type, is_het]);
		}
	}
	buf.destroy(buf);
	file.close();

	a.sort(function(x,y) {return y[0]-x[0]});
	if (het_num != null || tstv != null) {
		var i, s = 0, ts = 0, tv = 0;
		for (i = 0; i < a.length; ++i) {
			if (a[i][1] == 2) continue; // ignore MNPs as they are decomposed to SNPs
			else if (a[i][1] == 0) ++ts;
			else if (a[i][1] == 1) ++tv;
			if (a[i][2]) ++s;
			if (het_num != null && s > het_num) break;
			if (tstv != null && ts/tv < tstv && i > a.length>>1) break;
		}
		print(a[i-1][0]);
	} else {
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

	while (file.readline(buf) >= 0) {
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
			print(buf);
			continue;
		}
		var u = s2[i].split(",");
		var n_alleles = t[4].split(",").length + 1;
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
	var file = args.length > 0? new File(args[0]) : new File();
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
		var alt = [], dp = [0, 0], q = [], q_ref = 0, qsum = 0;
		for (var i = 0; i < t[4]; ++i) {
			var a = t[5*(i+1)], match;
			dp[0] += parseInt(t[5*(i+1) + 1]);
			dp[1] += parseInt(t[5*(i+1) + 2]);
			var qs = parseInt(t[5*(i+1) + 3]) + parseInt(t[5*(i+1) + 4]);
			qsum += qs;
			if (a == t[2]) {
				q_ref = qs;
				continue; // identical to the reference
			}
			if ((match = /^[A-Z]\+(\d+)([A-Z]+)/.exec(a)) != null) { // insertion
				alt.push(t[2] + match[2] + max_del);
			} else if ((match = /^[A-Z]-(\d+)([A-Z]+)/.exec(a)) != null) { // deletion
				alt.push(t[2] + max_del.substr(parseInt(match[1])));
			} else { // SNP
				alt.push(a);
			}
			q.push(qs);
		}
		if (alt.length == 0) continue; // not a variant
		var alt_sum = 0;
		for (var i = 0; i < q.length; ++i) alt_sum += q[i];
		q.unshift(q_ref);
		var gt;
		if (alt.length == 1 && q_ref == 0) {
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
		var out = [t[0], t[1], '.', t[2] + max_del, alt.join(","), alt_sum, '.', 'DP2='+dp[0]+','+dp[1]+';QSUM='+qsum, 'GT', gt];
		print(out.join("\t"));
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
		print("");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'upd1gt') b8_upd1gt(args);
	else if (cmd == 'qst1') b8_qst1(args);
	else if (cmd == 'deovlp') b8_deovlp(args);
	else if (cmd == 'plp2var') b8_plp2var(args);
	else if (cmd == 'var2vcf') b8_var2vcf(args);
	else warn("Unrecognized command");
}

main(arguments);
