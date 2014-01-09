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

var c, min_alt_cnt = 0, min_alt_frac = 0., show_depth = false;
while ((c = getopt(arguments, 'c:df:')) != null)
	if (c == 'c') min_alt_cnt = parseInt(getopt.arg);
	else if (c == 'd') show_depth = true;
	else if (c == 'f') min_alt_frac = parseFloat(getopt.arg);

var file = arguments.length > getopt.ind? new File(arguments[getopt.ind]) : new File();
var buf = new Bytes();

while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (!/[ACGTacgt]/.test(t[4])) continue;
	t[2] = t[2].toUpperCase();
	var i = 0, j = 0;
	var alleles = {}, cnt = [], sum = 0, alt_cnt = 0, depth = 0;
	while (i < t[4].length && j < t[5].length) {
		var b = t[4].charAt(i);
		if (b == '$') { ++i; continue; }
		if (b == '^') { i += 2; continue; }
		if (b == '*') { ++i, ++j; continue; }

		++depth;
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

		var ci;
		if (alleles[a] == null) alleles[a] = cnt.length, cnt.push([a, 0, 0, 0, 0]);
		ci = alleles[a];
		++cnt[ci][forward? 1 : 2];
		cnt[ci][forward? 3 : 4] += q;
		sum += q;
		if (a != t[2]) ++alt_cnt;
	}
	if (min_alt_cnt > 0 && alt_cnt < min_alt_cnt) continue;
	if (min_alt_frac > 0 && alt_cnt < depth * min_alt_frac) continue;
	
	var out = [t[0], t[1], t[2], show_depth? depth : sum, cnt.length];
	for (var i = 0; i < cnt.length; ++i)
		for (var j = 0; j < 5; ++j)
			out.push(cnt[i][j]);
	print(out.join("\t"));
}

buf.destroy();
file.close();
