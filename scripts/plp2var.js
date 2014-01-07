var file = arguments.length > 0? new File(arguments[0]) : new File();
var buf = new Bytes();

while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (!/[ACGTacgt]/.test(t[4])) continue;
	var i = 0, j = 0;
	var alleles = {}, cnt = [], sum = 0;
	while (i < t[4].length && j < t[5].length) {
		var b = t[4].charAt(i);
		if (b == '$') { ++i; continue; }
		if (b == '^') { i += 2; continue; }
		if (b == '*') { ++i, ++j; continue; }

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
	}
	
	var out = [t[0], t[1], t[2], sum, cnt.length];
	for (var i = 0; i < cnt.length; ++i)
		for (var j = 0; j < 5; ++j)
			out.push(cnt[i][j]);
	print(out.join("\t"));
}

buf.destroy();
file.close();
