var file = arguments.length? new File(arguments[0]) : new File();
var buf = new Bytes();

while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (t.length >= 5) { // indel
		if (t[4] != 1) continue;
		t.length = 4;
		var s = t[3].split(",");
		for (var i = 0; i < s.length; ++i) {
			var l = parseInt(s[i]);
			if (l != 1 && l != -1)
			{
				print(t.join("\t"));
				break;
			}
		}
	} else { // snp
		if (t[2] != 1) continue;
		t.length = 2;
		print(t.join("\t"));
	}
}
